"""Create cohort-wide alignment & pre-processing HTML report (`report_alignment`)."""
from __future__ import annotations

import itertools
from pathlib import Path

import dominate.tags as html_tags
import dominate.util as dom_util
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import load_histogram
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.stats import Stats
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots.util import choose_palette
import pandas as pd
import pysam

from . import report_util as util  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


class AlignmentSampleData:
    """Holds pre-processing, alignment, and target depth data for one sample."""

    def __init__(self, data_dir: Path):
        """Parse available report files in `data_dir`."""
        self.data_dir = data_dir
        self.sample_alias = data_dir.name

        self.fastcat_per_file = self._read_tsv(data_dir / "per-file-stats.tsv")
        self.length_hist = (
            load_histogram(data_dir, "length")
            if (data_dir / "length.hist").exists()
            else pd.DataFrame()
        )
        self.qual_hist = (
            load_histogram(data_dir, "quality")
            if (data_dir / "quality.hist").exists()
            else pd.DataFrame()
        )

        post_trim_path = data_dir / "post-trim-per-file-stats.tsv"
        if not post_trim_path.exists():
            post_trim_path = data_dir / "porechopped-per-file-stats.tsv"
        self.post_trim_per_file = self._read_tsv(post_trim_path)

        self.has_fastcat = (
            self.fastcat_per_file is not None
            and not self.fastcat_per_file.empty
            and not self.length_hist.empty
        )

        self.bamstats = self._read_tsv(data_dir / "bamstats.tsv")
        self.flagstat = self._read_tsv(data_dir / "bamstats-flagstat.tsv")
        self.depth = self._read_tsv(data_dir / "per-window-depth.tsv.gz")

        self.has_alignment = (
            self.bamstats is not None
            and not self.bamstats.empty
            and self.flagstat is not None
            and not self.flagstat.empty
        )

    @staticmethod
    def _read_tsv(path: Path) -> pd.DataFrame | None:
        if not path.exists():
            return None
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            return None
        return df if not df.empty else None

    def summary_metrics(self) -> dict:
        """Compute summary statistics for cards and cohort tables."""
        raw_reads = 0
        filtered_reads = 0
        trimmed_reads = 0
        mean_len = 0.0
        mean_qual = 0.0

        if self.fastcat_per_file is not None:
            raw_reads = int(self.fastcat_per_file["n_seqs"].sum())
        if not self.length_hist.empty:
            filtered_reads = int(self.length_hist["count"].sum())
        if self.post_trim_per_file is not None:
            trimmed_reads = int(self.post_trim_per_file["n_seqs"].sum())
            n_bases = float(self.post_trim_per_file["n_bases"].sum())
            if trimmed_reads > 0:
                mean_len = n_bases / trimmed_reads
                mean_qual = float(
                    self.post_trim_per_file.eval("n_seqs * mean_quality").sum()
                    / trimmed_reads
                )

        mapped_reads = 0
        unmapped_reads = 0
        mean_acc = 0.0
        if self.flagstat is not None:
            mapped_df = self.flagstat[self.flagstat["ref"] != "*"]
            unmapped_df = self.flagstat[self.flagstat["ref"] == "*"]
            mapped_reads = int(mapped_df["primary"].sum()) if not mapped_df.empty else 0
            unmapped_reads = (
                int(unmapped_df["unmapped"].sum()) if not unmapped_df.empty else 0
            )
        if self.bamstats is not None and not self.bamstats.empty:
            mapped_bs = self.bamstats[self.bamstats["ref"] != "*"]
            if not mapped_bs.empty and "acc" in mapped_bs.columns:
                mean_acc = float(mapped_bs["acc"].mean())
            if trimmed_reads == 0:
                trimmed_reads = int(len(self.bamstats))
                if "read_length" in self.bamstats.columns:
                    mean_len = float(self.bamstats["read_length"].mean())

        mean_target_depth = 0.0
        n_targets_detected = 0
        if self.depth is not None and not self.depth.empty:
            group_col = "target" if "target" in self.depth.columns else "ref"
            spans = (self.depth["end"] - self.depth["start"]).clip(lower=1)
            total_span = float(spans.sum())
            if total_span > 0:
                mean_target_depth = float((spans * self.depth["depth"]).sum() / total_span)
            for _, grp in self.depth.groupby(group_col):
                if grp["depth"].sum() > 0:
                    n_targets_detected += 1

        return {
            "raw_reads": raw_reads,
            "filtered_reads": filtered_reads,
            "trimmed_reads": trimmed_reads,
            "mapped_reads": mapped_reads,
            "unmapped_reads": unmapped_reads,
            "mean_length": round(mean_len, 1),
            "mean_quality": round(mean_qual, 1),
            "mean_acc": round(mean_acc, 1),
            "mean_target_depth": round(mean_target_depth, 1),
            "n_targets_detected": n_targets_detected,
        }


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report_alignment")
    parser.add_argument(
        "--report-fname",
        required=True,
        help="Name of output alignment HTML report file.",
    )
    parser.add_argument(
        "--data",
        required=True,
        type=Path,
        help="Directory with one sub-directory per sample containing report inputs.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        help="Reference FASTA file.",
    )
    parser.add_argument(
        "--target-bed",
        type=Path,
        help="Target Region BED file.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        help="Sample sheet with metadata.",
    )
    parser.add_argument(
        "--downsampling-size",
        type=int,
        default=0,
        help="Downsampled number of reads (per sample).",
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="CSV with program versions.",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="JSON file containing the workflow parameters.",
    )
    parser.add_argument(
        "--wf-version",
        required=True,
        help="Version of the executed workflow.",
    )
    return parser


def _load_target_regions(
    target_bed: Path | None, reference: Path | None
) -> pd.DataFrame:
    """Load Target Region intervals from BED or synthesize from Reference Sequence."""
    rows = []
    if target_bed and target_bed.exists():
        with open(target_bed) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
                    rows.append(
                        {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "target": name,
                            "length_bp": end - start,
                        }
                    )
    elif reference and reference.exists():
        with pysam.FastxFile(str(reference)) as fx:
            for entry in fx:
                seq_len = len(entry.sequence)
                rows.append(
                    {
                        "chrom": entry.name,
                        "start": 0,
                        "end": seq_len,
                        "target": entry.name,
                        "length_bp": seq_len,
                    }
                )
    return pd.DataFrame(rows)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("report_alignment")

    metadata = None
    if args.sample_sheet and args.sample_sheet.exists():
        metadata = pd.read_csv(
            args.sample_sheet,
            index_col="alias",
            dtype={"alias": str, "barcode": str, "type": str},
        )

    datasets = sorted(
        [AlignmentSampleData(d) for d in args.data.glob("*") if d.is_dir()],
        key=lambda x: x.sample_alias,
    )
    if metadata is not None:
        by_alias = {d.sample_alias: d for d in datasets}
        ordered = [by_alias[a] for a in metadata.index if a in by_alias]
        remaining = [d for d in datasets if d.sample_alias not in metadata.index]
        datasets = ordered + remaining

    targets_df = _load_target_regions(args.target_bed, args.reference)

    report = labs.LabsReport(
        "Diploid Amplicon Alignment & Pre-processing Report",
        "nf-amplicon-diploid",
        args.params,
        args.versions,
        args.wf_version,
    )
    populate_alignment_report(
        report, datasets, targets_df, metadata, args.downsampling_size
    )
    report.write(args.report_fname)
    logger.info(f"Alignment report written to '{args.report_fname}'.")


def populate_alignment_report(
    report,
    datasets: list[AlignmentSampleData],
    targets_df: pd.DataFrame,
    metadata: pd.DataFrame | None,
    downsampling_size: int,
):
    """Populate sections of the alignment report."""
    if not datasets:
        with report.add_section("Warning: No samples found", "Warning"):
            html_tags.p("No valid sample data directories were found for reporting.")
        return

    with report.add_section("Overview & Per-Sample Summary", "Overview"):
        html_tags.p(
            "This report summarizes input read pre-processing, adapter/barcode trimming, "
            "Reference Sequence alignment, and Target Region depth of coverage across all "
            f"{len(datasets)} sample(s)."
        )
        if downsampling_size and downsampling_size > 0:
            with html_tags.div(cls="alert alert-info"):
                html_tags.p(
                    html_tags.b("Post-filter downsampling active: "),
                    f"Filtered reads were downsampled up to {downsampling_size:,d} reads "
                    "per sample prior to trimming and alignment.",
                )

        if not targets_df.empty:
            html_tags.h5("Configured Target Regions")
            display_targets = targets_df.rename(
                columns={
                    "chrom": "Reference contig",
                    "start": "Start (0-based)",
                    "end": "End",
                    "target": "Target Region",
                    "length_bp": "Length (bp)",
                }
            )
            DataTable.from_pandas(display_targets, use_index=False)

        html_tags.h5("At a Glance (Per-Sample Tabs)")
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in datasets:
                m = d.summary_metrics()
                with tabs.add_dropdown_tab(d.sample_alias):
                    if not d.has_fastcat:
                        html_tags.p(
                            html_tags.b("Aligned BAM Shortcut: "),
                            "Sample was ingressed from an already-aligned BAM matching "
                            "the Reference Sequence; FASTQ conversion and trimming were skipped.",
                        )
                    cards = [
                        (f"{m['trimmed_reads']:,d}", "Processed / Input Reads"),
                        (f"{m['mapped_reads']:,d}", "Mapped Primary Reads"),
                        (f"{m['unmapped_reads']:,d}", "Unmapped Reads"),
                        (f"{m['mean_length']:.1f} bp", "Mean Read Length"),
                        (f"{m['mean_acc']:.1f}%", "Mean Alignment Accuracy"),
                        (f"{m['mean_target_depth']:.1f}x", "Mean Target Region Depth"),
                    ]
                    Stats(columns=3, items=cards)

    # Pre-processing section
    with report.add_section("Pre-processing & Trimming", "Pre-processing"):
        html_tags.p(
            "Summary of read counts, base yield, length bounds, and quality scores "
            "across raw, length/quality-filtered, and downsampled/trimmed stages."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in datasets:
                with tabs.add_dropdown_tab(d.sample_alias):
                    if not d.has_fastcat:
                        html_tags.p(
                            f"Sample '{d.sample_alias}' used the Aligned BAM Shortcut "
                            "in Reference-Based Mode (pre-aligned BAM input)."
                        )
                        continue
                    stage_rows = []
                    if d.fastcat_per_file is not None:
                        n_raw = int(d.fastcat_per_file["n_seqs"].sum())
                        stage_rows.append(
                            {
                                "Stage": "Raw",
                                "Reads": f"{n_raw:,d}",
                                "Bases": util.si_format(
                                    int(d.fastcat_per_file["n_bases"].sum())
                                ),
                                "Min length": f"{int(d.fastcat_per_file['min_length'].min()):,d}",
                                "Max length": f"{int(d.fastcat_per_file['max_length'].max()):,d}",
                                "Mean quality": round(
                                    float(
                                        d.fastcat_per_file.eval(
                                            "n_seqs * mean_quality"
                                        ).sum()
                                        / n_raw
                                    )
                                    if n_raw > 0
                                    else 0.0,
                                    1,
                                ),
                            }
                        )
                    if not d.length_hist.empty:
                        n_filt = int(d.length_hist["count"].sum())
                        mean_q = (
                            float(
                                d.qual_hist.eval("(start + end) / 2 * count").sum()
                                / n_filt
                            )
                            if (n_filt > 0 and not d.qual_hist.empty)
                            else 0.0
                        )
                        stage_rows.append(
                            {
                                "Stage": "Filtered",
                                "Reads": f"{n_filt:,d}",
                                "Bases": util.si_format(
                                    int(d.length_hist.eval("start * count").sum())
                                ),
                                "Min length": f"{int(d.length_hist['start'].min()):,d}",
                                "Max length": f"{int(d.length_hist['start'].max()):,d}",
                                "Mean quality": round(mean_q, 1),
                            }
                        )
                    if d.post_trim_per_file is not None:
                        n_trim = int(d.post_trim_per_file["n_seqs"].sum())
                        stage_rows.append(
                            {
                                "Stage": "Downsampled & Trimmed",
                                "Reads": f"{n_trim:,d}",
                                "Bases": util.si_format(
                                    int(d.post_trim_per_file["n_bases"].sum())
                                ),
                                "Min length": f"{int(d.post_trim_per_file['min_length'].min()):,d}",
                                "Max length": f"{int(d.post_trim_per_file['max_length'].max()):,d}",
                                "Mean quality": round(
                                    float(
                                        d.post_trim_per_file.eval(
                                            "n_seqs * mean_quality"
                                        ).sum()
                                        / n_trim
                                    )
                                    if n_trim > 0
                                    else 0.0,
                                    1,
                                ),
                            }
                        )
                    DataTable.from_pandas(pd.DataFrame(stage_rows), use_index=False)

        fastcat_datasets = [d for d in datasets if d.has_fastcat]
        if fastcat_datasets:
            html_tags.h5("Filtered Read Length & Quality Distributions")
            fastcat.SeqSummary(
                tuple(d.data_dir for d in fastcat_datasets),
                sample_names=tuple(d.sample_alias for d in fastcat_datasets),
            )

    # Alignment summary section
    aligned_datasets = [d for d in datasets if d.has_alignment]
    if not aligned_datasets:
        return

    with report.add_section("Alignment Statistics", "Alignment"):
        cohort_rows = []
        for d in aligned_datasets:
            m = d.summary_metrics()
            total_aln_reads = m["mapped_reads"] + m["unmapped_reads"]
            map_rate = (
                100.0 * m["mapped_reads"] / total_aln_reads
                if total_aln_reads > 0
                else 0.0
            )
            cohort_rows.append(
                {
                    "Sample": d.sample_alias,
                    "Primary Mapped": m["mapped_reads"],
                    "Unmapped": m["unmapped_reads"],
                    "Mapping Rate (%)": round(map_rate, 1),
                    "Mean Read Length (bp)": m["mean_length"],
                    "Mean Accuracy (%)": m["mean_acc"],
                    "Mean Target Depth (x)": m["mean_target_depth"],
                }
            )
        html_tags.h5("Cohort Alignment Summary")
        DataTable.from_pandas(pd.DataFrame(cohort_rows), use_index=False)

        html_tags.h5("Per-Sample Alignment Flagstats")
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in aligned_datasets:
                with tabs.add_dropdown_tab(d.sample_alias):
                    DataTable.from_pandas(d.flagstat, use_index=False)

    # Target Region depth section
    depth_datasets = [
        d for d in aligned_datasets if d.depth is not None and not d.depth.empty
    ]
    if not depth_datasets:
        return

    with report.add_section("Target Region Depth of Coverage", "Depth"):
        html_tags.p(
            "Depth of coverage profiles across each Target Region (or Reference Sequence "
            "contig when no BED file was supplied)."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in depth_datasets:
                with tabs.add_dropdown_tab(d.sample_alias):
                    df = d.depth.copy()
                    label_col = "target" if "target" in df.columns else "ref"
                    df["pos"] = 0.5 * (df["start"] + df["end"])
                    df["Target Region"] = df[label_col].astype(str)
                    palette = dict(
                        zip(
                            df["Target Region"].unique(),
                            itertools.cycle(choose_palette()),
                        )
                    )
                    plt = ezc.lineplot(
                        data=df[["pos", "depth", "Target Region"]].round(2),
                        x="pos",
                        y="depth",
                        hue="Target Region",
                        palette=palette,
                        marker=False,
                    )
                    plt.title = {
                        "text": f"Target Region coverage — {d.sample_alias}"
                    }
                    plt._fig.xaxis.axis_label = "Genomic coordinate (bp)"
                    plt._fig.yaxis.axis_label = "Sequencing depth (x)"
                    EZChart(plt, "epi2melabs")
                    dom_util.raw("<br>")
                    summary_by_target = (
                        df.groupby("Target Region", as_index=False)
                        .agg(
                            ref=("ref", "first"),
                            start=("start", "min"),
                            end=("end", "max"),
                            mean_depth=("depth", "mean"),
                            min_depth=("depth", "min"),
                            max_depth=("depth", "max"),
                        )
                        .round(1)
                    )
                    DataTable.from_pandas(summary_by_target, use_index=False)
