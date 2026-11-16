"""Create cohort-wide diploid de novo assembly HTML report (`report_denovo`)."""
from __future__ import annotations

import itertools
from pathlib import Path
from typing import Any

import dominate.tags as html_tags
import dominate.util as dom_util
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.stats import Stats
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots.util import choose_palette
import pandas as pd

from .dotplot import create_ezcharts_dotplot  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


class DeNovoSampleData:
    """Holds Stage 1 QC, Stage 2 consensus, depth, and PAF dot plot data for one sample."""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.sample_alias = data_dir.name

        self.qc_summary = self._load_qc_summary(data_dir / "qc-summary.tsv")
        self.consensus_records = self._load_consensus_fasta()
        self.depth = self._load_depth()
        self.paf_path = self._find_file(["*.denovo_to_ref.paf", "*.paf"])
        self.dotplot_png = self._find_file(["*.dotplot.png", "dotplot.png"])

    @staticmethod
    def _load_qc_summary(path: Path) -> pd.DataFrame | None:
        if not path.exists():
            return None
        try:
            df = pd.read_csv(path, sep="\t")
            return df if not df.empty else None
        except Exception:
            return None

    def _find_file(self, patterns: list[str]) -> Path | None:
        for pat in patterns:
            matches = list(self.data_dir.glob(pat))
            if matches:
                return matches[0]
        return None

    def _load_consensus_fasta(self) -> list[dict[str, Any]]:
        """Parse contig headers and lengths from consensus FASTA files."""
        fasta_files = list(self.data_dir.glob("*.consensus.fasta"))
        if not fasta_files:
            fasta_files = list(self.data_dir.glob("*.fasta"))
        if not fasta_files:
            return []

        records = []
        for fpath in fasta_files:
            try:
                with open(fpath, "r", encoding="utf-8") as f:
                    current_header = None
                    seq_chunks = []
                    for line in f:
                        line = line.strip()
                        if line.startswith(">"):
                            if current_header:
                                full_seq = "".join(seq_chunks)
                                gc = (full_seq.count("G") + full_seq.count("C") + full_seq.count("g") + full_seq.count("c"))
                                gc_pct = round(100.0 * gc / max(1, len(full_seq)), 1)
                                records.append({
                                    "name": current_header,
                                    "length": len(full_seq),
                                    "gc_pct": gc_pct,
                                })
                            current_header = line[1:].split()[0]
                            seq_chunks = []
                        else:
                            seq_chunks.append(line)
                    if current_header:
                        full_seq = "".join(seq_chunks)
                        gc = (full_seq.count("G") + full_seq.count("C") + full_seq.count("g") + full_seq.count("c"))
                        gc_pct = round(100.0 * gc / max(1, len(full_seq)), 1)
                        records.append({
                            "name": current_header,
                            "length": len(full_seq),
                            "gc_pct": gc_pct,
                        })
            except Exception:
                continue
        return records

    def _load_depth(self) -> pd.DataFrame | None:
        for fname in ["per-window-depth.tsv.gz", "per-window-depth.tsv", "depth.tsv.gz", "depth.tsv"]:
            p = self.data_dir / fname
            if p.exists():
                try:
                    df = pd.read_csv(p, sep="\t")
                    if not df.empty:
                        return df
                except Exception:
                    continue
        return None

    def status(self) -> str:
        """Return 'Passed' if any contig passed Stage 1 QC, else 'Failed'."""
        if self.qc_summary is not None and not self.qc_summary.empty:
            if "status" in self.qc_summary.columns:
                if (self.qc_summary["status"] == "passed").any():
                    return "Passed"
                return "Failed"
        if self.consensus_records:
            return "Passed"
        return "Failed"

    def assembly_method(self) -> str:
        if self.qc_summary is not None and not self.qc_summary.empty and "method" in self.qc_summary.columns:
            methods = self.qc_summary["method"].dropna().unique()
            if len(methods) > 0:
                return str(methods[0])
        return "Unknown"

    def mean_depth(self) -> float:
        if self.depth is not None and not self.depth.empty and "depth" in self.depth.columns:
            return round(float(self.depth["depth"].mean()), 1)
        if self.qc_summary is not None and not self.qc_summary.empty and "mean_depth" in self.qc_summary.columns:
            passed = self.qc_summary.query("status == 'passed'")
            if not passed.empty:
                return round(float(passed["mean_depth"].mean()), 1)
            return round(float(self.qc_summary["mean_depth"].mean()), 1)
        return 0.0


def argparser():
    """Argument parser for report_denovo CLI."""
    parser = wf_parser("report_denovo")
    parser.add_argument(
        "--report-fname",
        required=True,
        help="Name of output de novo assembly HTML report file.",
    )
    parser.add_argument(
        "--data",
        required=True,
        type=Path,
        help="Directory with one sub-directory per sample containing de novo assembly outputs.",
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="CSV file containing program versions.",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="JSON file containing workflow parameters.",
    )
    parser.add_argument(
        "--wf-version",
        required=True,
        help="Version of the executed workflow.",
    )
    return parser


def main(args):
    """Main function to render cohort de novo assembly report."""
    logger = get_named_logger("report_denovo")
    logger.info("Generating de novo assembly report from %s", args.data)

    sample_dirs = sorted(p for p in args.data.glob("*") if p.is_dir())
    datasets = [DeNovoSampleData(d) for d in sample_dirs]

    report = labs.LabsReport(
        "Diploid Amplicon De Novo Assembly Report",
        "nf-amplicon-diploid",
        args.params,
        args.versions,
        args.wf_version,
    )

    # 1. Overview Section
    with report.add_section("De Novo Overview", "Overview"):
        html_tags.p(
            "Summary of two-stage diploid de novo amplicon assembly. Stage 1 constructs and evaluates "
            "length-heterozygous draft contigs (miniasm/spoa -> racon -> medaka -> trim_and_qc), "
            "and Stage 2 resolves haplotype consensus sequences (hap1 / hap2) with depth and dot plot validation."
        )

        if not datasets:
            html_tags.p("No samples with de novo assembly data found.")
        else:
            cohort_rows = []
            for d in datasets:
                status_str = d.status()
                method_str = d.assembly_method()
                n_drafts = len(d.qc_summary) if d.qc_summary is not None else 0
                hap_desc = ", ".join(f"{r['name']} ({r['length']:,d} bp)" for r in d.consensus_records) or "None"
                tot_len = sum(r["length"] for r in d.consensus_records)
                paf_status = "Available" if d.paf_path else "N/A"

                cohort_rows.append({
                    "Sample": d.sample_alias,
                    "Status": status_str,
                    "Method": method_str,
                    "Draft Contigs": n_drafts,
                    "Final Consensus Sequences": hap_desc,
                    "Total Length (bp)": f"{tot_len:,d}" if tot_len > 0 else "0",
                    "Mean Depth (x)": d.mean_depth(),
                    "Dot Plot vs Ref": paf_status,
                })

            html_tags.h4("Cohort De Novo Assembly Summary")
            DataTable.from_pandas(pd.DataFrame(cohort_rows), use_index=False)

    # 2. Sample Details Section
    if datasets:
        with report.add_section("Sample Details", "Samples"):
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for d in datasets:
                    with tabs.add_dropdown_tab(d.sample_alias):
                        html_tags.h4(f"Sample: {d.sample_alias}")

                        status = d.status()
                        n_contigs = len(d.consensus_records)
                        tot_len = sum(r["length"] for r in d.consensus_records)
                        m_depth = d.mean_depth()
                        method = d.assembly_method()

                        cards = [
                            (status, "Assembly Status"),
                            (str(n_contigs), "Haplotypes / Contigs"),
                            (f"{tot_len:,d} bp" if tot_len > 0 else "0 bp", "Total Consensus Length"),
                            (f"{m_depth:.1f}x", "Mean Assembly Depth"),
                            (method, "Assembly Method"),
                        ]
                        Stats(columns=5, items=cards)
                        dom_util.raw("<br>")

                        # Stage 1 QC Table
                        if d.qc_summary is not None and not d.qc_summary.empty:
                            html_tags.h5("Stage 1 Draft Assembly QC (trim_and_qc)")
                            DataTable.from_pandas(d.qc_summary, use_index=False)
                            dom_util.raw("<br>")

                        # Consensus Sequences Table
                        if d.consensus_records:
                            html_tags.h5("Final Consensus Sequences")
                            cons_df = pd.DataFrame(d.consensus_records).rename(columns={
                                "name": "Contig / Haplotype",
                                "length": "Length (bp)",
                                "gc_pct": "GC Content (%)",
                            })
                            DataTable.from_pandas(cons_df, use_index=False)
                            dom_util.raw("<br>")

                        # Coverage Depth Profile
                        if d.depth is not None and not d.depth.empty:
                            html_tags.h5("De Novo Consensus Depth of Coverage")
                            df_depth = d.depth.copy()
                            df_depth["pos"] = 0.5 * (df_depth["start"] + df_depth["end"])
                            df_depth["Contig"] = df_depth["ref"].astype(str)

                            palette = dict(zip(
                                df_depth["Contig"].unique(),
                                itertools.cycle(choose_palette()),
                            ))
                            plt = ezc.lineplot(
                                data=df_depth[["pos", "depth", "Contig"]].round(2),
                                x="pos",
                                y="depth",
                                hue="Contig",
                                palette=palette,
                                marker=False,
                            )
                            plt.title = {"text": f"Consensus coverage — {d.sample_alias}"}
                            plt._fig.xaxis.axis_label = "Consensus coordinate (bp)"
                            plt._fig.yaxis.axis_label = "Sequencing depth (x)"
                            EZChart(plt, "epi2melabs")
                            dom_util.raw("<br>")

                        # Consensus Dot Plot against Reference
                        if d.paf_path:
                            html_tags.h5("Consensus Dot Plot against Reference")
                            html_tags.p(
                                "Alignment of de novo consensus sequences against reference sequence "
                                "(minimap2 -x asm20 -c --eqx). Forward alignments are solid lines; reverse alignments are dashed."
                            )
                            ez_dotplot = create_ezcharts_dotplot(d.paf_path, alias=d.sample_alias)
                            if ez_dotplot is not None:
                                EZChart(ez_dotplot, "epi2melabs")
                                dom_util.raw("<br>")

    report.write(args.report_fname)
    logger.info("De novo assembly report written to '%s'.", args.report_fname)
