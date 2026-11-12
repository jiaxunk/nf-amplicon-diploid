"""Create the cohort-wide diploid small-variant HTML report."""
from __future__ import annotations

import gzip
from pathlib import Path

import dominate.tags as html_tags
import dominate.util as dom_util
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.stats import Stats
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from .util import get_named_logger, wf_parser


def argparser():
    parser = wf_parser("report_snp")
    parser.add_argument("--report-fname", required=True)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--params", required=True)
    parser.add_argument("--wf-version", required=True)
    return parser


def _load_variants(sample_dir: Path) -> pd.DataFrame:
    """Load the normalized variant table produced for the sample."""
    path = sample_dir / "variants.tsv"
    cols = [
        "target", "chrom", "position", "ref", "alt", "genotype",
        "filter", "type", "indel_length", "zygosity", "phase_set",
        "haplotype", "vaf", "depth", "qual",
    ]
    if not path.exists():
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(path, sep="\t")
        for c in cols:
            if c not in df.columns:
                df[c] = None
        return df
    except Exception:
        return pd.DataFrame(columns=cols)


def _load_depth(sample_dir: Path) -> pd.DataFrame:
    """Load windowed haplotype depth data."""
    path = sample_dir / "haplotype-depth.tsv.gz"
    cols = ["target", "ref", "start", "end", "pos", "depth_total", "depth_hp1", "depth_hp2"]
    if not path.exists():
        plain_path = sample_dir / "haplotype-depth.tsv"
        if plain_path.exists():
            path = plain_path
        else:
            return pd.DataFrame(columns=cols)
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame(columns=cols)


def _render_summary_charts(variants: pd.DataFrame):
    """Render variant type, het/hom, and indel length summary charts."""
    if variants.empty:
        return

    html_tags.h5("Summary Distributions")
    with html_tags.div(cls="row"):
        # 1. Variant types
        with html_tags.div(cls="col-md-4"):
            html_tags.h6("Variant Types")
            type_counts = (
                variants["type"]
                .value_counts()
                .reset_index()
            )
            type_counts.columns = ["Variant Type", "Count"]
            p_type = ezc.barplot(data=type_counts, x="Variant Type", y="Count")
            p_type.title = {"text": "Variant Types"}
            EZChart(p_type, "epi2melabs")

        # 2. Zygosity (Het / Hom)
        with html_tags.div(cls="col-md-4"):
            html_tags.h6("Zygosity")
            zyg_counts = (
                variants["zygosity"]
                .value_counts()
                .reset_index()
            )
            zyg_counts.columns = ["Zygosity", "Count"]
            p_zyg = ezc.barplot(data=zyg_counts, x="Zygosity", y="Count")
            p_zyg.title = {"text": "Zygosity (Het / Hom)"}
            EZChart(p_zyg, "epi2melabs")

        # 3. Indel length distribution
        with html_tags.div(cls="col-md-4"):
            html_tags.h6("Indel Length Distribution")
            indels = variants[variants["type"] == "INDEL"].copy()
            if not indels.empty and "indel_length" in indels:
                len_counts = (
                    indels["indel_length"]
                    .value_counts()
                    .reset_index()
                )
                len_counts.columns = ["Net Length (bp)", "Count"]
                len_counts["Net Length (bp)"] = len_counts["Net Length (bp)"].astype(str)
                p_len = ezc.barplot(data=len_counts, x="Net Length (bp)", y="Count")
                p_len.title = {"text": "Indel Net Lengths"}
                EZChart(p_len, "epi2melabs")
            else:
                html_tags.p("No indels detected.")


def _render_target_tracks(target_name: str, df_depth: pd.DataFrame, df_var: pd.DataFrame):
    """Render haplotype coverage and VAF/Phase lollipop track for one target."""
    html_tags.h5(f"Target Region: {target_name}")

    # Coverage track
    target_depth = df_depth[df_depth["target"] == target_name].copy()
    if not target_depth.empty:
        # Check if phased depths exist
        has_hp = (
            "depth_hp1" in target_depth.columns
            and "depth_hp2" in target_depth.columns
            and (target_depth["depth_hp1"].sum() > 0 or target_depth["depth_hp2"].sum() > 0)
        )
        if has_hp:
            long_depth = []
            for _, r in target_depth.iterrows():
                long_depth.append({"pos": r["pos"], "depth": r["depth_total"], "Track": "Total"})
                long_depth.append({"pos": r["pos"], "depth": r["depth_hp1"], "Track": "HP:1"})
                long_depth.append({"pos": r["pos"], "depth": r["depth_hp2"], "Track": "HP:2"})
            depth_plot_df = pd.DataFrame(long_depth)
        else:
            depth_plot_df = pd.DataFrame({
                "pos": target_depth["pos"],
                "depth": target_depth["depth_total"],
                "Track": "Total",
            })

        p_cov = ezc.lineplot(
            data=depth_plot_df,
            x="pos",
            y="depth",
            hue="Track",
            marker=False,
        )
        p_cov.title = {"text": f"Haplotype Coverage — {target_name}"}
        p_cov._fig.xaxis.axis_label = "Genomic coordinate (bp)"
        p_cov._fig.yaxis.axis_label = "Depth (x)"
        EZChart(p_cov, "epi2melabs")
        dom_util.raw("<br>")

    # VAF / Phase lollipop track
    target_vars = df_var[df_var["target"] == target_name].copy()
    if not target_vars.empty:
        # Scatterplot
        p_vaf = ezc.scatterplot(
            data=target_vars,
            x="position",
            y="vaf",
            hue="haplotype",
        )
        p_vaf.title = {"text": f"Variant Allele Frequency & Phase — {target_name}"}
        p_vaf._fig.xaxis.axis_label = "Genomic coordinate (bp)"
        p_vaf._fig.yaxis.axis_label = "Variant Allele Frequency (VAF)"
        # Add lollipop stem segments from y0=0 to y1=vaf
        p_vaf._fig.segment(
            x0=target_vars["position"],
            y0=0,
            x1=target_vars["position"],
            y1=target_vars["vaf"],
            color="#888888",
            line_width=1.5,
            line_dash="dotted",
        )
        EZChart(p_vaf, "epi2melabs")
        dom_util.raw("<br>")


def main(args):
    logger = get_named_logger("report_snp")
    samples = sorted(path for path in args.data.glob("*") if path.is_dir())
    report = labs.LabsReport(
        "Diploid Amplicon Small Variant Report",
        "nf-amplicon-diploid",
        args.params,
        args.versions,
        args.wf_version,
    )

    with report.add_section("Variant Overview", "Overview"):
        if not samples:
            html_tags.p("No samples with small variant data were found.")
        else:
            summary = []
            for sample in samples:
                variants = _load_variants(sample)
                n_vars = len(variants)
                n_snv = int((variants["type"] == "SNV").sum()) if "type" in variants else 0
                n_indel = int((variants["type"] == "INDEL").sum()) if "type" in variants else 0
                n_het = int((variants["zygosity"] == "Heterozygous").sum()) if "zygosity" in variants else 0
                n_hom = int((variants["zygosity"] == "Homozygous").sum()) if "zygosity" in variants else 0
                n_phased = int((variants["haplotype"].isin(["HP:1", "HP:2", "Phased"])).sum()) if "haplotype" in variants else 0

                summary.append({
                    "Sample": sample.name,
                    "Total Variants": n_vars,
                    "SNVs": n_snv,
                    "Indels": n_indel,
                    "Heterozygous": n_het,
                    "Homozygous": n_hom,
                    "Phased (%)": round(100.0 * n_phased / n_vars, 1) if n_vars > 0 else 0.0,
                })

            html_tags.h4("Cohort Variant Calling Summary")
            DataTable.from_pandas(pd.DataFrame(summary), use_index=False)

    if samples:
        with report.add_section("Sample Details", "Samples"):
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in samples:
                    variants = _load_variants(sample)
                    depth = _load_depth(sample)

                    with tabs.add_dropdown_tab(sample.name):
                        html_tags.h4(f"Sample: {sample.name}")
                        n_vars = len(variants)
                        n_snv = int((variants["type"] == "SNV").sum()) if "type" in variants else 0
                        n_indel = int((variants["type"] == "INDEL").sum()) if "type" in variants else 0
                        n_het = int((variants["zygosity"] == "Heterozygous").sum()) if "zygosity" in variants else 0
                        n_hom = int((variants["zygosity"] == "Homozygous").sum()) if "zygosity" in variants else 0

                        cards = [
                            (f"{n_vars:,d}", "Total Small Variants"),
                            (f"{n_snv:,d}", "SNVs"),
                            (f"{n_indel:,d}", "Indels"),
                            (f"{n_het:,d}", "Heterozygous Sites"),
                            (f"{n_hom:,d}", "Homozygous Sites"),
                            (f"{n_het / max(1, n_hom):.2f}", "Het / Hom Ratio"),
                        ]
                        Stats(columns=3, items=cards)
                        dom_util.raw("<br>")

                        if variants.empty and depth.empty:
                            html_tags.p("No small variants were called.")
                        else:
                            # Haplotype Coverage & Lollipop Tracks
                            targets = []
                            if not depth.empty and "target" in depth.columns:
                                targets = list(depth["target"].unique())
                            elif not variants.empty and "target" in variants.columns:
                                targets = list(variants["target"].unique())

                            if targets:
                                html_tags.h4("Haplotype Coverage & Variant Tracks")
                                for t in targets:
                                    _render_target_tracks(t, depth, variants)

                            # Summary distributions
                            if not variants.empty:
                                _render_summary_charts(variants)
                                dom_util.raw("<br>")

                            # Phased variant call table
                            html_tags.h4("Variant Calls")
                            if variants.empty:
                                html_tags.p("No small variants were called.")
                            else:
                                display_cols = [
                                    "target", "chrom", "position", "ref", "alt",
                                    "genotype", "filter", "type", "phase_set",
                                    "haplotype", "vaf", "depth", "qual",
                                ]
                                present_cols = [c for c in display_cols if c in variants.columns]
                                display_df = variants[present_cols].rename(columns={
                                    "target": "Target",
                                    "chrom": "Chrom",
                                    "position": "Position",
                                    "ref": "Ref",
                                    "alt": "Alt",
                                    "genotype": "Genotype",
                                    "filter": "Filter",
                                    "type": "Type",
                                    "phase_set": "Phase Set",
                                    "haplotype": "Haplotype",
                                    "vaf": "VAF",
                                    "depth": "Depth",
                                    "qual": "Quality",
                                })
                                DataTable.from_pandas(display_df, use_index=False)

    report.write(args.report_fname)
    logger.info("SNP report written to '%s'.", args.report_fname)
