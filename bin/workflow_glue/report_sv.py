"""Create the cohort-wide diploid structural variant HTML report."""
from __future__ import annotations

import base64
from pathlib import Path
from typing import List, Tuple

import dominate.tags as html_tags
import dominate.util as dom_util
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.stats import Stats
from ezcharts.layout.snippets.table import DataTable
import pandas as pd
import pysam

from .util import get_named_logger, wf_parser


def argparser():
    parser = wf_parser("report_sv")
    parser.add_argument("--report-fname", required=True)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--params", required=True)
    parser.add_argument("--wf-version", required=True)
    return parser


def _load_sv_variants(sample_dir: Path) -> pd.DataFrame:
    """Load structural variant calls from TSV or VCF in sample_dir."""
    cols = [
        "id", "chrom", "pos", "end", "sv_type", "sv_len",
        "filter", "gt", "support", "vaf", "qual", "phase_set"
    ]

    tsv_files = list(sample_dir.glob("**/*sv*.tsv")) + list(sample_dir.glob("*.tsv"))
    for tsv_path in tsv_files:
        if tsv_path.is_file():
            try:
                df = pd.read_csv(tsv_path, sep="\t")
                # Normalize column names
                rename_map = {
                    "sv_id": "id",
                    "svtype": "sv_type",
                    "svlen": "sv_len",
                    "genotype": "gt",
                    "start": "pos",
                }
                df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
                for c in cols:
                    if c not in df.columns:
                        df[c] = None
                return df[cols]
            except Exception:
                pass

    # Look for any .nf_sv.vcf.gz or .vcf.gz
    vcf_files = list(sample_dir.glob("**/*.nf_sv.vcf.gz")) or list(sample_dir.glob("**/*.vcf.gz"))
    if not vcf_files:
        return pd.DataFrame(columns=cols)

    vcf_file = vcf_files[0]
    rows = []
    try:
        with pysam.VariantFile(str(vcf_file)) as vcf:
            for rec in vcf:
                sv_type = rec.info.get("SVTYPE", "")
                if not sv_type and rec.alts:
                    alt = rec.alts[0]
                    if alt.startswith("<") and alt.endswith(">"):
                        sv_type = alt[1:-1].split(":")[0]
                end = rec.stop if rec.stop is not None else rec.info.get("END", rec.pos + 1)
                sv_len = rec.info.get("SVLEN", None)
                if isinstance(sv_len, (list, tuple)):
                    sv_len = sv_len[0]
                support = rec.info.get("SUPPORT", None)
                if support is None and rec.samples:
                    support = rec.samples[0].get("DV", None)
                vaf = rec.info.get("VAF", None)
                if vaf is None and rec.samples:
                    vaf = rec.samples[0].get("VAF", rec.samples[0].get("AF", None))
                if isinstance(vaf, (list, tuple)):
                    vaf = vaf[0]
                gt = rec.samples[0].get("GT", (None, None)) if rec.samples else (None, None)
                gt_str = "/".join(str(g) if g is not None else "." for g in gt)
                ps = rec.samples[0].get("PS", None) if rec.samples else None
                filters = list(rec.filter.keys())
                rows.append({
                    "id": rec.id or f"{rec.chrom}_{rec.pos}_{sv_type}",
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "end": end,
                    "sv_type": sv_type,
                    "sv_len": sv_len,
                    "filter": ";".join(filters) if filters else "PASS",
                    "gt": gt_str,
                    "support": support if support is not None else 0,
                    "vaf": round(float(vaf), 4) if vaf is not None else 0.0,
                    "qual": round(float(rec.qual), 1) if rec.qual is not None else None,
                    "phase_set": ps,
                })
        return pd.DataFrame(rows, columns=cols)
    except Exception:
        return pd.DataFrame(columns=cols)


def _find_sv_plots(sample_dir: Path) -> dict[str, Path]:
    """Find all rendered Cuban plot PNGs for the sample mapped by SV ID."""
    plots: dict[str, Path] = {}
    search_dirs = [sample_dir / "sv_plots", sample_dir]
    for d in search_dirs:
        if d.is_dir():
            for p in d.glob("*.png"):
                plots[p.stem] = p
    return plots


def _render_image_base64(img_path: Path):
    """Embed an image as base64 in HTML."""
    try:
        with open(img_path, "rb") as f:
            b64_data = base64.b64encode(f.read()).decode("utf-8")
        html_tags.img(
            src=f"data:image/png;base64,{b64_data}",
            style="max-width: 100%; height: auto; border: 1px solid #ced4da; border-radius: 6px; box-shadow: 0 2px 4px rgba(0,0,0,0.08); margin-top: 8px; margin-bottom: 24px;",
            alt=img_path.name,
        )
    except Exception as e:
        html_tags.p(f"Could not load image {img_path.name}: {e}")


def main(args):
    logger = get_named_logger("report_sv")
    samples = sorted(path for path in args.data.glob("*") if path.is_dir())
    report = labs.LabsReport(
        "Diploid Amplicon Structural Variant Report",
        "nf-amplicon-diploid",
        args.params,
        args.versions,
        args.wf_version,
    )

    with report.add_section("Structural Variant Overview", "Overview"):
        if not samples:
            html_tags.p("No samples with structural variant data were found.")
        else:
            summary = []
            for sample in samples:
                sv_df = _load_sv_variants(sample)
                n_total = len(sv_df)
                n_del = int((sv_df["sv_type"] == "DEL").sum()) if "sv_type" in sv_df else 0
                n_ins = int((sv_df["sv_type"] == "INS").sum()) if "sv_type" in sv_df else 0
                n_dup = int((sv_df["sv_type"] == "DUP").sum()) if "sv_type" in sv_df else 0
                n_inv = int((sv_df["sv_type"] == "INV").sum()) if "sv_type" in sv_df else 0
                n_bnd = int((sv_df["sv_type"] == "BND").sum()) if "sv_type" in sv_df else 0
                n_pass = int((sv_df["filter"].str.contains("PASS", na=False)).sum()) if "filter" in sv_df else 0
                n_gt = int((sv_df["filter"].str.contains("GT", na=False)).sum()) if "filter" in sv_df else 0

                summary.append({
                    "Sample": sample.name,
                    "Total SVs": n_total,
                    "Deletions": n_del,
                    "Insertions": n_ins,
                    "Duplications": n_dup,
                    "Inversions": n_inv,
                    "Breakends": n_bnd,
                    "PASS Filter": n_pass,
                    "GT Filter": n_gt,
                })

            html_tags.h4("Cohort Structural Variant Calling Summary")
            DataTable.from_pandas(pd.DataFrame(summary), use_index=False)

    if samples:
        with report.add_section("Sample Details", "Samples"):
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for sample in samples:
                    sv_df = _load_sv_variants(sample)
                    plots = _find_sv_plots(sample)

                    with tabs.add_dropdown_tab(sample.name):
                        html_tags.h4(f"Sample: {sample.name}")
                        n_total = len(sv_df)
                        n_del = int((sv_df["sv_type"] == "DEL").sum()) if "sv_type" in sv_df else 0
                        n_ins = int((sv_df["sv_type"] == "INS").sum()) if "sv_type" in sv_df else 0
                        n_dup = int((sv_df["sv_type"] == "DUP").sum()) if "sv_type" in sv_df else 0
                        n_inv = int((sv_df["sv_type"] == "INV").sum()) if "sv_type" in sv_df else 0
                        n_bnd = int((sv_df["sv_type"] == "BND").sum()) if "sv_type" in sv_df else 0

                        cards = [
                            (f"{n_total:,d}", "Total Retained SVs"),
                            (f"{n_del:,d}", "Deletions"),
                            (f"{n_ins:,d}", "Insertions"),
                            (f"{n_dup:,d}", "Duplications"),
                            (f"{n_inv:,d}", "Inversions"),
                            (f"{n_bnd:,d}", "Breakends"),
                        ]
                        Stats(columns=3, items=cards)
                        dom_util.raw("<br>")

                        if sv_df.empty:
                            html_tags.p("No structural variants met the filtering criteria.")
                        else:
                            # Table of calls
                            html_tags.h4("Structural Variant Calls")
                            display_df = sv_df.rename(columns={
                                "id": "SV ID",
                                "sv_type": "Type",
                                "chrom": "Chrom",
                                "pos": "Start",
                                "end": "End",
                                "sv_len": "Length (bp)",
                                "gt": "Genotype",
                                "filter": "Filter",
                                "support": "Support",
                                "vaf": "VAF",
                                "qual": "Quality",
                                "phase_set": "Phase Set",
                            })
                            DataTable.from_pandas(display_df, use_index=False)
                            dom_util.raw("<br>")

                            # Cuban Breakpoint Evidence Plots
                            html_tags.h4("Cuban Breakpoint Evidence Plots")
                            for _, r in sv_df.iterrows():
                                sv_id = str(r["id"])
                                sv_type = str(r["sv_type"])
                                chrom = str(r["chrom"])
                                pos = int(r["pos"])
                                end = int(r["end"])
                                filter_val = str(r["filter"])
                                vaf_val = r["vaf"]
                                sup_val = r["support"]
                                gt_val = str(r["gt"])

                                with html_tags.div(cls="card my-3 p-3"):
                                    html_tags.h5(
                                        f"{sv_id} — {sv_type} on {chrom}:{pos:,}-{end:,} "
                                        f"(Genotype: {gt_val}, Filter: {filter_val}, Support: {sup_val}, VAF: {vaf_val})"
                                    )
                                    plot_path = plots.get(sv_id)
                                    if not plot_path:
                                        # Try matching by chrom_pos_type
                                        alt_name = f"{chrom}_{pos}_{sv_type}"
                                        plot_path = plots.get(alt_name)

                                    if plot_path and plot_path.exists():
                                        _render_image_base64(plot_path)
                                    else:
                                        html_tags.p("Breakpoint evidence plot not available for this variant.")

    report.write(args.report_fname)
    logger.info("Structural variant report written to '%s'.", args.report_fname)
