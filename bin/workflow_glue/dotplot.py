"""Parse PAF alignment records and generate consensus dot plots (static PNG & interactive ezcharts)."""
from __future__ import annotations

import itertools
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ezcharts.plots import BokehPlot
from ezcharts.plots.util import choose_palette

from .util import get_named_logger, wf_parser  # noqa: ABS101


def parse_paf(paf_path: Path | str, min_len: int = 0) -> list[dict[str, Any]]:
    """Parse PAF alignment records into structured dictionaries."""
    path = Path(paf_path)
    if not path.exists():
        return []

    records = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue

            q_name = parts[0]
            q_len = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            strand = parts[4]
            t_name = parts[5]
            t_len = int(parts[6])
            t_start = int(parts[7])
            t_end = int(parts[8])
            n_match = int(parts[9])
            block_len = int(parts[10])
            mapq = int(parts[11])

            if block_len < min_len:
                continue

            records.append({
                "query_name": q_name,
                "query_len": q_len,
                "q_start": q_start,
                "q_end": q_end,
                "strand": strand,
                "target_name": t_name,
                "target_len": t_len,
                "t_start": t_start,
                "t_end": t_end,
                "num_matches": n_match,
                "block_len": block_len,
                "mapq": mapq,
            })
    return records


def compute_segments(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compute (x0, y0, x1, y1) plot coordinates for alignment records.
    
    Target (reference) coordinates are plotted along the X-axis (t_start -> t_end).
    Query (consensus) coordinates are plotted along the Y-axis:
      - '+' strand: y0 = q_start, y1 = q_end (forward diagonal)
      - '-' strand: y0 = q_end, y1 = q_start (reverse diagonal)
    """
    segments = []
    for r in records:
        x0 = r["t_start"]
        x1 = r["t_end"]
        if r["strand"] == "+":
            y0 = r["q_start"]
            y1 = r["q_end"]
        else:
            y0 = r["q_end"]
            y1 = r["q_start"]

        segments.append({
            "x0": x0,
            "y0": y0,
            "x1": x1,
            "y1": y1,
            "query_name": r["query_name"],
            "query_len": r["query_len"],
            "target_name": r["target_name"],
            "target_len": r["target_len"],
            "strand": r["strand"],
            "block_len": r["block_len"],
            "mapq": r["mapq"],
        })
    return segments


def generate_dotplot_png(
    paf_path: Path | str,
    output_png: Path | str,
    alias: str = "",
    title: str | None = None,
    min_len: int = 0,
) -> Path:
    """Generate a static PNG dot plot using matplotlib."""
    out_path = Path(output_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    records = parse_paf(paf_path, min_len=min_len)
    segments = compute_segments(records)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
    plot_title = title or (f"Consensus Dot Plot — {alias}" if alias else "Consensus Dot Plot")

    if not segments:
        ax.text(
            0.5,
            0.5,
            "No alignments found",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=14,
            color="gray",
        )
        ax.set_title(plot_title, fontsize=14, fontweight="bold")
        ax.set_xlabel("Reference position (bp)")
        ax.set_ylabel("De Novo Consensus position (bp)")
    else:
        # Group segments by query name and strand for distinct styling
        palette = itertools.cycle(["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"])
        query_colors = {}
        for s in segments:
            q = s["query_name"]
            if q not in query_colors:
                query_colors[q] = next(palette)

        plotted_labels = set()
        max_x = max(max(s["x0"], s["x1"], s["target_len"]) for s in segments)
        max_y = max(max(s["y0"], s["y1"], s["query_len"]) for s in segments)
        target_names = sorted(list({s["target_name"] for s in segments}))

        for s in segments:
            q = s["query_name"]
            c = query_colors[q]
            linestyle = "-" if s["strand"] == "+" else "--"
            label_key = f"{q} ({s['strand']})"
            label = label_key if label_key not in plotted_labels else None
            if label:
                plotted_labels.add(label_key)

            ax.plot([s["x0"], s["x1"]], [s["y0"], s["y1"]], color=c, linestyle=linestyle, linewidth=1.5, label=label)

        ax.set_xlim(0, max_x)
        ax.set_ylim(0, max_y)
        ax.set_title(plot_title, fontsize=14, fontweight="bold")
        t_label = ", ".join(target_names) if target_names else "Reference"
        ax.set_xlabel(f"Reference position (bp) [{t_label}]")
        ax.set_ylabel("De Novo Consensus position (bp)")
        ax.grid(True, linestyle=":", alpha=0.6)
        if plotted_labels:
            ax.legend(loc="upper left", frameon=True)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def create_ezcharts_dotplot(
    paf_path: Path | str,
    alias: str = "",
    title: str | None = None,
    min_len: int = 0,
) -> BokehPlot | None:
    """Generate an interactive BokehPlot compatible with ezcharts for HTML embedding."""
    records = parse_paf(paf_path, min_len=min_len)
    segments = compute_segments(records)
    if not segments:
        return None

    plot_title = title or (f"Consensus Dot Plot — {alias}" if alias else "Consensus Dot Plot")
    max_x = max(max(s["x0"], s["x1"], s["target_len"]) for s in segments)
    max_y = max(max(s["y0"], s["y1"], s["query_len"]) for s in segments)
    target_names = sorted(list({s["target_name"] for s in segments}))
    t_label = ", ".join(target_names) if target_names else "Reference"

    bp = BokehPlot(
        width=650,
        height=600,
        x_range=(0, max_x),
        y_range=(0, max_y),
    )
    bp._fig.title.text = plot_title
    bp._fig.xaxis.axis_label = f"Reference coordinate (bp) [{t_label}]"
    bp._fig.yaxis.axis_label = "De Novo Consensus coordinate (bp)"

    palette = itertools.cycle(choose_palette())
    query_colors = {}
    for s in segments:
        q = s["query_name"]
        if q not in query_colors:
            query_colors[q] = next(palette)

    # Group segments for batch rendering in Bokeh
    grouped: dict[tuple[str, str], dict[str, list[Any]]] = {}
    for s in segments:
        key = (s["query_name"], s["strand"])
        if key not in grouped:
            grouped[key] = {"x0": [], "y0": [], "x1": [], "y1": []}
        grouped[key]["x0"].append(s["x0"])
        grouped[key]["y0"].append(s["y0"])
        grouped[key]["x1"].append(s["x1"])
        grouped[key]["y1"].append(s["y1"])

    for (q_name, strand), coords in grouped.items():
        color = query_colors[q_name]
        line_dash = "solid" if strand == "+" else "dashed"
        legend = f"{q_name} ({strand})"
        bp._fig.segment(
            x0=coords["x0"],
            y0=coords["y0"],
            x1=coords["x1"],
            y1=coords["y1"],
            color=color,
            line_width=2.5,
            line_dash=line_dash,
            legend_label=legend,
        )

    bp._fig.legend.location = "top_left"
    bp._fig.legend.click_policy = "hide"
    return bp


def argparser():
    """Argument parser for CLI entrypoint."""
    parser = wf_parser("dotplot")
    parser.add_argument(
        "--paf",
        required=True,
        type=Path,
        help="Path to PAF alignment file.",
    )
    parser.add_argument(
        "--output-png",
        required=True,
        type=Path,
        help="Path to output static PNG file.",
    )
    parser.add_argument(
        "--alias",
        default="",
        help="Sample alias for plot title.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom title for the dot plot.",
    )
    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=0,
        help="Minimum alignment block length to include.",
    )
    return parser


def main(args):
    """Main CLI entrypoint."""
    logger = get_named_logger("dotplot")
    logger.info(f"Generating dot plot for {args.paf} -> {args.output_png}")
    generate_dotplot_png(
        paf_path=args.paf,
        output_png=args.output_png,
        alias=args.alias,
        title=args.title,
        min_len=args.min_alignment_length,
    )
    logger.info("Done.")
