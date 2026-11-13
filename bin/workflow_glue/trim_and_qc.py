"""Select contig(s) if multiple were created, do some QC, modify depth / stats files."""
from pathlib import Path
import sys

import pandas as pd
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("trim_and_qc")

    logger.info("Read input files.")

    flagstat = pd.read_csv(args.flagstat, sep="\t", index_col=0)

    depths = pd.read_csv(
        args.depth, sep="\t", header=None, names=["ref", "start", "end", "depth"]
    )

    with pysam.FastxFile(args.consensus) as f:
        fastq = {entry.name: entry for entry in f}

    # make sure we got the same seq IDs in flagstat and the FASTQ file
    seq_ids = [x for x in flagstat.index if x != "*"]
    if set(fastq.keys()) != set(seq_ids):
        logger.warning("Sequence IDs in FASTQ and in flagstat file don't agree.")

    # if all reads in the BAM are unmapped, `mosdepth` annoyingly outputs an empty
    # per-base depth file; we need to account for that
    if depths.empty:
        mean_depths = pd.Series(0, index=seq_ids)
    else:
        # calculate the mean depths of the candidate consensus contigs
        mean_depths = depths.groupby("ref").apply(
            lambda df: df.eval("(end - start) * depth").sum() / df["end"].iloc[-1]
        )

    # perform QC checks on each contig
    qc_stats = pd.DataFrame(
        columns=[
            "method",
            "length",
            "mean_depth",
            "primary",
            "secondary",
            "supplementary",
            "primary_ratio",
            "status",
            "fail_reason",
        ]
    )
    qc_stats.index.name = "contig"
    for contig in seq_ids:
        mean_depth = mean_depths.get(contig, 0.0)
        status = "passed"
        contig_len = len(fastq[contig].sequence) if contig in fastq else 0
        fail_reasons = []
        if mean_depth < args.minimum_depth:
            status = "failed"
            fail_reasons.append("low depth")
        primary = flagstat.loc[contig, "primary"] if contig in flagstat.index else 0
        total = flagstat.loc[contig, "total"] if contig in flagstat.index else 0
        secondary = flagstat.loc[contig, "secondary"] if contig in flagstat.index else 0
        supplementary = flagstat.loc[contig, "supplementary"] if contig in flagstat.index else 0
        primary_ratio = (primary / total) if total > 0 else 0.0

        if primary_ratio < args.primary_threshold:
            status = "failed"
            fail_reasons.append("low primary ratio")

        qc_stats.loc[contig] = [
            args.asm_method,
            contig_len,
            round(float(mean_depth), 3),
            primary,
            secondary,
            supplementary,
            round(float(primary_ratio), 3),
            status,
            ", ".join(fail_reasons),
        ]

    qc_passed = qc_stats.query("status == 'passed'")
    failed = qc_passed.empty

    # create the output directory
    outdir = Path(args.outdir_fail if failed else args.outdir_pass)
    outdir.mkdir(parents=True, exist_ok=True)

    qc_stats.to_csv(outdir / args.qc_summary_tsv, sep="\t")

    if failed:
        logger.warning(
            "Consensus generation failed for this sample. No consensus candidate "
            "passed the quality checks."
        )
        sys.exit(0)

    # select up to max_contigs contigs with the largest mean depth
    max_c = args.max_contigs if args.max_contigs is not None else 2
    selected_contigs = (
        qc_passed["mean_depth"].sort_values(ascending=False).index[:max_c].tolist()
    )
    logger.info(
        f"Selected {len(selected_contigs)} contig(s): {selected_contigs}. "
        f"Write modified depth and stats files to {outdir}."
    )

    out_records = []
    for idx, ctg in enumerate(selected_contigs):
        entry = fastq[ctg]
        ctg_depths = depths.query("ref == @ctg")
        if not ctg_depths.empty and ctg_depths["depth"].max() > 0:
            passing = ctg_depths.query(
                "depth > depth.max() * @args.relative_depth_trim_threshold"
            )
            if not passing.empty:
                trim_start = int(passing.iloc[0]["start"])
                trim_end = int(passing.iloc[-1]["end"])
            else:
                trim_start = 0
                trim_end = len(entry.sequence)
        else:
            trim_start = 0
            trim_end = len(entry.sequence)

        if len(selected_contigs) == 1:
            entry.name = args.alias
        else:
            entry.name = f"{args.alias}_{idx + 1}"

        entry.sequence = entry.sequence[trim_start:trim_end]
        if entry.quality is not None:
            entry.quality = entry.quality[trim_start:trim_end]
        out_records.append(str(entry))

    with open(outdir / args.consensus.name, "w") as f:
        for rec in out_records:
            f.write(rec + "\n")

    logger.info("Done")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("trim_and_qc")
    parser.add_argument(
        "--alias",
        type=str,
        help="Sample alias.",
        required=True,
    )
    parser.add_argument(
        "--asm-method",
        type=str,
        help="Assembly method (is added to QC summary TSV).",
        required=True,
    )
    parser.add_argument(
        "--depth",
        type=Path,
        help="Path to mosdepth per-base BED file.",
        required=True,
    )
    parser.add_argument(
        "--flagstat",
        type=Path,
        help="Path to bamstats flagstat TSV file.",
        required=True,
    )
    parser.add_argument(
        "--consensus",
        type=Path,
        help="Path to consensus FASTQ/FASTA file.",
        required=True,
    )
    parser.add_argument(
        "--outdir-pass",
        type=Path,
        help="Name of output directory if passing QC.",
        required=True,
    )
    parser.add_argument(
        "--outdir-fail",
        type=Path,
        help="Name of output directory if failing QC.",
        required=True,
    )
    parser.add_argument(
        "--relative-depth-trim-threshold",
        type=float,
        help=(
            "Trim ends of consensus sequence if the relative coverage drops below "
            "this value."
        ),
        required=True,
    )
    parser.add_argument(
        "--minimum-depth",
        "--minimum-mean-depth",
        dest="minimum_depth",
        type=int,
        help="Minimum required depth to trust a consensus sequence.",
        required=True,
    )
    parser.add_argument(
        "--primary-threshold",
        "--primary-alignments-threshold",
        dest="primary_threshold",
        type=float,
        help=(
            "Threshold for proportion of primary alignments (will write a "
            "warning for samples with lower fraction)."
        ),
        required=True,
    )
    parser.add_argument(
        "--max-contigs",
        "--max-de-novo-contigs",
        dest="max_contigs",
        type=int,
        default=2,
        help="Maximum number of passing contigs to retain (default: 2).",
    )
    parser.add_argument(
        "--qc-summary-tsv",
        type=str,
        help="Name of TSV file into which to write QC summary.",
        required=True,
    )
    return parser
