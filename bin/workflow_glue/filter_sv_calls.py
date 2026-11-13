"""Target-scoped amplicon SV filtering per ADR-0002."""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import pysam

from .util import get_named_logger, wf_parser

_BND_MATE_RE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def argparser():
    parser = wf_parser("filter_sv_calls")
    parser.add_argument("--vcf", required=True, type=Path, help="Input Sniffles2 VCF file.")
    parser.add_argument("--target", required=True, type=Path, help="Target regions BED file.")
    parser.add_argument("--min-support", type=int, default=5, help="Minimum read support.")
    parser.add_argument("--min-vaf", type=float, default=0.1, help="Minimum variant allele frequency.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Filtered VCF output path (.vcf.gz).")
    parser.add_argument("--output-tsv", type=Path, default=None, help="Optional summary TSV output path.")
    return parser


def load_targets(target_bed: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Parse BED file into chromosome -> list of (start, end) intervals (1-based inclusive)."""
    targets: Dict[str, List[Tuple[int, int]]] = {}
    with open(target_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1]) + 1  # 0-based to 1-based
            end = int(parts[2])
            if chrom not in targets:
                targets[chrom] = []
            targets[chrom].append((start, end))
    return targets


def overlaps_target(
    chrom: str,
    pos: int,
    end: int,
    sv_type: str,
    alts: Tuple[str, ...],
    targets: Dict[str, List[Tuple[int, int]]],
) -> bool:
    """Check if SV span or either breakend locus for BND overlaps any target region."""
    if sv_type == "BND":
        # Check primary locus
        for t_start, t_end in targets.get(chrom, []):
            if t_start <= pos <= t_end:
                return True
        # Check mate locus
        if alts:
            m = _BND_MATE_RE.search(alts[0])
            if m:
                mate_chrom, mate_pos = m.group(1), int(m.group(2))
                for t_start, t_end in targets.get(mate_chrom, []):
                    if t_start <= mate_pos <= t_end:
                        return True
        return False

    # For non-BND SVs: interval overlap
    s = min(pos, end)
    e = max(pos, end)
    for t_start, t_end in targets.get(chrom, []):
        if max(s, t_start) <= min(e, t_end):
            return True
    return False


def main(args):
    logger = get_named_logger("filter_sv_calls")
    targets = load_targets(args.target)
    vcf_in = pysam.VariantFile(str(args.vcf))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Use wz mode for bgzipped VCF
    mode = "wz" if output_path.name.endswith(".gz") else "w"
    vcf_out = pysam.VariantFile(str(output_path), mode=mode, header=vcf_in.header)

    retained_records = []
    tsv_rows = []

    for record in vcf_in:
        # 1. Filter status: retain PASS or GT
        filters = list(record.filter.keys())
        # Sniffles2 records may have ['PASS'], ['GT'], or be unfiltered
        pass_or_gt = (
            len(filters) == 0
            or "PASS" in filters
            or "GT" in filters
        )
        if not pass_or_gt:
            continue

        # 2. Support
        support = record.info.get("SUPPORT", None)
        if support is None and record.samples:
            sample_data = record.samples[0]
            support = sample_data.get("DV", None)
        if support is None:
            support = 0
        if support < args.min_support:
            continue

        # 3. VAF
        vaf = record.info.get("VAF", None)
        if vaf is None and record.samples:
            sample_data = record.samples[0]
            vaf = sample_data.get("VAF", sample_data.get("AF", None))
        if isinstance(vaf, (list, tuple)):
            vaf = vaf[0]
        if vaf is None:
            vaf = 0.0
        if vaf < args.min_vaf:
            continue

        # 4. SV Type & End position
        sv_type = record.info.get("SVTYPE", "")
        if not sv_type and record.alts:
            alt = record.alts[0]
            if alt.startswith("<") and alt.endswith(">"):
                sv_type = alt[1:-1].split(":")[0]

        end = record.stop if record.stop is not None else record.info.get("END", record.pos + 1)
        if end is None:
            end = record.pos + 1

        # 5. Overlap with Target Region
        if not overlaps_target(record.chrom, record.pos, end, sv_type, record.alts, targets):
            continue

        vcf_out.write(record)
        retained_records.append(record)

        # Collect row for TSV summary
        gt = record.samples[0].get("GT", (None, None)) if record.samples else (None, None)
        gt_str = "/".join(str(g) if g is not None else "." for g in gt)
        sv_len = record.info.get("SVLEN", None)
        if isinstance(sv_len, (list, tuple)):
            sv_len = sv_len[0]
        if sv_len is None and sv_type != "BND":
            sv_len = end - record.pos

        qual = record.qual
        ps = record.samples[0].get("PS", None) if record.samples else None

        tsv_rows.append({
            "id": record.id or f"{record.chrom}_{record.pos}_{sv_type}",
            "chrom": record.chrom,
            "pos": record.pos,
            "end": end,
            "sv_type": sv_type,
            "sv_len": sv_len,
            "filter": ";".join(filters) if filters else "PASS",
            "gt": gt_str,
            "support": support,
            "vaf": round(float(vaf), 4),
            "qual": round(float(qual), 1) if qual is not None else None,
            "phase_set": ps,
        })

    vcf_out.close()
    vcf_in.close()

    # Index bgzipped VCF if applicable
    if str(output_path).endswith(".gz"):
        try:
            pysam.tabix_index(str(output_path), preset="vcf", force=True)
        except Exception as e:
            logger.warning("pysam.tabix_index failed (%s); falling back to tabix CLI", e)
            os.system(f"tabix -f -p vcf '{output_path}'")

    if args.output_tsv:
        df = pd.DataFrame(tsv_rows, columns=[
            "id", "chrom", "pos", "end", "sv_type", "sv_len",
            "filter", "gt", "support", "vaf", "qual", "phase_set"
        ])
        df.to_csv(args.output_tsv, sep="\t", index=False)

    logger.info("Retained %d SVs matching criteria out of input.", len(retained_records))
