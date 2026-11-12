"""Extract normalized small variants and windowed haplotype depth for reporting."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pysam

from .util import get_named_logger, wf_parser


def argparser():
    parser = wf_parser("extract_snp_data")
    parser.add_argument("--bam", required=True, type=Path, help="BAM file (aligned or haplotagged)")
    parser.add_argument("--vcf", required=True, type=Path, help="VCF file with small variants")
    parser.add_argument("--target", required=True, type=Path, help="Target region BED file")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--number-depth-windows",
        default=100,
        type=int,
        help="Number of depth windows per target region",
    )
    return parser


def parse_target_bed(target_bed: Path) -> list[dict]:
    """Parse BED file into a list of targets with chrom, start, end, name."""
    targets = []
    with open(target_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            targets.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return targets


def find_target(chrom: str, pos: int, targets: list[dict]) -> str:
    """Find target name for a genomic coordinate."""
    for t in targets:
        if t["chrom"] == chrom and t["start"] <= pos <= t["end"]:
            return t["name"]
    return chrom


def extract_variants(vcf_path: Path, targets: list[dict]) -> list[dict]:
    """Extract variant records from VCF with phasing, VAF, depth, and type."""
    rows = []
    vcf = pysam.VariantFile(str(vcf_path))
    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alts = record.alts or ("",)
        alt = alts[0] if alts else ""
        qual = float(record.qual) if record.qual is not None else 0.0

        if not record.filter.keys():
            filt = "PASS"
        else:
            filt = ";".join(record.filter.keys())

        ref_len = len(ref)
        alt_len = len(alt)
        if ref_len == 1 and alt_len == 1:
            var_type = "SNV"
            indel_len = 0
        else:
            var_type = "INDEL"
            indel_len = alt_len - ref_len

        samples = list(record.samples.values())
        if samples:
            sample = samples[0]
            gt_tuple = sample.get("GT", (None, None))
            is_phased = getattr(sample, "phased", False)
            ps = sample.get("PS", None)
            dp = sample.get("DP", 0) or 0

            af = sample.get("AF", None)
            ad = sample.get("AD", None)
            if af is not None:
                if isinstance(af, (list, tuple)):
                    vaf = float(af[0]) if af else 0.0
                else:
                    vaf = float(af)
            elif ad is not None and sum(ad) > 0:
                vaf = float(ad[1]) / float(sum(ad))
            else:
                vaf = 0.0

            if None in gt_tuple:
                gt_str = "./."
                zygosity = "Unknown"
                haplotype = "Unphased"
            else:
                sep = "|" if is_phased else "/"
                gt_str = f"{gt_tuple[0]}{sep}{gt_tuple[1]}"
                if gt_tuple[0] == gt_tuple[1]:
                    if gt_tuple[0] == 0:
                        zygosity = "Homozygous Reference"
                        haplotype = "Hom REF"
                    else:
                        zygosity = "Homozygous"
                        haplotype = "Hom ALT"
                else:
                    zygosity = "Heterozygous"
                    if is_phased:
                        if gt_tuple == (1, 0):
                            haplotype = "HP:1"
                        elif gt_tuple == (0, 1):
                            haplotype = "HP:2"
                        else:
                            haplotype = "Phased"
                    else:
                        haplotype = "Unphased"

            phase_set = str(ps) if ps is not None else "-"
        else:
            gt_str = "-/-"
            zygosity = "Unknown"
            haplotype = "Unknown"
            phase_set = "-"
            vaf = 0.0
            dp = 0

        target_name = find_target(chrom, pos, targets)

        rows.append({
            "target": target_name,
            "chrom": chrom,
            "position": pos,
            "ref": ref,
            "alt": alt,
            "genotype": gt_str,
            "filter": filt,
            "type": var_type,
            "indel_length": indel_len,
            "zygosity": zygosity,
            "phase_set": phase_set,
            "haplotype": haplotype,
            "vaf": round(float(vaf), 3),
            "depth": int(dp),
            "qual": round(float(qual), 1),
        })

    return rows


def extract_haplotype_depth(
    bam_path: Path,
    targets: list[dict],
    n_windows: int,
) -> list[dict]:
    """Calculate windowed depth for total, HP:1, and HP:2 across each target."""
    rows = []
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    for target in targets:
        chrom = target["chrom"]
        start = target["start"]
        end = target["end"]
        t_len = max(1, end - start)
        k = min(n_windows, t_len)
        step = t_len / k

        for i in range(k):
            w_start = int(start + i * step)
            w_end = int(start + (i + 1) * step)
            w_len = max(1, w_end - w_start)
            mid = 0.5 * (w_start + w_end)

            tot = 0
            hp1 = 0
            hp2 = 0
            try:
                for col in bam.pileup(chrom, w_start, w_end, truncate=True):
                    for read in col.pileups:
                        if not (read.is_del or read.is_refskip):
                            tot += 1
                            try:
                                hp = read.alignment.get_tag("HP")
                                if hp == 1:
                                    hp1 += 1
                                elif hp == 2:
                                    hp2 += 1
                            except KeyError:
                                pass
            except ValueError:
                pass

            rows.append({
                "target": target["name"],
                "ref": chrom,
                "start": w_start,
                "end": w_end,
                "pos": round(mid, 1),
                "depth_total": round(tot / w_len, 2),
                "depth_hp1": round(hp1 / w_len, 2),
                "depth_hp2": round(hp2 / w_len, 2),
            })

    return rows


def main(args):
    logger = get_named_logger("extract_snp_data")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    targets = parse_target_bed(args.target)
    logger.info("Parsed %d target regions from %s", len(targets), args.target)

    var_cols = [
        "target", "chrom", "position", "ref", "alt", "genotype",
        "filter", "type", "indel_length", "zygosity", "phase_set",
        "haplotype", "vaf", "depth", "qual",
    ]
    rows_variants = extract_variants(args.vcf, targets)
    var_tsv = args.out_dir / "variants.tsv"
    with open(var_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=var_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_variants)
    logger.info("Extracted %d variants to %s", len(rows_variants), var_tsv)

    depth_cols = ["target", "ref", "start", "end", "pos", "depth_total", "depth_hp1", "depth_hp2"]
    rows_depth = extract_haplotype_depth(args.bam, targets, args.number_depth_windows)
    depth_tsv = args.out_dir / "haplotype-depth.tsv.gz"
    with gzip.open(depth_tsv, "wt", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=depth_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_depth)
    logger.info("Wrote windowed haplotype depth to %s", depth_tsv)
