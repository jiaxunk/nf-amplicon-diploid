"""Unit and integration tests for filter_sv_calls CLI."""
from __future__ import annotations

import subprocess
from pathlib import Path
import pysam
import pytest


def _create_mock_vcf(path: Path, records: list[str]) -> Path:
    header = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FILTER=<ID=GT,Description=\"Genotype filter (e.g. amplicon boundary dropoff)\">",
        "##FILTER=<ID=LowQual,Description=\"Low quality call\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the variant\">",
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of variant reads\">",
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reference reads\">",
        "##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">",
        "##contig=<ID=chr16,length=90338345>",
        "##contig=<ID=chr17,length=83257441>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
    ]
    lines = header + records
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _create_mock_bed(path: Path, intervals: list[tuple[str, int, int, str]]) -> Path:
    lines = [f"{chrom}\t{start}\t{end}\t{name}" for chrom, start, end, name in intervals]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def test_filter_sv_calls_retention_and_filtering(tmp_path: Path):
    target_bed = _create_mock_bed(
        tmp_path / "target.bed",
        [("chr16", 1000, 50000, "amplicon_target")],
    )

    records = [
        # 1. PASS inside target, good support and VAF -> RETAINED
        "chr16\t2000\tDEL_PASS_1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-10000;END=12000;SUPPORT=20;VAF=0.50\tGT:DV:DR:VAF\t0/1:20:20:0.50",
        # 2. GT inside target, good support and VAF -> RETAINED
        "chr16\t15000\tDEL_GT_1\tN\t<DEL>\t45\tGT\tSVTYPE=DEL;SVLEN=-5000;END=20000;SUPPORT=15;VAF=0.45\tGT:DV:DR:VAF\t0/1:15:18:0.45",
        # 3. LowQual inside target -> DISCARDED (filter != PASS and != GT)
        "chr16\t25000\tDEL_LOWQUAL\tN\t<DEL>\t10\tLowQual\tSVTYPE=DEL;SVLEN=-2000;END=27000;SUPPORT=10;VAF=0.30\tGT:DV:DR:VAF\t0/1:10:20:0.30",
        # 4. PASS inside target, low support (< 5) -> DISCARDED
        "chr16\t28000\tDEL_LOW_SUP\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;SVLEN=-3000;END=31000;SUPPORT=3;VAF=0.40\tGT:DV:DR:VAF\t0/1:3:5:0.40",
        # 5. PASS inside target, low VAF (< 0.2) -> DISCARDED
        "chr16\t32000\tDEL_LOW_VAF\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;SVLEN=-3000;END=35000;SUPPORT=10;VAF=0.08\tGT:DV:DR:VAF\t0/1:2:23:0.08",
        # 6. PASS outside target (chr17) -> DISCARDED
        "chr17\t2000\tDEL_OUTSIDE\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-4000;END=6000;SUPPORT=25;VAF=0.52\tGT:DV:DR:VAF\t0/1:25:23:0.52",
        # 7. PASS BND overlapping target at POS -> RETAINED
        "chr16\t40000\tBND_PASS_1\tN\tN[chr16:60000[\t60\tPASS\tSVTYPE=BND;SUPPORT=12;VAF=0.48\tGT:DV:DR:VAF\t0/1:12:13:0.48",
    ]
    raw_vcf = _create_mock_vcf(tmp_path / "raw.vcf", records)
    out_vcf = tmp_path / "filtered.vcf.gz"
    out_tsv = tmp_path / "filtered.tsv"

    cmd = [
        "workflow-glue", "filter_sv_calls",
        "--vcf", str(raw_vcf),
        "--target", str(target_bed),
        "--min-support", "5",
        "--min-vaf", "0.20",
        "-o", str(out_vcf),
        "--output-tsv", str(out_tsv),
    ]
    subprocess.run(cmd, check=True)

    # Check output VCF exists and has index
    assert out_vcf.exists()
    assert Path(str(out_vcf) + ".tbi").exists()

    # Read records from output VCF
    with pysam.VariantFile(str(out_vcf)) as vcf_in:
        kept_ids = [rec.id for rec in vcf_in]

    assert kept_ids == ["DEL_PASS_1", "DEL_GT_1", "BND_PASS_1"]

    # Check TSV output
    assert out_tsv.exists()
    with open(out_tsv) as f:
        lines = [line.strip().split("\t") for line in f if line.strip()]
    header = lines[0]
    assert "id" in header
    assert "filter" in header
    assert len(lines) == 4  # 1 header + 3 variants


def test_filter_sv_calls_empty_when_no_pass(tmp_path: Path):
    target_bed = _create_mock_bed(
        tmp_path / "target.bed",
        [("chr16", 1000, 50000, "amplicon_target")],
    )
    records = [
        "chr16\t25000\tDEL_LOWQUAL\tN\t<DEL>\t10\tLowQual\tSVTYPE=DEL;SVLEN=-2000;END=27000;SUPPORT=10;VAF=0.30\tGT:DV:DR:VAF\t0/1:10:20:0.30",
    ]
    raw_vcf = _create_mock_vcf(tmp_path / "raw.vcf", records)
    out_vcf = tmp_path / "filtered_empty.vcf.gz"

    cmd = [
        "workflow-glue", "filter_sv_calls",
        "--vcf", str(raw_vcf),
        "--target", str(target_bed),
        "-o", str(out_vcf),
    ]
    subprocess.run(cmd, check=True)

    assert out_vcf.exists()
    assert Path(str(out_vcf) + ".tbi").exists()
    with pysam.VariantFile(str(out_vcf)) as vcf_in:
        kept = list(vcf_in)
        assert len(kept) == 0


def test_filter_sv_calls_hba_dataset(tmp_path: Path):
    """Verify retention of both the 38.6 kb PASS deletion and 19.3 kb GT deletion in the HBA dataset."""
    hba_vcf = Path("example/barcode08_sniffles.vcf.gz")
    target_bed = Path("example/target.bed")
    if not hba_vcf.exists() or not target_bed.exists():
        pytest.skip("example data not found")

    out_vcf = tmp_path / "barcode08.nf_sv.vcf.gz"
    out_tsv = tmp_path / "barcode08.nf_sv.tsv"

    cmd = [
        "workflow-glue", "filter_sv_calls",
        "--vcf", str(hba_vcf),
        "--target", str(target_bed),
        "--min-support", "5",
        "--min-vaf", "0.1",
        "-o", str(out_vcf),
        "--output-tsv", str(out_tsv),
    ]
    subprocess.run(cmd, check=True)

    assert out_vcf.exists()
    assert Path(str(out_vcf) + ".tbi").exists()
    assert out_tsv.exists()

    with pysam.VariantFile(str(out_vcf)) as vcf_in:
        kept_records = {
            rec.id: (rec.chrom, rec.pos, rec.stop, list(rec.filter.keys()), rec.info.get("SVLEN"))
            for rec in vcf_in
        }

    # Verify Sniffles2.DEL.3S0 (38.6 kb PASS deletion)
    assert "Sniffles2.DEL.3S0" in kept_records
    del_3s0 = kept_records["Sniffles2.DEL.3S0"]
    assert "PASS" in del_3s0[3]
    assert abs(del_3s0[4]) > 38000

    # Verify Sniffles2.DEL.2S0 (19.3 kb GT deletion)
    assert "Sniffles2.DEL.2S0" in kept_records
    del_2s0 = kept_records["Sniffles2.DEL.2S0"]
    assert "GT" in del_2s0[3]
    assert abs(del_2s0[4]) > 19000

