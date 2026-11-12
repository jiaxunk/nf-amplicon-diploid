"""Behavioral tests for the cohort SNP report CLI and SNP data extraction."""
from __future__ import annotations

import gzip
from pathlib import Path
import subprocess

import pandas as pd
import pysam


def _create_mock_snp_data(
    base_dir: Path,
    alias: str,
    phased: bool = True,
    empty_variants: bool = False,
) -> Path:
    sdir = base_dir / alias
    sdir.mkdir(parents=True, exist_ok=True)

    if empty_variants:
        var_cols = [
            "target", "chrom", "position", "ref", "alt", "genotype",
            "filter", "type", "indel_length", "zygosity", "phase_set",
            "haplotype", "vaf", "depth", "qual",
        ]
        pd.DataFrame(columns=var_cols).to_csv(sdir / "variants.tsv", sep="\t", index=False)
    else:
        if phased:
            variants = [
                {
                    "target": "amplicon_1",
                    "chrom": "chr16",
                    "position": 1000,
                    "ref": "A",
                    "alt": "G",
                    "genotype": "1|0",
                    "filter": "PASS",
                    "type": "SNV",
                    "indel_length": 0,
                    "zygosity": "Heterozygous",
                    "phase_set": "1000",
                    "haplotype": "HP:1",
                    "vaf": 0.48,
                    "depth": 120,
                    "qual": 45.0,
                },
                {
                    "target": "amplicon_1",
                    "chrom": "chr16",
                    "position": 1050,
                    "ref": "CT",
                    "alt": "C",
                    "genotype": "0|1",
                    "filter": "PASS",
                    "type": "INDEL",
                    "indel_length": -1,
                    "zygosity": "Heterozygous",
                    "phase_set": "1000",
                    "haplotype": "HP:2",
                    "vaf": 0.52,
                    "depth": 115,
                    "qual": 50.0,
                },
                {
                    "target": "amplicon_1",
                    "chrom": "chr16",
                    "position": 1200,
                    "ref": "G",
                    "alt": "GA",
                    "genotype": "1|1",
                    "filter": "PASS",
                    "type": "INDEL",
                    "indel_length": 1,
                    "zygosity": "Homozygous",
                    "phase_set": "-",
                    "haplotype": "Hom ALT",
                    "vaf": 0.99,
                    "depth": 130,
                    "qual": 60.0,
                },
            ]
        else:
            variants = [
                {
                    "target": "amplicon_1",
                    "chrom": "chr16",
                    "position": 1000,
                    "ref": "A",
                    "alt": "G",
                    "genotype": "0/1",
                    "filter": "PASS",
                    "type": "SNV",
                    "indel_length": 0,
                    "zygosity": "Heterozygous",
                    "phase_set": "-",
                    "haplotype": "Unphased",
                    "vaf": 0.50,
                    "depth": 100,
                    "qual": 40.0,
                },
            ]
        pd.DataFrame(variants).to_csv(sdir / "variants.tsv", sep="\t", index=False)

    # Mock haplotype depth windows
    depth_rows = []
    for i in range(10):
        start = 900 + i * 50
        end = start + 50
        pos = 0.5 * (start + end)
        if phased:
            depth_rows.append({
                "target": "amplicon_1",
                "ref": "chr16",
                "start": start,
                "end": end,
                "pos": pos,
                "depth_total": 100.0,
                "depth_hp1": 52.0,
                "depth_hp2": 48.0,
            })
        else:
            depth_rows.append({
                "target": "amplicon_1",
                "ref": "chr16",
                "start": start,
                "end": end,
                "pos": pos,
                "depth_total": 100.0,
                "depth_hp1": 0.0,
                "depth_hp2": 0.0,
            })

    with gzip.open(sdir / "haplotype-depth.tsv.gz", "wt") as f:
        pd.DataFrame(depth_rows).to_csv(f, sep="\t", index=False)

    return sdir


def test_report_snp_handles_empty_cohort_and_writes_report(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    report = tmp_path / "nf-amplicon-snp-report.html"
    versions = tmp_path / "versions.txt"
    versions.write_text("clair3,1.0\n")
    params = tmp_path / "params.json"
    params.write_text("{}")

    result = subprocess.run(
        [
            "workflow-glue", "report_snp", "--report-fname", str(report),
            "--data", str(data), "--versions", str(versions),
            "--params", str(params), "--wf-version", "v1.0.0",
        ],
        capture_output=True, text=True, check=True,
    )

    assert result.returncode == 0
    assert report.is_file()
    assert "No samples" in report.read_text()


def test_report_snp_renders_phased_sample_with_plots_and_tables(tmp_path):
    data = tmp_path / "data"
    _create_mock_snp_data(data, "sample_phased", phased=True)
    report = tmp_path / "nf-amplicon-snp-report.html"
    versions = tmp_path / "versions.txt"
    versions.write_text("clair3,1.0\nwhatshap,1.7\n")
    params = tmp_path / "params.json"
    params.write_text('{"phased": true, "snp": true}')

    result = subprocess.run(
        [
            "workflow-glue", "report_snp", "--report-fname", str(report),
            "--data", str(data), "--versions", str(versions),
            "--params", str(params), "--wf-version", "v1.0.0",
        ],
        capture_output=True, text=True, check=True,
    )

    assert result.returncode == 0
    assert report.is_file()
    html_content = report.read_text()
    # Sample name in tab
    assert "sample_phased" in html_content
    # Target region name
    assert "amplicon_1" in html_content
    # Phased variant table columns/content
    assert "HP:1" in html_content
    assert "HP:2" in html_content
    # Summary sections / charts
    assert "Haplotype Coverage" in html_content
    assert "Variant Calls" in html_content


def test_report_snp_renders_unphased_sample(tmp_path):
    data = tmp_path / "data"
    _create_mock_snp_data(data, "sample_unphased", phased=False)
    report = tmp_path / "nf-amplicon-snp-report.html"
    versions = tmp_path / "versions.txt"
    versions.write_text("clair3,1.0\n")
    params = tmp_path / "params.json"
    params.write_text('{"phased": false, "snp": true}')

    result = subprocess.run(
        [
            "workflow-glue", "report_snp", "--report-fname", str(report),
            "--data", str(data), "--versions", str(versions),
            "--params", str(params), "--wf-version", "v1.0.0",
        ],
        capture_output=True, text=True, check=True,
    )

    assert result.returncode == 0
    assert report.is_file()
    html_content = report.read_text()
    assert "sample_unphased" in html_content
    assert "Unphased" in html_content


def test_report_snp_handles_zero_variants(tmp_path):
    data = tmp_path / "data"
    _create_mock_snp_data(data, "sample_zero", empty_variants=True)
    report = tmp_path / "nf-amplicon-snp-report.html"
    versions = tmp_path / "versions.txt"
    versions.write_text("clair3,1.0\n")
    params = tmp_path / "params.json"
    params.write_text('{"phased": true, "snp": true}')

    result = subprocess.run(
        [
            "workflow-glue", "report_snp", "--report-fname", str(report),
            "--data", str(data), "--versions", str(versions),
            "--params", str(params), "--wf-version", "v1.0.0",
        ],
        capture_output=True, text=True, check=True,
    )

    assert result.returncode == 0
    assert report.is_file()
    html_content = report.read_text()
    assert "sample_zero" in html_content
    assert "No small variants were called." in html_content


def test_extract_snp_data_cli(tmp_path):
    # Create small BAM with an aligned read
    bam_path = tmp_path / "aligned.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr16", "LN": 5000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out_bam:
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 100
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1000
        a.mapping_quality = 60
        a.cigar = ((0, 100),)
        a.set_tag("HP", 1)
        a.set_tag("PS", 1000)
        out_bam.write(a)
    pysam.index(str(bam_path))

    # Create target BED
    target_bed = tmp_path / "target.bed"
    target_bed.write_text("chr16\t900\t1200\tamplicon_test\n")

    # Create VCF
    vcf_text = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr16\t1020\t.\tA\tC\t50\tPASS\t.\tGT:DP:AF:PS\t1|0:30:0.45:1000
chr16\t1050\t.\tTC\tT\t60\tPASS\t.\tGT:DP:AF:PS\t0|1:28:0.55:1000
"""
    vcf_file = tmp_path / "variants.vcf"
    vcf_file.write_text(vcf_text)
    vcf_gz = tmp_path / "variants.vcf.gz"
    pysam.tabix_compress(str(vcf_file), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

    out_dir = tmp_path / "extracted_data"
    out_dir.mkdir()

    result = subprocess.run(
        [
            "workflow-glue", "extract_snp_data",
            "--bam", str(bam_path),
            "--vcf", str(vcf_gz),
            "--target", str(target_bed),
            "--out-dir", str(out_dir),
            "--number-depth-windows", "10",
        ],
        capture_output=True, text=True, check=True,
    )
    assert result.returncode == 0
    var_file = out_dir / "variants.tsv"
    assert var_file.is_file()
    df_var = pd.read_csv(var_file, sep="\t")
    assert len(df_var) == 2
    assert "HP:1" in df_var["haplotype"].values
    assert "HP:2" in df_var["haplotype"].values
    assert "SNV" in df_var["type"].values
    assert "INDEL" in df_var["type"].values
    assert -1 in df_var["indel_length"].values

    depth_file = out_dir / "haplotype-depth.tsv.gz"
    assert depth_file.is_file()
    df_depth = pd.read_csv(depth_file, sep="\t")
    assert len(df_depth) > 0
    assert "depth_total" in df_depth.columns
    assert "depth_hp1" in df_depth.columns
