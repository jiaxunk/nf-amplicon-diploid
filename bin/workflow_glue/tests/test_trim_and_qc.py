"""Tests for trim_and_qc multi-contig retention and quality control."""
from pathlib import Path
import subprocess
import pysam
import pandas as pd
import pytest


def _create_mock_consensus(fasta_path: Path, contigs: dict[str, str]):
    """Write contigs to a FASTA file."""
    with open(fasta_path, "w") as f:
        for name, seq in contigs.items():
            f.write(f">{name}\n{seq}\n")


def _create_mock_depth(depth_path: Path, contig_depths: dict[str, list[tuple[int, int, int]]]):
    """Write mosdepth per-base bed: ref, start, end, depth."""
    with open(depth_path, "w") as f:
        for ref, intervals in contig_depths.items():
            for start, end, depth in intervals:
                f.write(f"{ref}\t{start}\t{end}\t{depth}\n")


def _create_mock_flagstat(flagstat_path: Path, stats: dict[str, dict[str, int]]):
    """Write bamstats flagstat TSV file."""
    rows = []
    for contig, data in stats.items():
        rows.append({
            "sample": contig,
            "total": data.get("total", 100),
            "primary": data.get("primary", 90),
            "secondary": data.get("secondary", 5),
            "supplementary": data.get("supplementary", 5),
            "unmapped": data.get("unmapped", 0),
        })
    df = pd.DataFrame(rows).set_index("sample")
    df.to_csv(flagstat_path, sep="\t")


def test_trim_and_qc_retains_two_passing_contigs(tmp_path: Path):
    """When two contigs meet depth and primary ratio thresholds, both are retained and trimmed."""
    seq1 = "A" * 100
    seq2 = "C" * 120
    consensus_path = tmp_path / "draft.fasta"
    _create_mock_consensus(consensus_path, {"ctg1": seq1, "ctg2": seq2})

    depth_path = tmp_path / "depth.bed"
    _create_mock_depth(depth_path, {
        "ctg1": [(0, 10, 5), (10, 90, 100), (90, 100, 5)],
        "ctg2": [(0, 15, 5), (15, 105, 80), (105, 120, 5)],
    })

    flagstat_path = tmp_path / "flagstat.tsv"
    _create_mock_flagstat(flagstat_path, {
        "ctg1": {"total": 100, "primary": 95, "secondary": 5, "supplementary": 0},
        "ctg2": {"total": 100, "primary": 90, "secondary": 10, "supplementary": 0},
    })

    outdir_pass = tmp_path / "pass"
    outdir_fail = tmp_path / "fail"

    cmd = [
        "workflow-glue", "trim_and_qc",
        "--alias", "sampleA",
        "--asm-method", "miniasm",
        "--depth", str(depth_path),
        "--flagstat", str(flagstat_path),
        "--consensus", str(consensus_path),
        "--outdir-pass", str(outdir_pass),
        "--outdir-fail", str(outdir_fail),
        "--relative-depth-trim-threshold", "0.1",
        "--minimum-depth", "30",
        "--primary-threshold", "0.7",
        "--max-contigs", "2",
        "--qc-summary-tsv", "qc-summary.tsv",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}"

    assert outdir_pass.exists()
    assert (outdir_pass / "qc-summary.tsv").exists()
    qc_df = pd.read_csv(outdir_pass / "qc-summary.tsv", sep="\t", index_col=0)
    assert len(qc_df) == 2
    assert (qc_df["status"] == "passed").all()

    out_fasta = outdir_pass / "draft.fasta"
    assert out_fasta.exists()
    with pysam.FastxFile(out_fasta) as f:
        retained = list(f)

    assert len(retained) == 2
    retained_names = [r.name for r in retained]
    assert "sampleA_1" in retained_names
    assert "sampleA_2" in retained_names
    lengths = {r.name: len(r.sequence) for r in retained}
    assert lengths["sampleA_1"] == 80
    assert lengths["sampleA_2"] == 90


def test_trim_and_qc_single_passing_contig_named_alias(tmp_path: Path):
    """When only a single contig passes QC, it is retained and named exactly as the sample alias."""
    seq1 = "G" * 150
    consensus_path = tmp_path / "draft.fasta"
    _create_mock_consensus(consensus_path, {"ctg1": seq1})

    depth_path = tmp_path / "depth.bed"
    _create_mock_depth(depth_path, {"ctg1": [(0, 150, 50)]})

    flagstat_path = tmp_path / "flagstat.tsv"
    _create_mock_flagstat(flagstat_path, {
        "ctg1": {"total": 100, "primary": 92, "secondary": 8, "supplementary": 0},
    })

    outdir_pass = tmp_path / "pass"
    outdir_fail = tmp_path / "fail"

    cmd = [
        "workflow-glue", "trim_and_qc",
        "--alias", "barcode01",
        "--asm-method", "spoa",
        "--depth", str(depth_path),
        "--flagstat", str(flagstat_path),
        "--consensus", str(consensus_path),
        "--outdir-pass", str(outdir_pass),
        "--outdir-fail", str(outdir_fail),
        "--relative-depth-trim-threshold", "0.1",
        "--minimum-mean-depth", "30",
        "--primary-alignments-threshold", "0.7",
        "--max-de-novo-contigs", "2",
        "--qc-summary-tsv", "qc-summary.tsv",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}"

    assert outdir_pass.exists()
    with pysam.FastxFile(outdir_pass / "draft.fasta") as f:
        retained = list(f)
    assert len(retained) == 1
    assert retained[0].name == "barcode01"
    assert len(retained[0].sequence) == 150


def test_trim_and_qc_all_failed_outputs_to_fail_dir(tmp_path: Path):
    """When all contigs fail QC, output is written to outdir-fail and no consensus FASTA is written."""
    seq1 = "T" * 80
    consensus_path = tmp_path / "draft.fasta"
    _create_mock_consensus(consensus_path, {"ctg1": seq1})

    depth_path = tmp_path / "depth.bed"
    # Low depth 5 < 30
    _create_mock_depth(depth_path, {"ctg1": [(0, 80, 5)]})

    flagstat_path = tmp_path / "flagstat.tsv"
    _create_mock_flagstat(flagstat_path, {
        "ctg1": {"total": 100, "primary": 30, "secondary": 70, "supplementary": 0},
    })

    outdir_pass = tmp_path / "pass"
    outdir_fail = tmp_path / "fail"

    cmd = [
        "workflow-glue", "trim_and_qc",
        "--alias", "sampleFail",
        "--asm-method", "miniasm",
        "--depth", str(depth_path),
        "--flagstat", str(flagstat_path),
        "--consensus", str(consensus_path),
        "--outdir-pass", str(outdir_pass),
        "--outdir-fail", str(outdir_fail),
        "--relative-depth-trim-threshold", "0.1",
        "--minimum-depth", "30",
        "--primary-threshold", "0.7",
        "--max-contigs", "2",
        "--qc-summary-tsv", "qc-summary.tsv",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}"

    assert not outdir_pass.exists()
    assert outdir_fail.exists()
    assert (outdir_fail / "qc-summary.tsv").exists()
    assert not (outdir_fail / "draft.fasta").exists()
    qc_df = pd.read_csv(outdir_fail / "qc-summary.tsv", sep="\t", index_col=0)
    assert len(qc_df) == 1
    assert qc_df.iloc[0]["status"] == "failed"
    assert "low depth" in qc_df.iloc[0]["fail_reason"]
    assert "low primary ratio" in qc_df.iloc[0]["fail_reason"]
