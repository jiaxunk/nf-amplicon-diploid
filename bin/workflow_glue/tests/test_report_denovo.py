"""Tests for cohort-wide de novo assembly HTML report generation."""
from pathlib import Path
import json
import pytest

MOCK_QC_SUMMARY = (
    "contig\tmethod\tlength\tmean_depth\tprimary\tsecondary\tsupplementary\tprimary_ratio\tstatus\tfail_reason\n"
    "draft_1\tminiasm\t1200\t150.5\t950\t20\t5\t0.974\tpassed\t\n"
    "draft_2\tminiasm\t1180\t142.0\t890\t15\t3\t0.980\tpassed\t\n"
)

MOCK_FAILED_QC_SUMMARY = (
    "contig\tmethod\tlength\tmean_depth\tprimary\tsecondary\tsupplementary\tprimary_ratio\tstatus\tfail_reason\n"
    "draft_1\tspoa\t400\t12.5\t50\t100\t0\t0.333\tfailed\tlow depth, low primary ratio\n"
)

MOCK_DEPTH_TSV = (
    "ref\tstart\tend\tdepth\n"
    "hap1\t0\t100\t140.0\n"
    "hap1\t100\t200\t160.0\n"
    "hap2\t0\t100\t135.0\n"
    "hap2\t100\t200\t150.0\n"
)

MOCK_PAF = (
    "hap1\t1200\t50\t1150\t+\tchr16\t5000\t1000\t2100\t1050\t1100\t60\ttp:A:P\n"
    "hap2\t1180\t40\t1140\t+\tchr16\t5000\t1000\t2100\t1040\t1100\t60\ttp:A:P\n"
)


def _setup_params_and_versions(tmp_path: Path):
    versions_file = tmp_path / "versions.txt"
    versions_file.write_text("workflow-glue,0.0.1\nminimap2,2.28\nmedaka,2.0.1\n")
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps({"analysis_mode": "denovo", "minimum_mean_depth": 30}))
    return versions_file, params_file


def test_report_denovo_with_paf_and_dotplot(tmp_path: Path):
    """Test generating de novo report with consensus depth, PAF dot plots, and QC summary."""
    from workflow_glue import report_denovo

    data_dir = tmp_path / "data"
    sample_dir = data_dir / "SAMPLE_A"
    sample_dir.mkdir(parents=True)

    (sample_dir / "qc-summary.tsv").write_text(MOCK_QC_SUMMARY)
    (sample_dir / "SAMPLE_A.consensus.fasta").write_text(">hap1\n" + "A" * 1200 + "\n>hap2\n" + "C" * 1180 + "\n")
    (sample_dir / "per-window-depth.tsv").write_text(MOCK_DEPTH_TSV)
    (sample_dir / "SAMPLE_A.denovo_to_ref.paf").write_text(MOCK_PAF)

    versions_file, params_file = _setup_params_and_versions(tmp_path)
    report_html = tmp_path / "nf-amplicon-denovo-report.html"

    parser = report_denovo.argparser()
    args = parser.parse_args([
        "--report-fname", str(report_html),
        "--data", str(data_dir),
        "--versions", str(versions_file),
        "--params", str(params_file),
        "--wf-version", "0.1.0",
    ])

    report_denovo.main(args)

    assert report_html.exists()
    content = report_html.read_text()
    assert "SAMPLE_A" in content
    assert "hap1" in content
    assert "hap2" in content
    assert "Consensus Dot Plot" in content


def test_report_denovo_multi_sample_pass_and_fail(tmp_path: Path):
    """Test generating report for cohort containing both passing and failing assemblies without PAF."""
    from workflow_glue import report_denovo

    data_dir = tmp_path / "data"

    s1_dir = data_dir / "SAMPLE_PASS"
    s1_dir.mkdir(parents=True)
    (s1_dir / "qc-summary.tsv").write_text(MOCK_QC_SUMMARY)
    (s1_dir / "SAMPLE_PASS.consensus.fasta").write_text(">consensus\n" + "G" * 1200 + "\n")

    s2_dir = data_dir / "SAMPLE_FAIL"
    s2_dir.mkdir(parents=True)
    (s2_dir / "qc-summary.tsv").write_text(MOCK_FAILED_QC_SUMMARY)

    versions_file, params_file = _setup_params_and_versions(tmp_path)
    report_html = tmp_path / "nf-amplicon-denovo-report.html"

    parser = report_denovo.argparser()
    args = parser.parse_args([
        "--report-fname", str(report_html),
        "--data", str(data_dir),
        "--versions", str(versions_file),
        "--params", str(params_file),
        "--wf-version", "0.1.0",
    ])

    report_denovo.main(args)

    assert report_html.exists()
    content = report_html.read_text()
    assert "SAMPLE_PASS" in content
    assert "SAMPLE_FAIL" in content
    assert "passed" in content.lower()
    assert "failed" in content.lower()
