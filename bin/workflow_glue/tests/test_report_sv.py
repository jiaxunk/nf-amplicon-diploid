"""Unit and integration tests for report_sv CLI."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import pytest


TINY_PNG = (
    b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
    b"\x08\x06\x00\x00\x00\x1f\x15c4\x00\x00\x00\nIDATx\x9cc\x00\x01\x00"
    b"\x00\x05\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82"
)


def _setup_mock_sv_cohort(root: Path):
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Sample 1: barcode01 has 2 SVs and 2 plots
    s1_var_dir = data_dir / "barcode01" / "variants"
    s1_var_dir.mkdir(parents=True, exist_ok=True)
    s1_plot_dir = data_dir / "barcode01" / "sv_plots"
    s1_plot_dir.mkdir(parents=True, exist_ok=True)

    tsv_content_1 = (
        "sv_id\tchrom\tstart\tend\tsvlen\tsvtype\tfilter\tgenotype\tsupport\tdepth\tvaf\n"
        "Sniffles2.DEL.1\tchr16\t21000\t59600\t-38600\tDEL\tPASS\t0/1\t18\t22\t0.45\n"
        "Sniffles2.DEL.2\tchr16\t40300\t59600\t-19300\tDEL\tGT\t0/1\t12\t20\t0.38\n"
    )
    (s1_var_dir / "barcode01.nf_sv.tsv").write_text(tsv_content_1)
    (s1_plot_dir / "Sniffles2.DEL.1.png").write_bytes(TINY_PNG)
    (s1_plot_dir / "Sniffles2.DEL.2.png").write_bytes(TINY_PNG)

    # Sample 2: barcode02 has 0 SVs
    s2_var_dir = data_dir / "barcode02" / "variants"
    s2_var_dir.mkdir(parents=True, exist_ok=True)
    s2_plot_dir = data_dir / "barcode02" / "sv_plots"
    s2_plot_dir.mkdir(parents=True, exist_ok=True)
    tsv_content_2 = "sv_id\tchrom\tstart\tend\tsvlen\tsvtype\tfilter\tgenotype\tsupport\tdepth\tvaf\n"
    (s2_var_dir / "barcode02.nf_sv.tsv").write_text(tsv_content_2)

    # params and versions
    params_path = root / "params.json"
    params_path.write_text(json.dumps({
        "target": "target.bed",
        "reference": "ref.fasta",
        "min_sv_length": 30,
        "min_sv_support": 5,
        "min_sv_vaf": 0.20,
    }))
    versions_path = root / "versions.txt"
    versions_path.write_text("workflow,0.1.0\nsniffles,2.6.3\ncuban,0.1.0\n")

    return data_dir, params_path, versions_path


def test_report_sv_generation(tmp_path: Path):
    data_dir, params_path, versions_path = _setup_mock_sv_cohort(tmp_path)
    out_html = tmp_path / "nf-amplicon-sv-report.html"

    cmd = [
        "workflow-glue", "report_sv",
        "--report-fname", str(out_html),
        "--data", str(data_dir),
        "--versions", str(versions_path),
        "--params", str(params_path),
        "--wf-version", "0.1.0",
    ]
    subprocess.run(cmd, check=True)

    assert out_html.exists()
    content = out_html.read_text(encoding="utf-8")
    assert "Amplicon Structural Variant Report" in content
    assert "barcode01" in content
    assert "barcode02" in content
    assert "Sniffles2.DEL.1" in content
    assert "Sniffles2.DEL.2" in content
    assert "data:image/png;base64," in content
    assert "No structural variants met the filtering criteria." in content
