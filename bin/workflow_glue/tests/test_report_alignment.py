"""Tests for cohort-wide alignment report (`workflow-glue report_alignment`)."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import subprocess

import pandas as pd


def _create_sample_report_dir(
    base_dir: Path, alias: str, with_fastcat: bool = True
) -> Path:
    sdir = base_dir / alias
    sdir.mkdir(parents=True, exist_ok=True)

    if with_fastcat:
        pd.DataFrame(
            [
                {
                    "filename": f"{alias}.fastq.gz",
                    "sample_name": alias,
                    "n_seqs": 1200,
                    "n_bases": 1800000,
                    "min_length": 400,
                    "max_length": 2500,
                    "mean_quality": 14.5,
                }
            ]
        ).to_csv(sdir / "per-file-stats.tsv", sep="\t", index=False)

        pd.DataFrame(
            [
                {"start": 500, "end": 1000, "count": 400},
                {"start": 1500, "end": 2000, "count": 600},
            ]
        ).to_csv(sdir / "length.hist", sep="\t", index=False, header=False)

        pd.DataFrame(
            [
                {"start": 12.0, "end": 14.0, "count": 500},
                {"start": 14.0, "end": 16.0, "count": 500},
            ]
        ).to_csv(sdir / "quality.hist", sep="\t", index=False, header=False)

        pd.DataFrame(
            [
                {
                    "filename": "trimmed.fastq.gz",
                    "sample_name": alias,
                    "n_seqs": 800,
                    "n_bases": 1150000,
                    "min_length": 380,
                    "max_length": 2400,
                    "mean_quality": 14.8,
                }
            ]
        ).to_csv(sdir / "post-trim-per-file-stats.tsv", sep="\t", index=False)

    # bamstats.tsv
    pd.DataFrame(
        [
            {
                "name": f"read_{i}",
                "ref": "chr16",
                "coverage": 95.0,
                "ref_coverage": 92.0,
                "read_length": 1500,
                "acc": 98.5,
            }
            for i in range(50)
        ]
    ).to_csv(sdir / "bamstats.tsv", sep="\t", index=False)

    # bamstats-flagstat.tsv
    pd.DataFrame(
        [
            {
                "ref": "chr16",
                "total": 50,
                "primary": 50,
                "secondary": 0,
                "supplementary": 0,
                "unmapped": 0,
                "qcfail": 0,
                "duplicate": 0,
            },
            {
                "ref": "*",
                "total": 2,
                "primary": 0,
                "secondary": 0,
                "supplementary": 0,
                "unmapped": 2,
                "qcfail": 0,
                "duplicate": 0,
            },
        ]
    ).to_csv(sdir / "bamstats-flagstat.tsv", sep="\t", index=False)

    # per-window-depth.tsv.gz
    depth_df = pd.DataFrame(
        [
            {
                "ref": "chr16",
                "start": 165000 + i * 200,
                "end": 165000 + (i + 1) * 200,
                "target": "HBA_region",
                "depth": 120.0 + (i % 5) * 10.0,
            }
            for i in range(20)
        ]
    )
    with gzip.open(sdir / "per-window-depth.tsv.gz", "wt") as gz:
        depth_df.to_csv(gz, sep="\t", index=False)

    return sdir


def test_report_alignment_generates_cohort_html_with_sample_tabs(tmp_path):
    """report_alignment produces nf-amplicon-alignment-report.html for FASTQ and BAM shortcut samples."""
    data_dir = tmp_path / "data"
    _create_sample_report_dir(data_dir, "barcode08", with_fastcat=True)
    _create_sample_report_dir(data_dir, "barcode09", with_fastcat=False)

    versions_txt = tmp_path / "versions.txt"
    versions_txt.write_text("minimap2,2.28\nsamtools,1.21\nmosdepth,0.3.8\n")

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps(
            {
                "analysis_mode": "reference",
                "reads_downsampling_size": 2000,
                "trimmer": "porechop",
            }
        )
    )

    ref_fa = tmp_path / "ref.fasta"
    ref_fa.write_text(">chr16\n" + "A" * 200000 + "\n")

    target_bed = tmp_path / "target.bed"
    target_bed.write_text("chr16\t165000\t185000\tHBA_region\n")

    report_html = tmp_path / "nf-amplicon-alignment-report.html"

    res = subprocess.run(
        [
            "workflow-glue",
            "report_alignment",
            "--report-fname",
            str(report_html),
            "--data",
            str(data_dir),
            "--reference",
            str(ref_fa),
            "--target-bed",
            str(target_bed),
            "--downsampling-size",
            "2000",
            "--versions",
            str(versions_txt),
            "--params",
            str(params_json),
            "--wf-version",
            "v1.0.0",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert res.returncode == 0
    assert report_html.exists()
    html_text = report_html.read_text()
    assert "barcode08" in html_text
    assert "barcode09" in html_text
    assert "HBA_region" in html_text
    assert "Target Region" in html_text or "Depth of coverage" in html_text
