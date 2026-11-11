"""Tests for length-mode-preserving read subsetting (`workflow-glue subset_reads`)."""
from __future__ import annotations

from pathlib import Path
import subprocess

import numpy as np
import pandas as pd


def _write_per_read_stats(path: Path, lengths: dict[str, int]) -> Path:
    df = pd.DataFrame(
        {
            "read_id": list(lengths.keys()),
            "read_length": list(lengths.values()),
        }
    )
    df.to_csv(path, sep="\t", index=False)
    return path


def test_subset_reads_preserves_bimodal_diploid_alleles(tmp_path):
    """When reads have a shorter deletion allele peak and a longer WT peak, both are preserved."""
    rng = np.random.default_rng(42)
    lengths = {}
    # 475 reads from shorter deletion allele (~800 bp)
    for i in range(475):
        lengths[f"del_read_{i}"] = int(rng.normal(800, 30))
    # 475 reads from longer wild-type allele (~2800 bp)
    for i in range(475):
        lengths[f"wt_read_{i}"] = int(rng.normal(2800, 50))
    # 50 chimeric concatemers (~5600 bp, top ~5%)
    for i in range(50):
        lengths[f"chimera_{i}"] = int(rng.normal(5600, 60))

    stats_tsv = _write_per_read_stats(tmp_path / "per-read-stats.tsv", lengths)

    res = subprocess.run(
        [
            "workflow-glue",
            "subset_reads",
            str(stats_tsv),
            "0.05",
            "200",
            "--take-longest-remaining",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    selected_ids = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    assert len(selected_ids) == 200

    # Chimeras should be dropped by the 0.05 longest fraction filter
    chimera_count = sum(1 for r in selected_ids if r.startswith("chimera_"))
    assert chimera_count == 0

    # Both the shorter deletion allele and longer WT allele peaks must be well-represented
    del_count = sum(1 for r in selected_ids if r.startswith("del_read_"))
    wt_count = sum(1 for r in selected_ids if r.startswith("wt_read_"))
    assert del_count >= 70, f"Expected shorter deletion allele to be preserved, got {del_count}"
    assert wt_count >= 70, f"Expected longer WT allele to be preserved, got {wt_count}"


def test_subset_reads_unimodal_takes_longest_after_drop(tmp_path):
    """For a single unimodal amplicon, subset_reads drops the top fraction and keeps the longest remaining."""
    lengths = {f"read_{i}": 1000 + i for i in range(100)}
    stats_tsv = _write_per_read_stats(tmp_path / "per-read-stats.tsv", lengths)

    res = subprocess.run(
        [
            "workflow-glue",
            "subset_reads",
            str(stats_tsv),
            "0.10",
            "20",
            "--take-longest-remaining",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    selected_ids = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    assert len(selected_ids) == 20
    # Top 10 reads (read_90..read_99) dropped; next 20 longest (read_89..read_70) kept
    expected_ids = [f"read_{i}" for i in range(89, 69, -1)]
    assert selected_ids == expected_ids


def test_subset_reads_zero_n_reads_keeps_all_after_drop(tmp_path):
    """When n_reads is 0, all reads remaining after drop_longest_frac are returned."""
    lengths = {f"read_{i}": 500 + i * 10 for i in range(20)}
    stats_tsv = _write_per_read_stats(tmp_path / "per-read-stats.tsv", lengths)

    res = subprocess.run(
        [
            "workflow-glue",
            "subset_reads",
            str(stats_tsv),
            "0.10",
            "0",
            "--take-longest-remaining",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    selected_ids = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    assert len(selected_ids) == 18
    assert "read_19" not in selected_ids
    assert "read_18" not in selected_ids
