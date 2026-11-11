"""Subset reads based on length (from fastcat per-read stats)."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from .util import get_named_logger, wf_parser  # noqa: ABS101


def _detect_valley_splits(lengths: np.ndarray) -> list[float]:
    """Detect length thresholds (valleys) separating distinct allele length modes.

    Returns a sorted list of valley split points when clearly separated modes exist,
    or an empty list if the read-length distribution is unimodal.
    """
    n_total = len(lengths)
    if n_total < 20:
        return []
    min_len = float(lengths.min())
    max_len = float(lengths.max())
    span = max_len - min_len
    min_peak_sep_bp = max(200.0, 0.15 * max_len)
    if span < min_peak_sep_bp:
        return []

    n_bins = min(40, max(12, int(np.sqrt(n_total))))
    counts, edges = np.histogram(lengths, bins=n_bins, range=(min_len, max_len))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = (max_len - min_len) / n_bins

    # Pad with zeros at boundaries so edge modes (bin 0 or bin -1) are detected by find_peaks
    smoothed_core = np.convolve(counts.astype(float), [0.25, 0.5, 0.25], mode="same")
    padded = np.concatenate([[0.0], smoothed_core, [0.0]])

    min_distance_bins = max(2, int(round(min_peak_sep_bp / bin_width)))
    min_height = max(3.0, 0.05 * n_total)
    raw_peaks, _ = find_peaks(
        padded,
        height=min_height,
        prominence=min_height * 0.5,
        distance=min_distance_bins,
    )
    peaks = [int(p) - 1 for p in raw_peaks if 1 <= int(p) <= len(smoothed_core)]
    smoothed = smoothed_core
    if len(peaks) <= 1:
        return []

    # Verify that adjacent peaks are separated by a genuine valley (< 50% of smaller peak)
    valid_peaks = [int(peaks[0])]
    splits = []
    for pk in peaks[1:]:
        pk = int(pk)
        prev_pk = valid_peaks[-1]
        valley_slice = smoothed[prev_pk : pk + 1]
        valley_rel_idx = int(np.argmin(valley_slice))
        valley_idx = prev_pk + valley_rel_idx
        valley_val = smoothed[valley_idx]
        smaller_peak = min(smoothed[prev_pk], smoothed[pk])
        if (
            centers[pk] - centers[prev_pk] >= min_peak_sep_bp
            and valley_val <= 0.5 * smaller_peak
        ):
            splits.append(float(centers[valley_idx]))
            valid_peaks.append(pk)
        elif smoothed[pk] > smoothed[prev_pk]:
            valid_peaks[-1] = pk

    return splits


def select_length_mode_preserving_reads(
    sorted_lengths: pd.Series, n_reads: int, logger=None
) -> pd.Index:
    """Select `n_reads` while preserving distinct diploid read-length modes."""
    if n_reads == 0 or n_reads >= len(sorted_lengths):
        if logger:
            logger.info("Take all remaining reads (sorted by length).")
        return sorted_lengths.index

    lengths = sorted_lengths.to_numpy(dtype=float)
    splits = _detect_valley_splits(lengths)
    if not splits:
        if logger:
            logger.info(f"Unimodal distribution: take {n_reads} longest remaining reads.")
        return sorted_lengths.iloc[:n_reads].index

    # Partition reads into mode clusters using the valley split thresholds
    boundaries = [-np.inf] + sorted(splits) + [np.inf]
    clusters = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        cluster = sorted_lengths[(sorted_lengths > lo) & (sorted_lengths <= hi)]
        if not cluster.empty:
            clusters.append(cluster)

    if len(clusters) <= 1:
        return sorted_lengths.iloc[:n_reads].index

    if logger:
        logger.info(
            f"Detected {len(clusters)} read-length modes; sampling proportionally "
            f"across modes up to {n_reads} reads."
        )

    # Proportional allocation across clusters
    total = sum(len(c) for c in clusters)
    raw_alloc = [n_reads * (len(c) / total) for c in clusters]
    alloc = [min(len(c), int(np.floor(x))) for c, x in zip(clusters, raw_alloc)]

    # Distribute any remainder by largest fractional part (and available capacity)
    remainder = n_reads - sum(alloc)
    frac_order = sorted(
        range(len(clusters)),
        key=lambda i: (raw_alloc[i] - alloc[i], len(clusters[i])),
        reverse=True,
    )
    for idx in frac_order:
        if remainder <= 0:
            break
        capacity = len(clusters[idx]) - alloc[idx]
        if capacity > 0:
            take_extra = min(remainder, capacity)
            alloc[idx] += take_extra
            remainder -= take_extra

    selected_series = pd.concat(
        [cluster.iloc[:k] for cluster, k in zip(clusters, alloc) if k > 0]
    ).sort_values(ascending=False)
    return selected_series.index


def main(args):
    """Run the entry point."""
    logger = get_named_logger("subsetReads")

    logger.info("Read per-read stats and sort lengths.")
    sorted_lengths = pd.read_csv(
        args.per_read_stats,
        sep="\t",
        index_col="read_id",
        usecols=["read_id", "read_length"],
    )["read_length"].sort_values(ascending=False)

    drop_longest_n = 0
    if args.drop_longest_frac:
        drop_longest_n = int(round(len(sorted_lengths) * args.drop_longest_frac))
        logger.info(f"Drop {drop_longest_n} longest reads.")
        sorted_lengths = sorted_lengths.iloc[drop_longest_n:]

    if args.take_longest_remaining:
        target_ids = select_length_mode_preserving_reads(
            sorted_lengths, args.n_reads, logger=logger
        )
    else:
        logger.info("Randomly select reads and write out IDs.")
        # select a subset and sort by length again (but only if we still got enough
        # reads)
        if args.n_reads > len(sorted_lengths) or args.n_reads == 0:
            target_ids = sorted_lengths.index
        else:
            target_ids = (
                sorted_lengths.sample(args.n_reads).sort_values(ascending=False).index
            )

    sys.stdout.write("\n".join(target_ids) + "\n")

    logger.info(
        f"Finished extracting {len(target_ids)} reads from '{args.per_read_stats}'."
    )


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("subset_reads")
    parser.add_argument(
        "per_read_stats", type=Path, help="Path to per-read stats TSV file"
    )
    parser.add_argument(
        "drop_longest_frac", type=float, help="Fraction of longest reads to drop"
    )
    parser.add_argument(
        "n_reads", type=int, help="Number of reads to select. If 0, take all."
    )
    parser.add_argument(
        "--take-longest-remaining",
        action="store_true",
        help=(
            "Whether to select the longest remaining reads (preserving distinct "
            "read-length modes) or randomly sample after dropping the longest reads"
        ),
    )
    return parser
