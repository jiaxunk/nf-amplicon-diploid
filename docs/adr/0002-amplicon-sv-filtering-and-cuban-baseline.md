# Amplicon SV filtering and target-scoped cuban baseline coverage

In amplicon sequencing, coverage drops sharply to zero at primer boundaries. Because Sniffles2's genotyping model assumes uniform genome-wide coverage, genuine large deletions anchored near amplicon ends are frequently flagged with `FILTER=GT` despite having thousands of supporting reads and high VAF. Similarly, `cuban` defaults to computing baseline coverage across the entire chromosome, which yields a near-zero baseline when a 20 kb amplicon is mapped to a whole chromosome.

We decided to retain both `PASS` and `GT` Sniffles2 calls that overlap the Target Region and pass `--min_sv_support` and `--min_sv_vaf` thresholds, and to explicitly pass the Target Region's mean depth (`--baseline-cov <sample>:<target_mean_depth>`) to `cuban`.
