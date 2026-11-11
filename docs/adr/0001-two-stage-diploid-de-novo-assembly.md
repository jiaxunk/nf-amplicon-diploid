# Two-stage diploid de novo amplicon assembly

Whole-genome assemblers (such as Flye or Hifiasm) struggle with short, ultra-high-depth PCR amplicons, while stock `wf-amplicon` (`miniasm`/`spoa` -> `racon` -> `medaka`) is strictly haploid: its read subsetting keeps only the longest reads (discarding shorter length-heterozygous deletion alleles) and `trim_and_qc` retains only the single highest-depth contig (collapsing SNV-heterozygous alleles).

We decided to extend the `miniasm`/`spoa` -> `racon` -> `medaka` backbone into a two-stage diploid de novo pipeline: first, sample reads across read-length peaks and allow `trim_and_qc` to retain up to two distinct allele contigs (`max_contigs = 2`) to preserve structural/length heterozygotes; second, align reads back to each draft contig and phase heterozygous variants (Clair3 + WhatsHap) to polish and emit haplotype-resolved consensus sequences (`hap1` and `hap2`).
