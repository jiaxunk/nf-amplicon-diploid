# Diploid Amplicon Analysis

An EPI2ME-compatible Nextflow workflow for analyzing Oxford Nanopore diploid amplicon sequencing data via reference-based variant calling and/or de novo haplotype-aware assembly.

## Language

**Analysis Mode**:
The top-level operational branch selected for a workflow run: Reference-Based Mode, De Novo Mode, or Both.
_Avoid_: Pipeline type, workflow mode

**Reference-Based Mode**:
An analysis path that aligns reads to a known reference sequence to call small variants, phase haplotypes, call structural variants, and generate variant-informed consensus sequences.
_Avoid_: Variant calling mode, mapping mode

**De Novo Mode**:
An analysis path that assembles amplicon consensus sequences directly from sequencing reads without prior reference alignment; if a Reference Sequence (and optional Target Region) is supplied, the assembled sequences are aligned against the reference (or target slice) for comparison and dot-plot visualization.
_Avoid_: Assembly-only mode

**Reference Sequence**:
The known genome, chromosome, or amplicon FASTA sequence used for read mapping in Reference-Based Mode or post-assembly comparison in De Novo Mode.
_Avoid_: Template, backbone

**Target Region**:
A genomic interval (`chrom`, `start`, `end`, `name`) defined in a BED file that scopes variant calling, coverage, and consensus generation in Reference-Based Mode, or slices the Reference Sequence for post-assembly dot-plot comparison in De Novo Mode.
_Avoid_: ROI, amplicon bed, interval list

**Haplotagged Alignment**:
A coordinate-sorted BAM file in which reads are annotated with haplotype (`HP`) and phase-set (`PS`) tags derived from phased heterozygous variants.
_Avoid_: Phased BAM, tagged reads

**Haplotype Consensus**:
A FASTA sequence representing one resolved allele (`hap1` or `hap2`) of a diploid Target Region (or a single IUPAC-ambiguity sequence when unphased), generated either by applying phased variant calls to the Reference Sequence or by de novo haplotype assembly.
_Avoid_: Draft contig, polished sequence

**Breakpoint Evidence Plot**:
A multi-panel figure rendered by `cuban` for a structural variant call, stacking target-normalized read depth, repeat annotations, and haplotype-grouped (`HP:1` / `HP:2`) split/clipped read alignments across the SV breakpoints.
_Avoid_: SV snapshot, pileup image

**Consensus Dot Plot**:
A 2D coordinate comparison between an assembled Haplotype Consensus sequence and the Reference Sequence (or Target Region slice) that reveals large insertions, deletions, and inversions as diagonal offsets or inversions.
_Avoid_: Alignment matrix, synteny plot
