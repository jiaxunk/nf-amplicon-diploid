# nf-amplicon-diploid

An Oxford Nanopore EPI2ME-compatible Nextflow workflow for targeted diploid amplicon sequencing. `nf-amplicon-diploid` provides comprehensive, end-to-end analysis for long-read amplicon resequencing, including reference-based diploid variant discovery, structural variant characterization, and reference-free *de novo* diploid amplicon assembly.

---

## Key Features

- **Three Analysis Modes (`--analysis_mode`)**:
  - `reference`: Reference-guided alignment, depth profiling, Clair3 SNP/indel calling, WhatsHap diploid phasing and haplotagging, Sniffles2 structural variant calling, Cuban breakpoint read support plots, and optional Medaka haplotyped consensus polishing.
  - `denovo`: Reference-free assembly via Miniasm or SPOA, iterative Racon polishing, Medaka neural network consensus polishing, Clair3 diploid variant calling on draft contigs, WhatsHap diploid phasing, hap1/hap2 splitting, read realignment, and PAF dot plots.
  - `both`: Executes both reference-based and *de novo* workflows in a single run, cross-maps *de novo* consensus contigs against the reference sequence (`<alias>.denovo_to_ref.bam`), and compiles all four aggregated HTML reports.
- **Rich HTML Reporting**:
  - `nf-amplicon-alignment-report.html`: Quality scores, length distributions, mapping metrics, and depth per target region.
  - `nf-amplicon-snp-report.html`: Cohort SNP/indel counts, Ti/Tv ratios, VAF distributions, and phasing block summaries.
  - `nf-amplicon-sv-report.html`: Structural variant classifications (DEL, INS, DUP, INV, BND), size profiles, and Cuban breakpoint evidence plots.
  - `nf-amplicon-denovo-report.html`: Assembly contig statistics, read-to-contig alignments, and interactive/static PAF dot plots.
- **Interactive IGV Session Generation (`--igv`)**: Generates top-level `igv.json` configured for multi-track visual inspection across reference, de novo, and both modes.
- **EPI2ME Desktop Ready**: Fully compliant schema with top-down visual ordering, bundled demo data for one-click **"Run demo analysis"**, and air-gapped `.2me` package export support.

---

## Import into EPI2ME Desktop

EPI2ME Desktop imports community workflows directly from GitHub:

1. Open **Workflows** in EPI2ME Desktop and choose **Import workflow**.
2. Enter the repository URL: `https://github.com/jiaxunk/nf-amplicon-diploid`.
3. Install the workflow, then launch it from **Installed**.
4. To test your installation, click **"Run demo analysis"** to run the pre-configured amplicon dataset.

---

## Quickstart (Command Line with Docker)

EPI2ME Desktop runs workflows in its managed Docker environment. For direct command-line execution, run Nextflow with the `standard` profile to use official containerized ONT images:

### 1. Reference-Based Mode

```bash
nextflow run . -profile standard \
  --fastq /path/to/fastq \
  --reference /path/to/reference.fasta \
  --target /path/to/target.bed \
  --sample_sheet /path/to/sample_sheet.csv \
  --analysis_mode reference \
  --igv \
  --out_dir results_reference
```

### 2. De Novo Assembly Mode

```bash
nextflow run . -profile standard \
  --fastq /path/to/fastq \
  --sample_sheet /path/to/sample_sheet.csv \
  --analysis_mode denovo \
  --reference /path/to/reference.fasta \
  --igv \
  --out_dir results_denovo
```
*(Note: `--reference` is optional in `denovo` mode; when provided, consensus contigs are also compared and cross-mapped against it).*

### 3. Combined "Both" Mode

```bash
nextflow run . -profile standard \
  --fastq /path/to/fastq \
  --reference /path/to/reference.fasta \
  --target /path/to/target.bed \
  --sample_sheet /path/to/sample_sheet.csv \
  --analysis_mode both \
  --igv \
  --out_dir results_both
```

---

## Inputs

- **Reads**: Provide either `--fastq` (FASTQ directory or file) or `--bam` (BAM/uBAM directory or file).
- **Reference**: `--reference` (FASTA). Required for `reference` and `both` modes; optional for `denovo`.
- **Target Regions**: `--target` (optional BED file). When omitted in `reference` mode, all reference contigs are treated as targets.
- **Sample Sheet**: `--sample_sheet` (optional CSV). Maps barcodes to aliases and metadata:
  ```csv
  barcode,alias,type
  barcode01,sample01,test_sample
  barcode02,sample02,test_sample
  ```
  For single unbarcoded runs, use `--sample <name>` instead.

---

## Key Parameters

| Category | Parameter | Default | Description |
| :--- | :--- | :--- | :--- |
| **Workflow** | `--analysis_mode` | `reference` | Analysis pipeline: `reference`, `denovo`, or `both`. |
| | `--igv` | `false` | Generate `igv.json` visual session file. |
| **Input / Sample** | `--fastq` / `--bam` | `null` | Mutually exclusive FASTQ or BAM input path. |
| | `--reference` | `null` | Reference FASTA file. |
| | `--target` | `null` | Optional BED file defining amplicon target regions. |
| | `--sample_sheet` | `null` | CSV mapping barcodes to sample aliases. |
| | `--sample` | `null` | Sample alias for single-sample runs. |
| **Preprocessing** | `--min_read_length` | `300` | Minimum read length retained. |
| | `--max_read_length` | `null` | Maximum read length retained. |
| | `--min_read_qual` | `10` | Minimum mean read quality (Q-score). |
| | `--reads_downsampling_size` | `2000` | Target read count per sample (`0` disables downsampling). |
| | `--trimmer` | `porechop` | Adapter trimmer: `porechop` or `none`. |
| | `--min_n_reads` | `40` | Minimum surviving reads required to proceed. |
| **Alignment** | `--force_remap` | `false` | Force remapping of already-aligned BAM inputs. |
| | `--minimap2_extra_args` | `null` | Additional arguments passed to minimap2. |
| **Variant Calling** | `--snp` | `true` | Enable Clair3 small variant calling. |
| | `--sv` | `true` | Enable Sniffles2 structural variant calling. |
| | `--phased` | `true` | Enable WhatsHap diploid phasing and haplotagging. |
| | `--consensus` | `false` | Generate polished reference-mode consensus sequences. |
| | `--min_coverage` | `20` | Minimum coverage threshold for variant calling. |
| | `--min_alignment_length` | `50` | Sniffles minimum alignment length threshold. |
| | `--cuban_max_reads` | `100` | Maximum read count plotted per Cuban SV breakpoint. |
| **De Novo Assembly**| `--max_de_novo_contigs` | `2` | Maximum contigs retained per amplicon (diploid = 2). |
| | `--force_spoa_length_threshold` | `2000` | Length threshold to use SPOA instead of Miniasm. |
| | `--minimum_mean_depth` | `30` | Minimum mean depth required for de novo assembly. |
| | `--medaka_target_depth_per_strand` | `150` | Medaka consensus target depth per strand. |

---

## Output Structure

Upon completion, `--out_dir` contains:

```text
output/
├── nf-amplicon-alignment-report.html  # Quality, mapping & target coverage report
├── nf-amplicon-snp-report.html        # Clair3 SNP/indel & phasing summary (reference/both)
├── nf-amplicon-sv-report.html         # Sniffles2 SV & Cuban breakpoint report (reference/both)
├── nf-amplicon-denovo-report.html     # Assembly statistics & dot plots (denovo/both)
├── igv.json                           # IGV session file (when --igv is enabled)
├── params.json                        # Parameter snapshot for reproducibility
├── versions.txt                       # Software tool version ledger
├── reference/                         # Reference / de novo consensus FASTA files
│   ├── reference.fa (.fai)
│   └── all-consensus-seqs.fasta (.fai)
├── execution/                         # Nextflow trace, timeline, and HTML report
└── <alias>/                           # Per-sample output directory
    ├── alignments/
    │   ├── <alias>.aligned.bam (.bai)
    │   └── <alias>.haplotagged.bam (.bai)
    ├── snp/
    │   ├── <alias>.nf_snp.vcf.gz (.tbi)
    │   └── <alias>.phased.vcf.gz (.tbi)
    ├── sv/
    │   ├── <alias>.nf_sv.vcf.gz (.tbi)
    │   └── cuban/                     # Cuban breakpoint PNG evidence plots
    └── denovo/
        ├── <alias>.consensus.fasta
        ├── <alias>.aligned.bam (.bai) # Reads aligned to de novo consensus
        ├── <alias>.denovo_to_ref.bam  # De novo consensus aligned to reference (both mode)
        ├── dotplots/                  # PAF dot plots (interactive HTML & static PNG)
        └── phases/                    # Phased hap1 / hap2 consensus contigs
```

---

## Packaging as Standalone `.2me` Bundle

For offline or air-gapped deployment in EPI2ME Desktop:

```bash
python3 bin/build_2me_package.py \
  --workflow-dir . \
  --output ./nf-amplicon-diploid-v1.0.0.2me
```

---

## License & Attribution

`nf-amplicon-diploid` is an Oxford Nanopore EPI2ME community workflow developed for diploid amplicon sequencing analysis. Distributed under the terms of the project license.
