# nf-amplicon-diploid

An Oxford Nanopore Nextflow workflow for diploid amplicon sequencing. The current release implements **reference-based read preprocessing, alignment, target-depth summaries, and an alignment report**. The planned de novo, combined, small-variant, and structural-variant analyses are not implemented yet.

## Import into EPI2ME Desktop

EPI2ME Desktop requires the repository README to be named exactly `README.md` for community workflow import. The workflow manifest must also match the GitHub owner/repository path and use a semantic version. See the [EPI2ME Desktop community workflow requirements](https://epi2me.nanoporetech.com/epi2me-docs/help/faq/).

1. Open **Workflows** in EPI2ME Desktop and choose **Import workflow**.
2. Enter `https://github.com/jiaxunk/nf-amplicon-diploid`.
3. Install the workflow, then launch it from **Installed**.

Use the repository URL directly; do not append `/tree/main` or a file path. Import instructions are also available in the [EPI2ME quick start guide](https://epi2me.nanoporetech.com/epi2me-docs/quickstart/).

## Current capabilities

- Accepts either FASTQ or BAM/uBAM input, as a file or directory tree.
- Supports single samples and multiplexed barcode directories, with an optional CSV sample sheet.
- Filters reads by minimum/maximum length and minimum mean quality.
- Downsamples after filtering to 2,000 reads per sample by default, preserving distinct read-length modes when detected. Set `--reads_downsampling_size 0` to disable downsampling.
- Trims reads with `porechop`, `dorado`, or `none`.
- Reuses aligned BAMs whose reference sequence names and lengths match the supplied FASTA. Set `--force_remap` to convert and remap them.
- Maps reads to the full reference, validates an optional target BED, and reports target-window depth and alignment metrics.
- Publishes sorted/indexed alignments per sample and a cohort alignment report.

Only `--analysis_mode reference` is currently available. Selecting `denovo` or `both` exits with an unsupported-mode message. Although the configuration exposes future variant and assembly parameters, those analyses are not implemented in this release. A packaged EPI2ME demo dataset is not yet configured.

## Inputs

Provide exactly one of `--fastq` or `--bam`, and provide `--reference` as a FASTA file. `--target` is an optional BED file; without it, every reference contig is treated as a target. For multiplexed inputs, provide a sample sheet with barcode and alias columns, for example:

```csv
barcode,alias,type
barcode01,sample01,test_sample
barcode02,sample02,test_sample
```

For a single input sample without a sample sheet, set `--sample` to the sample name.

## Run with Docker

EPI2ME Desktop runs the workflow in its managed environment. For command-line use, run Nextflow in the project's configured WSL/Docker environment with the `standard` profile; workflow tools run in the configured ONT containers.

```bash
nextflow run . -profile standard \
  --fastq /path/to/fastq \
  --reference /path/to/reference.fasta \
  --target /path/to/targets.bed \
  --sample_sheet /path/to/sample_sheet.csv \
  --out_dir results
```

For BAM input, replace `--fastq` with `--bam`. Omit `--target` when the entire reference is the target. The minimum Nextflow version in the workflow manifest is `23.04.2`.

## Main parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--fastq` / `--bam` | — | Mutually exclusive FASTQ or BAM/uBAM input. |
| `--reference` | required | Reference FASTA used for alignment. |
| `--target` | `null` | Optional target BED; otherwise all reference contigs are used. |
| `--sample` | `null` | Name for a single input sample without a sample sheet. |
| `--sample_sheet` | `null` | CSV mapping barcodes to aliases and sample metadata. |
| `--min_read_length` | `300` | Minimum read length retained. |
| `--max_read_length` | `null` | Optional maximum read length. |
| `--min_read_qual` | `10` | Minimum mean read quality. |
| `--reads_downsampling_size` | `2000` | Reads retained per sample after filtering; `0` disables downsampling. |
| `--trimmer` | `porechop` | One of `porechop`, `dorado`, or `none`. |
| `--min_n_reads` | `40` | Minimum reads retained after preprocessing. |
| `--force_remap` | `false` | Convert and remap BAMs even when they already match the reference. |
| `--minimap2_extra_args` | `null` | Additional arguments passed to minimap2. |
| `--out_dir` | `output` | Workflow output directory. |

## Outputs

- `output/<sample>/alignments/aligned.bam` and its index.
- `output/nf-amplicon-alignment-report.html`, with per-sample read, alignment, and target-depth summaries.
- `output/params.json` and `output/versions.txt`.
- Nextflow execution timeline, report, and trace files under `output/execution/`.

## Test data

Small input fixtures are under [`test_data/`](test_data/). The `diploid_amplicon/` fixture includes FASTQ, BAM/uBAM, a reference, a target BED, and a sample sheet. The reads are genomic data; generic filenames do not anonymize their sequence or embedded alignment metadata.
