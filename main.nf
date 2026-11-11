#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import Pinguscript
import WorkflowMain

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

workflow {
    WorkflowMain.initialise(workflow, params, log)
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.analysis_mode != 'reference') {
        error "Ticket 01 implements reference alignment; '${params.analysis_mode}' mode is not available yet."
    }
    if (params.analysis_mode in ['reference', 'both'] && !params.reference) {
        error "Analysis mode '${params.analysis_mode}' requires --reference."
    }
    if (!!params.fastq == !!params.bam) {
        error 'Provide exactly one of --fastq or --bam.'
    }
    if (!(params.trimmer in ['porechop', 'dorado', 'none'])) {
        error "Unsupported --trimmer '${params.trimmer}'; choose porechop, dorado, or none."
    }
    if (params.reads_downsampling_size != null && params.reads_downsampling_size < 0) {
        error '--reads_downsampling_size must be zero (disabled) or a positive integer.'
    }
    if (params.min_n_reads < 0) {
        error '--min_n_reads must be zero or a positive integer.'
    }

    ch_params = getParams()
    ch_reference = build_reference(Channel.fromPath(params.reference, checkIfExists: true))
    ch_target = params.target
        ? Channel.fromPath(params.target, checkIfExists: true)
        : make_full_target(ch_reference.map { ref, mmi, fai -> fai })

    if (params.fastq) {
        ch_ingress = fastq_ingress([
            input: params.fastq,
            sample: params.sample,
            sample_sheet: params.sample_sheet,
            analyse_unclassified: params.analyse_unclassified,
            analyse_fail: params.analyse_fail,
            stats: true,
            per_read_stats: true,
            allow_multiple_basecall_models: true,
        ])
        ch_fastq = ch_ingress
            .filter { meta, reads, stats -> reads != null }
            .map { meta, reads, stats -> [meta, reads instanceof List ? reads[0] : reads, stats] }
            | filter_downsample
    } else {
        ch_xam = xam_ingress([
            input: params.bam,
            sample: params.sample,
            sample_sheet: params.sample_sheet,
            analyse_unclassified: params.analyse_unclassified,
            analyse_fail: params.analyse_fail,
            stats: true,
            per_read_stats: true,
            allow_multiple_basecall_models: true,
        ])
        ch_routes = ch_xam.filter { meta, bam, bai, stats -> bam != null }
        ch_classified = classify_aligned_bam(
            ch_routes.combine(ch_reference)
                .map { meta, bam, bai, stats, ref, mmi, fai -> [meta, bam, bai, stats, ref, fai] }
        )
        ch_bam_routes = ch_classified.branch { meta, bam, bai, stats, match ->
            already_aligned: match == '1' && !params.force_remap
            remap: true
        }

        ch_bam_fastq = xam_ingress([
            input: params.bam,
            sample: params.sample,
            sample_sheet: params.sample_sheet,
            analyse_unclassified: params.analyse_unclassified,
            analyse_fail: params.analyse_fail,
            keep_unaligned: true,
            return_fastq: true,
            stats: true,
            per_read_stats: true,
            allow_multiple_basecall_models: true,
        ])
        ch_remap_aliases = ch_bam_routes.remap.map { meta, bam, bai, stats, match -> meta.alias }.collect()
        ch_fastq = ch_bam_fastq
            .combine(ch_remap_aliases)
            .filter { meta, reads, stats, remap_aliases ->
                meta.is_unaligned || meta.alias in remap_aliases
            }
            .map { meta, reads, stats, remap_aliases -> [meta, reads, stats] }
            | filter_downsample
        ch_already_aligned = ch_bam_routes.already_aligned
            .map { meta, bam, bai, stats, match -> [meta, bam, bai] }
    }

    ch_pretrimmed = ch_fastq
        .filter { meta, reads, report_data, n_filtered ->
            if (n_filtered.toInteger() == 0 || n_filtered.toInteger() < params.min_n_reads) {
                log.warn "Discarding sample '${meta.alias}': ${n_filtered} reads remain after read filtering; minimum is ${params.min_n_reads}."
                return false
            }
            true
        }
        .map { meta, reads, report_data, n_filtered -> [meta, reads, report_data] }

    if (params.trimmer == 'porechop') {
        ch_trimmed = porechop_trim(ch_pretrimmed)
    } else if (params.trimmer == 'dorado') {
        ch_trimmed = post_trim_stats(dorado_trim(ch_pretrimmed))
    } else {
        ch_trimmed = no_trim(ch_pretrimmed)
    }

    ch_usable_reads = ch_trimmed
        .filter { meta, reads, report_data, n_reads ->
            if (n_reads.toInteger() < params.min_n_reads) {
                log.warn "Discarding sample '${meta.alias}': ${n_reads} reads remain after trimming; minimum is ${params.min_n_reads}."
                return false
            }
            true
        }
        .map { meta, reads, report_data, n_reads -> [meta, reads, report_data] }

    ch_map_inputs = ch_usable_reads.combine(ch_reference).combine(ch_target)
        .map { meta, reads, report_data, ref, mmi, fai, target ->
            [meta, reads, report_data, ref, mmi, fai, target]
        }
    ch_mapped = map_reads(ch_map_inputs)

    if (params.bam) {
        ch_existing_inputs = ch_already_aligned.combine(ch_reference).combine(ch_target)
            .map { meta, bam, bai, ref, mmi, fai, target -> [meta, bam, bai, ref, fai, target] }
        ch_existing = keep_aligned_bam(ch_existing_inputs)
        ch_results = ch_mapped.mix(ch_existing)
    } else {
        ch_results = ch_mapped
    }

    ch_report_dirs = ch_results.map { meta, bam, bai, report_data -> report_data }.collect()
    ch_target_for_report = ch_target
    ch_versions = software_versions()
    make_alignment_report(
        ch_report_dirs,
        Channel.value(file(params.reference)),
        ch_target_for_report,
        ch_params,
        ch_versions
    )

    if (params.snp) {
        ch_snp = ch_results
            .combine(ch_reference)
            .combine(ch_target)
            .map { meta, bam, bai, report_data, ref, mmi, fai, target ->
                [meta, bam, bai, ref, fai, target, report_data]
            }
        call_clair3(ch_snp)
        if (params.phased) {
            phase_snp(call_clair3.out.calls)
            make_haplotype_consensus(phase_snp.out.consensus_inputs)
        } else {
            make_haplotype_consensus(call_clair3.out.calls)
        }
        make_snp_report(
            params.phased ? phase_snp.out.report_dirs.collect() : call_clair3.out.report_dirs.collect(),
            ch_params,
            ch_versions
        )
    }
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}

process build_reference {
    label 'wfamplicon'
    cpus 2
    memory '2 GB'
    input:
        path reference
    output:
        tuple path('reference.fa'), path('reference.mmi'), path('reference.fa.fai')
    script:
    """
    cp '$reference' reference.fa
    minimap2 -d reference.mmi reference.fa
    samtools faidx reference.fa
    """
}

process make_full_target {
    label 'wfamplicon'
    cpus 1
    memory '1 GB'
    input:
        path reference_fai
    output:
        path 'targets.bed'
    script:
    """
    awk 'BEGIN { OFS="\\t" } { print \$1, 0, \$2, \$1 }' '$reference_fai' > targets.bed
    """
}

process filter_downsample {
    label 'wfamplicon'
    cpus 2
    memory '4 GB'
    input:
        tuple val(meta), path(reads), path(ingress_stats)
    output:
        tuple val(meta), path('downsampled.fastq.gz'), path('report_data'), env(N_FILTERED)
    script:
    def max_len_arg = params.max_read_length ? "-M ${params.max_read_length}" : ''
    def n_reads = params.reads_downsampling_size ?: 0
    def take_longest_arg = params.take_longest_remaining_reads ? '--take-longest-remaining' : ''
    """
    mkdir -p report_data
    cp -a '$ingress_stats/.' report_data/ 2>/dev/null || true
    seqkit seq -m ${params.min_read_length} $max_len_arg -Q ${params.min_read_qual} '$reads' \
        | gzip -c > filtered.fastq.gz

    N_FILTERED=\$(zcat filtered.fastq.gz | awk 'NR % 4 == 1 { n++ } END { print n+0 }')
    if [ "\$N_FILTERED" -gt 0 ]; then
        fastcat --histograms filtered_histograms \
            -r >(bgzip -c > filtered-per-read-stats.tsv.gz) \
            filtered.fastq.gz > /dev/null
        mv filtered_histograms/* report_data/ 2>/dev/null || true
    else
        printf 'read_id\\tread_length\\n' > filtered-per-read-stats.tsv
        : > filtered-per-read-stats.tsv.gz
    fi

    if [ "$n_reads" -gt 0 ] && [ "\$N_FILTERED" -gt 0 ]; then
        workflow-glue subset_reads filtered-per-read-stats.tsv.gz \
            ${params.drop_frac_longest_reads} "$n_reads" $take_longest_arg > selected-read-ids.txt
        seqkit grep -f selected-read-ids.txt filtered.fastq.gz \
            | gzip -c > downsampled.fastq.gz
    else
        cp filtered.fastq.gz downsampled.fastq.gz
    fi
    export N_FILTERED
    """
}

process classify_aligned_bam {
    label 'wf_common'
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path(bam), path(bai), path(stats), path(reference), path(reference_fai)
    output:
        tuple val(meta), path(bam), path(bai), path(stats), env(REFERENCE_MATCH)
    script:
    """
    status=0
    workflow-glue check_sq_ref --xam '$bam' --ref '$reference' >/dev/null 2>check-ref.log || status=\$?
    if [ "\$status" -eq 0 ]; then
        export REFERENCE_MATCH=1
    elif [ "\$status" -eq 65 ]; then
        export REFERENCE_MATCH=0
    else
        cat check-ref.log >&2
        exit "\$status"
    fi
    """
}

process porechop_trim {
    label 'wfamplicon'
    cpus 4
    memory '4 GB'
    input:
        tuple val(meta), path(reads), path(report_data, stageAs: 'preprocess_data')
    output:
        tuple val(meta), path('trimmed.fastq.gz'), path('report_data'), env(N_READS)
    script:
    """
    porechop --threads $task.cpus ${params.trimmer_extra_args ?: '--discard_middle'} \
        -i '$reads' -o trimmed.fastq
    bgzip -@ $task.cpus trimmed.fastq
    mkdir -p report_data
    cp -a preprocess_data/. report_data/
    fastcat -f report_data/post-trim-per-file-stats.tsv \
        -i report_data/post-trim-runids.tsv -l report_data/post-trim-basecallers.tsv \
        trimmed.fastq.gz > /dev/null
    N_READS=\$(zcat trimmed.fastq.gz | awk 'NR % 4 == 1 { n++ } END { print n+0 }')
    export N_READS
    """
}

process dorado_trim {
    label 'wf_basecalling'
    cpus 4
    memory '4 GB'
    input:
        tuple val(meta), path(reads), path(report_data, stageAs: 'report_data')
    output:
        tuple val(meta), path('trimmed.fastq.gz'), path('report_data')
    script:
    """
    dorado --version > dorado-version.txt
    if ! dorado trim ${params.trimmer_extra_args ?: ''} --emit-fastq '$reads' \
        > trimmed.fastq 2>dorado-trim.log; then
        cat dorado-trim.log >&2
        echo 'The selected Dorado container could not trim this FASTQ input. Check its Dorado version/options, or choose porechop or none.' >&2
        exit 2
    fi
    bgzip -@ $task.cpus trimmed.fastq
    mv dorado-version.txt report_data/
    """
}

process post_trim_stats {
    label 'wfamplicon'
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path(reads), path(report_data, stageAs: 'preprocess_data')
    output:
        tuple val(meta), path('trimmed.fastq.gz'), path('report_data'), env(N_READS)
    script:
    """
    mkdir -p report_data
    cp -a preprocess_data/. report_data/
    fastcat -f report_data/post-trim-per-file-stats.tsv \
        -i report_data/post-trim-runids.tsv -l report_data/post-trim-basecallers.tsv \
        '$reads' > /dev/null
    N_READS=\$(zcat '$reads' | awk 'NR % 4 == 1 { n++ } END { print n+0 }')
    export N_READS
    """
}

process no_trim {
    label 'wfamplicon'
    cpus 1
    memory '1 GB'
    input:
        tuple val(meta), path(reads), path(report_data, stageAs: 'preprocess_data')
    output:
        tuple val(meta), path('trimmed.fastq.gz'), path('report_data'), env(N_READS)
    script:
    """
    cp '$reads' trimmed.fastq.gz
    mkdir -p report_data
    cp -a preprocess_data/. report_data/
    fastcat -f report_data/post-trim-per-file-stats.tsv \
        -i report_data/post-trim-runids.tsv -l report_data/post-trim-basecallers.tsv \
        trimmed.fastq.gz > /dev/null
    N_READS=\$(zcat trimmed.fastq.gz | awk 'NR % 4 == 1 { n++ } END { print n+0 }')
    export N_READS
    """
}

process map_reads {
    label 'wfamplicon'
    cpus 4
    memory '8 GB'
    publishDir "${params.out_dir}/${meta.alias}/alignments", mode: 'copy', pattern: 'aligned.bam*'
    input:
        tuple val(meta), path(reads), path(report_data, stageAs: 'preprocess_data'), path(reference), path(reference_index), path(reference_fai), path(target_bed)
    output:
        tuple val(meta), path('aligned.bam'), path('aligned.bam.bai'), path("${meta.alias}")
    script:
    """
    mkdir -p report_data
    cp -a preprocess_data/. report_data/
    awk -F '\\t' 'FNR==NR { lengths[\$1]=\$2; next } /^#/ || /^track/ || /^browser/ || NF==0 { next } ! (\$1 in lengths) || \$2 < 0 || \$3 <= \$2 || \$3 > lengths[\$1] { print "Invalid target BED record at line " FNR > "/dev/stderr"; exit 1 }' \
        '$reference_fai' '$target_bed'
    minimap2 -ax map-ont -t $task.cpus ${params.minimap2_extra_args ?: ''} \
        '$reference_index' '$reads' \
        | samtools sort -@ $task.cpus -o aligned.bam -
    samtools index -@ $task.cpus aligned.bam
    bamstats aligned.bam -s '${meta.alias}' -u -f report_data/bamstats-flagstat.tsv \
        -t $task.cpus > report_data/bamstats.tsv
    awk -F '\\t' -v n=${params.number_depth_windows} 'BEGIN { OFS="\\t" } /^#/ || /^track/ || /^browser/ || NF < 3 { next } { len=\$3-\$2; k=(len<n ? len : n); name=(NF>=4 ? \$4 : \$1 ":" \$2 "-" \$3); for (i=0; i<k; i++) { s=\$2+int(i*len/k); e=\$2+int((i+1)*len/k); print \$1,s,e,name "_window" i } }' \
        '$target_bed' > report_data/depth-windows.bed
    mosdepth --by report_data/depth-windows.bed --threads $task.cpus report_data/target aligned.bam
    { printf 'ref\\tstart\\tend\\tdepth\\n'; zcat report_data/target.regions.bed.gz | awk 'BEGIN { OFS="\\t" } { print \$1, \$2, \$3, \$NF }'; } \
        | bgzip -c > report_data/per-window-depth.tsv.gz
    mv report_data '${meta.alias}'
    """
}

process keep_aligned_bam {
    label 'wfamplicon'
    cpus 4
    memory '8 GB'
    publishDir "${params.out_dir}/${meta.alias}/alignments", mode: 'copy', pattern: 'aligned.bam*'
    input:
        tuple val(meta), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed)
    output:
        tuple val(meta), path('aligned.bam'), path('aligned.bam.bai'), path("${meta.alias}")
    script:
    """
    mkdir -p report_data
    awk -F '\\t' 'FNR==NR { lengths[\$1]=\$2; next } /^#/ || /^track/ || /^browser/ || NF==0 { next } ! (\$1 in lengths) || \$2 < 0 || \$3 <= \$2 || \$3 > lengths[\$1] { print "Invalid target BED record at line " FNR > "/dev/stderr"; exit 1 }' \
        '$reference_fai' '$target_bed'
    cp '$bam' aligned.bam
    cp '$bai' aligned.bam.bai
    printf 'Already aligned to the supplied reference; FASTQ conversion and trimming were skipped.\\n' \
        > report_data/already-aligned.txt
    bamstats aligned.bam -s '${meta.alias}' -u -f report_data/bamstats-flagstat.tsv \
        -t $task.cpus > report_data/bamstats.tsv
    awk -F '\\t' -v n=${params.number_depth_windows} 'BEGIN { OFS="\\t" } /^#/ || /^track/ || /^browser/ || NF < 3 { next } { len=\$3-\$2; k=(len<n ? len : n); name=(NF>=4 ? \$4 : \$1 ":" \$2 "-" \$3); for (i=0; i<k; i++) { s=\$2+int(i*len/k); e=\$2+int((i+1)*len/k); print \$1,s,e,name "_window" i } }' \
        '$target_bed' > report_data/depth-windows.bed
    mosdepth --by report_data/depth-windows.bed --threads $task.cpus report_data/target aligned.bam
    { printf 'ref\\tstart\\tend\\tdepth\\n'; zcat report_data/target.regions.bed.gz | awk 'BEGIN { OFS="\\t" } { print \$1, \$2, \$3, \$NF }'; } \
        | bgzip -c > report_data/per-window-depth.tsv.gz
    mv report_data '${meta.alias}'
    """
}

process software_versions {
    label 'wfamplicon'
    cpus 1
    memory '1 GB'
    publishDir params.out_dir, mode: 'copy', pattern: 'versions.txt'
    output:
        path 'versions.txt'
    script:
    """
    {
        echo "minimap2,\$(minimap2 --version | head -n 1)"
        echo "samtools,\$(samtools --version | head -n 1 | awk '{print \$2}')"
        echo "fastcat,\$(fastcat --version 2>&1 | head -n 1)"
        if command -v porechop >/dev/null 2>&1; then echo "porechop,\$(porechop --version 2>&1 | head -n 1)"; fi
        if command -v dorado >/dev/null 2>&1; then echo "dorado,\$(dorado --version 2>&1 | head -n 1)"; fi
    } > versions.txt
    """
}

process make_alignment_report {
    label 'wf_common'
    cpus 1
    memory '4 GB'
    publishDir params.out_dir, mode: 'copy', pattern: 'nf-amplicon-alignment-report.html'
    input:
        path sample_reports
        path reference
        path target_bed
        path params_json
        path versions
    output:
        path 'nf-amplicon-alignment-report.html'
    script:
    """
    workflow-glue report_alignment \
        --report-fname nf-amplicon-alignment-report.html \
        --data . \
        --reference '$reference' \
        --target-bed '$target_bed' \
        --downsampling-size ${params.reads_downsampling_size ?: 0} \
        --versions '$versions' \
        --params '$params_json' \
        --wf-version '${workflow.manifest.version}'
    """
}

process call_clair3 {
    label 'wf_human_snp'
    cpus 8
    memory '16 GB'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_snp.vcf.gz*'
    input:
        tuple val(meta), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed), path(report_data)
    output:
        tuple val(meta), path("${meta.alias}"), emit: report_dirs
        tuple val(meta), path(bam), path(reference), path(target_bed), path('clair3_out/merge_output.vcf.gz'), path('clair3_out/merge_output.vcf.gz.tbi'), path("${meta.alias}"), emit: calls
    script:
    def model = params.override_basecaller_cfg ?: params.clair3_model_path ?: ''
    def extra = params.clair3_extra_args ?: ''
    """
    mkdir -p '${meta.alias}'
    samtools faidx '$reference'
    model='${model ?: ''}'
    if [ -z "\$model" ] && [ -s '$report_data/post-trim-basecallers.tsv' ]; then
        model=\$(awk 'NR==2 { print \$2; exit }' '$report_data/post-trim-basecallers.tsv' | sed -E 's/^dna_r10\\.4\\.1_e8\\.2_/r1041_e82_/; s/@v([0-9]+)\\.([0-9]+)\\.([0-9]+)/_v\\1\\2\\3/')
    fi
    model=\$(printf '%s' "\$model" | sed -E 's/^dna_r10\\.4\\.1_e8\\.2_/r1041_e82_/; s/^r10\\.4\\.1_e8\\.2_/r1041_e82_/; s/@v([0-9]+)\\.([0-9]+)\\.([0-9]+)/_v\\1\\2\\3/')
    model="\${model:-r1041_e82_400bps_hac_v420}"
    export CONDA_PREFIX=\$(dirname "\$(dirname "\$(which clair3.py)")")
    find /opt -type d -name '*r1041*' 2>/dev/null | head -20 > available-models.txt
    if [ -d "/opt/custflow/epi2meuser/models/\$model" ]; then model="/opt/custflow/epi2meuser/models/\$model"; fi
    run_clair3.sh \\
        --bam_fn='$bam' \\
        --ref_fn='$reference' \\
        --threads=$task.cpus \\
        --platform='ont' \\
        --model_path="\$model" \\
        --bed_fn='$target_bed' \\
        --output='clair3_out' \\
        $extra
    tabix -f -p vcf clair3_out/merge_output.vcf.gz
    cp clair3_out/merge_output.vcf.gz '${meta.alias}/${meta.alias}.nf_snp.vcf.gz'
    tabix -f -p vcf '${meta.alias}/${meta.alias}.nf_snp.vcf.gz'
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\t%FILTER\\n' '${meta.alias}/${meta.alias}.nf_snp.vcf.gz' \\
        | awk 'BEGIN { OFS="\\t"; print "target","position","ref","alt","genotype","filter","type" } { t=(length(\$3)==1 && length(\$4)==1 ? "SNV" : "INDEL"); print \$1,\$2,\$3,\$4,\$5,\$6,t }' \\
        > '${meta.alias}/variants.tsv'
    """
}

process phase_snp {
    label 'wf_human_snp'
    cpus 4
    memory '8 GB'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_snp*'
    publishDir "${params.out_dir}/${meta.alias}/alignments", mode: 'copy', pattern: '*.haplotagged.bam*'
    input:
        tuple val(meta), path(bam), path(reference), path(target_bed), path(vcf), path(vcf_index), path(report_dir, stageAs: 'clair3-report')
    output:
        tuple val(meta), path('phased.bam'), path('phased.vcf.gz'), path('phased.vcf.gz.tbi'), path(reference), path(target_bed), path(report_dir), emit: consensus_inputs
        tuple val(meta), path("${meta.alias}"), emit: report_dirs
    script:
    """
    samtools faidx '$reference'
    samtools index '$bam'
    whatshap phase --ignore-read-groups --reference '$reference' --output phased.vcf.gz '$vcf' '$bam'
    tabix -f -p vcf phased.vcf.gz
    whatshap haplotag --ignore-read-groups --reference '$reference' --output phased.bam phased.vcf.gz '$bam'
    samtools index phased.bam
    whatshap haplotagphase --ignore-read-groups --only-indels --reference '$reference' --output phased-with-indels.vcf phased.vcf.gz phased.bam
    bgzip -c phased-with-indels.vcf > phased.vcf.gz
    tabix -f -p vcf phased.vcf.gz
    whatshap stats --gtf haploblocks.gtf phased.vcf.gz
    cp phased.bam '${meta.alias}.haplotagged.bam'
    cp phased.bam.bai '${meta.alias}.haplotagged.bam.bai'
    cp phased.vcf.gz '${meta.alias}.nf_snp.vcf.gz'
    cp phased.vcf.gz.tbi '${meta.alias}.nf_snp.vcf.gz.tbi'
    cp haploblocks.gtf '${meta.alias}.nf_snp.haploblocks.gtf'
    mkdir -p '${meta.alias}'
    cp clair3-report/variants.tsv '${meta.alias}/variants.tsv'
    """
}

process make_haplotype_consensus {
    label 'wf_human_snp'
    cpus 2
    memory '4 GB'
    publishDir "${params.out_dir}/${meta.alias}/consensus", mode: 'copy', pattern: '*.fasta'
    input:
        tuple val(meta), path(bam), path(reference), path(target_bed), path(vcf), path(vcf_index), path(report_dir)
    output:
        tuple val(meta), path('consensus')
    script:
    def phasedArg = params.phased ? 'true' : 'false'
    """
    mkdir -p consensus
    samtools index '$bam'
    while IFS=\$'\\t' read -r chrom start end target; do
        [[ -z "$chrom" || "$chrom" == \\#* ]] && continue
        target_ref="target.fa"
        samtools faidx '$reference' "\$chrom:\$((start + 1))-\$end" > "\$target_ref"
        samtools depth -aa -r "\$chrom:\$((start + 1))-\$end" '$bam' \
            | awk -v c="\$chrom" -v s="\$start" -v cutoff=${params.min_coverage} 'BEGIN { OFS="\\t" } \$3 < cutoff { print c,\$2-s-1,\$2-s }' > low-depth.bed
        bedtools maskfasta -fi "\$target_ref" -bed low-depth.bed -fo masked.fa
        if [ '$phasedArg' = true ]; then
            bcftools consensus -r "\$chrom:\$((start + 1))-\$end" -f masked.fa -H 1 '$vcf' \
                > "consensus/${meta.alias}_\${target}_hap1.fasta"
            bcftools consensus -r "\$chrom:\$((start + 1))-\$end" -f masked.fa -H 2 '$vcf' \
                > "consensus/${meta.alias}_\${target}_hap2.fasta"
        else
            bcftools consensus --iupac-codes -r "\$chrom:\$((start + 1))-\$end" -f masked.fa '$vcf' \
                > "consensus/${meta.alias}_\${target}.consensus.fasta"
        fi
    done < '$target_bed'
    """
}

process make_snp_report {
    label 'wf_common'
    cpus 1
    memory '4 GB'
    publishDir params.out_dir, mode: 'copy', pattern: 'nf-amplicon-snp-report.html'
    input:
        path sample_reports
        path params_json
        path versions
    output:
        path 'nf-amplicon-snp-report.html'
    script:
    """
    workflow-glue report_snp \\
        --report-fname nf-amplicon-snp-report.html \\
        --data . \\
        --versions '$versions' \\
        --params '$params_json' \\
        --wf-version '${workflow.manifest.version}'
    """
}
