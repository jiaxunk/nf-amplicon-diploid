#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import WorkflowMain

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

workflow {
    WorkflowMain.initialise(workflow, params, log)

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
            keep_unaligned: true,
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

        bam_to_fastq(ch_bam_routes.remap)
        ch_fastq = bam_to_fastq.out | filter_downsample
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

    ch_repeats = params.cuban_repeats
        ? Channel.fromPath(params.cuban_repeats, checkIfExists: true)
        : Channel.value([])

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
            ch_raw_consensus_inputs = phase_snp.out.consensus_inputs
            ch_raw_snp_reports = phase_snp.out.report_dirs
            ch_sniffles_bam = phase_snp.out.consensus_inputs.map { meta, bam, ref, target, vcf, tbi, rep ->
                [meta, bam, file("${bam}.bai"), true]
            }
        } else {
            ch_raw_consensus_inputs = call_clair3.out.calls
            ch_raw_snp_reports = call_clair3.out.report_dirs
            ch_sniffles_bam = ch_results.map { meta, bam, bai, report_data ->
                [meta, bam, bai, false]
            }
        }

        if (params.sv) {
            ch_sniffles_in = ch_sniffles_bam
                .combine(ch_reference)
                .map { meta, bam, bai, is_phased, ref, mmi, fai ->
                    [meta, bam, bai, ref, fai, is_phased]
                }
            call_sniffles(ch_sniffles_in)

            ch_filter_in = call_sniffles.out.raw_calls
                .combine(ch_target)
                .map { meta, raw_vcf, raw_tbi, target ->
                    [meta, raw_vcf, raw_tbi, target]
                }
            filter_sv_calls(ch_filter_in)

            ch_refine_in = ch_raw_consensus_inputs
                .map { meta, bam, ref, target, vcf, tbi, rep ->
                    [meta, bam, file("${bam}.bai"), ref, target, vcf, tbi]
                }
                .join(filter_sv_calls.out.filtered.map { meta, sv_vcf, sv_tbi, sv_tsv -> [meta, sv_vcf, sv_tbi] })
            refine_snp_zygosity(ch_refine_in)

            if (params.consensus) {
                make_haplotype_consensus(refine_snp_zygosity.out.consensus_inputs)
            }
            ch_snp_reports = refine_snp_zygosity.out.report_dirs.map { meta, dir -> dir }.collect()

            ch_cuban_in = ch_sniffles_bam
                .map { meta, bam, bai, is_phased -> [meta, bam, bai] }
                .join(filter_sv_calls.out.filtered)
                .join(ch_results.map { meta, bam, bai, rep -> [meta, rep] })
                .combine(ch_reference)
                .map { meta, bam, bai, sv_vcf, sv_tbi, sv_tsv, rep, ref, mmi, fai ->
                    [meta, bam, bai, sv_vcf, sv_tbi, sv_tsv, ref, fai, rep]
                }
            plot_cuban(ch_cuban_in, ch_repeats)

            ch_sv_reports = plot_cuban.out.report_dir.map { meta, dir -> dir }.collect()
            make_sv_report(
                ch_sv_reports,
                ch_params,
                ch_versions
            )
        } else {
            if (params.consensus) {
                make_haplotype_consensus(ch_raw_consensus_inputs)
            }
            ch_snp_reports = ch_raw_snp_reports.map { meta, dir -> dir }.collect()
        }

        make_snp_report(
            ch_snp_reports,
            ch_params,
            ch_versions
        )
    } else if (params.sv) {
        ch_sniffles_bam = ch_results.map { meta, bam, bai, report_data ->
            [meta, bam, bai, false]
        }
        ch_sniffles_in = ch_sniffles_bam
            .combine(ch_reference)
            .map { meta, bam, bai, is_phased, ref, mmi, fai ->
                [meta, bam, bai, ref, fai, is_phased]
            }
        call_sniffles(ch_sniffles_in)

        ch_filter_in = call_sniffles.out.raw_calls
            .combine(ch_target)
            .map { meta, raw_vcf, raw_tbi, target ->
                [meta, raw_vcf, raw_tbi, target]
            }
        filter_sv_calls(ch_filter_in)

        ch_cuban_in = ch_sniffles_bam
            .map { meta, bam, bai, is_phased -> [meta, bam, bai] }
            .join(filter_sv_calls.out.filtered)
            .join(ch_results.map { meta, bam, bai, rep -> [meta, rep] })
            .combine(ch_reference)
            .map { meta, bam, bai, sv_vcf, sv_tbi, sv_tsv, rep, ref, mmi, fai ->
                [meta, bam, bai, sv_vcf, sv_tbi, sv_tsv, ref, fai, rep]
            }
        plot_cuban(ch_cuban_in, ch_repeats)

        ch_sv_reports = plot_cuban.out.report_dir.map { meta, dir -> dir }.collect()
        make_sv_report(
            ch_sv_reports,
            ch_params,
            ch_versions
        )
    }
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

process bam_to_fastq {
    label 'wfamplicon'
    cpus 2
    memory '2 GB'
    input:
        tuple val(meta), path(bam), path(bai), path(stats), val(match)
    output:
        tuple val(meta), path('reads.fastq.gz'), path(stats)
    script:
    """
    samtools fastq -@ $task.cpus -T '*' '$bam' | gzip -c > reads.fastq.gz
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
        echo "sniffles,2.6.3"
        echo "cuban,0.1.0"
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
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_snp.vcf.gz*', enabled: !params.phased && !params.sv
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: 'variants.tsv', enabled: !params.phased && !params.sv
    input:
        tuple val(meta), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed), path(report_data)
    output:
        tuple val(meta), path("${meta.alias}"), emit: report_dirs
        tuple val(meta), path(bam), path(reference), path(target_bed), path('clair3_out/merge_output.vcf.gz'), path('clair3_out/merge_output.vcf.gz.tbi'), path("${meta.alias}"), emit: calls
        tuple val(meta), path("${meta.alias}.nf_snp.vcf.gz*"), optional: true, emit: unphased_vcf
        tuple val(meta), path("variants.tsv"), optional: true, emit: unphased_tsv
    script:
    def model = params.override_basecaller_cfg ?: (params.containsKey('clair3_model') && params.clair3_model ? params.clair3_model : params.clair3_model_path) ?: ''
    def extra = params.clair3_extra_args ?: ''
    def phasedArg = params.phased ? 'true' : 'false'
    """
    mkdir -p '${meta.alias}'
    samtools faidx '$reference'
    model='${model ?: ''}'
    if [ -z "\$model" ] && [ -s '$report_data/post-trim-basecallers.tsv' ]; then
        model=\$(awk 'NR==2 { print \$2; exit }' '$report_data/post-trim-basecallers.tsv')
    fi
    if [ -n "\$model" ]; then
        model=\$(printf '%s' "\$model" | sed -E 's/^dna_//; s/^r10\\.4\\.1_e8\\.2_/r1041_e82_/; s/@v([0-9]+)\\.([0-9]+)\\.([0-9]+)/_v\\1\\2\\3/; s/\\.([0-9]+)\$/\\1/')
    fi
    model="\${model:-r1041_e82_400bps_hac_v420}"
    export CONDA_PREFIX=\$(dirname "\$(dirname "\$(which clair3.py)")")
    if [ -d "/opt/custflow/epi2meuser/models/\$model" ]; then
        model="/opt/custflow/epi2meuser/models/\$model"
    fi
    run_clair3.sh \\
        --bam_fn='$bam' \\
        --ref_fn='$reference' \\
        --threads=$task.cpus \\
        --platform='ont' \\
        --model_path="\$model" \\
        --bed_fn='$target_bed' \\
        --output='clair3_out' \\
        --include_all_ctgs \\
        $extra
    tabix -f -p vcf clair3_out/merge_output.vcf.gz
    if [ '$phasedArg' = false ]; then
        cp clair3_out/merge_output.vcf.gz '${meta.alias}.nf_snp.vcf.gz'
        cp clair3_out/merge_output.vcf.gz.tbi '${meta.alias}.nf_snp.vcf.gz.tbi'
        workflow-glue extract_snp_data \\
            --bam '$bam' \\
            --vcf clair3_out/merge_output.vcf.gz \\
            --target '$target_bed' \\
            --out-dir '${meta.alias}'
        cp '${meta.alias}/variants.tsv' 'variants.tsv'
    fi
    """
}

process phase_snp {
    label 'wf_human_snp'
    cpus 4
    memory '8 GB'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_snp.vcf.gz*', enabled: !params.sv
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: 'variants.tsv', enabled: !params.sv
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.haploblocks.gtf'
    publishDir "${params.out_dir}/${meta.alias}/alignments", mode: 'copy', pattern: '*.haplotagged.bam*'
    input:
        tuple val(meta), path(bam), path(reference), path(target_bed), path(vcf), path(vcf_index), path(report_dir, stageAs: 'clair3-report')
    output:
        tuple val(meta), path('phased.bam'), path(reference), path(target_bed), path('phased.vcf.gz'), path('phased.vcf.gz.tbi'), path("${meta.alias}"), emit: consensus_inputs
        tuple val(meta), path("${meta.alias}"), emit: report_dirs
        tuple val(meta), path("${meta.alias}.haplotagged.bam*"), emit: haplotagged_bam
        tuple val(meta), path("${meta.alias}.nf_snp*"), emit: snp_vcf
        tuple val(meta), path("variants.tsv"), emit: variants_tsv
    script:
    """
    samtools faidx '$reference'
    samtools index '$bam'
    whatshap phase --indels --ignore-read-groups --reference '$reference' --output phased.vcf.gz '$vcf' '$bam'
    tabix -f -p vcf phased.vcf.gz
    whatshap haplotag --ignore-read-groups --skip-missing-contigs --reference '$reference' --output phased.bam phased.vcf.gz '$bam'
    samtools index phased.bam
    touch haploblocks.gtf
    whatshap stats --gtf haploblocks.gtf phased.vcf.gz || true
    cp phased.bam '${meta.alias}.haplotagged.bam'
    cp phased.bam.bai '${meta.alias}.haplotagged.bam.bai'
    cp phased.vcf.gz '${meta.alias}.nf_snp.vcf.gz'
    cp phased.vcf.gz.tbi '${meta.alias}.nf_snp.vcf.gz.tbi'
    cp haploblocks.gtf '${meta.alias}.nf_snp.haploblocks.gtf'
    mkdir -p '${meta.alias}'
    workflow-glue extract_snp_data \\
        --bam phased.bam \\
        --vcf phased.vcf.gz \\
        --target '$target_bed' \\
        --out-dir '${meta.alias}'
    cp '${meta.alias}/variants.tsv' 'variants.tsv'
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
        tuple val(meta), path('*.fasta'), emit: consensus
    script:
    def phasedArg = params.phased ? 'true' : 'false'
    """
    samtools faidx '$reference'
    samtools index '$bam'
    while IFS=\$'\\t' read -r chrom start end target; do
        [[ -z "\$chrom" || "\$chrom" == \\#* || "\$chrom" == track* || "\$chrom" == browser* ]] && continue
        target="\${target:-\${chrom}_\${start}_\${end}}"

        samtools depth -aa -r "\$chrom:\$((start + 1))-\$end" '$bam' \
            | awk -v c="\$chrom" -v cutoff=${params.min_coverage} 'BEGIN { OFS="\\t" } \$3 < cutoff { print c, \$2-1, \$2 }' > low-depth.bed

        if [ '$phasedArg' = true ]; then
            samtools faidx '$reference' "\$chrom:\$((start + 1))-\$end" \
                | bcftools consensus -i 'FILTER="PASS"' -H 1 -m low-depth.bed -f - '$vcf' \
                | sed "1s/.*/>${meta.alias}_\${target}_hap1/" \
                > "${meta.alias}_\${target}_hap1.fasta"
            samtools faidx '$reference' "\$chrom:\$((start + 1))-\$end" \
                | bcftools consensus -i 'FILTER="PASS"' -H 2 -m low-depth.bed -f - '$vcf' \
                | sed "1s/.*/>${meta.alias}_\${target}_hap2/" \
                > "${meta.alias}_\${target}_hap2.fasta"
        else
            samtools faidx '$reference' "\$chrom:\$((start + 1))-\$end" \
                | bcftools consensus -i 'FILTER="PASS"' --iupac-codes -m low-depth.bed -f - '$vcf' \
                | sed "1s/.*/>${meta.alias}_\${target}.consensus/" \
                > "${meta.alias}_\${target}.consensus.fasta"
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

process call_sniffles {
    label 'wf_human_sv'
    cpus 4
    memory '8 GB'
    input:
        tuple val(meta), path(bam), path(bai), path(reference), path(reference_fai), val(is_phased)
    output:
        tuple val(meta), path('raw_sv.vcf.gz'), path('raw_sv.vcf.gz.tbi'), emit: raw_calls
    script:
    def phase_arg = is_phased ? '--phase' : ''
    def extra = params.sniffles_extra_args ?: ''
    """
    samtools faidx '$reference'
    samtools index '$bam'
    sniffles \\
        --input '$bam' \\
        --vcf raw_sv.vcf.gz \\
        --reference '$reference' \\
        --minsvlen ${params.min_sv_length} \\
        --threads $task.cpus \\
        $phase_arg \\
        $extra
    tabix -f -p vcf raw_sv.vcf.gz
    """
}

process filter_sv_calls {
    label 'wfamplicon'
    cpus 1
    memory '2 GB'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_sv*'
    input:
        tuple val(meta), path(raw_vcf), path(raw_vcf_index), path(target_bed)
    output:
        tuple val(meta), path("${meta.alias}.nf_sv.vcf.gz"), path("${meta.alias}.nf_sv.vcf.gz.tbi"), path("${meta.alias}.nf_sv.tsv"), emit: filtered
    script:
    """
    workflow-glue filter_sv_calls \\
        --vcf '$raw_vcf' \\
        --target '$target_bed' \\
        --min-support ${params.min_sv_support} \\
        --min-vaf ${params.min_sv_vaf} \\
        -o '${meta.alias}.nf_sv.vcf.gz' \\
        --output-tsv '${meta.alias}.nf_sv.tsv'
    """
}

process refine_snp_zygosity {
    label 'wf_human_snp'
    cpus 2
    memory '4 GB'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: '*.nf_snp*'
    publishDir "${params.out_dir}/${meta.alias}/variants", mode: 'copy', pattern: 'variants.tsv'
    input:
        tuple val(meta), path(bam), path(bai), path(reference), path(target_bed), path(snp_vcf), path(snp_vcf_index), path(sv_vcf), path(sv_vcf_index)
    output:
        tuple val(meta), path(bam), path(reference), path(target_bed), path("${meta.alias}.nf_snp.vcf.gz"), path("${meta.alias}.nf_snp.vcf.gz.tbi"), path("${meta.alias}"), emit: consensus_inputs
        tuple val(meta), path("${meta.alias}"), emit: report_dirs
        tuple val(meta), path("${meta.alias}.nf_snp*"), emit: snp_vcf
        tuple val(meta), path("variants.tsv"), emit: variants_tsv
    script:
    """
    samtools faidx '$reference'
    samtools index '$bam'
    export PYTHONPATH=/opt/custflow/epi2meuser/conda/bin
    python3 /opt/custflow/epi2meuser/conda/bin/postprocess/SwitchZygosityBasedOnSVCalls.py \\
        --bam_fn '$bam' \\
        --clair3_vcf_input '$snp_vcf' \\
        --sv_vcf_input '$sv_vcf' \\
        --sv_filter_tag PASS,GT \\
        --vcf_output refined.vcf \\
        --threads $task.cpus

    cp refined.vcf.gz '${meta.alias}.nf_snp.vcf.gz'
    cp refined.vcf.gz.tbi '${meta.alias}.nf_snp.vcf.gz.tbi'
    mkdir -p '${meta.alias}'
    workflow-glue extract_snp_data \\
        --bam '$bam' \\
        --vcf '${meta.alias}.nf_snp.vcf.gz' \\
        --target '$target_bed' \\
        --out-dir '${meta.alias}'
    cp '${meta.alias}/variants.tsv' 'variants.tsv'
    """
}

process plot_cuban {
    label 'wfamplicon'
    cpus 2
    memory '4 GB'
    publishDir "${params.out_dir}/${meta.alias}/sv_plots", mode: 'copy', pattern: "${meta.alias}/sv_plots/*.png", saveAs: { fn -> file(fn).name }
    input:
        tuple val(meta), path(bam), path(bai), path(sv_vcf), path(sv_vcf_index), path(sv_tsv), path(reference), path(reference_fai), path(report_data)
        path repeats_bed
    output:
        tuple val(meta), path("${meta.alias}"), emit: report_dir
        tuple val(meta), path("${meta.alias}/sv_plots/*.png"), optional: true, emit: plots
    script:
    def repeats_arg = params.cuban_repeats ? "--repeats '$repeats_bed'" : "--no-repeats"
    def padding_arg = params.cuban_padding ? "--padding ${params.cuban_padding}" : ""
    def extra = params.cuban_extra_args ?: ""
    """
    mkdir -p '${meta.alias}/variants' '${meta.alias}/sv_plots'
    cp '$sv_vcf' '${meta.alias}/variants/${meta.alias}.nf_sv.vcf.gz'
    cp '$sv_vcf_index' '${meta.alias}/variants/${meta.alias}.nf_sv.vcf.gz.tbi'
    cp '$sv_tsv' '${meta.alias}/variants/${meta.alias}.nf_sv.tsv'

    cov=\$(awk '\$1=="total_region" { if (\$4 > 0.1) print \$4; else print "20.0"; exit }' '$report_data/target.mosdepth.summary.txt' 2>/dev/null || echo "20.0")
    if [ -z "\$cov" ]; then cov="20.0"; fi

    n_vars=\$(zgrep -v '^#' '$sv_vcf' | wc -l || echo 0)
    if [ "\$n_vars" -gt 0 ]; then
        cuban \\
            --vcf '$sv_vcf' \\
            --sample '${meta.alias}:$bam' \\
            --tech '${meta.alias}:lr' \\
            --baseline-cov "${meta.alias}:\$cov" \\
            --ref '$reference' \\
            --outdir '${meta.alias}/sv_plots' \\
            --window ${params.cuban_window} \\
            $padding_arg \\
            $repeats_arg \\
            $extra
    fi
    """
}

process make_sv_report {
    label 'wf_common'
    cpus 1
    memory '4 GB'
    publishDir params.out_dir, mode: 'copy', pattern: 'nf-amplicon-sv-report.html'
    input:
        path sample_reports
        path params_json
        path versions
    output:
        path 'nf-amplicon-sv-report.html'
    script:
    """
    workflow-glue report_sv \\
        --report-fname nf-amplicon-sv-report.html \\
        --data . \\
        --versions '$versions' \\
        --params '$params_json' \\
        --wf-version '${workflow.manifest.version}'
    """
}
