process Modkit_pileup {
    label 'Modkit_pileup'

    publishDir "${params.outdir}/DMR/bed",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)

    output:
    path("${sample_id}_modkit_pileup.bed"), emit: modkit_pileup_bed
    path("${sample_id}_modkit_pileup.log"), emit: modkit_pileup_log
    path("${sample_id}_modkit_pileup.bed.gz"), emit: modkit_pileup_bed_gz
    path("${sample_id}_modkit_pileup.bed.gz.tbi"), emit: modkit_pileup_bed_gz_tbi


    script:
    //nextflow run main.nf --input input_sheet.tsv --modkit_input results/align/demo_aligned_sorted.bam --region chr1 --filter_threshold 0.9 --mod_thresholds h:0.8 --modkit_motif "CGCG 0"

    def region = params.region ? "--region $params.region" : ''
    def filter_threshold = params.filter_threshold ? "--filter-threshold $params.filter_threshold" : ''
    def mod_thresholds = params.mod_thresholds ? "--mod-thresholds $params.mod_thresholds" : ''
    def ignore = params.ignore ? "--ignore $params.ignore" : '--ignore h'
    def motif = params.modkit_motif ? "--motif $params.modkit_motif" : ''

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    modkit pileup ${aligned_sorted_bam} ${sample_id}_modkit_pileup.bed \
      --ref ${genome} \
      --cpg --combine-strands \
      ${ignore} \
      ${region} ${filter_threshold} ${mod_thresholds} ${motif} \
      --log-filepath ${sample_id}_modkit_pileup.log
    
    bgzip -k ${sample_id}_modkit_pileup.bed
    tabix -p bed ${sample_id}_modkit_pileup.bed.gz
      
    """
}

process Modkit_DMR {
    label 'Modkit_DMR'

    publishDir "${params.outdir}/DMR/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(norm_pileup_bed)
    path(tumor_pileup_bed)

    output:
    path("${norm_id}_${tumor_id}_modkit_dmr.bed"), emit: modkit_pileup_bed
    path("${norm_id}_${tumor_id}_modkit_dmr.log"), emit: modkit_pileup_log


    script:
    //nextflow run main.nf --input input_sheet.tsv --norm_pileup_bed results/DMR/bed/demo_modkit_pileup.bed.gz --tumor_pileup_bed results/DMR/bed/demo1_modkit_pileup.bed.gz
    norm_index = norm_pileup_bed.toString().indexOf("_modkit")
    norm_id = norm_pileup_bed.toString().substring(0, norm_index)
    tumor_index = tumor_pileup_bed.toString().indexOf("_modkit")
    tumor_id = tumor_pileup_bed.toString().substring(0, tumor_index)

    def regions = params.regions ? "--region " + file(params.regions) : ''

    """
    tabix -p bed ${norm_pileup_bed}
    tabix -p bed ${tumor_pileup_bed}

    modkit dmr pair \
      -a ${norm_pileup_bed} \
      -b ${tumor_pileup_bed} \
      -o ${norm_id}_${tumor_id}_modkit_dmr.bed \
      --ref ${genome} \
      ${regions} \
      --base C \
      --threads 4 \
      -f --header \
      --log-filepath ${norm_id}_${tumor_id}_modkit_dmr.log
      
    """
}