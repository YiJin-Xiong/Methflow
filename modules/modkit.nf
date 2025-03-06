process Modkit_pileup {
    label 'Modkit_pileup'

    publishDir "${params.outdir}/DMR/haplotype/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(haplotagged_bam)

    output:
    env 'bed_num', emit: bed_num
    path "*.bed", emit: bed_files

    script:
    //nextflow run main.nf --input input_sheet.tsv --modkit_input results/whatshap/demo_whatshap_haplotagged.bam --region chr1 --filter_threshold 0.9 --mod_thresholds h:0.8 --modkit_motif "CGCG 0"

    def region = params.region ? "--region $params.region" : ''
    def filter_threshold = params.filter_threshold ? "--filter-threshold $params.filter_threshold" : ''
    def mod_thresholds = params.mod_thresholds ? "--mod-thresholds $params.mod_thresholds" : ''
    def ignore = params.ignore ? "--ignore $params.ignore" : '--ignore h'
    def motif = params.modkit_motif ? "--motif $params.modkit_motif" : ''

    """
    samtools index ${haplotagged_bam}

    samtools faidx ${genome}

    modkit pileup ${haplotagged_bam} . \
      --ref ${genome} \
      --cpg --combine-strands \
      --partition-tag HP \
      --prefix ${sample_id}_haplotype \
      ${ignore} \
      ${region} ${filter_threshold} ${mod_thresholds} ${motif}
    
    bed_num=\$(ls -1 *.bed | wc -l)
    echo \$bed_num
    
    """
}

process Modkit_DMR {
    label 'Modkit_DMR'

    publishDir "${params.outdir}/DMR/",
        mode: "copy",
        pattern: "*"

    input:
    val bed_num
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(bed_hp1)
    path(bed_hp2)

    output:
    path("${sample_id}_1&2_modkit_dmr.bed"), emit: modkit_dmr_bed
    path("${sample_id}_1&2_modkit_dmr.log"), emit: modkit_dmr_log


    script:
    //nextflow run main.nf --input input_sheet.tsv --norm_pileup_bed results/DMR/bed/demo_modkit_pileup.bed.gz --tumor_pileup_bed results/DMR/bed/demo1_modkit_pileup.bed.gz
    def dmr_regions = params.dmr_regions ? "--region " + file(params.dmr_regions) : ''

    """
    tabix -p bed ${bed_hp1}
    tabix -p bed ${bed_hp2}

    modkit dmr pair \
      -a ${bed_hp1} \
      -b ${bed_hp2} \
      -o "${sample_id}_1&2_modkit_dmr.bed" \
      --ref ${genome} \
      ${dmr_regions} \
      --base C \
      --threads 4 \
      -f --header \
      --log-filepath "${sample_id}_1&2_modkit_dmr.log"
      
    """
}