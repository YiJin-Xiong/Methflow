process MODKIT_PILEUP {
    label 'Modkit_pileup'

    conda "${params.project_dir}/environment.yml"

    publishDir "${params.outdir}/${sample_id}/DMR/haplotype/",
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

    def pileup_region    = params.pileup_region    ? "--region $params.pileup_region" : ''
    def filter_threshold = params.filter_threshold ? "--filter-threshold $params.filter_threshold" : ''
    def mod_thresholds   = params.mod_thresholds   ? "--mod-thresholds $params.mod_thresholds" : ''
    def ignore           = params.ignore           ? "--ignore $params.ignore" : '--ignore h'
    def motif            = params.modkit_motif     ? "--motif $params.modkit_motif" : ''

    """
    samtools index ${haplotagged_bam}

    samtools faidx ${genome}

    modkit pileup ${haplotagged_bam} . \
      --ref ${genome} \
      --cpg --combine-strands \
      --partition-tag HP \
      --prefix ${sample_id}_haplotype \
      ${ignore} \
      ${pileup_region} ${filter_threshold} ${mod_thresholds} ${motif}
    
    bed_num=\$(ls -1 *.bed | wc -l)
    echo \$bed_num
    
    """
}