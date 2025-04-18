process DSS_PREPROCESS {
    label 'DSS_preprocess'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/DSS/",
        mode: "copy",
        pattern: "*"

    input:
    val bed_num
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(bed_hp1)
    path(bed_hp2)

    output:
    path("${sample_id}_haplotype_preprocessed_1.bed"), emit: bed_preprocessed_hp1
    path("${sample_id}_haplotype_preprocessed_2.bed"), emit: bed_preprocessed_hp2
    
    script:
    def coverage = ( params.DSS_coverage ) ? "--DSS_coverage "+ params.DSS_coverage : 1

    """
    gunzip -c ${bed_hp1} \
      | awk 'BEGIN{FS="\t"; OFS="\t"} {print \$1, \$3, \$10, \$12}' \
      | awk -F '\t' '\$3 >= ${coverage}' \
      | awk 'BEGIN {print "chro\tpos\tN\tX"} {print}' \
      > ${sample_id}_haplotype_preprocessed_1.bed
    
    gunzip -c ${bed_hp2} \
      | awk 'BEGIN{FS="\t"; OFS="\t"} {print \$1, \$3, \$10, \$12}' \
      | awk -F '\t' '\$3 >= ${coverage}' \
      | awk 'BEGIN {print "chro\tpos\tN\tX"} {print}' \
      > ${sample_id}_haplotype_preprocessed_2.bed

    """
}