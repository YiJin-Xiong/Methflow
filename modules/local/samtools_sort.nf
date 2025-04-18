process SAMTOOLS_SORT {
    label 'samtools_sort'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/alignment/",
        mode: "copy",
        pattern: "${sample_id}_aligned_sorted.bam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_bam)

    output:
    path("${sample_id}_aligned_sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools sort -m 4G \
        -@ 20 ${aligned_bam} \
        -o ${sample_id}_aligned_sorted.bam \
    
    """
}