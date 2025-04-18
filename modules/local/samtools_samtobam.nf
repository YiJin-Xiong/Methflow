process SAMTOOLS_SAMTOBAM {
    label 'samtools_sam2bam'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/alignment/",
        mode: "copy",
        pattern: "${sample_id}_aligned.bam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sam)

    output:
    path("${sample_id}_aligned.bam"), emit: aligned_bam

    script:
    """
    samtools view -@20 -bh \
        ${aligned_sam} \
        > ${sample_id}_aligned.bam
    
    """
}