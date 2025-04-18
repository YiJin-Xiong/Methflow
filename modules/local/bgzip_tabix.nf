process BGZIP_TABIX {
    label 'bgzip_tabix'
    tag "${bed_file}"

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/DMR/haplotype/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(sample_id), path(bed_file)

    output:
    path("${bed_file}.gz", emit: bed_gz)
    path("${bed_file}.gz.tbi", emit: bed_gz_tbi)

    script:

    """
    bgzip -k ${bed_file}
    tabix -p bed ${bed_file}.gz
    
    """
}
