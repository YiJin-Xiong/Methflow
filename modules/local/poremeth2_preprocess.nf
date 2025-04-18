process POREMETH2_PREPROCESS {
    label 'poremeth2_preprocess'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/DMR/PoreMeth2/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(haplotagged_bam)

    output:
    //path("${sample_name}_extracted.tsv"), emit: extracted_tsv
    path("${sample_name}_extracted_sorted.entropy.file.tsv"), emit: extracted_entroyp_tsv
    
    script:
    sample_name = haplotagged_bam.name.split('_whatshap')[0]

    """
    modkit extract full \
        --mapped-only --cpg --force \
        --threads 10 \
        --reference ${genome} \
        ${haplotagged_bam} \
        ${sample_name}_extracted.tsv
    
    sh ${params.project_dir}/bin/ModkitResorter.sh ${sample_name}_extracted.tsv

    rm ${sample_name}_extracted.tsv

    """
}