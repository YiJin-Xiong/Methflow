process WHATSHAP_PHASE {
    label 'WhatsHap_phase'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/phase/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)
    path(clair3_vcf)

    output:
    path("${sample_id}_SNV_PASS.vcf"), emit: SNV_pass_vcf
    path("${sample_id}_phased.vcf.gz"), emit: phased_vcf_gz
    
    script:
    def phase_sample             = params.phase_sample     ? "--sample "     + params.phase_sample : ''
    def phase_chromosome         = params.phase_chromosome ? "--chromosome " + params.phase_chromosome : ''
    def phase_exclude_chromosome = params.phase_exclude_chromosome ? "--exclude-chromosome " + params.phase_exclude_chromosome : ''
    def phase_mapq               = params.phase_mapq       ? "--mapq "       + params.phase_mapq : ''
    def only_snvs                = params.only_snvs           ? "--only-snvs" : ''

    //samtools index ${aligned_sorted_bam}
    //samtools faidx ${genome}

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    gunzip -c ${clair3_vcf} | \
        awk '/^#/ || (\$4 != "." && \$5 != "." && length(\$4) == 1 && length(\$5) == 1 && \$7 =="PASS")' - \
        > ${sample_id}_SNV_PASS.vcf
    
    whatshap phase --ignore-read-groups \
        -o ${sample_id}_phased.vcf \
        --reference=${genome} \
        --tag=HP \
        ${phase_sample} ${phase_chromosome} ${phase_mapq} \
        ${only_snvs} ${phase_exclude_chromosome} \
        ${sample_id}_SNV_PASS.vcf ${aligned_sorted_bam}

    bgzip -@ 4 ${sample_id}_phased.vcf

    """
}
