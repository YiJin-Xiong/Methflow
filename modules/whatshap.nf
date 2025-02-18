process WhatsHap {
    label 'WhatsHap'

    publishDir "${params.outdir}/whatshap/",
        mode: "copy",
        pattern: "*.vcf"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)
    path(clair3_vcf)

    output:
    path("${sample_id}_SNV_PASS.vcf"), emit: SNV_pass_vcf
    path("${sample_id}_phased.vcf"), emit: SNV_phased_vcf

    script:

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    gunzip -c ${clair3_vcf} | \
        awk '/^#/ || (\$4 != "." && \$5 != "." && length(\$4) == 1 && length(\$5) == 1 && \$7 =="PASS")' - \
        > ${sample_id}_SNV_PASS.vcf
    
    whatshap phase --ignore-read-groups \
        -o ${sample_id}_phased.vcf \
        --reference=${genome} \
        ${sample_id}_SNV_PASS.vcf ${aligned_sorted_bam}


    """
}