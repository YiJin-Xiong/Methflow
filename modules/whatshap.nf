process WhatsHap {
    label 'WhatsHap'

    publishDir "${params.outdir}/whatshap/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)
    path(clair3_vcf)

    output:
    path("${sample_id}_SNV_PASS.vcf"), emit: SNV_pass_vcf
    path("${sample_id}_phased.vcf.gz"), emit: phased_vcf_gz
    path("${sample_id}_phased.vcf.gz.tbi"), emit: phased_vcf_gz_tbi
    path("${sample_id}_whatshap_haplotagged.bam"), emit: whatshap_haplotagged_bam
    path("${sample_id}_whatshap_haplotagged.bam.bai"), emit: whatshap_haplotagged_bam_bai
    path("${sample_id}_whatshap_haplotagged.readlist"), emit: whatshap_haplotagged_readlist
    
    script:
    def ignore_linked_read = params.ignore_linked_read ? "--ignore-linked-read" : ''
    def tag_supplementary = params.tag_supplementary ? "--tag-supplementary" : ''
    def haplotag_sample = params.haplotag_sample ? "--sample" + params.haplotag_sample : ''

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

    bgzip -@ 4 ${sample_id}_phased.vcf && \
        tabix -p vcf ${sample_id}_phased.vcf.gz
    
    whatshap haplotag --ignore-read-groups \
        --output-haplotag-list ${sample_id}_whatshap_haplotagged.readlist \
        -o ${sample_id}_whatshap_haplotagged.bam \
        --reference ${genome} \
        ${ignore_linked_read} \
        ${tag_supplementary} \
        ${haplotag_sample} \
        ${sample_id}_phased.vcf.gz ${aligned_sorted_bam}
    
    samtools index -@ 4 ${sample_id}_whatshap_haplotagged.bam

    """
}
