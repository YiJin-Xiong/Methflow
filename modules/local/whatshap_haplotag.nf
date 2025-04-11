process WHATSHAP_HAPLOTAG {
    label 'WhatsHap_haplotag'

    conda "${params.project_dir}/environment.yml"

    publishDir "${params.outdir}/${sample_id}/whatshap/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)
    path(phased_vcf_gz)

    output:
    path("${sample_id}_phased.vcf.gz.tbi"), emit: phased_vcf_gz_tbi
    path("${sample_id}_whatshap_haplotagged.bam"), emit: whatshap_haplotagged_bam
    path("${sample_id}_whatshap_haplotagged.bam.bai"), emit: whatshap_haplotagged_bam_bai
    path("${sample_id}_whatshap_haplotagged.readlist"), emit: whatshap_haplotagged_readlist
    
    script:
    def ignore_linked_read    = params.ignore_linked_read    ? "--ignore-linked-read" : ''
    def tag_supplementary     = params.tag_supplementary     ? "--tag-supplementary" : ''
    def skip_missing_contigs  = params.skip_missing_contigs  ? "--skip-missing-contigs" : ''
    def haplotag_sample       = params.haplotag_sample       ? "--sample " + params.haplotag_sample : ''
    def ploidy                = params.ploidy                ? "--ploidy " + params.ploidy : ''

    """
    samtools faidx ${genome}

    samtools index ${aligned_sorted_bam}

    tabix -p vcf ${phased_vcf_gz}

    whatshap haplotag --ignore-read-groups \
        --output-haplotag-list ${sample_id}_whatshap_haplotagged.readlist \
        -o ${sample_id}_whatshap_haplotagged.bam \
        --reference ${genome} \
        ${ignore_linked_read} \
        ${tag_supplementary} \
        ${skip_missing_contigs} \
        ${haplotag_sample} ${ploidy} \
        ${phased_vcf_gz} ${aligned_sorted_bam}
    
    samtools index -@ 4 ${sample_id}_whatshap_haplotagged.bam

    """
}
