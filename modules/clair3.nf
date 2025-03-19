process Clair3 {
    label 'Clair3'

    publishDir "${params.outdir}/vcf/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)

    output:
    //path("temp"), emit: clair3_output_dir
    path("${sample_id}_clair3_merge.vcf.gz"), emit: clair3_merge_vcf
    path("${sample_id}_clair3_pileup.vcf.gz.tbi"), emit: clair3_merge_vcf_tbi

    script:
    def include_all_ctgs = ( species == 'Human') ? "" : "--include_all_ctgs"
    def clar3_model_path = file(params.clar3_model_path)
    //--whatshap="${params.whatshap_path}" \

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    ~/tool/Clair3/run_clair3.sh \
        --bam_fn=${aligned_sorted_bam} \
        --ref_fn=${genome} \
        --threads=8 \
        --platform="ont" \
        --model_path=${clar3_model_path} \
        --output="temp" \
        ${include_all_ctgs}
    
    cp temp/merge_output.vcf.gz ${sample_id}_clair3_merge.vcf.gz
    cp temp/pileup.vcf.gz.tbi ${sample_id}_clair3_pileup.vcf.gz.tbi

    rm -r temp

    """
}