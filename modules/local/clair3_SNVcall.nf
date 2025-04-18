process CLAIR3_SNVCALL {
    label 'Clair3'

    conda     (params.enable_conda ? "${params.project_dir}/environment-clair3.yml" : null)
    container (params.use_docker ? "${params.clair3_docker_name}" : "${params.clair3_singularity_name}")
    //container "hkubal/clair3:latest"

    publishDir "${params.outdir}/${sample_id}/SNVcall/",
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
    def clar3_model_path = file(params.clar3_model_path)
    def include_all_ctgs = ( species == 'Human') ? "" : "--include_all_ctgs"
    def bed_fn           = params.bed_fn         ? "--bed_fn=" + file(params.bed_fn) : ''
    def vcf_fn           = params.vcf_fn         ? "--vcf_fn=" + file(params.vcf_fn) : ''
    def ctg_name         = params.ctg_name       ? "--ctg_name=" + params.ctg_name : ''
    def sample_name      = params.sample_name    ? "--sample_name=" + params.sample_name : ''
    def qual             = params.qual           ? "--qual=" + params.qual : ''

    //--whatshap="${params.whatshap_path}" \
    // /opt/bin/run_clair3.sh \

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    /opt/bin/run_clair3.sh \
        --bam_fn=${aligned_sorted_bam} \
        --ref_fn=${genome} \
        --threads=8 \
        --platform="ont" \
        --model_path=${clar3_model_path} \
        --output="temp" \
        ${bed_fn} ${vcf_fn} ${ctg_name} ${sample_name} ${qual} \
        ${include_all_ctgs}
    
    cp temp/merge_output.vcf.gz ${sample_id}_clair3_merge.vcf.gz
    cp temp/pileup.vcf.gz.tbi ${sample_id}_clair3_pileup.vcf.gz.tbi

    rm -r temp

    """
}