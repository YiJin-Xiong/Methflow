process Clair3 {
    label 'Clair3'

    publishDir "${params.outdir}/SNV/",
        mode: "copy",
        pattern: "*_aligned.sam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(aligned_sorted_bam)

    output:
    path("SNV"), emit: clair3_output

    script:
    def include_all_ctgs = ( species == 'Human') ? "" : "--include_all_ctgs"
    def clar3_model_path = file(params.clar3_model_path)

    """
    samtools index ${aligned_sorted_bam}

    samtools faidx ${genome}

    ~/tool/Clair3/run_clair3.sh \
        --bam_fn=${aligned_sorted_bam} \
        --ref_fn=${genome} \
        --threads=8 \
        --platform="ont" \
        --model_path=${clar3_model_path} \
        --output="SNV" \
        --whatshap="${params.whatshap_path}" \
        ${include_all_ctgs}

    """
}