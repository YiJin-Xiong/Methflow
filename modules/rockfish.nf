process Rockfish_inference {
    label 'rockfish'

    publishDir "${params.outdir}/modcall/",
        mode: "copy",
        pattern: "*_predictions.tsv"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir)
    path(basecall_bam_dir)

    output:
    path("${sample_id}_predictions.tsv"), emit: rockfish_predictions

    script:
    def rockfish_model_dir = file(params.rockfish_model)
    
    """
    rockfish inference \
        -i ${pod5_dir} \
        --bam_path ${basecall_bam_dir} \
        --model_path ${rockfish_model_dir} -r \
        -o ${sample_id}_predictions.tsv
    
    """
}
