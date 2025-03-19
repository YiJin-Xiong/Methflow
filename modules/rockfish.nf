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

    def motif = params.motif ? "--motif $params.motif" : ''
    def idx = params.idx ? "--idx $params.idx" : ''
    def mapq_filter = params.mapq_filter ? "--mapq_filter $params.mapq_filter" : ''

    def rockfish_window = params.rockfish_window ? "--window $params.rockfish_window" : ''
    println(rockfish_window)
    
    """
    rockfish inference \
        -i ${pod5_dir} \
        --bam_path ${basecall_bam_dir} \
        --model_path ${rockfish_model_dir} -r \
        ${motif} \
        ${idx} \
        ${mapq_filter} \
        ${rockfish_window} \
        -o ${sample_id}_predictions.tsv
    
    """
}
