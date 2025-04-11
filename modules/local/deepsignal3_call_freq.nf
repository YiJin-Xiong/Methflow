process DEEPSIGNAL3_CALL_FREQ {
    label 'deepsignal3'

    queue true

    publishDir "${params.outdir}/${sample_id}/modcall/",
        mode: "copy",
        pattern: "*_pod5.CG.call_mods.frequency.tsv"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(call_mods_tsv_dir)

    output:
    path("${sample_id}_pod5.CG.call_mods.frequency.tsv"), emit: deepsignal3_call_freq

    script:
    
    """
    deepsignal3 call_freq \
        --input_path ${call_mods_tsv_dir} \
        --result_file ${sample_id}_pod5.CG.call_mods.frequency.tsv \
        --bed --sort
    
    """
}
