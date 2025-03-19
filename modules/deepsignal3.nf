process Deepsignal3_call_mods {
    label 'deepsignal3'

    publishDir "${params.outdir}/modcall/",
        mode: "copy",
        pattern: "*_pod5.CG.call_mods.tsv"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), val(species), val(sample_rate)
    path(basecall_bam_dir)

    output:
    path("${sample_id}_pod5.CG.call_mods.tsv"), emit: deepsignal3_callmods

    script:
    def deepsignal3_model = ''

    if ( analyte_type.length()!=0 ){
        if ( species.toLowerCase() == 'human' ){
            if (sample_rate.toLowerCase() == '4k'){
                deepsignal3_model = "human_r1041_4khz_CG_epoch7.ckpt"
            }
            else if (sample_rate.toLowerCase() == '5k'){
                deepsignal3_model = "human_r1041_5khz_CG_epoch5.ckpt"
            }
            else{
                deepsignal3_model = "human_r1041_4khz_CG_epoch7.ckpt"
            }
        }
        else if ( species.toLowerCase() == 'plant' ){
            if (sample_rate.toLowerCase() == '4k'){
                deepsignal3_model = "plant_r1041_4khz_C_epoch7.ckpt"
            }
            else if (sample_rate.toLowerCase() == '5k'){
                deepsignal3_model = "plant_r1041_5khz_C_epoch4.ckpt"
            }
            else{
                deepsignal3_model = "plant_r1041_5khz_C_epoch4.ckpt"
            }
        } else{
            deepsignal3_model = "human_r1041_4khz_CG_epoch7.ckpt"
        }
    } else{
        deepsignal3_model = "human_r1041_4khz_CG_epoch7.ckpt"
    }

    def deepsignal3_model_dir = file("models/deepsignal3_models/" + deepsignal3_model)
    println(deepsignal3_model_dir)
    
    """
    deepsignal3 --pod5 call_mods \
        --input_path ${pod5_dir} \
        --bam ${basecall_bam_dir} \
        --model_path ${deepsignal3_model_dir} \
        --result_file ${sample_id}_pod5.CG.call_mods.tsv \
        --nproc 16 --nproc_gpu 4 \
        --seq_len 21 --signal_len 15 -b 8192
    
    """
}

process Deepsignal3_call_freq {
    label 'deepsignal3'

    publishDir "${params.outdir}/modcall/",
        mode: "copy",
        pattern: "*_pod5.CG.call_mods.frequency.tsv"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), val(species), val(sample_rate)
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
