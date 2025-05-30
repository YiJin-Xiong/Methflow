process DORADO_BASECALL {
    // tag "$meta.id"
    label 'dorado_basecall'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/basecall/",
        mode: "copy",
        pattern: "*_basecall.bam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)

    output:
    path("${sample_id}_basecall.bam"), emit: basecall_bam

    script:
    def model            = ''
    def modified_bases   = ''
    if ( analyte_type.length()!=0 ){
        if ( analyte_type.toLowerCase() == 'dna' ){
            modified_bases = "5mCG_5hmCG"
            if ( pore_type == 'r9' || pore_type == 'R9' ){
                model = "dna_r9.4.1_e8_hac@v3.3"
            }
            else if ( pore_type == 'r10' || pore_type == 'R10' ){
                model = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
            }
        }
        else if ( analyte_type.substring(0,6).toLowerCase() == 'rna002' ){
            modified_bases = "m6A,m5C"
            model = "rna002_70bps_hac@v3"
        }
        else if ( analyte_type.substring(0,6).toLowerCase() == 'rna004' ){
            modified_bases = "m6A,m5C"
            model = "rna004_130bps_hac@v5.1.0"
        } else{
            model = "hac"
        }
    } else{
        model = "hac"
    }

    // modified_bases = modified_bases ? "--modified-bases $modified_bases" : ''
    modified_bases       = params.modified_bases ? "--modified-bases $params.modified_bases" : "--modified-bases $modified_bases"
    def dorado_model_dir = model=="hac"          ? model : file("models/dorado_models/" + model)
    def max_reads        = params.max_reads      ? "--max-reads $params.max_reads" : ''
    def read_ids         = params.read_ids       ? "--read-ids $params.read_ids" : ''
    def min_qscore       = params.min_qscore     ? "--min-qscore $params.min_qscore" : ''
    def trim             = params.trim           ? "$params.trim" : "all"
    def kit_name         = params.kit_name       ? "--kit-name $params.kit_name" : ''
    def gpu_device       = params.gpu            ? "--device cuda:all" : ''
    //println(dorado_model_dir + modified_bases + pod5_dir + max_reads + read_ids + trim + "\n" + task.cpus)

    ///dorado download \
    ///    --model ${model} --models-directory models/dorado_models
    //--batchsize 128 \
    """
    dorado basecaller ${dorado_model_dir} ${pod5_dir} \
        ${modified_bases} \
        ${max_reads} \
        ${read_ids} \
        ${min_qscore} \
        ${kit_name} \
        --trim ${trim} \
        ${gpu_device} \
        --emit-moves \
        --batchsize 128 \
        > ${sample_id}_basecall.bam
    """
}