process DORADO_ALIGN {
    label 'dorado_align'

    conda "${params.project_dir}/environment.yml"

    publishDir "${params.outdir}/${sample_id}/basecall/",
        mode: "copy",
        pattern: "*_aligned.bam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(genome)
    path(basecall_bam_dir)

    output:
    path("${sample_id}_aligned.bam"), emit: alignment

    script:
    //如果用户想自己设置align的参数，如-k 15 -w 10，可以输入参数--mm2_opts
    def mm2_opts = params.mm2_opts ? "--mm2_opts \"$params.mm2_opts\"" : ''
    //println( basecall_bam_dir + genome + mm2_opts)

    """
    dorado aligner ${genome} ${basecall_bam_dir} \
        ${mm2_opts} \
        > ${sample_id}_aligned.bam
    """
}