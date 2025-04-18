process MINIMAP2_ALIGN {
    label 'minimap2'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/alignment/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(basecalled_fastq)

    output:
    path("${sample_id}_aligned.sam"), emit: minimap2_aligned_sam

    script:
    //DNA -ax map-ont
    //RNA -ax splice -uf 
    def preset       = ( analyte_type == 'DNA' )  ? "-ax map-ont" : "-ax splice"
    def stranded     = ( analyte_type == 'RNA')   ? "-uf" : ""
    def kmer         = ( params.minimap2_kmer )   ? "-k " + params.minimap2_kmer : ""
    def window       = ( params.minimap2_window ) ? "-w " + params.minimap2_window : ""
    def minimap2_f   = ( params.minimap2_f )      ? "-f " + params.minimap2_f : ""
    def minimap2_g   = ( params.minimap2_g )      ? "-g " + params.minimap2_g : ""
    def minimap2_G   = ( params.minimap2_G )      ? "-G " + params.minimap2_G : ""
    def minimap2_n   = ( params.minimap2_n )      ? "-n " + params.minimap2_n : ""
    def minimap2_m   = ( params.minimap2_m )      ? "-m " + params.minimap2_m : ""
    def minimap2_p   = ( params.minimap2_p )      ? "-p " + params.minimap2_p : ""
    def minimap2_N   = ( params.minimap2_N )      ? "-N " + params.minimap2_N : ""
    def minimap2_A   = ( params.minimap2_A )      ? "-A " + params.minimap2_A : ""
    def minimap2_B   = ( params.minimap2_B )      ? "-B " + params.minimap2_B : ""

    // println("preset:" + preset + " stranded:" + stranded + " kmer:"  + kmer + " window:" + window)
    // println("minimap2_f:" + minimap2_f + " minimap2_g:" + minimap2_g + " minimap2_G:"  + minimap2_G + " minimap2_n:" + minimap2_n + " minimap2_m:" + minimap2_m)
    // println("minimap2_p:" + minimap2_p + " minimap2_N:" + minimap2_N + " minimap2_A:"  + minimap2_A + " minimap2_B:" + minimap2_B )
    // println( genome + " fastq: " + basecalled_fastq)
    //-t $task.cpus \\ 
    """
    minimap2 \
        ${preset} \
        ${stranded} \
        ${kmer} \
        ${window} \
        ${minimap2_f} ${minimap2_g} ${minimap2_G} ${minimap2_n} ${minimap2_m} \
        ${minimap2_p} ${minimap2_N} ${minimap2_A} ${minimap2_B} \
        -L -y \
        -t 40 \
        ${genome} ${basecalled_fastq} \
        > ${sample_id}_aligned.sam
    
    """
}