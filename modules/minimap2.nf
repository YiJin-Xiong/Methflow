process Minimap2_align {
    label 'minimap2'

    publishDir "${params.outdir}/align/",
        mode: "copy",
        pattern: "*_aln.sam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(basecalled_fastq)

    output:
    path("${sample_id}_aln.sam"), emit: minimap2_align_sam

    script:
    //DNA -ax map-ont
    //RNA -ax splice -uf 
    def preset    = ( analyte_type == 'DNA' ) ? "-ax map-ont" : "-ax splice"
    def stranded  = ( analyte_type == 'RNA') ? "-uf" : ""
    def kmer      = ( params.minimap2_kmer ) ? "-k"+params.minimap2_kmer : ""
    def window    = ( params.minimap2_window ) ? "-w"+params.minimap2_window : ""
    //println(preset + stranded + kmer + window)
    //-t $task.cpus \\
    """
    minimap2 \
        ${preset} \
        ${stranded} \
        ${kmer} \
        ${window} \
        -L -y \
        -t 40 \
        ${genome} ${basecalled_fastq} \
        > ${sample_id}_aln.sam
    
    """
}