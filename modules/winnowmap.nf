process Winnowmap_align {
    label 'winnowmap'

    publishDir "${params.outdir}/align/",
        mode: "copy",
        pattern: "*_aligned.sam"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(basecalled_fastq)

    output:
    path("${sample_id}_merylDB"), emit: winnowmap_merylDB
    path("${sample_id}_repetitive.txt"), emit: winnowmap_txt
    path("${sample_id}_aligned.sam"), emit: winnowmap_align_sam

    script:
    def meryl_kmer = ( params.meryl_kmer ) ? "k="+params.meryl_kmer : "k=15"
    def meryl_filter = ( params.meryl_filter ) ? params.meryl_filter : "greater-than"
    def meryl_distinct = ( params.meryl_distinct ) ? "distinct="+params.meryl_distinct : "distinct=0.9998"
    def preset    = ( analyte_type == 'DNA' ) ? "-ax map-ont" : "-ax splice"
    def stranded  = ( analyte_type == 'RNA') ? "-uf" : ""
    def kmer      = ( params.winnowmap_kmer ) ? "-k"+params.mwinnowmap_kmer : ""
    def window    = ( params.winnowmap_window ) ? "-w"+params.winnowmap_window : ""

    """
    meryl count ${meryl_kmer} \
        output ${sample_id}_merylDB ${genome}

    meryl print \
        ${meryl_filter} ${meryl_distinct} ${sample_id}_merylDB \
        > ${sample_id}_repetitive.txt
    
    winnowmap \
        -W ${sample_id}_repetitive.txt \
        ${preset} \
        ${stranded} \
        ${kmer} \
        ${window} \
        -L -y \
        -t 40 \
        ${genome} ${basecalled_fastq} \
        > ${sample_id}_aligned.sam

    """
}