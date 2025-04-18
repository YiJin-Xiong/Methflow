process WINNOWMAP_ALIGN {
    label 'winnowmap'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/alignment/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(basecalled_fastq)

    output:
    path("${sample_id}_merylDB"), emit: winnowmap_merylDB
    path("${sample_id}_repetitive.txt"), emit: winnowmap_txt
    path("${sample_id}_aligned.sam"), emit: winnowmap_aligned_sam

    script:
    def meryl_kmer    = ( params.meryl_kmer )       ? "k=" + params.meryl_kmer : "k=15"
    def meryl_filter  = ( params.meryl_filter )     ? params.meryl_filter : "greater-than distinct=0.9998"
    def preset        = ( analyte_type == 'DNA' )   ? "-ax map-ont" : "-ax splice"
    def stranded      = ( analyte_type == 'RNA')    ? "-uf" : ""
    def kmer          = ( params.winnowmap_kmer )   ? "-k " + params.winnowmap_kmer : ""
    def window        = ( params.winnowmap_window ) ? "-w " + params.winnowmap_window : ""
    def winnowmap_f   = ( params.winnowmap_f )      ? "-f " + params.winnowmap_f : ""
    def winnowmap_g   = ( params.winnowmap_g )      ? "-g " + params.winnowmap_g : ""
    def winnowmap_G   = ( params.winnowmap_G )      ? "-G " + params.winnowmap_G : ""
    def winnowmap_n   = ( params.winnowmap_n )      ? "-n " + params.winnowmap_n : ""
    def winnowmap_m   = ( params.winnowmap_m )      ? "-m " + params.winnowmap_m : ""
    def winnowmap_p   = ( params.winnowmap_p )      ? "-p " + params.winnowmap_p : ""
    def winnowmap_A   = ( params.winnowmap_A )      ? "-A " + params.winnowmap_A : ""
    def winnowmap_B   = ( params.winnowmap_B )      ? "-B " + params.winnowmap_B : ""

    """
    meryl count ${meryl_kmer} \
        output ${sample_id}_merylDB ${genome}

    meryl print \
        ${meryl_filter} ${sample_id}_merylDB \
        > ${sample_id}_repetitive.txt
    
    winnowmap \
        -W ${sample_id}_repetitive.txt \
        ${preset} \
        ${stranded} \
        ${kmer} \
        ${window} \
        ${genome} ${basecalled_fastq} \
        ${winnowmap_f} ${winnowmap_g} ${winnowmap_G} ${winnowmap_n} \
        ${winnowmap_m} ${winnowmap_p} ${winnowmap_A} ${winnowmap_B} \
        -L -y \
        -t 40 \
        > ${sample_id}_aligned.sam

    """
}