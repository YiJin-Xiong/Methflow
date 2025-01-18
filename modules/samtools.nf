process Samtools_bam2fastq {
    label 'samtools_bam2fastq'

    publishDir "${params.outdir}/align/",
        mode: "copy",
        pattern: "${sample_id}_basecall.fastq"

    input:
    tuple val(group_id), val(sample_id)
    path(basecall_bam_dir)

    output:
    path("${sample_id}_basecall.fastq"), emit: basecalled_fastq

    script:
    """
    samtools fastq -@20 -T qs,du,ns,ts,mx,ch,st,rn,fn,sm,sd,sv,dx,RG,mv,MN,MM,ML \
        ${basecall_bam_dir} \
        > ${sample_id}_basecall.fastq
    
    """
}

process Samtools_sam2bam {
    label 'samtools_sam2bam'

    publishDir "${params.outdir}/align/",
        mode: "copy",
        pattern: "${sample_id}_aln.bam"

    input:
    tuple val(group_id), val(sample_id)
    path(minimap2_align_sam)

    output:
    path("${sample_id}_aln.bam"), emit: aligned_bam

    script:
    """
    samtools view -@20 -bh \
        ${minimap2_align_sam} \
        > ${sample_id}_aln.bam
    
    """
}

process Samtools_sort {
    label 'samtools_sort'

    publishDir "${params.outdir}/align/",
        mode: "copy",
        pattern: "${sample_id}_basecall.fastq"

    input:
    tuple val(group_id), val(sample_id)
    path(aligned_bam)

    output:
    path("${sample_id}_aligned_sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools sort -m 4G \
        -@ 20 ${aligned_bam} \
        -o ${sample_id}_aligned_sorted.bam \
    
    """
}