/*
 * samtobam and sort with SAMTOOLS
 */

 include { SAMTOOLS_SAMTOBAM } from '../../modules/local/samtools_samtobam'
 include { SAMTOOLS_SORT } from '../../modules/local/samtools_sort'


workflow SAMTOBAM_SORT_SAMTOOLS {
    take:
    ch_input // channel: [val]
    ch_aligned_sam

    main:
    
    SAMTOOLS_SAMTOBAM( ch_input,ch_aligned_sam )
    ch_aligned_bam = SAMTOOLS_SAMTOBAM.out.aligned_bam
    SAMTOOLS_SORT( ch_input,ch_aligned_bam )
    ch_aligned_sorted_bam = SAMTOOLS_SORT.out.aligned_sorted_bam
    
    emit:
    ch_aligned_bam
    ch_aligned_sorted_bam
}