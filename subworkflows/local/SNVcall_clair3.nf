/*
 * SNVcall with CLAIR3
 */

 include { CLAIR3_SNVCALL } from '../../modules/local/clair3_SNVcall'

workflow SNVCALL_CLAIR3 {
    take:
    ch_input // channel: [val]
    ch_aligned_sorted_bam

    main:
    CLAIR3_SNVCALL( ch_input,ch_aligned_sorted_bam )
    ch_clair3_merge_vcf = CLAIR3_SNVCALL.out.clair3_merge_vcf

    emit:
    ch_clair3_merge_vcf
}