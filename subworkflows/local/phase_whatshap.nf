/*
 * Phase with WHATSHAP
 */

 include { WHATSHAP_PHASE } from '../../modules/local/whatshap_phase'
 include { WHATSHAP_HAPLOTAG } from '../../modules/local/whatshap_haplotag'

workflow PHASE_WHATSHAP {
    take:
    ch_input // channel: [val]
    ch_aligned_sorted_bam
    ch_clair3_merge_vcf

    main:
    WHATSHAP_PHASE( ch_input,ch_aligned_sorted_bam,ch_clair3_merge_vcf )
    ch_phased_vcf_gz = WHATSHAP_PHASE.out.phased_vcf_gz
    WHATSHAP_HAPLOTAG( ch_input,ch_aligned_sorted_bam,ch_phased_vcf_gz )
    ch_whatshap_haplotagged_bam = WHATSHAP_HAPLOTAG.out.whatshap_haplotagged_bam

    emit:
    ch_whatshap_haplotagged_bam
}