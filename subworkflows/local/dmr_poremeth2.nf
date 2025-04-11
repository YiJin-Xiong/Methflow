/*
 * DMR with POREMETH2
 */
include { POREMETH2_PREPROCESS as PoreMeth2_preprocess_bam1 } from '../../modules/local/poremeth2_preprocess'
include { POREMETH2_PREPROCESS as PoreMeth2_preprocess_bam2} from '../../modules/local/poremeth2_preprocess'

 include { POREMETH2_DMR } from '../../modules/local/poremeth2_dmr'

workflow DMR_POREMETH2 {
    take:
    ch_input // channel: [val]
    ch_haplotagged_bam1
    ch_haplotagged_bam2

    main:
    PoreMeth2_preprocess_bam1( ch_input,ch_haplotagged_bam1 )
    PoreMeth2_preprocess_bam2( ch_input,ch_haplotagged_bam2 )
    ch_extracted_entroyp_tsv1 = PoreMeth2_preprocess_bam1.out.extracted_entroyp_tsv
    ch_extracted_entroyp_tsv2 = PoreMeth2_preprocess_bam2.out.extracted_entroyp_tsv
    POREMETH2_DMR( ch_input,ch_extracted_entroyp_tsv1,ch_extracted_entroyp_tsv2 )

    emit:
    ch_input
}