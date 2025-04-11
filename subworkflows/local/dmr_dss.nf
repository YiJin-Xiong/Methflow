/*
 * DMR with DSS
 */
 include { MODKIT_PILEUP } from '../../modules/local/modkit_pileup'
 include { DSS_PREPROCESS } from '../../modules/local/dss_preprocess'
 include { DSS_DMR } from '../../modules/local/dss_dmr'

workflow DMR_DSS {
    take:
    ch_input // channel: [val]
    ch_modkit_dmr_trigger
    ch_bed_hp1
    ch_bed_hp2

    main:
    DSS_PREPROCESS( ch_modkit_dmr_trigger,ch_input,ch_bed_hp1,ch_bed_hp2 )
    ch_bed_preprocessed_hp1 = DSS_PREPROCESS.out.bed_preprocessed_hp1
    ch_bed_preprocessed_hp2 = DSS_PREPROCESS.out.bed_preprocessed_hp2
    DSS_DMR( ch_modkit_dmr_trigger,ch_input,ch_bed_preprocessed_hp1,ch_bed_preprocessed_hp2)

    emit:
    ch_input
}