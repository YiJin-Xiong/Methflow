/*
 * Modcall with ROCKFISH
 */

 include { ROCKFISH_MODCALL } from '../../modules/local/rockfish_modcall'

workflow MODCALL_ROCKFISH {
    take:
    ch_input // channel: [val]
    ch_basecall_bam

    main:

    ROCKFISH_MODCALL( ch_input , ch_basecall_bam )
    ch_rockfish_predictions = ROCKFISH_MODCALL.out.rockfish_predictions

    emit:
    ch_rockfish_predictions
}