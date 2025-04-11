/*
 * Basecall with DORADO
 */

 include { DORADO_BASECALL } from '../../modules/local/dorado_basecall'

workflow BASECALL_DORADO {
    take:
    ch_input // channel: [val]

    main:

    DORADO_BASECALL( ch_input )
    ch_basecall_bam = DORADO_BASECALL.out.basecall_bam
    
    emit:
    ch_basecall_bam
}