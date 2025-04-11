/*
 * Modcall with DEEPSIGNAL3
 */

 include { DEEPSIGNAL3_CALL_MODS } from '../../modules/local/deepsignal3_call_mods'
 include { DEEPSIGNAL3_CALL_FREQ } from '../../modules/local/deepsignal3_call_freq'


workflow MODCALL_DEEPSIGNAL3 {
    take:
    ch_input // channel: [val]
    ch_basecall_bam

    main:
    ch_input
    .map { row ->
        def path = row[4].toString()
        def prefixDir = '/' + path.tokenize('/')[0..-2].join('/') + '/'
        [row[0], row[1], row[2], row[3], prefixDir, row[5], row[6], row[7]]
    }
    .unique { it[4] }
    .set { ch_input_dir }

    ch_input_dir.view()

    DEEPSIGNAL3_CALL_MODS( ch_input_dir , ch_basecall_bam )
    ch_deepsignal3_call_mods = DEEPSIGNAL3_CALL_MODS.out.deepsignal3_callmods
    DEEPSIGNAL3_CALL_FREQ( ch_input_dir , ch_deepsignal3_call_mods )
    ch_deepsignal3_call_freq = DEEPSIGNAL3_CALL_FREQ.out.deepsignal3_call_freq
    
    emit:
    ch_deepsignal3_call_mods
    ch_deepsignal3_call_freq
}