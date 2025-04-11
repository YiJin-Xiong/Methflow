/*
 * Pileup with MODKIT
 */
 include { MODKIT_PILEUP } from '../../modules/local/modkit_pileup'
 include { BGZIP_TABIX } from '../../modules/local/bgzip_tabix'

workflow PILEUP_MODKIT {
    take:
    ch_input // channel: [val]
    ch_whatshap_haplotagged_bam

    main:
    ch_sample_id = ch_input.map { it[1] }.first()

    MODKIT_PILEUP( ch_input,ch_whatshap_haplotagged_bam )

    MODKIT_PILEUP.out.bed_files.flatten()
        .map { bed_file -> tuple(ch_sample_id.get(), bed_file) } 
        .set { ch_modkit_output_beds }

    //ch_modkit_output_beds.view()
    //ch_sample_id.view()
    
    BGZIP_TABIX( ch_modkit_output_beds )
    ch_all_bedgz = BGZIP_TABIX.out.bed_gz.collect()

    ch_bed_hp1 = ch_all_bedgz.flatten().filter{ it.name.contains('_1.bed.gz') }
    ch_bed_hp2 = ch_all_bedgz.flatten().filter{ it.name.contains('_2.bed.gz') }
    ch_modkit_dmr_trigger = MODKIT_PILEUP.out.bed_num.filter { it.toInteger() >= 2 }

    ch_modkit_dmr_trigger.view()
    
    emit:
    ch_modkit_dmr_trigger
    ch_bed_hp1
    ch_bed_hp2
}