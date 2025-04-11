/*
 * Alignment with WINNOWMAP
 */

 include { WINNOWMAP_ALIGN } from '../../modules/local/winnowmap_align'


workflow ALIGN_WINNOWMAP {
    take:
    ch_input // channel: [val]
    ch_basecalled_fastq

    main:
    
    WINNOWMAP_ALIGN( ch_input,ch_basecalled_fastq )
    ch_winnowmap_aligned_sam = WINNOWMAP_ALIGN.out.winnowmap_aligned_sam
    
    emit:
    ch_winnowmap_aligned_sam
}