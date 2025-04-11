/*
 * Alignment with MINIMAP2
 */

 include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'


workflow ALIGN_MINIMAP2 {
    take:
    ch_input // channel: [val]
    ch_basecalled_fastq

    main:
    
    MINIMAP2_ALIGN( ch_input,ch_basecalled_fastq )
    ch_minimap2_aligned_sam = MINIMAP2_ALIGN.out.minimap2_aligned_sam
    
    emit:
    ch_minimap2_aligned_sam
}