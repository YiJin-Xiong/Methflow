/*
 * BamtoFastq with SAMTOOLS
 */

 include { SAMTOOLS_BAMTOFASTQ } from '../../modules/local/samtools_bamtofastq'


workflow BAMTOFASTQ_SAMTOOLS {
    take:
    ch_input // channel: [val]
    ch_basecall_bam

    main:
    
    SAMTOOLS_BAMTOFASTQ( ch_input,ch_basecall_bam )
    ch_basecalled_fastq = SAMTOOLS_BAMTOFASTQ.out.basecalled_fastq
    
    emit:
    ch_basecalled_fastq

}