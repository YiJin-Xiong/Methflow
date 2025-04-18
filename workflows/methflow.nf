/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist

def checkOptions(){
    if ( ! params.input ){
        exit 1, 'No input sheet.'
    }
    
    if ( params.protocol != 'DNA' && params.protocol != 'RNA' ) {
        exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'RNA' "
    }

    if ( params.modcall_tool ){
        if (params.modcall_tool != 'rockfish' && params.modcall_tool != 'deepsignal3' && params.modcall_tool != 'all' ) {
            exit 1, "Invalid modcall_tool option: ${params.modcall_tool}. Valid options: 'rockfish', 'deepsignal3', 'all' "
        }
    }

    if ( params.align_tool ){
        if (params.align_tool != 'minimap2' && params.align_tool != 'winnowmap' ) {
            exit 1, "Invalid align_tool option: ${params.align_tool}. Valid options: 'minimap2', 'winnowmap' "
        }
    }

    if ( params.DMR_tool ){
        if (params.DMR_tool != 'modkit' && params.DMR_tool != 'DSS' && params.DMR_tool != 'poremeth2' && params.DMR_tool != 'all' ) {
            exit 1, "Invalid DMR_tool option: ${params.DMR_tool}. Valid options: 'modkit', 'DSS', 'poremeth2', 'all' "
        }
    }
    
}


////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////





/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
  */

include { SHOW_HELP                      } from '../subworkflows/local/show_help'
include { CHECK_INPUTSHEET               } from '../subworkflows/local/check_inputsheet'
include { BASECALL_DORADO                } from '../subworkflows/local/basecall_dorado'
include { MODCALL_ROCKFISH               } from '../subworkflows/local/modcall_rockfish'
include { MODCALL_DEEPSIGNAL3            } from '../subworkflows/local/modcall_deepsignal3'
include { BAMTOFASTQ_SAMTOOLS            } from '../subworkflows/local/bamtofastq_samtools'
include { ALIGN_MINIMAP2                 } from '../subworkflows/local/align_minimap2'
include { ALIGN_WINNOWMAP                } from '../subworkflows/local/align_winnowmap'
include { SAMTOBAM_SORT_SAMTOOLS         } from '../subworkflows/local/samtobam_sort_samtools'
include { SNVCALL_CLAIR3                 } from '../subworkflows/local/SNVcall_clair3'
include { PHASE_WHATSHAP                 } from '../subworkflows/local/phase_whatshap'
include { PILEUP_MODKIT                  } from '../subworkflows/local/pileup_modkit'
include { DMR_DSS                        } from '../subworkflows/local/dmr_dss'
include { DMR_POREMETH2                  } from '../subworkflows/local/dmr_poremeth2'

include { MODKIT_DMR } from '../modules/local/modkit_DMR'

workflow METHFLOW {

    if ( params.help ){
        SHOW_HELP()
        exit 0
    }

    // Check mandatory options (missing input or protocol will exit the run.)
    // Check tools options (invalid tool options will exit the run.)
    checkOptions()

    // Check input_sheet (invalid format will exit the run.)
    if ( params.input.endsWith(".txt") || params.input.endsWith(".tsv") || params.input.endsWith(".csv")) {
        ch_input_path = Channel.fromPath( params.input, checkIfExists: true )
    } else {
        exit 1, "Invalid input format. Valid format: '.txt', '.tsv', 'csv' "
    }

    // =============================================================================
    /*
     * SUBWORKFLOW: Read in inputsheet, validate and stage input files
     */
    CHECK_INPUTSHEET( ch_input_path ).set{ ch_input }
    ch_input.view()

    // =============================================================================
    /*
     * SUBWORKFLOW: Basecall with dorado
     */
    BASECALL_DORADO( ch_input )
    ch_basecall_bam = BASECALL_DORADO.out.ch_basecall_bam

    // =============================================================================
    /*
     * SUBWORKFLOW: Convert Bam files to Fastq files with samtools
     */
    BAMTOFASTQ_SAMTOOLS( ch_input,ch_basecall_bam )
    ch_basecalled_fastq = BAMTOFASTQ_SAMTOOLS.out.ch_basecalled_fastq

    if ( params.align_tool == 'winnowmap' ){
        /*
         * SUBWORKFLOW: Alignment with winnowmap
         */
        ALIGN_WINNOWMAP( ch_input,ch_basecalled_fastq )
        ch_winnowmap_aligned_sam = ALIGN_WINNOWMAP.out.ch_winnowmap_aligned_sam
        SAMTOBAM_SORT_SAMTOOLS( ch_input,ch_winnowmap_aligned_sam )
        ch_aligned_sorted_bam = SAMTOBAM_SORT_SAMTOOLS.out.ch_aligned_sorted_bam

    } else {
        /*
         * SUBWORKFLOW: Alignment with minimap2 [default]
         */
        ALIGN_MINIMAP2( ch_input,ch_basecalled_fastq )
        ch_minimap2_aligned_sam = ALIGN_MINIMAP2.out.ch_minimap2_aligned_sam
        SAMTOBAM_SORT_SAMTOOLS( ch_input,ch_minimap2_aligned_sam )
        ch_aligned_sorted_bam = SAMTOBAM_SORT_SAMTOOLS.out.ch_aligned_sorted_bam
        
    }

    // =============================================================================
    if ( params.modcall_tool == 'deepsignal3' ){
        /*
         * SUBWORKFLOW: Modcall with deepsignal3
         */
        MODCALL_DEEPSIGNAL3( ch_input,ch_aligned_sorted_bam )
        
    } else if ( params.modcall_tool == 'all' ){
        /*
         * SUBWORKFLOW: Modcall with rockfish and deepsignal3
         */
        MODCALL_DEEPSIGNAL3( ch_input,ch_aligned_sorted_bam )
        MODCALL_ROCKFISH( ch_input,ch_aligned_sorted_bam )
        
    } else {
        /*
         * SUBWORKFLOW: Modcall with rockfish [default]
         */
        MODCALL_ROCKFISH( ch_input,ch_aligned_sorted_bam )
    }


    /*
     * SUBWORKFLOW: SNV_calling with clair3
     */
    SNVCALL_CLAIR3( ch_input,ch_aligned_sorted_bam )
    ch_clair3_merge_vcf = SNVCALL_CLAIR3.out.ch_clair3_merge_vcf


    /*
     * SUBWORKFLOW: phase with whatshap
     */
    PHASE_WHATSHAP( ch_input,ch_aligned_sorted_bam,ch_clair3_merge_vcf)
    ch_whatshap_haplotagged_bam = PHASE_WHATSHAP.out.ch_whatshap_haplotagged_bam



    if ( params.DMR_tool == 'DSS' ){
        /*
         * SUBWORKFLOW: DMR with DSS
         */
        PILEUP_MODKIT( ch_input,ch_whatshap_haplotagged_bam )
        ch_modkit_dmr_trigger = PILEUP_MODKIT.out.ch_modkit_dmr_trigger
        ch_bed_hp1 = PILEUP_MODKIT.out.ch_bed_hp1
        ch_bed_hp2 = PILEUP_MODKIT.out.ch_bed_hp2
        DMR_DSS( ch_input,ch_modkit_dmr_trigger,ch_bed_hp1,ch_bed_hp2 )
        
    } else if ( params.DMR_tool == 'poremeth2' ){
        /*
         * SUBWORKFLOW: DMR with poremeth2
         */

        // def bam1 = file(params.bam1)
        // def bam2 = file(params.bam2)
        // DMR_POREMETH2( ch_input, bam1, bam2)
        
        
    } else if ( params.DMR_tool == 'all' ){
        /*
         * SUBWORKFLOW: DMR with modkit , DSS and poremeth2
         */
        
        
    } else {
        /*
         * SUBWORKFLOW: DMR with modkit [default]
         */
        PILEUP_MODKIT( ch_input,ch_whatshap_haplotagged_bam )
        ch_modkit_dmr_trigger = PILEUP_MODKIT.out.ch_modkit_dmr_trigger
        ch_bed_hp1 = PILEUP_MODKIT.out.ch_bed_hp1
        ch_bed_hp2 = PILEUP_MODKIT.out.ch_bed_hp2
        MODKIT_DMR( ch_modkit_dmr_trigger,ch_input,ch_bed_hp1,ch_bed_hp2 )
    }



    
    
}