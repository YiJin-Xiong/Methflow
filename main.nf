//import modules
include { Dorado_basecall } from './modules/dorado'
include { Dorado_align } from './modules/dorado'
include { Rockfish_inference } from './modules/rockfish'
include { Deepsignal3_call_mods } from './modules/deepsignal3'
include { Deepsignal3_call_freq } from './modules/deepsignal3'
include { Minimap2_align } from './modules/minimap2'
include { Winnowmap_align } from './modules/winnowmap'
include { Samtools_bam2fastq } from './modules/samtools'
include { Samtools_sam2bam } from './modules/samtools'
include { Samtools_sort } from './modules/samtools'
include { Clair3 } from './modules/clair3'
include { WhatsHap } from './modules/whatshap.nf'
include { Modkit_pileup } from './modules/modkit.nf'
include { Modkit_DMR } from './modules/modkit.nf'
include { bgzip_tabix } from './modules/bgzip_tabix.nf'
include { DSS_preprocess } from './modules/dss.nf'
include { DSS } from './modules/dss.nf'
include { PoreMeth2_preprocess as PoreMeth2_preprocess_bam1 } from './modules/poremeth2.nf'
include { PoreMeth2_preprocess as PoreMeth2_preprocess_bam2} from './modules/poremeth2.nf'
include { PoreMeth2 } from './modules/poremeth2.nf'



workflow {
    //nextflow run main.nf --input input_sheet.tsv 
    //--mm2_opts "-k 15 -w 10"
    //--rockfish_window 12
    //--align results/bam/demo_aligned.bam
    //--input_path pod5原始文件，要求必须是文件夹/
    //--minimap2_kmer 14 --minimap2_window 10
    //--clair3_input results/align/demo_aligned_sorted.bam
    //--clair3_output_dir results/vcf/clair3_output_dir
    //nextflow run main.nf --input input_sheet.tsv --whatshap_bam results/align/demo_aligned_sorted.bam --whatshap_vcf results/vcf/demo_clair3_merge.vcf.gz --ignored_linked_read --tag_supplementary
    //nextflow run main.nf --input input_sheet.tsv --modkit_input results/whatshap/demo_whatshap_haplotagged.bam
    //nextflow run main.nf --input input_sheet.tsv --modkit_input results/temp/subset_1/subset_1_whatshap_haplotagged.bam
    ////nextflow run main.nf --input input_sheet.tsv --norm_pileup_bed results/DMR/bed/demo_modkit_pileup.bed.gz --tumor_pileup_bed results/DMR/bed/demo1_modkit_pileup.bed.gz

    //nextflow run main.nf --input input_sheet.tsv --bam1 results/temp/subset_1/subset_1_whatshap_haplotagged.bam --bam2 results/temp/subset_2/subset_2_whatshap_haplotagged.bam

    // generate input files, and send into Channels for pipelines
    if ( params.input.endsWith(".txt") || params.input.endsWith(".tsv") ) {
        Channel.fromPath( params.input, checkIfExists: true )
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{ it -> [it.Group_ID, it.Sample_ID, it.Analyte_Type, it.Pore_Type, file(it.Path), file(it.Genome), it.Species, it.Sample_Rate ] }
            .filter{ it[4].exists() && it[5].exists() }
            .set{ input_dorado }
    } else {
        exit 1, "--input must be in tsv format, see ./demo/input_sheet.tsv for more information!"
    }

    input_dorado.map { row ->
        def speciesValue = row[6]
        if (speciesValue == null || speciesValue == "") {
            exit 1, "Please speficy \"Species\" in the input_sheet.tsv!"
        }
        return row
    }

    //查看通道里的内容 
    input_dorado.view()

    //传入Group_ID, Sample_ID, Analyte_Type, Pore_Type, Path
    //Dorado_basecall( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4] ] } )

    //传入Sample_ID, Genome, 和Dorado_basecall的输出文件basecall_bam
    //Dorado_align( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[5] ] } , Dorado_basecall.out.basecall)

    //Rockfish
    //Rockfish_inference( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4] ] } , Dorado_basecall.out.basecall)

    //DeepSignal3
    //def alignm = file(params.align)
    //Deepsignal3_call_mods( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[6], it[7] ] }, alignm)
    //Deepsignal3_call_mods( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[6], it[7] ] }, Dorado_align.out.alignment)

    //Deepsignal3_call_freq( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[6], it[7] ] }, Deepsignal3_call_mods.out.deepsignal3_callmods)


    //Minimap2_workflow 1. bam2fastq 2. align 3. sam2bam 4. sort
    //1. bam2fastq
    //Samtools_bam2fastq( input_dorado.map {it -> [ it[0], it[1] ] } , Dorado_basecall.out.basecall )
    //2. align
    //Minimap2_align( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , Samtools_bam2fastq.out.basecalled_fastq )
    //3. sam2bam
    //Samtools_sam2bam( input_dorado.map {it -> [ it[0], it[1] ] } , Minimap2_align.out.minimap2_align_sam)
    //4. sort
    //Samtools_sort( input_dorado.map {it -> [ it[0], it[1] ] } , Samtools_sam2bam.out.aligned_bam)


    //Winnowmap_workflow 1. bam2fastq 2. align 3. sam2bam 4. sort
    //1. bam2fastq
    //Samtools_bam2fastq( input_dorado.map {it -> [ it[0], it[1] ] } , Dorado_basecall.out.basecall )
    //2. align
    //Winnowmap_align( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , Samtools_bam2fastq.out.basecalled_fastq )
    //3. sam2bam
    //Samtools_sam2bam( input_dorado.map {it -> [ it[0], it[1] ] } , Winnowmap_align.out.winnowmap_align_sam)
    //4. sort
    //Samtools_sort( input_dorado.map {it -> [ it[0], it[1] ] } , Samtools_sam2bam.out.aligned_bam)


    //Clair3
    //def clair3input = file(params.clair3_input)
    //Clair3( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , clair3input )
    //////////Clair3( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , Samtools_sort.out.aligned_sorted_bam )

    //whatshap
    //  def whatshapbam = file(params.whatshap_bam)
    //  def whatshapvcf = file(params.whatshap_vcf)
    //  WhatsHap( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , whatshapbam , whatshapvcf)
    // WhatsHap( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , clair3input , Clair3.out.clair3_merge_vcf)
    ///////////////WhatsHap( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , Samtools_sort.out.aligned_sorted_bam , Clair3.out.clair3_merge_vcf)

    //modkit
    /////////def modkitinput = file(params.modkit_input)
    /////////Modkit_pileup( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } ,  modkitinput )
    /////Modkit_pileup( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , Samtools_sort.out.aligned_sorted_bam )

    /////////modkit_output_beds = Modkit_pileup.out.bed_files.flatten()
    //modkit_output_beds.view()
    /////////bgzip_tabix(modkit_output_beds )

    /////////all_bedgz = bgzip_tabix.out.bed_gz.collect()
    /////////bed_hp1 = all_bedgz.flatten().filter{ it.name.contains('_1.bed.gz') }
    /////////bed_hp2 = all_bedgz.flatten().filter{ it.name.contains('_2.bed.gz') }
    //bed_hp1.view()
    //bed_hp2.view()

    //Modkit_pileup.out.bed_num.view()
    /////////modkit_dmr_trigger = Modkit_pileup.out.bed_num.filter { it.toInteger() >= 2 }
    /////////Modkit_DMR( modkit_dmr_trigger, input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , bed_hp1, bed_hp2 )

    /////////DSS_preprocess( modkit_dmr_trigger, input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , bed_hp1, bed_hp2 )
    /////////DSS( modkit_dmr_trigger, input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] } , DSS_preprocess.out.bed_preprocessed_hp1, DSS_preprocess.out.bed_preprocessed_hp2 )

    def bam1 = file(params.bam1)
    def bam2 = file(params.bam2)
    PoreMeth2_preprocess_bam1(input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] }, bam1)
    PoreMeth2_preprocess_bam2(input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] }, bam2)
    PoreMeth2(PoreMeth2_preprocess_bam1.out.extracted_entroyp_tsv , PoreMeth2_preprocess_bam2.out.extracted_entroyp_tsv)


}
