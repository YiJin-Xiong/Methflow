//import modules
include { Dorado_basecall } from './modules/dorado'
include { Dorado_align } from './modules/dorado'
include { Rockfish_inference } from './modules/rockfish'

workflow {
    //nextflow run main.nf --input input_sheet.tsv
    //--mm2_opts "-k 15 -w 10"
    //--rockfish results/bam/demo_basecall.bam

    // generate input files, and send into Channels for pipelines
    if ( params.input.endsWith(".txt") || params.input.endsWith(".tsv") ) {
        Channel.fromPath( params.input, checkIfExists: true )
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{ it -> [it.Group_ID, it.Sample_ID, it.Analyte_Type, it.Pore_Type, file(it.Path), file(it.Genome) ] }
            .filter{ it[4].exists() && it[5].exists() }
            .set{ input_dorado }
    } else {
        exit 1, "--input must be in tsv format, see ./demo/input_sheet.tsv for more information!"
    }

    //查看通道里的内容 
    input_dorado.view()

    //传入Group_ID, Sample_ID, Analyte_Type, Pore_Type, Path
    Dorado_basecall( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4] ] } )

    //传入Sample_ID, Genome, 和Dorado_basecall的输出文件basecall_bam
    Dorado_align( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[5] ] } , Dorado_basecall.out.basecall)

    //Rockfish
    //Rockfish_inference( input_dorado.map {it -> [ it[0], it[1], it[2], it[3], it[4] ] } , Dorado_basecall.out.basecall)

}
