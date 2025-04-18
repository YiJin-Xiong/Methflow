process POREMETH2_DMR {
    label 'poremeth2_dmr'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/DMR/PoreMeth2/",
        mode: "copy",
        pattern: "*"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(test_file)
    path(control_file)

    output:
    path("*_poremeth2_log.log")
    path("*_ExpQualityPlot.png"), optional: true
    path("*_DMRTable.txt")
    path("*_DMRAnnotatedTable.txt"), emit: poremeth2
    
    script:
    def omega               = ( params.omega )               ? "--omega "+ params.omega : ""
    def eta                 = ( params.eta )                 ? "--eta "+ params.eta : ""
    def FW                  = ( params.fw )                  ? "--FW "+ params.fw : ""
    def AnnotationType      = ( params.annotation_type )     ? "--AnnotationType "+ params.annotation_type : ""
    def Annotate_Assembly   = ( params.annotate_assembly )   ? "--Annotate_Assembly "+ params.annotate_assembly : ""
    def Statistics_Assembly = ( params.statistics_assembly ) ? "--Statistics_Assembly "+ params.statistics_assembly : ""
    def BetaThr             = ( params.betathr )             ? "--BetaThr "+ params.betathr : ""
    def EntropyThr          = ( params.entropythr )          ? "--EntropyThr "+ params.entropythr : ""
    def PValueThr           = ( params.pvaluethr )           ? "--PValueThr "+ params.pvaluethr : ""
    def AnalysisClass       = ( params.analysis_class )      ? "--AnalysisClass "+ params.analysis_class : ""
    def betacovThr          = ( params.betacovthr )          ? "--betacovThr "+ params.betacovthr : ""

    test_name = test_file.name.split('_extracted')[0]
    control_name = control_file.name.split('_extracted')[0]

    """
    file_prefix="${test_name}_${control_name}_PoreMeth2"

    python ${params.project_dir}/bin/call_poremeth2.py \
        --test ${test_file} \
        --control ${control_file} \
        --test_name ${test_name} \
        --control_name ${control_name} \
        --out_dir ./ \
        --out_prefix \${file_prefix} \
        --overwrite \
        ${omega} ${eta} ${FW} \
        ${AnnotationType} ${Annotate_Assembly} \
        ${Statistics_Assembly} ${BetaThr} \
        ${EntropyThr} ${PValueThr} ${AnalysisClass} ${betacovThr} \
        > \${file_prefix}_poremeth2_log.log

    """
}