process DSS_preprocess {
    label 'DSS_preprocess'

    publishDir "${params.outdir}/DSS/",
        mode: "copy",
        pattern: "*"

    input:
    val bed_num
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(bed_hp1)
    path(bed_hp2)

    output:
    path("${sample_id}_haplotype_preprocessed_1.bed"), emit: bed_preprocessed_hp1
    path("${sample_id}_haplotype_preprocessed_2.bed"), emit: bed_preprocessed_hp2
    
    script:
    def coverage = ( params.DSS_coverage ) ? "--DSS_coverage "+ params.DSS_coverage : 1

    """
    gunzip -c ${bed_hp1} \
      | awk 'BEGIN{FS="\t"; OFS="\t"} {print \$1, \$3, \$10, \$12}' \
      | awk -F '\t' '\$3 >= ${coverage}' \
      | awk 'BEGIN {print "chro\tpos\tN\tX"} {print}' \
      > ${sample_id}_haplotype_preprocessed_1.bed
    
    gunzip -c ${bed_hp2} \
      | awk 'BEGIN{FS="\t"; OFS="\t"} {print \$1, \$3, \$10, \$12}' \
      | awk -F '\t' '\$3 >= ${coverage}' \
      | awk 'BEGIN {print "chro\tpos\tN\tX"} {print}' \
      > ${sample_id}_haplotype_preprocessed_2.bed

    """
}

process DSS {
    label 'DSS'

    publishDir "${params.outdir}/DSS/",
        mode: "copy",
        pattern: "*"

    input:
    val bed_num
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(bed_preprocessed_hp1)
    path(bed_preprocessed_hp2)

    output:
    path("*_logdmr.log")
    path("*_DMLtest.txt"), optional: true
    path("*_callDML.txt"), optional: true
    path("*_callDMR.txt"), optional: true, emit: dss_dmr
    
    script:
    def min_len = ( params.min_len ) ? "--minlen "+ params.min_len : ""
    def min_CG = ( params.min_CG ) ? "--minCG "+ params.min_CG : ""
    def smoothing_span = ( params.smoothing_span ) ? "--smoothing_span "+ params.smoothing_span : ""
    def smoothing_flag = ( params.smoothing_flag ) ? "--smoothing_flag "+ params.smoothing_flag : ""
    def equal_disp = ( params.equal_disp ) ? "--equal_disp "+ params.equal_disp : ""
    def pval_cutoff = ( params.pval_cutoff ) ? "--pval_cutoff "+ params.pval_cutoff : ""
    def delta_cutoff = ( params.delta_cutoff ) ? "--delta_cutoff "+ params.delta_cutoff : ""
    def pct_sig = ( params.pct_sig ) ? "--pct_sig "+ params.pct_sig : ""

    """
    file_prefix="${sample_id}_DSS"

    python ${projectDir}/bin/call_dss.py \
        --case ${bed_preprocessed_hp1} \
        --control ${bed_preprocessed_hp2} \
        --out_dir ./ \
        --out_prefix \${file_prefix} \
        --overwrite \
        ${min_len} ${min_CG} \
        ${smoothing_span} ${smoothing_flag} \
        ${equal_disp} ${pval_cutoff} \
        ${delta_cutoff} ${pct_sig} \
        > \${file_prefix}_logdmr.log 2>&1
    
    dmr_txt=\${file_prefix}_callDMR.txt

    if [[ ! -f \${dmr_txt} ]] ; then
        echo "### ERROR: no expected output file"
        exit -1
    fi

    """
}