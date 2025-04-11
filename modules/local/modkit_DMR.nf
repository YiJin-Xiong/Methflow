process MODKIT_DMR {
    label 'Modkit_DMR'

    conda "${params.project_dir}/environment.yml"

    publishDir "${params.outdir}/${sample_id}/DMR/",
        mode: "copy",
        pattern: "*"

    input:
    val bed_num
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(bed_hp1)
    path(bed_hp2)

    output:
    path("${sample_id}_1&2_modkit_dmr.bed"), emit: modkit_dmr_bed
    path("${sample_id}_1&2_modkit_dmr.log"), emit: modkit_dmr_log

    script:
    //nextflow run main.nf --input input_sheet.tsv --norm_pileup_bed results/DMR/bed/demo_modkit_pileup.bed.gz --tumor_pileup_bed results/DMR/bed/demo1_modkit_pileup.bed.gz
    def dmr_regions = params.dmr_regions ? "--region " + file(params.dmr_regions) : ''

    """
    tabix -p bed ${bed_hp1}
    tabix -p bed ${bed_hp2}

    modkit dmr pair \
      -a ${bed_hp1} \
      -b ${bed_hp2} \
      -o "${sample_id}_1&2_modkit_dmr.bed" \
      --ref ${genome} \
      ${dmr_regions} \
      --base C \
      --threads 4 \
      -f --header \
      --log-filepath "${sample_id}_1&2_modkit_dmr.log"
      
    """
}