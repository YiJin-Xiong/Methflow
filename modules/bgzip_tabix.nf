process bgzip_tabix {
    label 'bgzip_tabix'

    publishDir "${params.outdir}/DMR/haplotype/",
        mode: "copy",
        pattern: "*"

    input:
    path(bed_file)


    output:
    path("${bed_file}.gz", emit: bed_gz)
    path("${bed_file}.gz.tbi", emit: bed_gz_tbi)

    script:

    """
    bgzip -k ${bed_file}
    tabix -p bed ${bed_file}.gz
    
    """
}
