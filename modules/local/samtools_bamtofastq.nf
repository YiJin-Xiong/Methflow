process SAMTOOLS_BAMTOFASTQ {
    label 'samtools_bam2fastq'

    conda     (params.enable_conda ? "${params.project_dir}/environment.yml" : null)
    container (params.use_docker ? "${params.docker_name}" : "${params.singularity_name}")

    publishDir "${params.outdir}/${sample_id}/alignment/",
        mode: "copy",
        pattern: "${sample_id}_basecall.fastq"

    input:
    tuple val(group_id), val(sample_id), val(analyte_type), val(pore_type), path(pod5_dir), path(genome), val(species), val(sample_rate)
    path(basecall_bam_dir)

    output:
    path("${sample_id}_basecall.fastq"), emit: basecalled_fastq

    script:
    """
    samtools fastq -@20 -T qs,du,ns,ts,mx,ch,st,rn,fn,sm,sd,sv,dx,RG,mv,MN,MM,ML \
        ${basecall_bam_dir} \
        > ${sample_id}_basecall.fastq
    
    """
}