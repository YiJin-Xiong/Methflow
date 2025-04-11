/*
 * Show help Message
 */

workflow SHOW_HELP {
    log.info"""
    methflow - Nextflow PIPELINE (v$workflow.manifest.version)
    ------------------------------------------------------
    Usage:
    Typical pipeline command:

    nextflow run main.nf --input input_sheet.tsv --protocol DNA

    Input/output options:
      --input                       Path to comma-separated file containing information about the samples in the experiment.
      --protocol                    Input sample type. Valid options: 'DNA', 'RNA',  and 'PacBio HiFi'.
      --outdir                      The output directory where the results will be saved. [default: results]
    
    Tools options:
      --modcall_tool                Choose which tool to use for modcall, 'rockfish', 'deepsiganl3', 'all'. [default: rockfish]
      --align_tool                  Choose which tool to use for alignment, 'minimap2', 'winnowmap'. [default: minimap2]
      --DMR_tool                    Choose which tool to use for DMR, 'modkit', 'DSS', 'poremeth2', 'all'. [default: modkit]
    
    Basecall(dorado) options:
      --modified-bases              A space separated list of modified base codes. Choose from: pseU, 4mC_5mC, inosine_m6A, 5mCG, 5mCG_5hmCG, 5mC_5hmC, 5mC, m5C, 6mA, m6A, m6A_DRACH. [nargs: 1 or more] 
      --max_reads                   Limit the number of reads to be basecalled. [nargs=0..1] [default: 0]
      --read_ids                    A file with a newline-delimited list of reads to basecall. If not provided, all reads will be basecalled. [default: ""]
      --kit_name                    Enable barcoding with the provided kit name.
      --trim                        Specify what to trim. Options are 'none', 'all', 'adapters', and 'primers'. Default behaviour is to trim all detected adapters, primers, or barcodes. [default: ""]
      --min_qscore                  Discard reads with mean Q-score below this threshold. [nargs=0..1] [default: 0]
      --mm2_opts                    Optional minimap2 options string. For multiple arguments surround with double quotes.

    Modcall(rockfish) options:
      --motif                       Specify a sequence motif.
      --idx                         Specify an index for a specific part of the data or a sequence. 
      --mapq_filter                 Define a mapping quality filter for rockfish inference.
      --rockfish_window             Define a window size for rockfish inference.
      --unique_aln                  Only considers unique alignments.

    Modcall(deepsignal3) options:
      --seq_len                     Len of kmer. [default: 21]
      --signal_len                  Signal num of one base. [default: 15]
      --model_type                  Type of model to use, 'both_bilstm', 'seq_bilstm' or 'signal_bilstm'. [default: 'both_bilstm']

    Alignment(minimap2) options:
      --minimap2_kmer               k-mer size for minimap2 (no larger than 28) [default: 15]
      --minimap2_window             minimizer window size for minimap2 [default: 10]
      --minimap2_f                  filter out top FLOAT fraction of repetitive minimizers [default: 0.0002]
      --minimap2_g                  stop chain enlongation if there are no minimizers in INT-bp [default: 5000]
      --minimap2_G                  max intron length (effective with -xsplice; changing -r) [default: 200k]
      --minimap2_n                  minimal number of minimizers on a chain [default: 3]
      --minimap2_m                  minimal chaining score (matching bases minus log gap penalty) [default: 40]
      --minimap2_p                  min secondary-to-primary score ratio [default: 0.8]
      --minimap2_N                  retain at most INT secondary alignments [default: 5]
      --minimap2_A                  matching score [default: 2]
      --minimap2_B                  mismatch penalty (larger value for lower divergence) [default: 4]
    
    Alignment(winnowmap) options:
      --meryl_kmer                  create mers of size K bases (mandatory).
      --meryl_filter                Filtering criterion to use, 'less-than N', 'greater-than N', 
                                    'equal-to N', 'not-equal-to N'. [default: "greater-than distinct=0.9998"]
      --winnowmap_kmer              k-mer size for winnowmap (no larger than 28) [default: 15]
      --winnowmap_window            minimizer window size for winnowmap [default: 10]
      --winnowmap_f                 filter out top FLOAT fraction of repetitive minimizers [default: 0.0]
      --winnowmap_g                 stop chain enlongation if there are no minimizers in INT-bp [default: 5000]
      --winnowmap_G                 max intron length (effective with -xsplice; changing -r) [default: 200k]
      --winnowmap_n                 minimal number of minimizers on a chain [default: 3]
      --winnowmap_m                 minimal chaining score (matching bases minus log gap penalty) [default: 40]
      --winnowmap_p                 min secondary-to-primary score ratio [default: 0.8]
      --winnowmap_A                 matching score [default: 2]
      --winnowmap_B                 mismatch penalty (larger value for lower divergence) [default: 4]
    
    SNV_calling(clair3) options:
      --bed_fn                      Call variants only in the provided bed regions. [FILE]
      --vcf_fn                      Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided. [FILE]
      --ctg_name                    The name of the sequence to be processed.
      --sample_name                 Define the sample name to be shown in the VCF file.
      --qual                        If set, variants with >\$qual will be marked PASS, or LowQual otherwise.

    phasing(whatshap) options:
      --phase_sample                Name of a sample to phase. If not given, all samples in the input VCF are phased. Can be used multiple times.
      --phase_chromosome            Name of chromosome to phase. If not given, all chromosomes in the input VCF are phased. Can be used multiple times.
      --phase_exclude_chromosome    Name of chromosome not to phase.
      --phase_mapq                  Minimum mapping quality [default: 20].
      --only_snvs                   Phase only SNVs.
      --ignore_linked_read          Ignore linkage information stored in BX tags of the reads.
      --tag_supplementary           Also tag supplementary alignments. Supplementary alignments are assigned to the same haplotype as the primary alignment.
                                    [default: only tag primary alignments]
      --skip_missing_contigs        Skip reads that map to a contig that does not exist in the VCF.
      --haplotag_sample             Name of a sample to phase. If not given, all samples in the input VCF are phased. Can be used multiple times.
      --ploidy                      Ploidy [default: 2].

    DMR(modkit) options:
      --pileup_region               Process only the specified region of the BAM when performing pileup.
                                    Format should be <chrom_name>:<start>-<end> or <chrom_name>. Commas are allowed.
      --filter_threshold            Specify the filter threshold globally or per-base.
      --mod_thresholds              Specify a passing threshold to use for a base modification, independent of the threshold for the primary sequence base or the default.
      --ignore                      Ignore specific types of base modifications.
      --modkit_motif                Output pileup counts for only sequence motifs provided.
      --dmr_regions                 BED file of regions over which to compare methylation levels.
    
    DMR(DSS) options:
      --min_len                     Minimum length (in basepairs) required for DMR. [default: 100]
      --min_CG                      Minimum number of CpG sites required for DMR. [default: 3]
      --smoothing_span              The size of smoothing window, in basepairs. [default: 500]
      --smoothing_flag              A flag to indicate whether to appyly smoothing, 'TRUE/FALSE'. [default: TRUE]
      --equal_disp                  Specify whether to assume equal dispersion for all groups. [default: TRUE]
      --pval_cutoff                 A threshold of p-values for calling DMR. [default: 0.001]
      --delta_cutoff                A threshold for defining DMR. [default: 0.05]
      --pct_sig                     In all DMRs, the percentage of CG sites with significant p-values (less than p.threshold) must be greater than this threshold. [default: 0.5]

    DMR(poremeth2) options:
      --omega                       Optional parameter that modulates the relative weight between the experimental and the biological variance. [default: 0.1]
      --eta                         Optional parameter that represents the baseline probability the mean process (m_i) changes its value for the HSLM algorithm. [default: 1e-05]
      --fw                          The minimum number of datapoints for a DMR to be called (DMRs made of a number of CpGs smaller than FW are discarded). [default: 3]
      --annotation_type             Specify whether to annotate DMRs on genic elements only ('Genes') or genic elements 
                                    and regulatory features ('GenesReg'). [default: 'Genes']
      --annotate_assembly           Specify reference version to use for annotation, 'hg19' or 'hg38'. [default: 'hg19']
      --statistics_assembly         Specify reference version to use for statistics, 'hg19' or 'hg38'. [default: 'hg19']
      --betathr                     Delta beta threshold applied for DMRs' classification. [default: 0.2]
      --entropythr                  Delta S threshold applied for DMRs' classification. [default: 0.2]
      --pvaluethr                   The p.value threshold to consider a DMR. [default: 0.05]
      --analysis_class              Define the summary to return ('All', 'Beta', 'Entropy') [default: 'All']
      --betacovthr                  A threshold to filter input data based on beta_cov. [default: 0.1]


    Running environment options:
      --docker_name     Docker name used for pipeline, default is '/:latest'
      --singularity_name    Singularity name used for pipeline, default is 'docker:///:latest'
      --singularity_cache   Singularity cache dir, default is 'local_singularity_cache'
      --conda_name      Conda name used for pipeline, default is 'nanome'
      --conda_base_dir  Conda base directory, default is '/opt/conda'
      --gpu                     Use gpu.

    -profile options:
      Use this parameter to choose a predefined configuration profile. Profiles can give configuration presets for different compute environments.

      test      A test demo config
      docker    A generic configuration profile to be used with Docker, pulls software from Docker Hub: /:latest
      singularity   A generic configuration profile to be used with Singularity, pulls software from: docker:///:latest
      conda     Please only use conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity. Check our GitHub for how to install local conda environment

    """.stripIndent()
}