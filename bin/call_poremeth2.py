#! /usr/bin/env python3

import os
import glob
import argparse
import subprocess
import time

def openfile(file):
    '''
    Opens a file
    '''
    if file.endswith('.tsv'):
        opened_file = open(file,'r')
    else:
        raise ValueError("The file is not in TSV format.")
    return opened_file


def main_poremeth2(args):
    """
    This is the DMA module which does differential methylation analysis
    using PoreMeth2 package to find differentially methylated regions.
    """
    t_start = time.time()
    if os.path.isdir(os.path.abspath(args.test)):
        tests = []
        for (dirpath, filenames) in os.walk(os.path.abspath(args.test)):
            for filename in filenames:
                tests.append(dirpath + '/' + filename)
    else:
        tests = [os.path.abspath(args.test)]

    if os.path.isdir(os.path.abspath(args.control)):
        controls = []
        for (dirpath, filenames) in os.walk(os.path.abspath(args.control)):
            for filename in filenames:
                controls.append(dirpath + '/' + filename)
    else:
        controls = [os.path.abspath(args.control)]
    

    testname = args.test_name
    controlname = args.control_name
    out_dir = os.path.abspath(args.out_dir)
    out_prefix = out_dir + '/' + (args.out_prefix)
    Rscript = args.Rscript  # os.path.abspath(args.Rscript)
    script = os.path.abspath(args.script_file)

    omega = args.omega
    eta = args.eta
    FW = args.FW
    AnnotationType = args.AnnotationType
    Annotate_Assembly = args.Annotate_Assembly
    Statistics_Assembly = args.Statistics_Assembly
    BetaThr = args.BetaThr
    EntropyThr = args.EntropyThr
    PValueThr = args.PValueThr
    AnalysisClass = args.AnalysisClass
    betacovThr = args.betacovThr


    # check if outputs exist
    check_outs = [x for x in glob.glob("{}*DM*.txt".format(out_prefix))]
    if check_outs and not args.overwrite:
        raise FileExistsError("The selected output files {} already "
                                "exist. Select --overwrite option if you "
                                "want to overwrite them or use a different "
                                "prefix".format(check_outs))

    subprocess.call(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(Rscript,
                                                           script,
                                                           " ".join(tests),
                                                           " ".join(controls),
                                                           testname,
                                                           controlname,
                                                           out_prefix,
                                                           omega,
                                                           eta,
                                                           FW,
                                                           AnnotationType,
                                                           Annotate_Assembly,
                                                           Statistics_Assembly,
                                                           BetaThr,
                                                           EntropyThr,
                                                           PValueThr,
                                                           AnalysisClass,
                                                           betacovThr),
        shell=True)
    t_end1 = time.time()
    print("===call PoreMeth2 costs {:.1f} seconds".format(t_end1 - t_start))


def main():
    parser = argparse.ArgumentParser()
    poremeth2_input = parser.add_argument_group("required arguments")
    poremeth2_input.add_argument("--test", "-te",
                           action="store",
                           type=str,
                           required=True,
                           help=("Tables output by ModkitResorter.sh on the Test samples. "
                                 "If multiple files, files must be in the same directory "
                                 "and give the path to the directory."))
    poremeth2_input.add_argument("--control", "-co",
                           action="store",
                           type=str,
                           required=True,
                           help=("Tables output by ModkitResorter.sh on the Control samples. "
                                 "If multiple files, files must be in the same directory "
                                 "and give the path to the directory."))
    poremeth2_input.add_argument("--test_name",
                           action="store",
                           type=str,
                           required=True,
                           help=("Sample_id of test file."))
    poremeth2_input.add_argument("--control_name",
                           action="store",
                           type=str,
                           required=True,
                           help=("Sample_id of control file."))
    poremeth2_input.add_argument("--out_dir", "-o",
                           action="store",
                           type=str,
                           required=True,
                           help="The path to the output directory")
    poremeth2_input.add_argument("--out_prefix", "-op",
                           action="store",
                           type=str,
                           required=True,
                           help="The prefix for the output files")
    
    poremeth2_opt = parser.add_argument_group("optional arguments")
    poremeth2_opt.add_argument("--Rscript", "-rs",
                         action="store",
                         type=str,
                         required=False,
                         default="Rscript",
                         help="The path to a particular instance of "
                              "Rscript to use")
    poremeth2_opt.add_argument("--script_file", "-sf",
                         action="store",
                         type=str,
                         required=False,
                         default=os.path.join(os.path.dirname(
                             os.path.realpath(__file__)
                         ),
                             "PoreMeth2.R"),
                         help="The path to the PoreMeth2.R script file")
        
    poremeth2_opt.add_argument("--betacovThr",
                            action="store",
                            type=int,
                            default=0,
                            required=False,
                            help=("A threshold to filter input data based on beta_cov."
                                "Default is 0.1."))
    poremeth2_opt.add_argument("--omega",
                         action="store",
                         type=float,
                         default=0.1,
                         required=False,
                         help=("Optional parameter that modulates the relative weight "
                               "between the experimental and the biological variance. "
                               "When omega is close to 1, the biological variance is "
                               "much larger than the experimental one, while "
                               "for values of omega close to 0 the experimental noise "
                               "gives the leading contribution to the total variance. "
                               "We suggest to use omega in the range 0.1-0.5."
                               "Default is 0.1."))
    poremeth2_opt.add_argument("--eta",
                         action="store",
                         type=float,
                         default=1e-05,
                         required=False,
                         help=("Optional parameter that represents the baseline probability "
                               "the mean process (m_i) changes its value for the HSLM "
                               "algorithm. Suggested values are inside 10^-7-10^-3 range."
                               "Default is 1e-05."))
    poremeth2_opt.add_argument("--FW",
                         action="store",
                         type=int,
                         default=3,
                         required=False,
                         help=("The minimum number of datapoints for a DMR to be called "
                               "(DMRs made of a number of CpGs smaller than FW are discarded)."
                               "Default is 3."))
    poremeth2_opt.add_argument("--AnnotationType",
                         action="store",
                         type=str,
                         default="Genes",
                         required=False,
                         help=("Optional argument that specifies whether to annotate "
                               "DMRs on genic elements only ('Genes') or genic elements "
                               "and regulatory features ('GenesReg'). Default is Genes."))
    poremeth2_opt.add_argument("--Annotate_Assembly",
                         action="store",
                         type=str,
                         default="hg19",
                         required=False,
                         help=("Optional parameter that specifies the reference version to "
                               "use for annotation (‘hg19’ or ‘hg38’). Default is hg19."))
    poremeth2_opt.add_argument("--Statistics_Assembly",
                         action="store",
                         type=str,
                         default="hg19",
                         required=False,
                         help=("Optional parameter that specifies the reference version to "
                               "use for statistics (‘hg19’ or ‘hg38’). Default is hg19."))
    poremeth2_opt.add_argument("--BetaThr",
                         action="store",
                         type=float,
                         default=0.2,
                         required=False,
                         help=("Delta beta threshold applied for DMRs' classification. "
                               "Default is 0.2."))
    poremeth2_opt.add_argument("--EntropyThr",
                         action="store",
                         type=float,
                         default=0.1,
                         required=False,
                         help=("Delta S threshold applied for DMRs' classification. "
                               "Default is 0.2."))
    poremeth2_opt.add_argument("--PValueThr",
                         action="store",
                         type=float,
                         default=0.05,
                         required=False,
                         help=("The p.value threshold to consider a DMR. "
                               "Default is 0.05."))
    poremeth2_opt.add_argument("--AnalysisClass",
                         action="store",
                         type=str,
                         default="All",
                         required=False,
                         help=("Define the summary to return ('All', 'Beta', 'Entropy'). "
                               "Default is 'All'"))
    poremeth2_opt.add_argument("--overwrite", "-ow",
                         action="store_true",
                         required=False,
                         help="If output files exist overwrite them")

    args = parser.parse_args()

    main_poremeth2(args)


if __name__ == '__main__':
    main()
