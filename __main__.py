import argparse
import idr_pipeline
import extra_funcs
import os
import sys


def main(prefix, w_dir, t_docs, t_names, c_docs, c_prefix, genome_size,
         macs2_file_type, self_threshold, true_threshold,
         pool_threshold, chrom_sizes, filter_chr, filter_neg,
         peak_type, other_macs2_params, results_file, idr_path,
         macs2_path, bedtools_path, samtools_path,
         filter_using_idr):
    idrtype="PYTHON"

    if len(t_docs) < 2:
        sys.exit("There should be at least 2 true replicates to perform IDR")

    t_names_list=[]
    if not t_names:
        for i in range(0, len(t_docs)):
            t_names_list.append("treatment%i" % i)
    else:
        if len(t_names) != len(t_docs):
            sys.exit("--t_names and --t_docs must be of equal length!")
        else:
            t_names_list=t_names

    macs2_file_type = macs2_file_type.upper()
    if macs2_file_type != "BAM":
        if macs2_file_type != "BED":
            if macs2_file_type != "BEDPE":
                if macs2_file_type != "BAMPE":
                    sys.exit("--macs2_file_type must be BAM, BED, BAMPE, or BEDPE")

    peak_type = peak_type.upper()
    if peak_type != "NARROW":
        if peak_type != "BROAD":
            sys.exit("--peak_type must be NARROW or BROAD")

    if not results_file:
        results_file = ("%s/%s_results.txt" % (w_dir,prefix))

    if not extra_funcs.between0and1(self_threshold):
        sys.exit("selfPseud_threshold must be a float between 0 and 1")

    if not extra_funcs.between0and1(true_threshold):
        sys.exit("trueRep_threshold must be a float between 0 and 1")

    if not extra_funcs.between0and1(pool_threshold):
        sys.exit("poolPseud_threshold must be a float between 0 and 1")

    for t_file in t_docs:
        if not os.path.isfile(t_file):
            sys.exit("Treatment file %s does not exist" % t_file)

    if len(c_docs) > 0:
        for c_file in c_docs:
            if not os.path.isfile(c_file):
                sys.exit("Control file %s does not exist" % c_file)

    if not os.path.isfile(chrom_sizes):
        sys.exit("Chromosome size file %s does not exist" % chrom_sizes)

    idr_pipeline.run(idrtype, prefix, w_dir, t_docs, t_names_list,
        c_docs, c_prefix, genome_size, macs2_file_type, self_threshold,
        true_threshold, pool_threshold, chrom_sizes, filter_chr,
        filter_neg, peak_type, other_macs2_params,
        results_file, idr_path, macs2_path, bedtools_path,
        samtools_path)

    if filter_using_idr:
        idr_pipeline.filter_with_idr(idrtype, prefix, w_dir, t_names,
                                     self_threshold, true_threshold,
                                     pool_threshold, peak_type,
                                     bedtools_path)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This script generates \
        self-pseudo replicates and pooled-pseudo replicates for a ChIP-seq \
        experiment. It then run MACS2 to call ChIP peaks on the true replicates, \
        pooled replicates, self-pseudo replicates, and pooled-pseudo replicates. \
        Next, IDR is run on all pairs of true replicates, self-pseudoreplicates, \
        ad pooled pseudo replicates. The results file gives the number of peaks \
        that pass IDR for each pairwise comparison. Using the maximum number of \
        peaks that pass IDR given two true replicates, a conservative set of peaks \
        is output from the pooled MACS2 results.")
    parser.add_argument('prefix', help='In this script we generate a few \
        files that will contain this prefix in front of it. I suggest using a prefix \
        that represents the experiment')
    parser.add_argument('working_dir', help='Directory to output files')
    parser.add_argument('selfPseud_threshold', help='IDR Threshold (p-value) for \
        self-pseudo replicates', type = float)
    parser.add_argument('trueRep_threshold', help='IDR Threshold (p-value) for \
        true replicates', type = float)
    parser.add_argument('poolPseud_threshold', help='IDR Threshold (p-value) for \
        pooled-pseudo replicates', type = float)
    parser.add_argument('chrom_sizes', help='file which contains the length of the \
        chromosomes')
    parser.add_argument('-td','--t_docs', help='bam or bed files for treatment \
        samples', nargs='+', required=True)
    parser.add_argument('-cd','--c_docs', help='bam or bed files for input \
        (control) samples', nargs='*', required=False, default=[])
    parser.add_argument('-tn','--t_names', help='names for treatment samples (list \
        must be equal to t_docs list)', nargs='+', required=False)
    parser.add_argument('-cp','--c_prefix', help='input (control samples) will be \
        pooled into one mapping file with this prefix and .bam or .bed as its suffix \
        restively (default:control)', default="control", required=False)
    parser.add_argument('-gs','--genome_size', help='effective genome size \
        (default: 101274395 which is the effective genome size of Arabidopsis)', \
        type=int, required=False, default=101274395)
    parser.add_argument('-ft','--macs2_file_type', help='MACS2 file type: BAM or \
        BED (default:BAM)', type=str, required=False, default="BAM")
    parser.add_argument('-fc','--filter_chr', help='list of chromosomes to filter out \
        after running MACS2', nargs='*', required=False)
    parser.add_argument('-fn','--filter_neg', help='negative control bed file to \
        filter out after running', nargs='*', required=False)
    parser.add_argument('-pt','--peak_type', help='NARROW or BROAD peak type \
        (default: NARROW)', required=False, default="NARROW")
    parser.add_argument('-omp', '--other_macs2_params', help="other parameters to use \
        in macs2 (see their github page for more info)")
    parser.add_argument('-O','--results_file', help='Name of IDR results file; \
        Default: $expname_results.txt', required=False)
    parser.add_argument('--filter_using_idr', help='filter using idr',
                        action='store_true')
    parser.add_argument('-mp','--macs2_path', help='path to macs2', default = "")
    parser.add_argument('-ip','--idr_path', help='path to idr scripts', default = "")
    parser.add_argument('-bp','--bedtools_path', help='path to bedtools', default = "")
    parser.add_argument('-sp','--samtools_path', help='path to samtools', default = "")

    args = parser.parse_args()

    main(args.prefix, args.working_dir, args.t_docs, args.t_names,
         args.c_docs, args.c_prefix, args.genome_size, args.macs2_file_type, args.selfPseud_threshold,
         args.trueRep_threshold, args.poolPseud_threshold, args.chrom_sizes, args.filter_chr,
         args.filter_neg, args.peak_type, args.other_macs2_params,
         args.results_file, args.idr_path, args.macs2_path, args.bedtools_path,
         args.samtools_path, args.filter_using_idr)
