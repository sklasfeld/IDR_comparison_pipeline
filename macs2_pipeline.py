#!/usr/bin/env python

# author: Samantha Klasfeld
# date: 10/2/2017

import os
import numpy as np
import sys
import pandas as pd

# Using MACS2 output filter out specific features
# narrowPeakFile = output from MACS2
# prefix = prefix for output files
# outputdir = output directory for MACS2
# peak_type = "narrow" or "broad"
# qval2 = cutoff to the peak summit region
# filter_chr = filter list of chrs
# filter_neg = filter negative ctrl
# bedtools_path = PATH to bedtools scripts
def filternarrowpeak(narrowPeakFile, prefix, outputdir, peak_type,
                     qval2, filter_chr, filter_neg, bedtools_path):
    np_table = pd.read_csv(narrowPeakFile, sep='\t', header=None,
                           names=["chr", "start", "stop", "name", "signal", "strand",
                                  "signalValue", "pValue", "qValue", "peak"],
                           dtype={'chr':'str','start':'int', 'stop':'int',
                                  'name':'str', 'signal':'float', 'strand':'str',
                                  'signalValue':'float', 'pValue':'float', 'qValue':'float',
                                  'peak':'int'})

    # remove unwanted chromosomes from output
    if filter_chr:
        if len(filter_chr) > 0:
            sys.stdout.write("REMOVE FROM OUTPUT CHROMOSOME(S): %s...\n" % (",".join(filter_chr)))
            np_table = np_table[~np_table["chr"].isin(filter_chr)]

    # filter via qvalue
    if qval2 >= 0:
        sys.stdout.write("Cutoff MACS2 results at %f...\n" % qval2)
        qval2_min = -np.log10(qval2)
        np_table = np_table[np_table["qValue"] >= qval2_min]


    out = ("%s/%s_postFilter.%sPeak" % (outputdir, prefix, peak_type))
    if filter_neg:
        if len(filter_neg) > 0:
            # output current narrowpeak file
            before_subtraction_f = ("%s/%s_prePostFilter.%sPeak" % (outputdir, prefix, peak_type))
            np_table.to_csv(before_subtraction_f, sep="\t", index=False, header=False)
            # run bedtools subtract
            subtract_neg = ("%sbedtools subtract -A -F %s -a %s -b %s > %s"
                % (bedtools_path, filter_neg[1], before_subtraction_f, filter_neg[0], out))
            sys.stdout.write("%s\n" % subtract_neg)
            sys.stdout.flush()
            os.system(subtract_neg)
    else:
        np_table.to_csv(out, sep="\t", index=False, header=False)
    return out

# run MACS2 Pipeline
# t_docs = treatment input files
# t_names = prefixes for each corresponding treatment input file
# c_docs = control input files
# genome_size = effective size of genome
# macs2_file_type = "BED" or "BAM"
# qval2 = cutoff to the peak summit region
# filter_chr = filter list of chrs
# filter_neg = filter negative ctrl
# peak_type = "narrow" or "broad"
# pool_reps = run each treatment file seperately (False) or together (True)
# other_param = other MACS2 parameters
# outputdir = output directory for MACS2
# macs2_path = PATH to macs2 script
# bedtools_path = PATH to bedtools scripts
def macs2_peaks(t_docs, t_names, c_docs, genome_size, macs2_file_type, qval2,
                filter_chr, filter_neg, peak_type, other_param,
                outputdir, macs2_path, bedtools_path):
    peak_type = peak_type.lower()
    params = ""
    if len(c_docs) > 0:
        c_docs_string = " ".join(c_docs)
        params = ("%s -c %s" % (params, c_docs_string))
    if peak_type == "broad":
            params = ("%s --broad " % params)

    params = ("%s -f %s %s" % (params, macs2_file_type, other_param))

    sys.stdout.write("run MACS2 on each of the replicates...\n")
    sys.stdout.flush()
    out_files = []
    for i in range(0, len(t_docs)):
        out1 = ("%s/%s_peaks.%sPeak" % (outputdir, t_names[i], peak_type))
        if not os.path.isfile(out1):
            run_macs2 = ("%smacs2 callpeak -t %s %s -g %s "
                "--outdir %s -n %s" % (macs2_path, t_docs[i],
                                       params.strip(), genome_size,
                                       outputdir, t_names[i]))
            sys.stdout.write("%s\n" % run_macs2)
            sys.stdout.flush()
            os.system(run_macs2)

            if qval2 >= 0 or len(filter_chr) > 0:
                out2 = filternarrowpeak(out1, t_names[i],
                                        outputdir, peak_type,
                                        qval2, filter_chr,
                                        filter_neg, bedtools_path)
                out_files.append(out2)
            else:
                out_files.append(out1)
    return out_files

def macs2_pooled_peaks(prefix, t_docs, t_names, c_docs, genome_size, macs2_file_type, qval2,
    filter_chr, filter_neg, peak_type, other_param,
    outputdir, macs2_path, bedtools_path):
        peak_type = peak_type.lower()
        params = ""
        if len(c_docs) > 0:
            c_docs_string = " ".join(c_docs)
            params = ("%s -c %s" % (params, c_docs_string))
        if peak_type == "broad":
                params = ("%s --broad " % params)

        params = ("%s -f %s %s" % (params, macs2_file_type, other_param))
        chipout = ("%s/%s_peaks.%sPeak" % (outputdir, prefix, peak_type))
        if not os.path.isfile(chipout):
            t_docs_string = " ".join(t_docs)
            sys.stdout.write("running MACS2 on pooled treatment files: %s\n" %
                            t_docs_string)
            sys.stdout.flush()
            # pool the treatment files
            params = ("-t %s %s" % (t_docs_string, params.strip()))
            run_macs2 = ("%smacs2 callpeak %s -g %s --outdir %s -n %s"
                     % (macs2_path, params, genome_size, outputdir, prefix))
            sys.stdout.write("%s\n" % run_macs2)
            sys.stdout.flush()
            os.system(run_macs2)
        if qval2 >= 0 or filter_chr or filter_neg:
            out = filternarrowpeak(chipout, prefix, outputdir, peak_type,
                                          qval2, filter_chr, filter_neg,
                                          bedtools_path)
            return out
        else:
            return chipout

