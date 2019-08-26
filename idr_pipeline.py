import os
import sys
from __builtin__ import float
from __builtin__ import int

import macs2_pipeline
import extra_funcs
import subprocess
import re
import numpy as np
def pass_idr(fname, threshold, idrtype):
    assert not extra_funcs.empty(fname), ("Missing File %s!" % fname)
    pass_idr_out = 0
    if idrtype == "PYTHON":
        cmd = ("""awk '(10**-$12) <= %f {print $0}' %s | wc -l"""
               % (threshold, fname))
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        pass_idr_out = int(output.strip())

    # THIS R SCRIPT IS DEPRECATED
    # if idrtype == "RSCRIPT":
    #     cmd = ("""awk '$11 <= %f {print $0}' %s | wc -l"""
    #        % (threshold, fname))
    #     ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
    #                       stderr=subprocess.STDOUT)
    #     output = ps.communicate()[0]
    #     pass_idr_out = int(output.strip()) - 1

    return  pass_idr_out

def run_idrscript(idrtype, idrpath, peak_type, peakfile1, \
                  peakfile2, threshold, idr_result_prefix, \
                  scriptpath, chromsizes):
    out = ("%s-overlapped-peaks.txt" % idr_result_prefix)
    if extra_funcs.wc(out) == 0:
        run_idr_script = (("%sidr --samples "
                          "%s %s "
                          "--rank p.value --plot "
                          "--input-file-type %sPeak "
                          "--soft-idr-threshold %s "
                          "-o %s")
                         % (idrpath, peakfile1, peakfile2,
                            peak_type.lower(), threshold,
                            out))
        # THIS RSCRIPT IS DEPRECATED
        # if idrtype == "RSCRIPT":
        #     bro = "F"
        #     if peak_type == "broad":
        #         bro = "T"
        #     run_idr_script = (("Rscript %sbatch-consistency-analysis.r "
        #                       "%s %s -1 %s 0 %s p.value %s %s")
        #                      % (idrpath, peakfile1, peakfile2,
        #                         idr_result_prefix, bro,
        #                         scriptpath, chromsizes))
        extra_funcs.runcmdline(run_idr_script)
    return out

# run IDR pipeline
#   idrtype = PYTHON or RSCRIPT
# 	expName = name of the experiment
# 	t_docs = treatment input files
# 	t_names = prefixes for each corresponding treatment input file
# 	c_docs = control input files
# 	genome_size = effective size of genome
# 	macs2_file_type = "BED" or "BAM"
# 	selfPseud_threshold = IDR Threshold (p-value) for self-pseudo replicates
# 	trueRep_threshold = IDR Threshold (p-value) for true replicates
#	poolPseud_threshold = IDR Threshold (p-value) for pooled-pseudo replicates
#	chrom_sizes = file which contains the length of the chromosomes
# 	filter_chr = list of chromosomes to filter out after running MACS2
# 	filter_neg = negative control to filter out after running
# 	peak_type = "narrow" or "broad"
# 	other_macs2_params = other parameters to run with macs2
#	results_file = Name of IDR results file; Default: $expName_results.txt
# 	idr_path = PATH to idr scripts
# 	macs2_path = PATH to macs2 script
# 	bedtools_path = PATH to bedtools scripts
#	samtools_path = PATH to samtools scripts
def run(idrtype, expName, working_dir, t_docs, t_names, c_docs, c_prefix, genome_size,
                 macs2_file_type, self_threshold, true_threshold,
                 pool_threshold, chrom_sizes, filter_chr, filter_neg,
                 peak_type, other_macs2_params, results_file, idr_path,
                 macs2_path, bedtools_path, samtools_path):
    assert len(t_docs) == len(t_names), "TDOCS and TNAMES must be equal"
    sys.stdout.write("\nSETTINGS\n")
    sys.stdout.write("----------------\n")
    sys.stdout.write("PREFIX: %s\n" % expName)
    sys.stdout.write("WORKING DIRECTORY: %s\n" % working_dir)    
    sys.stdout.write("TREATMENT FILES: %s\n" % ",".join(t_docs))
    sys.stdout.write("TREATMENT NAMES: %s\n" % ",".join(t_names))
    sys.stdout.write("CONTROL FILES: %s\n" % ",".join(c_docs))
    sys.stdout.write("GENOME SIZE: %i\n" % int(genome_size))
    sys.stdout.write("MACS2 FILE TYPE: %s\n" % macs2_file_type)
    sys.stdout.write("SELF PSEUDO REPS THRESHOLD: %f\n" % float(self_threshold))
    sys.stdout.write("TRUE REPS THRESHOLD: %f\n" % float(true_threshold))
    sys.stdout.write("POOLED PSEUDO REPS THRESHOLD: %f\n" % float(pool_threshold))
    sys.stdout.write("CHROM SIZE FILE: %s\n" % chrom_sizes)
    sys.stdout.write("PEAK TYPE: %s\n" % peak_type)
    sys.stdout.write("OTHER MACS2 PARAMS: %s\n" % other_macs2_params)
    sys.stdout.write("RESULTS FILE: %s\n\n" % results_file)


    # SET MACS2 PARAMETERS
    qval2 = -1
    peak_type = peak_type.lower()
    if len(c_docs) > 1:
        c_docs_str = " ".join(c_docs)
        cdocs_pooled = ("%s/%s.bam" % (working_dir,c_prefix))
        if not os.path.isfile(cdocs_pooled):
            make_ctrl = ("""%ssamtools merge %s %s""" % (samtools_path,
                                                         cdocs_pooled,
                                                         c_docs_str))
            sys.stdout.write("%s\n" % make_ctrl)
            sys.stdout.flush()
            os.system(make_ctrl)

        c_docs = [cdocs_pooled]

    # SET OUTPUT FILE
    if len(results_file) == 0:
        results_file = ("%s/%s_results.txt" % (working_dir, expName))

    file_r = open(results_file, 'a')
    
    # SET IDR PARAMETERS
    idrCode_p = idr_path.split("/")
    script_path = "/".join(idrCode_p[:-2])
    num_reps = len(t_docs)

    # 1. DEAL WITH TRUE REPLICATES
    file_r.write("# TRUE REPS (p <= %.3f)\n" % float(true_threshold))

    if not os.path.isdir(working_dir):
        mkdir_trueRep = ("""mkdir %s""" % working_dir)
        sys.stdout.write("%s\n" % mkdir_trueRep)
        sys.stdout.flush()
        os.system(mkdir_trueRep)

    ## 1A. MACS2 with TRUE REPLICATES
    sys.stdout.write("True replicates peak calling...\n")
    sys.stdout.flush()
    trueReps_macs_done = True
    for true_rep_prefix in t_names:
        true_rep_out = ("%s/%s.%sPeak" % (working_dir, true_rep_prefix, peak_type))
        if not os.path.isfile(true_rep_out):
            trueReps_macs_done = False
    if not trueReps_macs_done:
        macs2_pipeline.macs2_peaks(t_docs, t_names, c_docs, genome_size,
                                   macs2_file_type, qval2, filter_chr, filter_neg,
                                   peak_type, other_macs2_params, working_dir,
                                   macs2_path, bedtools_path)

    ## 1B. SORTING with TRUE REPLICATES (PREP FOR IDR)

    for t_i in range(0, len(t_docs)):
        sys.stdout.write("Sorting peaks for TrueRep: %s...\n" % t_names[t_i])
        sys.stdout.flush()
        out = ("%s/%s_sorted_peaks.%sPeak" % (working_dir, t_names[t_i],
                                              peak_type))
        if not os.path.isfile(out) or extra_funcs.empty(out):
            sys.stdout.write("sorted file: %s\n" % out)
            sys.stdout.flush()
            in_narrowPeak = (("%s/%s_postFilter.%sPeak") % (working_dir,
                                                            t_names[t_i],
                                                            peak_type))
            if not os.path.isfile(in_narrowPeak):
                in_narrowPeak = ("%s/%s_peaks.%sPeak" % (working_dir,
                                                         t_names[t_i],
                                                         peak_type))

            run_sort = ("""sort -k 8nr,8nr %s  | head -n 100000 > %s""" %
                        (in_narrowPeak, out))
            sys.stdout.write("%s\n" % run_sort)
            sys.stdout.flush()
            os.system(run_sort)

    ## 1C. IDR with TRUE REPLICATES

    sys.stdout.write("Running IDR for True Replicates...\n")
    sys.stdout.flush()
    for t_i in range(0, len(t_docs) - 1):
        for t_j in range(t_i + 1, len(t_docs)):
            idr_result_prefix = ("%s/%s_v_%s" % (working_dir, t_names[t_i], t_names[t_j]))
            peakf1 = ("%s/%s_sorted_peaks.%sPeak" % (working_dir, t_names[t_i], peak_type))
            peakf2 = ("%s/%s_sorted_peaks.%sPeak" % (working_dir, t_names[t_j], peak_type))
            out = run_idrscript(idrtype, idr_path, peak_type, peakf1, peakf2, \
                                true_threshold, idr_result_prefix, script_path, chrom_sizes)
            pass_threshold = pass_idr(out, true_threshold, idrtype)
            out_basename = os.path.splitext(os.path.basename(out))[0]
            file_r.write("%s\t%s\n" % (
                out_basename.replace("-overlapped-peaks",""), 
                pass_threshold))


    # 2. DEAL WITH SELF PSEUDO REPLICATES
    file_r.write("# SELF PSEUDO REPS (p <= %.3f)\n" % float(self_threshold))
    
    ## 2A. CREATING SELF PSEUDO REPLICATES
    idr_prefix = []
    for t_i in range(0, len(t_docs)):
        sys.stdout.write("PseudoReplicates of true replicate (%s)\n" \
                         % t_names[t_i])
        sys.stdout.flush()

        file_name = t_docs[t_i]

        output_stub = t_names[t_i]  # prefix for pseudoReplicate files

        out1 = ("%s/%s.pr1.tagAlign.bam" % (working_dir, output_stub))
        out2 = ("%s/%s.pr2.tagAlign.bam" % (working_dir, output_stub))
        if (not os.path.isfile(out1) and not os.path.isfile(out2)) or \
            extra_funcs.empty(out1, samtools_path) or extra_funcs.empty(out2, samtools_path):
            sys.stdout.write("counting lines in %s...\n" % file_name)
            sys.stdout.flush()
            nlines = extra_funcs.split_f(file_name, num_reps, samtools_path)
            # nlines = extra_funcs.split_f(file_name, 2, samtools_path)
            sys.stdout.write("shuffle the lines in the file and split it into two parts\n")
            sys.stdout.flush()
            extra_funcs.shufNsplitBam(file_name, nlines, working_dir, output_stub, samtools_path)
        if (not os.path.isfile(out1) and not os.path.isfile(out2)) or \
            extra_funcs.empty(out1, samtools_path) or extra_funcs.empty(out2, samtools_path):
            err_mes = (("ERROR: %s and/or %s were not built. " + \
                "There should be Error Messages Above.\n") % (out1,out2))
            sys.exit(err_mes)

        ## 2B. MACS2 with SELF PSEUDO REPLICATES

        sys.stdout.write("Calling peaks for self pseudo-replicate %s_pr1...\n"
                         % (output_stub))
        sys.stdout.flush()
        out1 = ("%s/%s_pr1_peaks.%sPeak" % (working_dir, output_stub, peak_type))
        self_pseudo_inputs = []
        self_pseudo_prefixs = []
        run_macs = 0
        if not os.path.isfile(out1):
            run_macs = 1
            ip_file = ("%s/%s.pr1.tagAlign.bam" % (working_dir, output_stub))
            macs2_prefix = ("%s_pr1" % output_stub)
            self_pseudo_inputs.append(ip_file)
            self_pseudo_prefixs.append(macs2_prefix)
        else:
            sys.stderr.write("Warning: Will not override previous "
                             "version of %s\n" % out1)

        sys.stdout.write("Calling peaks for self pseudo-replicate %s_pr2...\n"
                         % (output_stub))
        sys.stdout.flush()
        out2 = ("%s/%s_pr2_peaks.%sPeak" % (working_dir, output_stub, peak_type))
        if not os.path.isfile(out2):
            run_macs = 1
            ip_file = ("%s/%s.pr2.tagAlign.bam" % (working_dir, output_stub))
            macs2_prefix = ("%s_pr2" % output_stub)
            self_pseudo_inputs.append(ip_file)
            self_pseudo_prefixs.append(macs2_prefix)
        else:
            sys.stderr.write("Warning: Will not override previous "
                             "version of %s\n" % out2)
            sys.stdout.flush()

        if run_macs == 1:
            macs2_pipeline.macs2_peaks(self_pseudo_inputs, self_pseudo_prefixs,
                                       c_docs, genome_size, macs2_file_type,
                                       qval2, filter_chr, filter_neg,
                                       peak_type, other_macs2_params,
                                       working_dir, macs2_path, bedtools_path)

        ## 2C. SORTING with SELF PSEUDO REPLICATES (PREP FOR IDR)
        sys.stdout.write("Sorting peaks for self-pseudo replicate1 of %s...\n"
                         % t_names[t_i])
        sys.stdout.flush()

        out1 = ("%s/%s_pr1.sorted_peaks.%sPeak" % (working_dir, output_stub,
                                                   peak_type))
        if not os.path.isfile(out1) or extra_funcs.empty(out1):
            in_narrowPeak = ("%s/%s_pr1_postFilter.%sPeak" % (working_dir,
                                                              output_stub,
                                                              peak_type))
            if not os.path.isfile(in_narrowPeak):
                in_narrowPeak = ("%s/%s_pr1_peaks.%sPeak" % (working_dir,
                                                             output_stub,
                                                             peak_type))
            run_sort = ("""sort -k 8nr,8nr %s  | head -n 100000 > %s"""
                        % (in_narrowPeak, out1))
            sys.stdout.write("%s\n" % run_sort)
            sys.stdout.flush()
            os.system(run_sort)

        sys.stdout.write("Sorting peaks for self-pseudo replicate2 of %s...\n"
                         % t_names[t_i])
        sys.stdout.flush()

        out2 = ("%s/%s_pr2.sorted_peaks.%sPeak" % (working_dir, output_stub,
                                                   peak_type))
        if not os.path.isfile(out2) or extra_funcs.empty(out2):
            in_narrowPeak = ("%s/%s_pr2_postFilter.%sPeak" % (working_dir,
                                                              output_stub,
                                                              peak_type))
            if not os.path.isfile(in_narrowPeak):
                in_narrowPeak = (("%s/%s_pr2_peaks.%sPeak") % (working_dir,
                                                              output_stub,
                                                              peak_type))
            run_sort = ("""sort -k 8nr,8nr %s  | head -n 100000 > %s"""
                        % (in_narrowPeak, out2))
            sys.stdout.write("%s\n" % run_sort)
            sys.stdout.flush()
            os.system(run_sort)

        ## 2D. IDR with SELF PSEUDO REPLICATES
        sys.stdout.write("Running IDR for Self-Pseudo Replicates for %s...\n"
                         % t_names[t_i])
        sys.stdout.flush()
        idr_result_prefix = ("%s/%s" % (working_dir, output_stub))
        #output_stub = ('chipSampleRep%s' % str(t_i + 1))  # prefix for pseudoReplicate files
        peakf1 = ("%s/%s_pr1.sorted_peaks.%sPeak" % (working_dir, output_stub, peak_type))
        peakf2 = ("%s/%s_pr2.sorted_peaks.%sPeak" % (working_dir, output_stub, peak_type))
        out = run_idrscript(idrtype, idr_path, peak_type, peakf1, peakf2, \
                                self_threshold, idr_result_prefix, script_path, chrom_sizes)
        pass_threshold = pass_idr(out, self_threshold, idrtype)
        out_basename = os.path.splitext(os.path.basename(out))[0]
        file_r.write("%s\t%s\n" % (
            out_basename.replace("-overlapped-peaks",""), 
            pass_threshold))

    ## 	2E. CREATE PLOT FOR IDR FOR SELF PSEUDO REPLICATES
    #idr_prefixes = " ".join(idr_prefix)
    #out = ("%s/%s-plot.ps" % (working_dir, working_dir))
    #if not os.path.isfile(out):
    #    run_plot = ("Rscript %sbatch-consistency-plot.r %s %i %s/%s %s"
    #                % (idr_path, script_path, num_reps, working_dir,
    #                   working_dir, idr_prefixes))
    #    sys.stdout.write(("%s\n" % run_plot))
    #    sys.stdout.flush()
    #    os.system(run_plot)


    # 3. DEAL WITH POOLED REPLICATES

    ### 3A. POOL REPLICATES

    output_stub = ("%s_pooled" % expName)  # pooled prefix
    sys.stdout.write("pooling...\n")
    sys.stdout.flush()
    pooled_bam = ("%s/%s.bam" % (working_dir, output_stub))
    if not os.path.isfile(pooled_bam):
        all_tBAMfiles = " ".join(t_docs)
        pool_tFils = ("%ssamtools merge %s %s" % (samtools_path,
                                                  pooled_bam, all_tBAMfiles))
        sys.stdout.write("%s\n" % pool_tFils)
        sys.stdout.flush()
        os.system(pool_tFils)
    if not os.path.isfile(pooled_bam):
        sys.exit("ERROR: %s was not built. There should be Error Messages Above.\n" % pooled_bam)

    ## 3B. MACS2 with POOLED REPLICATES

    sys.stdout.write("Calling peaks for pooled data of all replicates\n")
    sys.stdout.flush()

    out = ("%s/%s_peaks.%sPeak" % (working_dir, output_stub, peak_type))
    if not os.path.isfile(out):
        macs2_pipeline.macs2_pooled_peaks(output_stub, t_docs, t_names, c_docs, genome_size,
                                   macs2_file_type, qval2, filter_chr, filter_neg,
                                   peak_type, other_macs2_params, working_dir,
                                   macs2_path, bedtools_path)
    else:
        sys.stderr.write("Warning: Will not override previous "
                         "version of %s\n" % out)





    # 4. DEAL WITH POOLED PSEUDO REPLICATES
    file_r.write("# POOLED PSEUDO REPS (p <= %.3f)\n" % float(pool_threshold))

    sys.stdout.write("PseudoReplicates of Pooled Replicates\n")
    sys.stdout.flush()

    ## 4A. CREATING POOLED PSEUDO REPLICATES



    ### 4Ai. SPLIT & SHUFFLE
    sys.stdout.write("shuffle the lines in the file and split it into X "
                     "parts where X is the number of IP replicates\n")
    sys.stdout.flush()
    out1 = ("%s/%s.pr1.tagAlign.bam" % (working_dir, output_stub))
    out2 = ("%s/%s.pr2.tagAlign.bam" % (working_dir, output_stub))
    if not (os.path.isfile(out1) and os.path.isfile(out2)):
        sys.stdout.write("counting...\n")
        sys.stdout.flush()
        nlines = extra_funcs.split_f(pooled_bam, num_reps, samtools_path)
        # nlines = extra_funcs.split_f(pooled_bam, 2, samtools_path)
        if nlines == 0:
            no_file_err = ("no contents in %s" % pooled_bam)
            sys.exit(no_file_err)
        extra_funcs.shufNsplitBam(pooled_bam, nlines, working_dir, output_stub, samtools_path)
    if not (os.path.isfile(out1) and os.path.isfile(out2)):
        err_mes = (("ERROR: %s and/or %s were not built. " + \
            "There should be Error Messages Above.\n") % (out1,out2))
        sys.exit(err_mes)

    ## 4B. MACS2 with POOLED PSEUDO REPLICATES
    pool_pseudo_inputs = []
    pool_pseudo_prefixs = []
    run_macs = 0

    out1 = ("%s/%s_pr1_peaks.%sPeak" % (working_dir, output_stub, peak_type))
    if not os.path.isfile(out1):
        run_macs = 1
        ip_file = ("%s/%s.pr1.tagAlign.bam" % (working_dir, output_stub))
        macs2_prefix = ("%s_pr1" % output_stub)
        pool_pseudo_inputs.append(ip_file)
        pool_pseudo_prefixs.append(macs2_prefix)
    else:
        sys.stderr.write("Warning: Will not override previous "
                         "version of %s\n" % out1)

    out2 = ("%s/%s_pr2_peaks.%sPeak" % (working_dir, output_stub, peak_type))
    if not os.path.isfile(out2):
        run_macs = 1
        ip_file = ("%s/%s.pr2.tagAlign.bam" % (working_dir, output_stub))
        macs2_prefix = ("%s_pr2" % output_stub)
        pool_pseudo_inputs.append(ip_file)
        pool_pseudo_prefixs.append(macs2_prefix)
    else:
        sys.stderr.write("Warning: Will not override previous "
                         "version of %s\n" % out2)
        sys.stdout.flush()

    sys.stdout.write("Calling peaks for pooled pseudo replicates...\n")
    sys.stdout.flush()

    if run_macs == 1:
        macs2_pipeline.macs2_peaks(pool_pseudo_inputs, pool_pseudo_prefixs,
                                   c_docs, genome_size, macs2_file_type,
                                   qval2, filter_chr, filter_neg,
                                   peak_type, other_macs2_params,
                                   working_dir, macs2_path, bedtools_path)

    ## 4C. SORTING with POOLED PSEUDO REPLICATES (PREP FOR IDR)
    sys.stdout.write("Sorting peaks of pooled pseudo replicate 1...\n")
    sys.stdout.flush()

    out1 = ("%s/%s_pr1.sorted_peaks.%sPeak" % (working_dir, output_stub, peak_type))

    if not os.path.isfile(out1) or extra_funcs.empty(out1):
        in_narrowPeak = ("%s/%s_pr1_postFilter.%sPeak" % (working_dir,
                                                          output_stub,
                                                          peak_type))
        if not os.path.isfile(in_narrowPeak):
            in_narrowPeak = ("%s/%s_pr1_peaks.%sPeak" % (working_dir,
                                                        output_stub,
                                                        peak_type))
        run_sort = ("""sort -k 8nr,8nr %s  | head -n 100000 > %s"""
                    % (in_narrowPeak, out1))
        sys.stdout.write("%s\n" % run_sort)
        sys.stdout.flush()
        os.system(run_sort)
    if not os.path.isfile(out1) or extra_funcs.empty(out1):
        sys.exit(("\nCANNOT RUN IDR because %s was not created!\n" % out1))

    sys.stdout.write("Sorting peaks of pooled pseudo replicate 2...\n")
    sys.stdout.flush()

    out2 = ("%s/%s_pr2.sorted_peaks.%sPeak" % (working_dir, output_stub,
                                               peak_type))
    if not os.path.isfile(out2) or extra_funcs.empty(out2):
        in_narrowPeak = ("%s/%s_pr2_postFilter.%sPeak" % (working_dir,
                                                          output_stub,
                                                          peak_type))
        if not os.path.isfile(in_narrowPeak):
            in_narrowPeak = ("%s/%s_pr2_peaks.%sPeak" % (working_dir,
                                                         output_stub,
                                                         peak_type))
        run_sort = ("""sort -k 8nr,8nr %s  | head -n 100000 > %s"""
                    % (in_narrowPeak, out2))
        sys.stdout.write("%s\n" % run_sort)
        sys.stdout.flush()
        os.system(run_sort)
    if not os.path.isfile(out2) or extra_funcs.empty(out2):
        sys.exit(("\nCANNOT RUN IDR because %s was not created!\n" % out1))

    ## 4D. IDR with POOLED PSEUDO REPLICATES
    sys.stdout.write("Running IDR for Pooled Pseudo Replicates...\n")
    sys.stdout.flush()

    idr_prefix = ("%s/%s" % (working_dir, output_stub))  # pooled pseudo prefix
    peakf1 = ("%s/%s_pr1.sorted_peaks.%sPeak" % (working_dir, output_stub, peak_type))
    peakf2 = ("%s/%s_pr2.sorted_peaks.%sPeak" % (working_dir, output_stub, peak_type))
    out = run_idrscript(idrtype, idr_path, peak_type, peakf1, peakf2, \
                        pool_threshold, idr_prefix, script_path, chrom_sizes)
    pass_threshold = pass_idr(out, pool_threshold, idrtype)
    out_basename = os.path.splitext(os.path.basename(out))[0]
    file_r.write("%s\t%s\n" % (
        out_basename.replace("-overlapped-peaks",""),
        pass_threshold))

    ## 	3E. CREATE PLOT FOR IDR FOR POOLED PSEUDO REPLICATES
    #out = ("%s/%s-plot.ps" % (working_dir, working_dir))
    #if not os.path.isfile(out):
    #    run_plot = ("Rscript %sbatch-consistency-plot.r %s %i %s/%s %s"
    #                % (idr_path, script_path, 1, working_dir,
    #                   working_dir, idr_result_prefix))
    #    sys.stdout.write(("%s\n" % run_plot))
    #    sys.stdout.flush()
    #    os.system(run_plot)

    file_r.close()


def filter_with_idr(idrtype, expName, working_dir, tnames, self_threshold, true_threshold, 
        pool_threshold, peak_type, bedtools_path):
    peak_type = peak_type.lower()
    numTreats = len(tnames)
    # check for necessary macs2 pooled file
    pool_macs = ("%s/%s_pooled_postFilter.%sPeak"
                 % (working_dir, expName, peak_type))
    sys.stdout.write("checking for %s...\n" % pool_macs)
    if not (os.path.isfile(pool_macs)):
        pool_macs = ("%s/%s_pooled_peaks.%sPeak"
                     % (working_dir, expName, peak_type))
        if not (os.path.isfile(pool_macs)):
            sys.exit("Missing %s. Try running pooling the files and running "
                     " the MACS2 command first." \
                     % pool_macs)

    # POOLED-CONSISTENCY FILE
    pool_file = (("%s/%s_pooled-overlapped-peaks.txt") % (working_dir, expName))
    sys.stdout.write("checking for %s...\n" % pool_file)
    sys.stdout.flush()
    if not (os.path.isfile(pool_file)):
        sys.exit("Missing %s. Try running IDR command first." \
                 % pool_file)

    self_files = []
    true_files = []
    for i in range(0, numTreats):
        # SELF-CONSISTENCY FILES
        self_file = ("%s/%s-overlapped-peaks.txt" %
                     (working_dir, tnames[i]))
        self_files.append(self_file)
        sys.stdout.write("checking for %s...\n" % self_file)
        sys.stdout.flush()
        if not (os.path.isfile(self_file)):
            sys.exit("Missing %s. Try running IDR command first."
                     % self_file)
        for j in range(i + 1, numTreats):
            # TRUE-REP FILES
            rep1 = tnames[i]
            rep2 = tnames[j]
            true_file = ("%s/%s_v_%s-overlapped-peaks.txt"
                         % (working_dir, rep1, rep2))
            true_files.append(true_file)
            sys.stdout.write("checking for %s...\n" % true_file)
            sys.stdout.flush()
            if not (os.path.isfile(true_file)):
                sys.exit("Missing %s. Try running IDR command first."
                         % true_file)

    # write idr.threshold file
    sys.stdout.write("write idr.threshold file...\n")
    sys.stdout.flush()
    theshold_file = ("%s/%s.idr.threshold" % (working_dir,expName))
    idr_tf = open(theshold_file, 'w')


    np_p = pass_idr(pool_file, float(pool_threshold), idrtype)
    numPeaks_pool = np_p

    s_num_list = []
    t_num_list = []
    for i in range(0, numTreats):
        rep_val = i + 1
        sys.stdout.write("rep %i...\n" % rep_val)
        sys.stdout.flush()
        nps_self = pass_idr(self_files[i], float(self_threshold), idrtype)

        numPeaks_true = []
        for j in range(0, len(true_files)):
            cur_rep = tnames[i]
            tf = true_files[j]
            if re.search(cur_rep, tf):
                # get list of numPeak for each true replicate idr file
                np_t = pass_idr(tf, float(true_threshold), idrtype)
                sys.stdout.write("%i\n" % int(np_t))
                sys.stdout.flush()
                numPeaks_true.append(int(np_t))

        max_nps_true = 0
        if (len(numPeaks_true) > 0):
            max_nps_true = np.max(numPeaks_true)
        idr_tf.write("%s\trep%i\t%s\t%i\t%i\n" %
                     (expName, rep_val, nps_self,
                      max_nps_true, numPeaks_pool))
        s_num_list.append(int(nps_self))
        t_num_list.append(max_nps_true)

    idr_tf.close()

    # create idr.optimal.threshold
    sys.stdout.write("write idr.optimal.threshold file...\n")
    sys.stdout.flush()
    optimal_file = ("%s/%s.idr.optimal.threshold" % (working_dir, expName))
    idr_optimal = open(optimal_file, 'w')
    s = "|".join(str(x) for x in s_num_list)
    p = numPeaks_pool
    t = np.max(t_num_list)
    conservative_threshold = t
    optimal_threshold = np.max([conservative_threshold, p])
    idr_optimal.write("%s\t%i\t%i\t%s\t%i\t%i\n" % (expName,
                                                    optimal_threshold,
                                                    numTreats, s, t, p))

    # get chipSampleRepAll.optimal.narrowPeak file using IDR results
    optimal_bed_idr = ("%s/%s.idr.optimal.%sPeak" % (working_dir, expName, peak_type))
    final_filter_1 = (("""sort -k8nr,8nr %s > %s.tmp""") % (pool_macs,
                                                            optimal_bed_idr))
    final_filter_2 = (("""head -n %i %s.tmp| %sbedtools sort -i > %s""") \
                      % (optimal_threshold, optimal_bed_idr, bedtools_path, optimal_bed_idr))
    final_filter_3 = ("rm %s.tmp " % optimal_bed_idr)
    sys.stdout.write("%s\n" % final_filter_1)
    sys.stdout.flush()
    os.system(final_filter_1)
    sys.stdout.write("%s\n" % final_filter_2)
    sys.stdout.flush()
    os.system(final_filter_2)
    sys.stdout.write("%s\n" % final_filter_3)
    sys.stdout.flush()
    os.system(final_filter_3)

    # get chipSampleRepAll.conservative.narrowPeak file using IDR results
    conservative_bed_idr = ("%s/%s.idr.conservative.%sPeak"
                            % (working_dir, expName, peak_type))
    final_filter_1 = ("""sort -k8nr,8nr %s > %s.tmp""" % (pool_macs,
                                                          conservative_bed_idr))
    final_filter_2 = ("""head -n %i %s.tmp| %sbedtools sort -i > %s"""
                      % (conservative_threshold, conservative_bed_idr,
                         bedtools_path, conservative_bed_idr))
    final_filter_3 = ("rm %s.tmp " % conservative_bed_idr)
    sys.stdout.write("%s\n" % final_filter_1)
    sys.stdout.flush()
    os.system(final_filter_1)
    sys.stdout.write("%s\n" % final_filter_2)
    sys.stdout.flush()
    os.system(final_filter_2)
    sys.stdout.write("%s\n" % final_filter_3)
    sys.stdout.flush()
    os.system(final_filter_3)
