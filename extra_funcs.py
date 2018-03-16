#!/usr/bin/env python

# author: Samantha Klasfeld
# date: 10/2/2017

import os
import sys
import subprocess
import math



# Global Definitions

# check if value is between 0 and 1
def between0and1(x):
    if x<=1 and x>=0:
        return True
    else:
        return False

def convert_delim(infile, outfile, old_delim, new_delim):
    with open(infile) as infile:
        out = open(outfile, 'w')
        for line in infile:
            fields = line.split(old_delim)
            out.write(new_delim.join(fields))
        out.close()


def file_len(fname):
    num_lines = sum(1 for line in open(fname))
    return num_lines


def runcmdline(cmd):
    sys.stdout.write("%s\n" % cmd)
    sys.stdout.flush()
    os.system(cmd)


def shufNsplitBam(bamfile, nlines, outputdir, output_stub, samtools_path=""):
    out = ("%s/%s" % (outputdir, output_stub))
    # get header of files
    header_cmd = ("%ssamtools view -H %s > %s.header" % (samtools_path, bamfile, out))
    # get lines of sam file without header
    shuf_sam1 = ("%ssamtools view %s > %s_preshuf.sam" % (samtools_path, bamfile, out))
    # shuffle the lines of sam file
    shuf_sam2 = ("shuf %s_preshuf.sam -o %s_shuf.sam" % (out, out))
    # split the lines of the same file into two files
    shuf_sam_split = ("split --numeric-suffixes=0 -l %i %s_shuf.sam %s" % (nlines, out, out))

    # add headers to split files to turn them into sam files
    get_pr1_1 = ("cat %s.header %s00 > %s00.sam" % (out, out, out))
    get_pr2_1 = ("cat %s.header %s01 > %s01.sam" % (out, out, out))
    # convert sam files to bam files
    get_pr1_2 = ("%ssamtools view -b %s00.sam -o %s00.bam" % (samtools_path, out, out))
    get_pr2_2 = ("%ssamtools view -b %s01.sam -o %s01.bam" % (samtools_path, out, out))
    # sort the bam files
    get_pr1_3 = ("%ssamtools sort %s00.bam -o %s.pr1.tagAlign.bam" % (samtools_path, out, out))
    get_pr2_3 = ("%ssamtools sort %s01.bam -o %s.pr2.tagAlign.bam" % (samtools_path, out, out))
    # remove all files that are not the final output
    rm_oldfiles = ("rm %s.header %s_preshuf.sam %s_shuf.sam %s00 %s01  %s00.sam %s01.sam %s00.bam %s01.bam "\
                   % (out, out, out, out, out ,out, out, out, out))
    runcmdline(header_cmd)
    runcmdline(shuf_sam1)
    runcmdline(shuf_sam2)
    runcmdline(shuf_sam_split)
    runcmdline(get_pr1_1)
    runcmdline(get_pr2_1)
    runcmdline(get_pr1_2)
    runcmdline(get_pr2_2)
    runcmdline(get_pr1_3)
    runcmdline(get_pr2_3)
    runcmdline(rm_oldfiles)



def split_f(file_name, num, samtools_path=""):
    nlines = wc(file_name, samtools_path)
    splitby = int(math.ceil((nlines + 1) / num))
    return splitby


def wc(file_name, samtools_path=""):
    if (not os.path.isfile(file_name)):
        sys.stdout.write("Count: 0\n")
        sys.stdout.flush()
        return 0
    else:
        cmd = ("wc -l %s" % file_name)
        if os.path.splitext(file_name)[1] == ".gz":
            cmd = ("gunzip -c %s | wc -l" % file_name)
        elif os.path.splitext(file_name)[1] == ".bam" or \
            os.path.splitext(file_name)[1] == ".sam" or \
            os.path.splitext(file_name)[1] == ".tagAlign":
            cmd = ("%ssamtools view -c %s" % (samtools_path, file_name))
        sys.stdout.write("%s\n" % cmd)
        sys.stdout.flush()
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        #print(output)
        break_out = output.split()
        sys.stdout.write("Count: %i\n" % int(break_out[0]))
        sys.stdout.flush()
        return int(break_out[0])

def empty(file_name, samtools_path=""):
    if wc(file_name, samtools_path) > 0:
            return False
    else:
            return True
