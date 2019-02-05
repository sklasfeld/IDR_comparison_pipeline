#! /bin/bash

# This Script was created by Samantha Klasfeld
# January 30, 2019

# Example Script: countDuplicates.sh bam_file peak_file out_file

### Warning: no suffixes needed
if [ $# -ne 4 ]; then
   echo 'Example Script: countDuplicates.sh <bam_file> <peak_file> <dup_out_file> <all_out_file>'
   return
fi

READS=${1%.bam} # bam with duplicates annotated already
PEAKS=${2%.narrowPeak} # narrowPeak file includes peak regions
OUT_COVERAGE1=$3 #contains coverage of duplicate reads
OUT_COVERAGE2=$4 #contains coverage of all reads
NTHREADS=2 ## number of threads


echo '# obtain duplicate reads...'
if [ -f "${READS}.dups.bam" ];
then
   echo "# File ${READS}.dups.bam exist."
else
   samtools view -@ ${NTHREADS} -f 1024 -F 772 -b ${READS}.bam -o ${READS}.dups.bam
fi

echo '# index duplicate reads...'
if [ -f "${READS}.dups.bam.bai" ];
then
   echo "# File ${READS}.dups.bam.bai exist."
else
   samtools index ${READS}.dups.bam
fi

echo '# count total duplicate reads...'
TOTAL_DUPS=`samtools view -c ${READS}.dups.bam`

echo "NUMBER OF DUPLICATES IN TOTAL: ${TOTAL_DUPS}"

echo '# running bedtools coverage...'
bedtools coverage -a ${PEAKS}.narrowPeak -b ${READS}.dups.bam > ${OUT_COVERAGE1}
echo '# running bedtools coverage...'
bedtools coverage -a ${PEAKS}.narrowPeak -b ${READS}.bam > ${OUT_COVERAGE2}

echo '# calculating duplicate coverage in peaks...'
DUP_COUNT=`awk '{x=x+$(NF-3)}END{print x}' ${OUT_COVERAGE1}`
READ_COUNT=`awk '{x=x+$(NF-3)}END{print x}' ${OUT_COVERAGE2}`
PEAK_AREA=`awk '{x=x+$(NF-1)}END{print x}' ${OUT_COVERAGE}`
echo "NUMBER OF DUPLICATES IN PEAKS: ${DUP_COUNT}"
echo "NUMBER OF READS IN PEAKS: ${DUP_COUNT}"
echo "NUMBER OF BP IN PEAKS: ${PEAK_AREA}"