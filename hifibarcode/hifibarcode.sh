#!/usr/bin/env bash

if [[ $# -ne 7 ]]; then
        echo "usage bash hifibarcode.sh <prefix> <fastq.1> <fastq.2> <primer.mapFile> <indexFile> <cluster minimum size |5> <cpu number>"
        exit 0
fi
pre="$1"
fqf="$2"
fqs="$3"
primer="$4"
idx="$5"
min="$6"
cpu="$7"
dir=$(dirname "$0")

echo "set minimum depth of end sequences of ${min} and cpu number of ${cpu}\n"
echo "set primer map file as ${primer} and index file as ${idx}"
echo "\n\nmake sure the above settings are correct"

#1 demultiplex the input fastq
perl ${dir}/bin/1_demultiplex.pl -fq1 ${fqf} -fq2 ${fqs} -pri ${primer} -ind ${idx} -out ${pre}
echo "done for demultiplex the input fastq"

#2 unique and cluster for each file
perl ${dir}/bin/2_uniqu_sort_cluster.Pro.pl ${pre}
echo "done for unique and cluster for each file"

#3 overlap for end files
perl ${dir}/bin/3_sep_extract_overlap.pl ${pre} ${fqf} ${fqs}
echo "done for overlap for end files"

#4 cluster for overlapped end file
perl ${dir}/bin/4_cluster_fromend.pl ${pre} {$min}
echo "done for cluster for overlapped end file"

#5 assemble
perl ${dir}/bin/5_forgap_filling.pl ${pre} ${cpu}
echo "done for assembly"

# rename output
perl ${dir}/bin/6_rename_kmer.pl ${pre} 658 1
echo "finished"


