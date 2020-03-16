#!/bin/sh
#
#  Reserve 1 CPUs for this job
#$ -pe parallel 1
#  Request 10G of RAM
#$ -l h_vmem=10G
#
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
# -N $name
#
#  The path used for the standard output stream of the job
# -o /ebio/abt6_projects7/bacterial_strain_analysis/data/processed_DATA
#
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
# -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
#  Run job from current working directory
#$ -cwd
#
#  Send email when the job begins, ends, aborts, or is suspended
##$ -m beas

merged=$1
usearch=/ebio/abt6_projects7/bacterial_strain_analysis/code/usearch10

#515F_GI502
rm "$merged"_515F_GI502F.fastq
#515
grep "^GAGTG[CT]CAGC[AC]GCCGCGGTAA[ATCG]*" $merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 21 -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_515F_GI502F.fastq
rm "$merged"_temp*

#GI
grep "^GTAAAGATAAATGGGTCATCTAA[ATCG]*" $merged -B 1 -A 2 --no-group-separator > "$merged"_temp.fq
$usearch -fastx_truncate "$merged"_temp.fq -stripleft 23 -fastqout "$merged"_temp.fq.stripped
cat "$merged"_temp.fq.stripped >> "$merged"_515F_GI502F.fastq
rm "$merged"_temp*
