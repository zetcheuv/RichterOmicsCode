#!/bin/bash
#@seb2016 sebastien.hergalant@inserm.fr
# 2) splicing-aware alignement step
#samples must each be placed in a dedicated subdirectory of 0-seq/, in fastq.gz format
#main dependencies: hisat2 picard samtools merge_logs.R hg38_reference
#parallel computation for multi-sample handling can be tweaked with nb_threads

if [[ "$#" -lt 1 ]] || [[ ! -d $1 ]]; then
        echo "usage: $(basename $0) working_directory"
        echo "working_dir must exist and should have already been processed by 0-/1-rnaseq...bash"
        exit 1
fi

# some changeable vars
export align_options="--score-min L,0,-0.20 --sp 10,3 --no-unal --dta" #for hisat2
export nb_threads=15
export ref_genome="hg38snptr"
export MD="MarkDuplicates"

# programs and files
export bowtie2="/usr/local/bin/hisat2"
export HISAT2_INDEXES="/home/chac/ngs/ref_librairies/"
export picard="/usr/local/bin/picard.jar"
export samtools="/usr/local/bin/samtools"
export plot_stats="/usr/local/bin/plot-bamstats"
#find merge_logs.R in the same dir as $0
merge_logs="$(dirname $(readlink -f $0))/2-merge_logs.R"
# absolute path with link resolution for starting_dir
export starting_dir=$(readlink -f $1)
export seq_dir="${starting_dir}/0-seq"
export align_dir="${starting_dir}/2-alignments"

# create working dir and make sure we are in
mkdir -p ${align_dir}
cd ${starting_dir}
# even with -p 20, we can do 20 alignments together, no sweat
for sample in $(ls -1 ${seq_dir}/ |sort -R) ; do echo $sample ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   sample={}
   [[ ! -d ${seq_dir}/${sample} ]] && exit 0 ;
   cd ${seq_dir}/${sample} ;

   # there are only links to seqfiles of interest in each dir, no need for complicated pattern check anymore
   R1=""; for r1 in *_R1_*.fastq.gz; do R1="$R1$r1,"; done ; R1=${R1:0:${#R1}-1} ; # there is a comma left at the end
   
   # align R1 on ref genome
   ${bowtie2} -p ${nb_threads} --mm ${align_options} -x ${ref_genome} -U ${R1} -S ${align_dir}/${sample}.sam &> ${align_dir}/${sample}.mapper.log ;
}'

# compute global alignment rates for all samples
${merge_logs} "${align_dir}/.*\\.mapper\\.log" > ${align_dir}/all_samples.mapper.merged.log
# move to align_dir
cd ${align_dir}
# compress, sort, index, stats on alignment files / careful with the use of -@ samtools sort option with xargs !!!
for sam in $(ls -1 *.sam |sort -R) ; do echo $sam ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   sam={}
   # no metrics -> M=/dev/null
   ${samtools} sort -m 2G -@ 3 -O bam $sam -T ${sam%.*} -o ${sam%.*}.so.bam &&
   rm -f $sam &&
   java -Xmx4G -jar ${picard} $MD I=${sam%.*}.so.bam O=${sam%.*}.somd.bam M=/dev/null &> ${sam%.*}.picard.log &&
   rm -f ${sam%.*}.so.bam &&
   ${samtools} index ${sam%.*}.somd.bam &&
   ${samtools} flagstat ${sam%.*}.somd.bam > ${sam%.*}.stats.log &&
   ${samtools} stats ${sam%.*}.somd.bam > ${sam%.*}.somd.stats &&
   ${plot_stats} -p ${sam%.*}_plots/ ${sam%.*}.somd.stats ;
}'

# finally, compute global stats for all bams
${merge_logs} ".*\\.stats\\.log" > all_samples.stats.merged.log
for somd in *.somd.stats ; do grep "^SN" $somd |cut -f 2- > ${somd}.sn ; done
${merge_logs} ".*\\.somd\\.stats\\.sn" > all_samples.somd.merged.log
rm -fr *.stats.sn

