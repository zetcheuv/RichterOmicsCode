#!/bin/bash
#@seb2016 sebastien.hergalant@inserm.fr
# 4) variant calling from RNAseq data
#main dependencies: bcftools samtools bedtools hg38_references

if [[ "$#" -lt 1 ]] || [[ ! -d $1 ]]; then
        echo "usage: $(basename $0) working_directory"
        echo "working_dir must exist and should contain 2-alignments subdirectectory with sorted bams inside"
        exit 1
fi

# some changeable vars
ref_genome="hg38"
export min_depth=10 # minimum base coverage for variant calling
export min_qual=10 # minimum reported variant quality 
export isec_num=2 # minimum number of samples presenting with the same variants
export nb_threads=20 # max 30
export mutation_rate="1e-3" # default 1e-3

# programs and files
if [ ${ref_genome} == "hg19" ]; then export ploidy="GRCh37"; fi
if [ ${ref_genome} == "hg38" ]; then export ploidy="GRCh38"; fi
export ref_fasta="/home/chac/ngs/ref_sequences/${ref_genome}.fa"
export bcftools="/usr/local/bin/bcftools"
export samtools="/usr/local/bin/samtools"
export plot_stats="/usr/local/bin/plot-vcfstats"
bedtools="/usr/local/bin/bedtools"
# absolute path with link resolution for starting_dir
starting_dir=$(readlink -f $1)
align_dir="${starting_dir}/2-alignments"
coverage_dir="${starting_dir}/3-coverage"
export variant_dir="${starting_dir}/4-variants"
export isec_region="${variant_dir}/isec${isec_num}_regions.txt"
# variants file is the same as in 5-rnaseq_var-annots.bash
variants_final="${variant_dir}/variants_final.bcf"

# create working dir
mkdir ${variant_dir} &> /dev/null
# and make sure we are in the right dir
cd ${align_dir}
# variant calling
   # [ -B is needed in mpileup for calling long indels with bwa ]
   # -A keeps orphan reads
   # -C 50 adjusts mapping quality (recommended)
   # -q with high value is valid here
   # -Q 30 keeps only bases with good sequencing quality 
   # -m 2 a gap must occur at least twice to be considered significant
   # -v in bcftools is needed for variant calling (otherwise the call is made on EVERY positions)
   # -m in bcftools call for multiallelic variants (useful for cancer clones but maybe not so much for gg) / -c for consensus call
   # filter read depth at base position: DP<10
for bam in $(ls -1 *.somd.bam |sort -R) ; do echo $bam ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   bam={}
   ${samtools} mpileup -uf ${ref_fasta} -t "DP,SP,AD,ADF,ADR" -A -C 50 -q 10 -Q 30 -m 2 $bam |${bcftools} call --ploidy ${ploidy} -v -m -P ${mutation_rate} - |${bcftools} filter -O b -o ${variant_dir}/${bam%%.*}.bcf -e "DP<${min_depth} || QUAL<${min_qual}" - ;
   ${bcftools} index ${variant_dir}/${bam%%.*}.bcf ;
}'

# intersecting all variants present at least in 3 samples from the whole
cd ${variant_dir}
${bcftools} isec -n +${isec_num} -o ${isec_region} *.bcf 
# now recompute all polymorphisms for each sample for every (no -v in bcftools call) position
# in the intersection file so as to keep track of reference genotypes instead of missing ones
# filter by read depth and qual
cd ${align_dir}
for bam in $(ls -1 *.somd.bam |sort -R) ; do echo $bam ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   bam={}
   ${samtools} mpileup -uf ${ref_fasta} -t "DP,SP,AD,ADF,ADR" -A -C 50 -l ${isec_region} -q 10 -Q 30 -m 2 $bam |${bcftools} call --ploidy ${ploidy} -m -P ${mutation_rate} - |${bcftools} filter -O b -o ${variant_dir}/${bam%%.*}_isec${isec_num}.bcf -e "DP<5" - ;
   ${bcftools} index ${variant_dir}/${bam%%.*}_isec${isec_num}.bcf ;
}'
# merge all samples into one single variant file for SNPs and INDELs only (weird results with -o and wildcard *.bcf so redirection only)
# normalise indels to reference and rejoin all sites then split multiallelic sites (step by step due to a bug in bcftools norm)
# it is ok by now to change all resulting ./1 to 0/1 genotypes
# [ bcftools1.3 merge: resulting QUAL is always the max QUAL from all samples, non negociable at the moment ]
cd ${variant_dir}
${bcftools} merge -m both -O b *_isec${isec_num}.bcf |${bcftools} norm -m +both -f ${ref_fasta} -c s -O b - |${bcftools} norm -m -both -f ${ref_fasta} - |sed 's/\t\.\/1:/\t0\/1:/g' |${bcftools} +fill-tags -O b -o ${variants_final}
${bcftools} index ${variants_final}
# remove refiltered single/merged sample files and indexes
rm -f *_isec${isec_num}.bcf*

# stats for all samples
for bcf in $(ls -1 *.bcf |sort -R) ; do echo $bcf ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   bcf={}
   ${bcftools} stats $bcf > ${bcf%%.*}.stats ;
   ${plot_stats} -p ${bcf%%.*}_plots/ ${bcf%%.*}.stats ;
}'

