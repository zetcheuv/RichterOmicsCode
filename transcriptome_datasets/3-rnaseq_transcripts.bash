#!/bin/bash
#@seb2016 sebastien.hergalant@inserm.fr
# 3) transcript reconstruction step
#main dependencies: stringtie hg38_gtf_file

if [[ "$#" -lt 1 ]] || [[ ! -d $1 ]]; then
        echo "usage: $(basename $0) working_directory"
        echo "working_dir must exist and should contain 2-alignments"
	     echo "subdirectory with sorted bams inside"
        exit 1
fi

# some vars
nb_threads=30
ref_gtf="/home/chac/ngs/ref_sequences/hg38.gtf"
transcript_label="trID"

# programs and files
stringtie="/usr/local/bin/stringtie"
stringtie2counts="/usr/local/bin/stringtie2counts.py"
sampledir_pattern="_stats"
# absolute path with link resolution for starting_dir
starting_dir=$(readlink -f $1)
align_dir="${starting_dir}/2-alignments"
transcript_dir="${starting_dir}/3-transcripts"
merged_gtf="${transcript_dir}/all-samples.merged.gtf"

# create working dir
mkdir ${transcript_dir} &> /dev/null
# and make sure we are in the right dir
cd ${transcript_dir}

# compute isoform reconstruction from read alignments
for bam in ${align_dir}/*.somd.bam
do
   # -f 0<>1: min isoform abundance
   # -j float: min number of spliced reads to allow splice junction
   # -c float: min read coverage
   # -M 0<>1: max fraction of multiple match reads at one place
   ${stringtie} $bam -o $(basename ${bam%.*.*}).gtf -G ${ref_gtf} -p ${nb_threads} -f 0.2 -j 3 -c 10 -M 0.5 -l ${transcript_label}
done
# merge process
# -m minimum input transcript length
${stringtie} --merge *.gtf -G ${ref_gtf} -o ${merged_gtf} -l ${transcript_label} -m 200
# final abundance estimation for each sample
for bam in ../2-alignments/*.somd.bam
do
   sample="$(basename ${bam%.*.*})"
   ${stringtie} $bam -B -e -o "${sample}${sampledir_pattern}/${sample}.gtf" -G ${merged_gtf} -p ${nb_threads} -A "${sample}.abundance" -C "${sample}.cover" -f 0.2 -j 3 -c 10 -M 0.5
done
# transform into raw read counts for genes and transcripts
# -l integer: read length is needed here
avlength=$(grep "average length:" ${align_dir}/all_samples.somd.merged.log |cut -f2)
${stringtie2counts} -i . -l ${avlength%.*} -p ".*${sampledir_pattern}" -s ${transcript_label}

# gene and transcript matrix annotations
mkdir all-samples_counts &> /dev/null
mv gene_count_matrix.csv transcript_count_matrix.csv all-samples_counts/
cd all-samples_counts/
awk 'BEGIN{FS=OFS="\t"}{if($3=="transcript") print $9}' ${merged_gtf} > transcripts-symbols.gtf
sed -n 's/^gene_id "\(.*\)"; transcript_id ".*"; gene_name "\(.*\)"; ref_gene_id "\(.*\)";.*$/\1\t\2\/\/\3\/\/\1/p' transcripts-symbols.gtf |awk -F '\t' '{if(substr($1,1,4)=="ENSG"){ print substr($0,1,length($0)-17) } else print }' |sort -u > ID2symbols.txt
awk -F'\t' 'NF>1 { if(a[$1]!="") {a[$1] = a[$1]"|"$2} else {a[$1]=$2} } END {for(i in a) {print i"\t"a[i]} }' ID2symbols.txt |awk -F'\t' '{ split($2,a,"|") ; if(length(a) == 1){ print } else { for(i=1;i<=length(a);i++){ if(match(a[i],/^[A-Za-z0-9]+\/\/ENSG[0-9]+\/\/trID\.[0-9]+$/) && !match(a[i],/^(MIR[0-9]+)|(SNOR[A-D][0-9]+)|(sno)|(KIAA[0-9]+)|(RNA5S)|(LINC[0-9]+).*\/\/ENSG[0-9]+\/\/trID\.[0-9]+$/)){ b=a[1]; a[1]=a[i]; a[i]=b; b=a[1] ; for(j in a){ if(j!=1) b=b"|"a[j] ; } print $1"\t"b ; break; } else if(length(a)==i) print ; } } }' > TR2IDs.txt
sed -n 's/^gene_id "\(.*\)"; transcript_id "\(.*\)"; gene_name "\(.*\)"; ref_gene_id "\(.*\)";.*$/\2\t\3\/\/\2\/\/\1/p' transcripts-symbols.gtf > ENST2IDs.txt
rm -f transcripts-symbols.gtf ID2symbols.txt

