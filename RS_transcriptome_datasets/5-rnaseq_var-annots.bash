#!/bin/bash
#@seb2016 sebastien.hergalant@inserm.fr
# 5) variant annotation step
#main dependencies: bcftools annovar hg38_reference annovar_annotation_databases

if [[ "$#" -lt 1 ]] || [[ ! -d $1 ]]; then
        echo "usage: $(basename $0) working_directory"
        echo "working_dir must exist and should contain 4-variants"
	echo "subdirectectory all_samples bcf file(s) inside"
        exit 1
fi

#some vars
export ref_genome="hg38"
export nb_threads=20 # max 30
export var_protocol=refGene,knownGene,cytoband,avsnp147,cosmic70,clinvar_20160302,dbscsnv11,dbnsfp30a
export var_operation=g,g,r,f,f,f,f,f
export var_argument='-neargene 1000 -exonicsplicing -splicing_threshold 10','-neargene 1000 -exonicsplicing -splicing_threshold 10',,,,,,

# programs and files
export bcftools="/usr/local/bin/bcftools"
export annovar="/usr/local/bin/table_annovar.pl"
export ref_annots="/home/chac/ngs/ref_annots/"
# absolute path with link resolution for starting_dir
starting_dir=$(readlink -f $1)
variant_dir="${starting_dir}/4-variants"
annots_dir="${starting_dir}/5-var_annots"
# variants file is the same as in 4-exome_variants.bash
variants_final="${variant_dir}/variants_final.bcf"
var_table="${annots_dir}/variants.table"
snips="${annots_dir}/snips"
indels="${annots_dir}/indels"


# create working dir
mkdir ${annots_dir} &> /dev/null
# and make sure we are in the right directory
cd ${annots_dir}
# extract all variants with sample genotypes
# -H for column header and sed for renaming samples and removing unwanted column infos from the header
${bcftools} query -H -f '%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%TGT]\n' ${variants_final} |sed "1s/\(# \)\|\(\[[0-9]\+\]\)\|\(\.somd\.bam\:GT\)//g" > ${var_table}
# separate SNPs from INDELs
head -1 ${var_table} |cut --complement -f3 |tee ${snips}.table ${indels}.table > /dev/null
awk 'BEGIN {FS=OFS="\t"}{if($3=="SNP") print}' ${var_table} |cut --complement -f3 >> ${snips}.table
awk 'BEGIN {FS=OFS="\t"}{if($3=="INDEL") print}' ${var_table} |cut --complement -f3 >> ${indels}.table

# annotate with annovar // beware: hg38 is not available nor it is "buildable" for ensGene and snp138
${bcftools} view -v snps ${variants_final} > ${snips}.vcf
${bcftools} view -v indels ${variants_final} > ${indels}.vcf
${annovar} ${snips}.vcf -vcfinput -dot2underline ${ref_annots} -buildver ${ref_genome} -outfile $(basename ${snips}) -remove -protocol "${var_protocol}" -operation "${var_operation}" -argument "${var_argument}" -nastring .
${annovar} ${indels}.vcf -vcfinput -dot2underline ${ref_annots} -buildver ${ref_genome} -outfile $(basename ${indels}) -remove -protocol "${var_protocol}" -operation "${var_operation}" -argument "${var_argument}" -nastring .
rm -f ${snips}.avinput ${indels}.avinput
mv -f ${snips}.${ref_genome}_multianno.vcf ${snips}.vcf
mv -f ${indels}.${ref_genome}_multianno.vcf ${indels}.vcf

# add index to column to all resulting files (they are already sorted)
awk 'BEGIN {FS=OFS="\t"}{if(NR==1) {i=1; printf "index\t"} else {printf "SNP%i\t",i++} ;print}' ${snips}.table > ${snips}.gt
awk 'BEGIN {FS=OFS="\t"}{if(NR==1) {i=1; printf "index\t"} else {printf "INDEL%i\t",i++} ;print}' ${indels}.table > ${indels}.gt
awk 'BEGIN {FS=OFS="\t"}{if(NR==1) {i=1; printf "index\t"} else {printf "SNP%i\t",i++} ;print}' ${snips}.${ref_genome}_multianno.txt > ${snips}.annots
awk 'BEGIN {FS=OFS="\t"}{if(NR==1) {i=1; printf "index\t"} else {printf "INDEL%i\t",i++} ;print}' ${indels}.${ref_genome}_multianno.txt > ${indels}.annots
rm -f ${snips}.${ref_genome}_multianno.txt ${indels}.${ref_genome}_multianno.txt ${snips}.table ${indels}.table
# indexing of the vcf files / cut \x3b \x3d \x2c weird chars from input
(${bcftools} view -h ${snips}.vcf |tac |sed "1s/\(\.somd\.bam\)//g" |tac; paste <(${bcftools} view -H ${snips}.vcf |cut -f1-2) <(awk 'BEGIN {FS=OFS="\t"}{if(NR!=1) print $1}' ${snips}.annots) <(${bcftools} view -H ${snips}.vcf |cut -f4- |sed 's/\(\\x3b\)\|\(\\x3d\)\|\(\\x2c\)/\//g')) > ${snips}_annots.vcf
(${bcftools} view -h ${indels}.vcf |tac |sed "1s/\(\.somd\.bam\)//g" |tac; paste <(${bcftools} view -H ${indels}.vcf |cut -f1-2) <(awk 'BEGIN {FS=OFS="\t"}{if(NR!=1) print $1}' ${indels}.annots) <(${bcftools} view -H ${indels}.vcf |cut -f4- |sed 's/\(\\x3b\)\|\(\\x3d\)\|\(\\x2c\)/\//g')) > ${indels}_annots.vcf
mv -f ${snips}_annots.vcf ${snips}.vcf
mv -f ${indels}_annots.vcf ${indels}.vcf

# now annotate all bcf in variant_dir and be done with it
cd ${variant_dir}
for bcf in $(ls -1 *.bcf |sort -R) ; do echo $bcf ; done |xargs -I{} --max-procs ${nb_threads} bash -c '{
   bcf={}
   [[ $(${bcftools} view -h $bcf |egrep "ANNOVAR_DATE" |wc -l) -ne 0 ]] && exit 0 ; #do not reannotate!!
   ${bcftools} view $bcf > annotating.${bcf%.*}.vcf ;
   ${annovar} annotating.${bcf%.*}.vcf -vcfinput -dot2underline ${ref_annots} -buildver ${ref_genome} -outfile annotating.${bcf%.*} -remove -protocol "${var_protocol}" -operation "${var_operation}" -argument "${var_argument}" -nastring . ;
   (${bcftools} view -h annotating.${bcf%.*}.${ref_genome}_multianno.vcf |tac |sed "1s/\(\.somd\.bam\)//g" |tac; ${bcftools} view -H annotating.${bcf%.*}.${ref_genome}_multianno.vcf |sed '\''s/\(\\x3b\)\|\(\\x3d\)\|\(\\x2c\)/\//g'\'') |${bcftools} view -O b -o $bcf - ;
   rm -f annotating.${bcf%.*}.* ;
   ${bcftools} index $bcf ;
}'

