#!/bin/bash
#@seb2016 sebastien.hergalant@inserm.fr
# 1) QC step
#samples must each be placed in a dedicated subdirectory of 0-seq/, in fastq.gz format
#main dependencies: fastqc qc_utils.R 

# casava option detects automatically different samples and separates R1/R2 files but put all the rest together
# when files are not from casava output (ex: F3/F5-DNA) -> complains to stderr but works anyway
# also works with colourspace -> conversions are made on the fly and do not seem to affect stats strongly

if [[ "$#" -lt 1 ]] || [[ ! -d $1 ]]; then
        echo "usage: $(basename $0) working_directory"
        echo "working_dir must exist and should have already been scanned (see 0-..._seqfind.bash)"
        exit 1
fi

# some changeable vars
nb_threads=30

# programs and files
fastqc="/usr/local/bin/fastqc"
# find utils_qc.R in the same dir as $0
utils_qc="$(dirname $(readlink -f $0))/1-qc_utils.R"
# absolute path with link resolution for starting_dir
starting_dir=$(readlink -f $1)
seq_dir="${starting_dir}/0-seq"
qc_dir="${starting_dir}/1-qc"
report_file="${qc_dir}/all_samples-report.txt"

# create working dir
mkdir ${qc_dir} &> /dev/null
# and make sure we are in the right dir
cd ${starting_dir}

# fastq.gz files are found in subfolders (one for each sample) in seq_dir
$fastqc -o ${qc_dir} -q -t ${nb_threads} --casava ${seq_dir}/*/*.fastq.gz

# some useful general reports...
cd ${qc_dir}
echo -e "QC modules failed:\n------------------" > ${report_file}
for file in *.zip ;do unzip -p $file */summary.txt |grep ^FAIL ; done >> ${report_file}
echo -e "\n###############################\nSequence info:\n--------------" >> ${report_file}
for file in *.zip ; do unzip -p $file */fastqc_data.txt |grep -E "Filename|Sequences" ; done >> ${report_file}
# ...and plots
${utils_qc} ${qc_dir}
