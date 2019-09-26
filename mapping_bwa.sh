#!/bin/bash


#Description
#This script map single or paired-end reads using BWA MEM algorithm.

usage="$(basename "$0") [-h] -1 file1.fastq [-2 file2.fastq] -i path/to/bwa/index/prefix_bwaidx -o /path/to/output/dir \n
This script uses BWA to map single or paired-end data reads using the MEM algorithm.\n
\n
Argument:\n
    -h  show this help text\n
    -1  first fastq file\n
    -2  second fastq file (ptional)\n
    -i  region to map in case of local mapping (either '45S' or 'chr3')\n
    -o  output directory\n
\n
Author:\n
    written by Johan Zicola (johan.zicola@gmail.com)\n
    2019-09-25\n
"

# Definition
# h: help display
# 1: fastq file input
# 2: second fastq file in case the data are paired-end
# i: path bwa index file (should contain the prefix of the index)
# o: name output directory

while getopts "h?:1:2:i:o:" arg; do
	case $arg in
		h|\?) 
			echo -e $usage | fold -s -w 80
			exit 0
			;;
		1) 
			if [ ! -e $OPTARG ]; then
				echo "Fastq file incorrect or not existing"
				exit 1
			elif [ -d $OPTARG ]; then
				echo "Fastq file argument is a directory. Provide a fastq file instead"
				exit 1
			else
				fastq=$OPTARG
				only_fastq=$(basename $fastq)
				path_fastq=$(dirname $fastq)
				name_file=$(echo $only_fastq | cut -d'.' -f1)
			fi
			;;			
		2)
			if [ ! -e $OPTARG ]; then
				echo "Fastq file #2 mates incorrect or not existing"
				exit 1
			elif [ -d $OPTARG ]; then
				echo "Fastq file #2 mates argument is a directory. Provide a fastq file instead"
				exit 1
			else
				fastq2=$OPTARG
			fi
			;;
		i)
			index=$OPTARG
			# Check if directory containing bwa index exists
			path_index=$(dirname $OPTARG)
			# Check if the prefix_bwaidx.amb file exists
			index_name=${OPTARG}.amb
			
			if [ ! -d $OPTARG ]; then
				echo "Path to bwa index file incorrect or not existing"
				exit 1
			elif [ ! -e $index_name ]; then
				echo "bwa index file $index_name not existing"
				exit 1
			fi
			;;
		o)
			if [ ! -d $OPTARG ]; then
				echo "Path to output directory incorrect or not existing"
				exit 1
			else
				# Remove final slash of the path if any
				path_output=${OPTARG%/}
			fi
			;;
	esac
done



# Return help message if no arguments given
if [ $OPTIND -eq 1 ]; then 
	echo -e "\n!!!!!!!!!!!!No options were passed!!!!!!!!!!!!!!!!\n"
	echo -e $usage | fold -s -w 80
	exit 0
fi

# Set output directory. Take path of fastq file if -o not 
# provided by user
if [ -z $path_output ]; then
	# Remove final slash if any
	path_output=${path_fastq%/}
fi

# Display arguments defined by user
echo "Output directory: $path_output"

echo "Fastq file: $fastq"

if [ ! -z ${fastq2} ]; then
	echo "Fastq file 2: $fastq2"
else
	echo "Fastq file 2: No 2nd fastq file provided, data are assumed to be single end"
fi

echo "Path to bwa index file: $path_index"



# Check installed softwares
for i in {bwa,samtools,zgrep}; do
	command -v $i >/dev/null 2>&1 || { echo >&2 "$i is not installed.  Aborting."; exit 1; }
done



#Create a function giving the date (prefer that to alias DATE="date +%Y-%m-%d_%H:%M:%S"
# as alias will require shopt -s expand_aliases to be able to be displayed when called $(DATE)
# function does not need any changed in shell options
give_date(){
    date=$(date +%Y-%m-%d_%H:%M:%S)
    echo $date
}


mappingBWA(){

    #Map using the MEM algorithm of BWA, keep only primary mapped reads and index them (use 4 threads)
    echo "Mapping at $(give_date)" >> ${path_output}/log.txt
   	
	# Use paired-end mode if second fastq file was provided
	if [ ! -z ${fastq2} ]; then 
		#Map, convert sam to bam, remove unmapped reads (-F 4), sort the bam file according to read position
		echo -e "bwa mem -M -t 4 $index ${fastq_file1} ${fastq_file2} | samtools view -bS -F 4 - | samtools sort - -o ${path_output}/${name_file}.sorted.bam\t" >> ${path_output}/log.txt
		bwa mem -M -t 4 $index ${fastq_file1} ${fastq_file2} | samtools view -bS -F 4 - | samtools sort - -o ${path_output}/${name_file}.sorted.bam
    # If only 1 fastq file provided, use SE mode
	else
		echo -e "bwa mem -M -t 4 $index ${fastq_file1} | samtools view -bS -F 4 - | samtools sort - -o ${path_output}/${name_file}.sorted.bam\t" >> ${path_output}/log.txt
		bwa mem -M -t 4 $index ${fastq_file1} | samtools view -bS -F 4 - | samtools sort - -o ${path_output}/${name_file}.sorted.bam
}

mappingBWA "$@"
