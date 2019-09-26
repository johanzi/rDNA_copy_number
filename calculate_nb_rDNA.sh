#!/usr/bin/env bash 

# This script takes as -1 file with number of read count for 18S and
# -2 file with read count for chr3 10Mb. It returns an estimated value of the
# the number of 45S rDNA gene copies


# Check if python3 is installed
command -v python3 >/dev/null 2>&1 || { echo >&2 "python3 is not installed.  Aborting."; exit 1; }

while getopts ":1:2:" arg; do
    case $arg in
        1)
            if [ ! -e $OPTARG ]; then
                echo "File 1 name incorrect or not existing"
                exit 1
            elif [ -d $OPTARG ]; then
                echo "File 1 name is a directory. Provide a file name"
                exit 1
            else
				file_18S=$OPTARG
            fi
            ;;
        2)
			if [ ! -e $OPTARG ]; then
                echo "File 2 name incorrect or not existing"
                exit 1
            elif [ -d $OPTARG ]; then
                echo "File 2 name is a directory. Provide a file name"
                exit 1
            else
            	file_chr3=$OPTARG
            fi
            ;;
   esac
done

# Return help message if no arguments given. 4 is the total number of
# arguments which should be provided (flags included)
if [ ! $# -eq 4 ]; then
    echo -e "\n!!!!!!!!!!!!Wrong number of argument!!!!!!!!!!!!!!!!\n"
	echo -e "\nExample: bash $(basename "$0") -1 file_count18S.txt -2 file_count_chr3.txt\n"
    exit 1 
fi

count_18S=$(cat $file_18S)
count_chr3=$(cat $file_chr3)

# Normalize size based on length of each region (in bp)
normalized_18S_count=$(python3 -c "print($count_18S/1808)")
normalized_chr3_count=$(python3 -c "print($count_chr3/1E7)")

# Divide the normalized values
nb_copies=$(python3 -c "print($normalized_18S_count/$normalized_chr3_count)")

# Return in STDOUT the value
echo $nb_copies



