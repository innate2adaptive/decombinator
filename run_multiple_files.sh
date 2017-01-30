 #!/bin/bash

function decombine {
 		if [ "$2" == "a" ]
 			then chain="alpha"
 		elif [ "$2" = "b" ]
 			then chain="beta"
 		fi
 		python SC_Decombinator.py -fq $1 -c $chain -o $3 -nbc
 		echo " "
 		echo "unzipping..."
		length=${#fastq_file}
		remove_ext=$(($length - 6))
		gzip -d "./dcr_"$chain"_"${fastq_file:0:remove_ext}".n12.gz"
		echo " "
		echo "moving to results folder..."
		mv "./dcr_"$chain"_"${fastq_file:0:remove_ext}".n12" $directory"/dcr_"$chain"_"${fastq_file:0:remove_ext}"_"$3".n12"
		echo " "
		echo "cleaning up Decombinator..."
}

directory="/home/tom/complex/MP1/bashtest"
files=()
filecount=0

while read file_name || [[ -n $file_name ]]; do
 	files[$filecount]="$file_name"
 	(( filecount++ ))	
done <list_of_files

echo " "
 for fastq_file in "${files[@]}"
 	do
 		decombine $fastq_file a forward
		rm dcr_*
		decombine $fastq_file a reverse
		rm dcr_*
		decombine $fastq_file b forward
		rm dcr_*
		decombine $fastq_file b reverse
		rm dcr_*
 	done
# for fastq_file in "${files[@]}"
# 	do
# 		length=${#fastq_file}
# 		remove_ext=$(($length - 6))
# 		echo "dcr_alpha_"${fastq_file:0:remove_ext}".n12.gz"
# 	done