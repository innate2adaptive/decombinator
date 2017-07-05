 #!/bin/bash

function removepath {
	b=$(basename $1)
	echo $b
}

function decombine {
 		if [ "$2" == "a" ]
 			then chain="alpha"
 		elif [ "$2" = "b" ]
 			then chain="beta"
 		fi
 		python SC_Decombinator.py -fq $1 -c $chain -o $3 -nbc
 		echo " "
 		echo "unzipping..."
		remove_path=$(removepath $fastq_file)
		length=${#remove_path}
		remove_ext=$(($length - 6))
		gzip -d "./dcr_"$chain"_"${remove_path:0:remove_ext}".n12.gz"
		echo " "
		echo "moving to results folder..."
		mv "./dcr_"$chain"_"${remove_path:0:remove_ext}".n12" $directory"/dcr_"$chain"_"${remove_path:0:remove_ext}"_"$3".n12"
		echo " "
		echo "cleaning up Decombinator..."
}

directory=$2
files=()
filecount=0
list_of_files=$1

while read file_name || [[ -n $file_name ]]; do
 	files[$filecount]="$file_name"
 	(( filecount++ ))	
done <$list_of_files

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