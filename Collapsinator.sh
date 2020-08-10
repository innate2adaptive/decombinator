#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
#running Collpasinator  through loop 
chmod u+r+x Collapsinator.sh
#for file in *.fq.gz
for file in *.n12.gz
do 
python Collapsinator.py -in $file &
done
wait
