#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
#running Decombinator though loop 
for file in LTX_*.fq
do 
python Decombinator.py -fq file -sp human 
done
