
#! /bin/bash

mkdir Fastqc

for i in *fq.gz
do
    fastqc -t 16 -o Fastqc/ $i > $i.log 
done

