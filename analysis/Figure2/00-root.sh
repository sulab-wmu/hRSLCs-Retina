#!/usr/bin/env bash

#PBS -N Drop-PBMC-rsem
#PBS -l mem=150gb
#PBS -l ncpus=16
#PBS -j oe

workspace=/share2/pub/zhouyj/zhouyj/Liu/20230810_bulk/bulk_fatsq
cd $workspace

for idx in  {S01T0001,S01T0002,S01T0003,S01T0004,S01T0005,S01T0006}
do
	bash -v ./03-mapping/01-RNA-seq.sh workspace=${workspace} , idx=${idx}
done


