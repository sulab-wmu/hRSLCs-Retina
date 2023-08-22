
idx=${1}
fq_file_pre=${2}
#/share/pub/xiongyc/project/scRNA/JiangFanChen/drop-seq/01-fq
output_file_pre=${3}
cpu=${4}

echo "#######"${output_file_pre}
#### mapping
ref=/share/pub/xiongyc/ref/hg38/ebi/star_bulk_149
fq_R1=${fq_file_pre}/${idx}_R1.fq.gz
fq_R2=${fq_file_pre}/${idx}_R2.fq.gz
bam_aligned=${output_file_pre}/01-bam_star_${idx}_
log=${output_file_pre}/01-log_align.log


STAR \
	--genomeDir ${ref} \
	--readFilesIn ${fq_R1} ${fq_R2} \
	--readFilesCommand zcat \
	--runMode alignReads \
	--runThreadN ${cpu} \
	--genomeLoad NoSharedMemory \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outSAMunmapped Within \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS nM NM ch \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--sjdbScore 2 \
	--twopassMode Basic \
	--twopass1readsN -1 \
	--outFileNamePrefix ${bam_aligned} \
	--outFilterIntronMotifs None \
	--alignSoftClipAtReferenceEnds Yes \
	--outFilterMultimapScoreRange 1 \
	--chimJunctionOverhangMin 15 \
	--sjdbOverhang 149 \
	--outFilterScoreMinOverLread 0.33 \
	--outFilterMatchNminOverLread 0.33 \
	--limitSjdbInsertNsj 1200000 \
	--outSAMstrandField intronMotif \
	--chimSegmentMin 15 \
	--chimSegmentReadGapMax 0 \
	--chimOutType Junctions WithinBAM SoftClip \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimOutJunctionFormat 1 \
	--outSAMattrRGline ID:rg1 SM:sm1 \
	--outWigType wiggle \
	--outWigStrand Unstranded  \
	--outWigNorm RPM


#### RSEM quantify
rsem_output=${output_file_pre}/03_${idx}.rsem

rsem-calculate-expression \
	--num-threads ${cpu} \
	--fragment-length-max 1000 \
	--no-bam-output \
	--paired-end \
	--strandedness none \
	--forward-prob 0.5 \
	--alignments \
	--append-names \
	${bam_aligned}Aligned.toTranscriptome.out.bam  \
	${ref}/../star_rsem_149/star_rsem_149  \
	${rsem_output}

<<COMMENT

rsem-calculate-expression --paired-end \
	--alignments \
	-p 8 \
	${bam_aligned}Aligned.toTranscriptome.out.bam  \
	${ref}/../star_rsem_149/star_rsem_149  \
	${rsem_output}

COMMENT



