workspace=/share2/pub/zhouyj/zhouyj/Liu/20230810_bulk/bulk_fatsq/
cd $workspace

source /share/pub/zhouyj/anaconda3/bin/activate SC

fq_file_pre=${workspace}/01-fq
output_file_pre=${workspace}/03-mapping

for idx in  {S01T0001,S01T0002,S01T0003,S01T0004,S01T0005,S01T0006}
do
echo "####"${idx}"#####"

dir_output=${output_file_pre}/${idx}
if [ ! -d ${dir_output} ]; then
   mkdir ${dir_output}
fi
sh ${workspace}/03-mapping/01-RNA-seq-pipeline-STAR.sh ${idx} ${fq_file_pre} ${dir_output} 16
#sh ${workspace}/03-mapping/01-RNA-seq-pipeline-rsem.sh ${idx} ${fq_file_pre} ${dir_output} 16
done



