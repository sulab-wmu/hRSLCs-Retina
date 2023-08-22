species='hg'
DATA_FOLDER="/share2/pub/xiongyc/xiongyc/project/scRNA/HuiLiu/20220107_snRNAseq/05-tScenic/fat"
RESOURCES_FOLDER=${DATA_FOLDER}/../cisTarget_databases

cd ${DATA_FOLDER}

echo "#####################grn"
#python 0611-00-file_format_transform.py ${DATA_FOLDER}/organoid-genes_x_cells.csv
dataIn=/share2/pub/zhouyj/zhouyj/Liu/analysis/Reg/drived_data/scenic_final2/organoid-genes_x_cells.csv.loom

tf="hs_hgnc_tfs.txt"
if [[ "$species" == "mm" ]];then
	tf="mm_mgi_tfs.txt"
fi

echo "############# TF:"${tf}
pyscenic grn \
        --num_workers 16 \
        -o  ${dataIn}.out.expr_mat.adjacencies.tsv  \
        ${dataIn}  \
        ${RESOURCES_FOLDER}/${tf}


echo "###################ctx"
fea1="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
fea2="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
anno="motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
if [[ "$species" == "mm" ]];then
	fea1="mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
	fea2="mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
	anno="motifs-v9-nr.mgi-m0.001-o0.0.tbl"
fi

echo "################# fea1 and anno: "${fea1}", "${fea2}", "${anno}

pyscenic ctx \
        ${dataIn}.out.expr_mat.adjacencies.tsv \
        ${RESOURCES_FOLDER}/v2/${fea1} \
        ${RESOURCES_FOLDER}/v2/${fea2} \
        --annotations_fname ${RESOURCES_FOLDER}/${anno} \
        --expression_mtx_fname ${dataIn} \
        --mode "dask_multiprocessing" \
        --output ${dataIn}.out.regulons.csv \
        --num_workers 16

echo "###################### aucell"

pyscenic aucell \
        ${dataIn} \
        ${dataIn}.out.regulons.csv \
        -o ${dataIn}.out.auc_mtx.loom \
        --num_workers 16
        
        
