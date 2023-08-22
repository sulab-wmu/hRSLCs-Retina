#PBS -N pyscenic
#PBS -q workq
#PBS -l mem=150gb
#PBS -l ncpus=30
source /share/pub/xiongyc/program/conda/install/bin/activate jinli
cd /share2/pub/chenchg/chenchg/SingleCell/LiHui/final/pyscenic/

##change the format of doucment
python file_format_transform.py 

##grn
pyscenic grn \
         --num_workers 16 \
         --output SingleCell_out.expr_mat.adjacencies.tsv  \
         --method grnboost2 \
          SingleCell.loom \
         /share2/pub/chenchg/chenchg/SingleCell/LiHui/final_version/pyscenic/cisTarget_databases/hs_hgnc_tfs.txt

##ctx
pyscenic ctx \
        SingleCell_out.expr_mat.adjacencies.tsv \
        /share2/pub/chenchg/chenchg/SingleCell/LiHui/final_version/pyscenic/cisTarget_databases/v2/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
        /share2/pub/chenchg/chenchg/SingleCell/LiHui/final_version/pyscenic/cisTarget_databases/v2/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
        --annotations_fname /share2/pub/chenchg/chenchg/SingleCell/LiHui/final_version/pyscenic/cisTarget_databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname SingleCell.loom \
        --mode "dask_multiprocessing" \
        --output SingleCell.out.regulons.csv \
        --num_workers 16

##aucell
pyscenic aucell \
         SingleCell.loom\
         SingleCell.out.regulons.csv \
        --output SingleCell.out.auc_mtx.loom \
        --num_workers 16
