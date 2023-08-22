
#! /bin/bash



#for i in /share/pub/xiongyc/project/scRNA/JiangFanChen/data/batch5/N2102608_ZW_80-611488948_SEQ/210309-X2B_L002/*.gz
for i in /share2/pub/zhouyj/zhouyj/Liu/20230810_bulk/Data/*.gz
do
    z=`echo $i | awk -F '/' '{print $NF}' | awk -F '[_-]' '{print $1"_"$NF}'  `
    echo $z
    cp $i $z
    #ln -s $i $z
done


