#! /bin/bash
# for i in `ls dir`; do ./cluster_mcl.sh $i gtex; done

name=`basename $1 _Spearman.txt`
i=5.0

if [ $2 = "gtex" ]; then
  inDir='/opt/data/GTEx/processed/fi/wel-spearman/'
  outDir='/opt/data/GTEx/processed/fi/mcl-spearman'
elif [ $2 = "tcga" ]; then
  inDir='/opt/data/TCGA/processed/fi/wel-spearman/'
  outDir='/opt/data/TCGA/processed/fi/mcl-spearman'
else
   echo "Please provide a consortia name - ex tcga or gtex"
fi

fileExt='_Spearman.txt'

echo $name

mcl $inDir$name$fileExt --abc -I $i -odir $outDir
