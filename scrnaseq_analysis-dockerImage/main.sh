SKIP_QC=0
while getopts s flag
do
    case "${flag}" in
        s) SKIP_QC=1;;
    esac
done
OUTPUT=/app/out
SRR_DATA_PATH="/app/data"
QUALITY_OUTPUT="$OUTPUT/quality"
mkdir $QUALITY_OUTPUT
KALLISTO_OUTPUT="$OUTPUT/kallisto"
SEURAT_INPUT="$KALLISTO_OUTPUT/counts_unfiltered"
SEURAT_OUTPUT="$OUTPUT/seurat"
mkdir $SEURAT_OUTPUT

#bash download.sh -n $SRR_ID -o $SRR_OUTPUT
# error control for download process
if [ $SKIP_QC -eq 0 ]; then
    bash quality_check.sh -d $SRR_DATA_PATH  -o $QUALITY_OUTPUT &
fi
# define inside docker env
# kb ref -i $KALLISTO_INDEX -g $KALLISTO_T2G -f1 $KALLISTO_TRANSCRIPTOME $HUMAN_GENOME_FASTA $HUMAN_GENOME_GTF
bash kallisto.sh -d $SRR_DATA_PATH -o $KALLISTO_OUTPUT
Rscript seurat.R -d $SEURAT_INPUT -t $KALLISTO_T2G -o $SEURAT_OUTPUT