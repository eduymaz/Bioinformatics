while getopts d:o: flag
do
    case "${flag}" in
        d) SRR_PATH=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

FILE_COUNT=$(ls -1 $SRR_PATH | wc -l)
FILES=($(ls $SRR_PATH))
if [ $FILE_COUNT -eq 1 ]; then
    kb count -i $KALLISTO_INDEX -g $KALLISTO_T2G -o $OUTPUT -x 10xv3 $SRR_PATH/${FILES[0]}
    return 0
fi
# kb ref -i $KALLISTO_INDEX -g $KALLISTO_T2G -f1 $KALLISTO_TRANSCRIPTOME $HUMAN_GENOME_FASTA $HUMAN_GENOME_GTF
kb count -i $KALLISTO_INDEX -g $KALLISTO_T2G -o $OUTPUT -x 10xv3 $SRR_PATH/${FILES[0]} $SRR_PATH/${FILES[1]}
