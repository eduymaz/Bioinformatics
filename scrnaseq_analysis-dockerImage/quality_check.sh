while getopts d:o: flag
do
    case "${flag}" in
        d) SRR_PATH=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

for file in $(ls $SRR_PATH)
do
    fastqc -t 5 $SRR_PATH/$file -o $OUTPUT &
done
wait

