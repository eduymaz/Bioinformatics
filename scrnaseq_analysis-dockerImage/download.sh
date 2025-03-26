while getopts n:o: flag
do
    case "${flag}" in
        n) SRR_ID=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

SRR_ERROR=0
fastq-dump --split-files $SRR_ID || SRR_ERROR=1

if [ $SRR_ERROR -eq 1 ]; then
    echo 'SRR error'
    exit 1
fi