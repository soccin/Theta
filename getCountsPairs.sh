#!/bin/bash



VCF=$1
NORMAL=$2
TUMOR=$3
SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [[ ! "$LD_LIBRARY_PATH" =~ htslib/htslib-1.9 ]]; then
    export LD_LIBRARY_PATH=/opt/common/CentOS_7-dev/htslib/htslib-1.9:$LD_LIBRARY_PATH
fi

getBAMSMTag() {
    samtools view -H $1 \
        | fgrep SM: \
        | tr '\t' '\n' \
        | fgrep SM: \
        | uniq \
        | sed 's/SM://' \
        | xargs \
        | tr ' ' ','
}

normalId=$(getBAMSMTag $NORMAL)
tumorId=$(getBAMSMTag $TUMOR)

D1=$(echo $VCF | md5sum - | awk '{print $1}' | perl -pe 's/^(.......).*/\1/')
D2=$(echo $normalId $tumorId | md5sum - | awk '{print $1}' | perl -pe 's/^(..).*/\1/')

ODIR=Counts/$D1/$D2
mkdir -p $ODIR
echo $ODIR

echo $normalId $tumorId
$SDIR/code/snp-pileup -g -P 500 \
    -v $VCF \
    $ODIR/counts_${tumorId}___${normalId}_.txt.gz $NORMAL $TUMOR

md5sum $ODIR/counts_${tumorId}___${normalId}_.txt.gz >$ODIR/counts_${tumorId}___${normalId}_.txt.gz.md5
