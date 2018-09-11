#!/bin/bash


if [ $# -lt 7 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName seqYName lenght similarity WL fixedL output.frags"
   echo ""
   exit -1
fi

MGDIR=WWJJIITEMPFILE
genome1=$(basename "$1")
genome2=$(basename "$2")
genome1="${genome1%.*}"
genome2="${genome2%.*}"

MYROUTE=$PWD
BINARIES=/home/galaxy/galaxy/tools/metagecko/binaries/gecko_gg/bin


cp $1 $MGDIR/${genome1}.fasta
cp $2 $MGDIR/${genome2}.fasta

cd $MGDIR

echo "workflow.sh $MGDIR/${genome1}.fasta $MGDIR/${genome2}.fasta $3 $4 $5 $6"


echo "my route is $PWD"

source /home/galaxy/galaxy/tools/metagecko/binaries/gecko_gg/bin/workflow.sh $MGDIR/${genome1}.fasta $MGDIR/${genome2}.fasta $3 $4 $5 $6

cd $MYROUTE

echo "mv $MGDIR/results/${genome1}-${genome2}.frags $7"
mv $MGDIR/results/${genome1}-${genome2}.frags $7

#rm -rf $MGDIR/*
#rm -rf $MGDIR/intermediateFiles
#rm -rf $MGDIR/results


