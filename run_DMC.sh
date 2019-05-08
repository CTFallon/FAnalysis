#!/bin/bash


NAME=$1
MACRO="datamc/$NAME.py"
OUT="$NAME.root"

echo "Submitting jobs for $NAME"


#for fileID in QCD16 TTJets16 WJets16 ZJets16 QCD17 TTJets17 WJets17 ZJets17
for fileID in Data16 Data17
do
	if [ ! -d "outputs/dmc/$fileID" ]; then
		mkdir outputs/dmc/$fileID
	fi
	python main.py $MACRO $fileID inputTree_skims.txt outputs/dmc/$fileID $OUT 2 >& outputs/dmc/$fileID/$NAME.out &
	echo "Submitted  $fileID"
	sleep 10
done
echo Working...

wait

echo Finished.
