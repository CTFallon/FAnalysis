#!/bin/bash


NAME=$1
MACRO="/noiseFilters/$NAME.py"
OUT="$NAME.root"

echo "Submitting jobs for $NAME"

for fileID in QCD16 TTJets16 WJets16 ZJets16 Data16 QCD17 TTJets17 WJets17 ZJets17 Data17
do
	if [ ! -d "outputs/noiseFilters/$fileID" ]; then
		mkdir outputs/noiseFilters/$fileID
	fi
	python main.py $MACRO $fileID inputTree_skims.txt outputs/noiseFilters/$fileID $OUT 2 >& outputs/noiseFilters/$fileID/$NAME.out &
	echo "Submitted  $fileID"
	sleep 30
done
echo Working...

wait

echo Finished.
