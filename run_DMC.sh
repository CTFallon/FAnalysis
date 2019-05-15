#!/bin/bash


NAME=$1
MACRO="datamc/$NAME.py"
OUT="$NAME.root"

echo "Submitting jobs for $NAME"


for fileID in Data16 Data17 Data18PRE Data18POST ST16 ST17 ST18PRE ST18POST QCD16 QCD17 QCD18PRE QCD18POST TTJets16 TTJets17 TTJets18PRE TTJets18POST WJets16 WJets17 WJets18PRE WJets18POST ZJets16 ZJets17 ZJets18PRE ZJets18POST
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
