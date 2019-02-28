#!/bin/bash


NAME=$1
MACRO="$NAME.py"
OUT="$NAME.root"

if [ $NAME = "makeFriend" ]
then
OUT="outFriend.root"
fi
echo "Submitting jobs for $NAME"

for fileID in QCD16 TTJets16 WJets16 ZJets16 QCD17 TTJets17 WJets17 ZJets17
do
	if [ ! -d "outputs/$fileID" ]; then
		mkdir outputs/$fileID
	fi
	python main.py $MACRO $fileID inputTree.txt outputs/$fileID $OUT 2 >& outputs/$fileID/$NAME.out &
	echo "Submitted  $fileID"
	sleep 30
done
echo Working...

wait

echo Finished.
