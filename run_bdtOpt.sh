#!/bin/bash


NAME=$1
NAMEEXT=${2:-""}
MACRO="bdtOpt/$NAME.py"
OUT="$NAME$NAMEEXT.root"

echo "Submitting jobs for $NAME"
date +%T

for fileID in QCD16 QCD17 QCD18PRE QCD18POST TTJets16 TTJets17 TTJets18PRE TTJets18POST WJets16 WJets17 WJets18PRE WJets18POST ZJets16 ZJets17 ZJets18PRE ZJets18POST base
#for fileID in base z15 z20 z25 z35 z40 z45 d10 d40 d60 d80 d00 r02 r05 r07 r08 adh adl
do
	if [ ! -d "outputs/bdtOpt/$fileID" ]; then
		mkdir outputs/bdtOpt/$fileID
	fi
	python main.py $MACRO $fileID inputTree_skims.txt outputs/bdtOpt/$fileID $OUT 2 >& outputs/bdtOpt/$fileID/$NAME$NAMEEXT.out &
	echo "Submitted  $fileID"
	echo $(date +%T) $(jobs | wc -l)
	sleep 10
done
echo Working...

sleep 60 #sleep for a minute in case macro is bugged
echo Time Jobs Left
while [[ $(jobs | wc -c) -ne 0 ]] # while jobs returns not zero letters:
do
	echo $(date +%T) $(jobs | wc -l) # print the time and number of lines jobs returns
	sleep 300
done

echo Finished.
date +%T
