#!/bin/bash


NAME=$1
NAMEEXT=${2:-""}
MACRO="sfABCD/$NAME.py"
OUT="$NAME$NAMEEXT.root"

echo "Submitting jobs for $NAME"
echo $(date +%T) $(jobs | wc -l)

for fileID in Data16 Data17 Data18PRE Data18POST QCD16 QCD17 QCD18PRE QCD18POST TTJets16 TTJets17 TTJets18PRE TTJets18POST WJets16 WJets17 WJets18PRE WJets18POST ZJets16 ZJets17 ZJets18PRE ZJets18POST
#for fileID in base
#for fileID in QCD17 TTJets17 WJets17 ZJets17 base
do
	if [ ! -d "outputs/sfABCD/$fileID" ]; then
		mkdir outputs/sfABCD/$fileID
	fi
	python main.py $MACRO $fileID inputTree_skims.txt outputs/sfABCD/$fileID $OUT 2 >& outputs/sfABCD/$fileID/$NAME$NAMEEXT.out &
	echo "Submitted  $fileID"
	echo $(date +%T) $(jobs | wc -l)
	sleep 10
done
echo Working...

sleep 60 #sleep for a minute in case macro is bugged
echo Time Jobs Left
while [[ $(jobs | wc -l) -ne 0 ]] # while jobs returns not zero lines:
do
	echo $(date +%T) $(jobs | wc -l) # print the time and number of lines jobs returns
	jobs
	sleep 120 # query every 2 mins
done

echo Finished.
echo $(date +%T) $(jobs | wc -l)
