#!/bin/bash


NAME=$1
MACRO="fsopt/$NAME.py"
OUT="$NAME.root"

echo "Submitting jobs for $NAME"


#for fileID in  QCD16 TTJets16 WJets16 ZJets16 QCD17 TTJets17 WJets17 ZJets17 #base-16 base-17 
#for fileID in  z10 z20 z40 d10 d50 d100 r01 r05 r08 ah al
for fileID in base-16 base-17 a2 z05 z06 z07 z08 z09 z10 z11 z12 z13 z14 z15 z16 z17 z18 z19 z20 z21 z22 z23 z24 z25 z26 z27 z28 z29 z31 z32 z33 z34 z35 z36 z37 z38 z39 z40 z41 z42 z43 z44 z45 d1 d5 d10 d30 d40 d50 d60 d70 d80 d90 d100 r00 r01 r02 r04 r05 r06 r07 r08 r09 r10 al ah
#for fileID in TTJets16 TTJets17
do
	if [ ! -d "outputs/fsopt/$fileID" ]; then
		mkdir outputs/fsopt/$fileID
	fi
	python main.py $MACRO $fileID inputTree_skims.txt outputs/fsopt/$fileID $OUT 2 >& outputs/fsopt/$fileID/$NAME.out &
	echo "Submitted  $fileID"
	sleep 10
done
echo Working...

wait

echo Finished.
