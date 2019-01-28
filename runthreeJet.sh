#!/bin/bash


NAME=$1
MACRO="$NAME.py"
OUT="$NAME.root"

if [ $NAME = "makeFriend" ]
then
OUT="outFriend.root"
fi
echo "Submitting jobs for $NAME"

for fileID in 3232-16 3232-17 323p mz1 mz15 mz2 mz25 mz35 mz4 alow ahigh rinv0 rinv1 rinv2 rinv4 rinv5 rinv6 rinv7 rinv8 rinv9 rinv10 md50 md5
do
	if [ ! -d "outputs/$fileID" ]; then
		mkdir outputs/$fileID
	fi
	python main.py $MACRO $fileID inputTree.txt outputs/$fileID $OUT 2 >& outputs/$fileID/$NAME.out &
	echo "Submitted  $fileID"
	sleep 3
done
echo Working...

wait

echo Finished.
