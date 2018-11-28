#!/bin/bash


NAME=$1
MACRO="$NAME.py"
OUT="$NAME.root"

if [ $NAME = "makeFriend" ]
then
OUT="outFriend.root"
fi
echo Submitting jobs for $NAME 
python main.py $MACRO inputTestISR1232.txt inputTree.txt cusNtup/1232 $OUT 2 >& cusNtup/1232/$NAME.out &
echo Submitted  1 / 13
sleep 2
python main.py $MACRO inputTestISR2232.txt inputTree.txt cusNtup/2232 $OUT 2 >& cusNtup/2232/$NAME.out &
echo Submitted  2 / 13
sleep 2
python main.py $MACRO inputTestISR4232.txt inputTree.txt cusNtup/4232 $OUT 2 >& cusNtup/4232/$NAME.out &
echo Submitted  3 / 13
sleep 2
python main.py $MACRO inputTestISR3132.txt inputTree.txt cusNtup/3132 $OUT 2 >& cusNtup/3132/$NAME.out &
echo Submitted  4 / 13
sleep 2
python main.py $MACRO inputTestISR3532.txt inputTree.txt cusNtup/3532 $OUT 2 >& cusNtup/3532/$NAME.out &
echo Submitted  5 / 13
sleep 2
python main.py $MACRO inputTestISR3h32.txt inputTree.txt cusNtup/3h32 $OUT 2 >& cusNtup/3h32/$NAME.out &
echo Submitted  6 / 13
sleep 2
python main.py $MACRO inputTestISR3212.txt inputTree.txt cusNtup/3212 $OUT 2 >& cusNtup/3212/$NAME.out &
echo Submitted  7 / 13
sleep 2
python main.py $MACRO inputTestISR3252.txt inputTree.txt cusNtup/3252 $OUT 2 >& cusNtup/3252/$NAME.out &
echo Submitted  8 / 13
sleep 2
python main.py $MACRO inputTestISR3272.txt inputTree.txt cusNtup/3272 $OUT 2 >& cusNtup/3272/$NAME.out &
echo Submitted  9 / 13
sleep 2
python main.py $MACRO inputTestISR3232.txt inputTree.txt cusNtup/3232 $OUT 2 >& cusNtup/3232/$NAME.out &
echo Submitted 10 / 13
sleep 2
python main.py $MACRO inputTestISR3231.txt inputTree.txt cusNtup/3231 $OUT 2 >& cusNtup/3231/$NAME.out &
echo Submitted 11 / 13
sleep 2
python main.py $MACRO inputTestISR3235.txt inputTree.txt cusNtup/3235 $OUT 2 >& cusNtup/3235/$NAME.out &
echo Submitted 12 / 13
sleep 2
python main.py $MACRO inputTestISR323p1.txt inputTree.txt cusNtup/323p1 $OUT 2 >& cusNtup/323p1/$NAME.out &
echo Submitted 13 / 13
echo Working...

wait

echo Finished.
