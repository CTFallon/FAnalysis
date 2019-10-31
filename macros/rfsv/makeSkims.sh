#!/bin/bash


NJOBS=0

for dSet in rf rD rH fD fH
do
	echo "Dataset $dSet"

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 20 0.3 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 2000 20 0.3 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 4000 20 0.3 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 50 0.3 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 100 0.3 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 20 0.5 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 20 0.7 peak

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 20 0.3 low

	NJOBS=$((NJOBS+1))
	echo "Runnig Job $NJOBS of 45"
	python rfsv_makeHistos.py $dSet 3000 20 0.3 high
done
