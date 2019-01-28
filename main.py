# handler for FAnalysis
# make sure that all the input files are in order
# maybe make this a function in the baseAnalysis class?
# calls the analysis code



import sys, os
import analysisBase as ab

if os.path.lexists("analysisClass.py"):
	os.remove("analysisClass.py")
os.symlink("{}".format("macros/"+sys.argv[1]), "analysisClass.py")

import analysisClass as tm

tm.addLoop()

# arguments
# [0] main.py
# [1] macro name/path that has the loop function
# [2] input files
# [3] unique identifier for subject dataset
# [4] list of tree names
# [5] directory name for extra root files https://stackoverflow.com/questions/17255737/importing-variables-from-another-file
# [6] outputFile name
# python main.py <macro>.py <fileID> <inputTrees>.txt <dir> <outputName>

#for iArg in range(len(sys.argv)):
#	print(str(iArg), sys.argv[iArg])


if not os.path.exists("macros/"+sys.argv[1]):
	exit("Macro path doesn't exist")

if not os.path.exists('input_conf/'+sys.argv[3]):
	exit("Tree List doesn't exist")


analysis = ab.baseClass(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
analysis.run()

# create instance ofobject analysisClass with string args (inputList, treeList, outFile)

# 
