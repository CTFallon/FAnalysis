# handler for FAnalysis
# make sure that all the input files are in order
# maybe make this a function in the baseAnalysis class?
# calls the analysis code



import sys, os
import analysisBase as ab

if os.path.lexists("analysisClass.py"):
	os.remove("analysisClass.py")
os.symlink("{}".format(sys.argv[1]), "analysisClass.py")

import analysisClass as tm

tm.addLoop()

# arguments
# [0] main.py
# [1] macro name/path that has the loop function
# [2] list of input files
# [3] list of tree names
# [4] outputFile name/path

for iArg in range(len(sys.argv)):
	print(str(iArg), sys.argv[iArg])


if not os.path.exists(sys.argv[1]):
	exit("Macro path doesn't exist")
print("----------")
print("Macro is:")
print(sys.argv[1])
print("----------")
if not os.path.exists(sys.argv[2]):
	exit("Input List doesn't exist")
iFile = open(sys.argv[2])
print("Input Files are:")
for line in iFile:
	print(line[:-1])# don't print new line symbol
print("----------")

if not os.path.exists(sys.argv[3]):
	exit("Tree List doesn't exist")
iTree = open(sys.argv[3])
print("Tree Names are:")
for line in iTree:
	print(line[:-1])# don't print new line symbol
print("----------")


analysis = ab.baseClass(sys.argv[2], sys.argv[3], sys.argv[4])
analysis.run()

# create instance ofobject analysisClass with string args (inputList, treeList, outFile)

# 
