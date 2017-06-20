# attempt to do modular analysis with pyroot

# loop over all events
# make one cut (BTags?)
# plot one variable (HT?)

import ROOT as rt
from sys import argv,modules
import FAparser
import importlib

# argv
# [1] - input list, .ini configParser file
# [2] - output file, full path name of the output file

if len(argv) != 4:
	print "Please use the following syntax:"
	print "python modTrial.py <name_of_loop> <input_parser_file.ini> <name_of_output.root>"
	exit()


importlib.import_module(argv[1])

config = FAparser.Fparser(argv[2])

analysis = modules[argv[1]].analysisClass()

# create/open output file
outFile = rt.TFile(argv[3],"RECREATE")

analysis.run(config)
