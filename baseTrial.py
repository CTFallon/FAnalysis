# first attempt to do analysis with pyroot

# loop over all events
# make one cut (BTags?)
# plot one variable (HT?)

import ROOT as rt
from sys import argv
import ConfigParser as cp
import json

# argv
# [1] - input list, .ini configParser file
# [2] - output file, full path name of the output file

config = cp.SafeConfigParser()
config.read(argv[1])

bkgList = config.sections()



# create/open output file
outFile = rt.TFile(argv[2],"RECREATE")

# next, create any stack histograms we want to fill
h_GenPart_Mag_stack = rt.THStack("h_GenPart_Mag_stack","GenPart_Mag_stack")
h_GenPart_3P_stack = rt.THStack("h_GenPart_3P_stack","GenPart_3P_stack")
h_GenPart_Pt_stack = rt.THStack("h_GenPart_Pt_stack","GenPart_Pt_stack")



for bkgType in bkgList:
	tempChain = rt.TChain(bkgType)
	for iFile in range(config.getint(bkgType,'nfile')):
		fileName = config.get(bkgType,'file'+str(iFile)+".Name")
		tempChain.Add(fileName)

	# from here to just after the cuts should be indpendant
	# to make easy changes and save mulitple analysises

	h_Mag_temp = rt.TH1F("h_GenPart_Mag_"+bkgType,"GenPart_Mag_"+bkgType,1000,-1,45)
	h_3P_temp = rt.TH1F("h_GenPart_3P_"+bkgType,"GenPart_3P_"+bkgType,1000,0,1000)
	h_Pt_temp = rt.TH1F("h_GenPart_Pt_"+bkgType,"GenPart_Pt_"+bkgType,1000,0,1000)
	h_Mag_temp.SetFillColor(config.getint(bkgType,'color'))
	h_3P_temp.SetFillColor(config.getint(bkgType,'color'))
	h_Pt_temp.SetFillColor(config.getint(bkgType,"color"))
	nEvents = tempChain.GetEntries()
	print "Numeber of Events is " + str(nEvents) + "."

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print "Processing Event Number:", iEvent, "in background", bkgType
		tempChain.GetEvent(iEvent)
		if tempChain.BTags == 3:
			for particle in tempChain.GenParticles:
				h_Mag_temp.Fill(particle.Mag(),tempChain.Weight)
				h_3P_temp.Fill(particle.P(),tempChain.Weight)
				h_Pt_temp.Fill(particle.Pt(),tempChain.Weight)
	h_GenPart_Mag_stack.Add(h_Mag_temp,"hist")
	h_GenPart_3P_stack.Add(h_3P_temp,"hist")
	h_GenPart_Pt_stack.Add(h_Pt_temp,"hist")
	h_Mag_temp.Write()
	h_3P_temp.Write()
	h_Pt_temp.Write()
	del tempChain
	del h_Mag_temp
	del h_3P_temp
	del h_Pt_temp
# This ends the "indpendant" code

# next part should write all objects, be sure to add all objects
# to an array/list to loop over and write
h_GenPart_Mag_stack.Write()
h_GenPart_3P_stack.Write()
h_GenPart_Pt_stack.Write()



