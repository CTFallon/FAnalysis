from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

def loop(self):
	# set up trees or chains
	#f = rt.TFile.Open(self.inputFileList[0])
	#tree = f.Get(self.treeNameList[0])
	tree = self.getChain(self.treeNameList[0])
	# added friend tree
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))


	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()

	nBins = 1000
	# list of branches to plot
	plotDict = {#key = var name, value = [varType, nBins, binLow, binHigh, title]
				# varType can be "s" - single value (ie 'MET')
				#				 "sF" - single value, but a function of (not sure if this exists, but just in case)
				#				 "vA" - vector, all values (ie 'JetsAK8_girth')
				#				 "vI" - vector, only index value (ie 'JetsAK8_girth[0]')
				#				 "vAF" - vector, all values but function (ie 'JetsAK8.Pt()' - Pt of all ak8 Jets)
				#				 "vIF" - vector, indexed function (ie 'JetsAK8[0].Pt()', only Pt of leading AK8 Jet)
				#				 "vR", "vRF
	'JetsAK8_bdtSVJtag[0]':["vI",nBins,0,1,self.fileID+"; SVJ BDT Output, lead jet; Events"],
	'JetsAK8_bdtSVJtag[1]':["vI",nBins,0,1,self.fileID+"; SVJ BDT Output, sub jet; Events"],
	'JetsAK8_bdtSVJtag':["vA",nBins,0,1,self.fileID+"; SVJ BDT Output, all; Events"]
	}

	hist_leadBDT = self.makeTH1F("hist_leadBDT",plotDict['JetsAK8_bdtSVJtag[0]'][4], plotDict['JetsAK8_bdtSVJtag[0]'][1],plotDict['JetsAK8_bdtSVJtag[0]'][2],plotDict['JetsAK8_bdtSVJtag[0]'][3])
	hist_subBDT = self.makeTH1F("hist_subBDT", plotDict['JetsAK8_bdtSVJtag[1]'][4],plotDict['JetsAK8_bdtSVJtag[1]'][1],plotDict['JetsAK8_bdtSVJtag[1]'][2],plotDict['JetsAK8_bdtSVJtag[1]'][3])
	hist_allBDT = self.makeTH1F("hist_allBDT", plotDict['JetsAK8_bdtSVJtag'][4],plotDict['JetsAK8_bdtSVJtag'][1],plotDict['JetsAK8_bdtSVJtag'][2],plotDict['JetsAK8_bdtSVJtag'][3])
	listOfBranchesOneMightNeed = ["RunNum","HEMOptVetoFilter" ,"madHT","GenElectrons","GenMuons","GenTaus","GenMET","MT_AK8","JetsAK8_bdtSVJtag","Weight","puWeight"]
	listOfBranches = tree.GetListOfBranches()
	for branch in listOfBranches:
		if branch.GetName() in listOfBranchesOneMightNeed:
			branch.SetStatus(1)
		else:
			branch.SetStatus(0)
#	tree.SetBranchStatus("*",0)
#	tree.SetBranchStatus("RunNum" ,1)
#	tree.SetBranchStatus("HEMOptVetoFilter" ,1)
#	tree.SetBranchStatus("madHT",1)
#	tree.SetBranchStatus("GenElectrons",1)
#	tree.SetBranchStatus("GenMuons",1)
#	tree.SetBranchStatus("GenTaus",1)
#	tree.SetBranchStatus("GenMET",1)
#	tree.SetBranchStatus("MT_AK8",1)
#	tree.SetBranchStatus("JetsAK8_bdtSVJtag",1)
#	tree.SetBranchStatus("Weight",1)
#	tree.SetBranchStatus("puWeight",1)

	if (("Jets" in self.fileID) or ("QCD" in self.fileID)): # only need to do this for MC bkg
		if "16" in self.fileID:
			lumi = 35921.036
			print("2016 Lumi")
		elif "17" in self.fileID:
			lumi = 41521.331
			print("2017 Lumi")
		elif "18PRE" in self.fileID:
			lumi = 21071.460
			print("2018 pre Lumi")
		elif "18POST" in self.fileID:
			lumi = 38621.232
			print("2018 post Lumi")
		elif "18" in self.fileID:
			lumi = 59692.692
			print("2018 full lumi")
		else:
			print("Dont know what total lumi to use. default to 40 fb-1")
			lumi = 40000.
	else:
		lumi = 1
	histDict = {}

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET, self.fileID):
				continue
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)):
			#print("Weights: {} {}".format(tree.Weight, tree.puWeight))
			weight = tree.Weight*tree.puWeight
			if weight == 0.0:
				print("Weights: {} {}".format(tree.Weight, tree.puWeight))
		else:
			weight = 1.
		if ("18" in self.fileID):
			# four cases:
			# Data18PRE only use events from runs < 319077
			# mc18PRE use eveything
			# Data18POST only use envets from runs >= 319077 and pass HEMO
			# mc18POST only use evnets that pass HEMO
			if "PRE" in self.fileID:
				if (("Data" in self.fileID) and (tree.RunNum >= 319077)):
					continue
				else:
					pass
			elif ("POST" in self.fileID):
				if tree.HEMOptVetoFilter == 0:
					continue
				else:
					if (("Data" in self.fileID) and (tree.RunNum < 319077)):
						continue
				
				# getattr is funky for methods. for a public variable of a class, getattr(obj, attr) works.
				# for a public function, needs getattr(obj,func)(args of func)
				# since JetsAK8[0].Pt() is a public function, the proper way to get that value using getattr is
				# getattr(JetsAK8[0],Pt)(), which we do here, but JetsAK8[0] is replaced by getattr(tree,JetsAK8)[0]
		#fill histograms
		#if tree.MT_AK8 > 4000 or tree.MT_AK8 < 2000:
		#	continue
		hist_leadBDT.Fill(tree.JetsAK8_bdtSVJtag[0], weight*lumi)
		hist_subBDT.Fill(tree.JetsAK8_bdtSVJtag[1], weight*lumi)
		for tagVal in tree.JetsAK8_bdtSVJtag[:2]:
			hist_allBDT.Fill(tagVal, weight*lumi)



def addLoop():
	baseClass.loop = loop


