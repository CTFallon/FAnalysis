from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

# Comparing N# to N'# (Annapaola's slides, July 31)

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
	#makeTH3F(self, name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax)
	# We want to compare everything to SVJ# and RT value, to cut on when we plot

	nBinsX = 3
	binLowX = 0.0
	binHighX = 3.0
	
	nBinsY = 100
	binLowY = 0.0
	binHighY = 1.0
	
	nBins = 1000

	# list of branches to plot
	plotDict = {
	'LeadJetPt':[nBins,200,3000,self.fileID+";Num. SVJ Tags;R_{T};Jet Pt, lead jet"],
	'SubLeadJetPt':[nBins,200,3000,self.fileID+";Num. SVJ Tags;R_{T}; Jet Pt, sub jet"],
	'BothJetPt':[nBins,200,3000,self.fileID+";Num. SVJ Tags;R_{T}; Jet Pt, both jets"],
	'MET':[nBins,200,3000,self.fileID+";Num. SVJ Tags;R_{T}; MET"],
	'METPhi':[nBins,-rt.TMath.Pi(),rt.TMath.Pi(),self.fileID+";Num. SVJ Tags;R_{T}; METPhi"],
	'DeltaEta':[60,1.0,2.2,self.fileID+";Num. SVJ Tags;R_{T};#Delta#eta"]
	}

	histDict_3d = {}
	# define the WP
	WP = 0.6
	for plotVar, histSpecs in plotDict.items():
		histDict_3d[plotVar] = self.makeTH3F(self.fileID+"_"+plotVar, histSpecs[3], nBinsX, binLowX, binHighX, nBinsY, binLowY, binHighY, histSpecs[0], histSpecs[1], histSpecs[2])
	
	listOfBranchesOneMightNeed = ["RunNum","HEMOptVetoFilter" ,"madHT","GenElectrons","GenMuons","GenTaus","GenMET","MT_AK8","JetsAK8_bdtSVJtag","Weight","puWeightNew"]
	listOfBranches = tree.GetListOfBranches()
	for branch in listOfBranches:
		if branch.GetName() in listOfBranchesOneMightNeed:
			branch.SetStatus(1)
		else:
			branch.SetStatus(0)

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

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET, self.fileID):
				continue
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)):
			weight = tree.Weight*tree.puWeightNew
			if weight == 0.0:
				print("Weights: {} {}".format(tree.Weight, tree.puWeightNew))
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
		numSVJ = int(bool(tree.JetsAK8_bdtSVJtag[0]>WP)) + int(bool(tree.JetsAK8_bdtSVJtag[1]>WP))
		RTval = tree.MET/tree.MT_AK8
		deltaEta = abs(tree.JetsAK8[0].Eta()-tree.JetsAK8[1].Eta())
		#Filling the Histograms:
		histDict_3d['LeadJetPt'].Fill(numSVJ, RTval, tree.JetsAK8[0].Pt(), lumi*weight)
		histDict_3d['SubLeadJetPt'].Fill(numSVJ, RTval, tree.JetsAK8[1].Pt(), lumi*weight)
		histDict_3d['BothJetPt'].Fill(numSVJ, RTval, tree.JetsAK8[0].Pt(), lumi*weight)
		histDict_3d['BothJetPt'].Fill(numSVJ, RTval, tree.JetsAK8[1].Pt(), lumi*weight)
		histDict_3d['MET'].Fill(numSVJ, RTval, tree.MET, lumi*weight)
		histDict_3d['METPhi'].Fill(numSVJ, RTval, tree.METPhi, lumi*weight)
		histDict_3d["DeltaEta"].Fill(numSVJ, RTval, deltaEta, lumi*weight)



def addLoop():
	baseClass.loop = loop


