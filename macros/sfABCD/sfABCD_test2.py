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
	
	nBins = 1000

	# list of branches to plot
	plotDict = {
	'LeadJetPt':[nBins,200,3000,self.fileID+";Jet Pt, lead jet;events"],
	'SubLeadJetPt':[nBins,200,3000,self.fileID+";Jet Pt, sub jet;events"],
	'BothJetPt':[nBins,200,3000,self.fileID+";Jet Pt, both jets;events"],
	'MT':[nBins,1500,5000,self.fileID+"; MT;events"],
	'MET':[nBins,200,3000,self.fileID+"; MET;events"],
	'METPhi':[nBins,-rt.TMath.Pi(),rt.TMath.Pi(),self.fileID+"; METPhi;events"],
	'DeltaEta':[44,0.0,2.2,self.fileID+";#Delta#eta;events"],
	'nSVJTag':[3,0,3,self.fileID+";Num. SVJ Tags;events"],
	'RT':[65,0.15,0.8,self.fileID+";R_{T};events"]
	}

	histDict = {}
	# define the WP
	WP = 0.6
	listOfRegions = ["A0","A1","A2","B0","B1","B2"]
	for plotVar, histSpecs in plotDict.items():
		for region in listOfRegions:
			histDict[plotVar+"_"+region] = self.makeTH1F(self.fileID+"_"+plotVar+"_"+region, region+"_"+histSpecs[3], histSpecs[0], histSpecs[1], histSpecs[2])
	
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
		if tree.MT_AK8 < 1850:
			continue
		numSVJ = int(bool(tree.JetsAK8_bdtSVJtag[0]>WP)) + int(bool(tree.JetsAK8_bdtSVJtag[1]>WP))
		RTval = tree.MET/tree.MT_AK8
		deltaEta = abs(tree.JetsAK8[0].Eta()-tree.JetsAK8[1].Eta())
		# dEta cut, dont worry about. highdeta and deta skims have that covered.
		#region A0 rt > .25, dEta < 1.5, 0 SVJ
		#region A1 rt > .25, dEta < 1.5, 1 SVJ
		#region A2 rt > .25, dEta < 1.5, 2 SVJ
		#region B0 0.15 < rt < .25, dEta < 1.5, 0 SVJ
		#region B1 0.15 < rt < .25, dEta < 1.5, 1 SVJ
		#region B2 0.15 < rt < .25, dEta < 1.5, 2 SVJ
		if ((RTval > 0.25) and (numSVJ == 0)):
			histDict['LeadJetPt_A0'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_A0'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_A0'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_A0'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_A0'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_A0'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_A0"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_A0"].Fill(numSVJ, lumi*weight)
			histDict["RT_A0"].Fill(RTval, lumi*weight)
			histDict['MT_A0'].Fill(tree.MT_AK8, lumi*weight)
		elif ((RTval > 0.25) and (numSVJ == 1)):
			histDict['LeadJetPt_A1'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_A1'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_A1'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_A1'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_A1'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_A1'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_A1"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_A1"].Fill(numSVJ, lumi*weight)
			histDict["RT_A1"].Fill(RTval, lumi*weight)
			histDict['MT_A1'].Fill(tree.MT_AK8, lumi*weight)
		elif ((RTval > 0.25) and (numSVJ == 2)):
			histDict['LeadJetPt_A2'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_A2'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_A2'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_A2'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_A2'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_A2'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_A2"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_A2"].Fill(numSVJ, lumi*weight)
			histDict["RT_A2"].Fill(RTval, lumi*weight)
			histDict['MT_A2'].Fill(tree.MT_AK8, lumi*weight)
		elif ((RTval > 0.15) and (RTval < 0.25) and (numSVJ == 0)):
			histDict['LeadJetPt_B0'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_B0'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_B0'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_B0'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_B0'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_B0'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_B0"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_B0"].Fill(numSVJ, lumi*weight)
			histDict["RT_B0"].Fill(RTval, lumi*weight)
			histDict['MT_B0'].Fill(tree.MT_AK8, lumi*weight)
		elif ((RTval > 0.15) and (RTval < 0.25) and (numSVJ == 1)):
			histDict['LeadJetPt_B1'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_B1'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_B1'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_B1'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_B1'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_B1'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_B1"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_B1"].Fill(numSVJ, lumi*weight)
			histDict["RT_B1"].Fill(RTval, lumi*weight)
			histDict['MT_B1'].Fill(tree.MT_AK8, lumi*weight)
		elif ((RTval > 0.15) and (RTval < 0.25) and (numSVJ == 2)):
			histDict['LeadJetPt_B2'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['SubLeadJetPt_B2'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['BothJetPt_B2'].Fill(tree.JetsAK8[0].Pt(), lumi*weight)
			histDict['BothJetPt_B2'].Fill(tree.JetsAK8[1].Pt(), lumi*weight)
			histDict['MET_B2'].Fill(tree.MET, lumi*weight)
			histDict['METPhi_B2'].Fill(tree.METPhi, lumi*weight)
			histDict["DeltaEta_B2"].Fill(deltaEta, lumi*weight)
			histDict["nSVJTag_B2"].Fill(numSVJ, lumi*weight)
			histDict["RT_B2"].Fill(RTval, lumi*weight)
			histDict['MT_B2'].Fill(tree.MT_AK8, lumi*weight)
		else:
			print("Yall dun goofed! {} {} {}".format(numSVJ, RTval, detlaEta))


def addLoop():
	baseClass.loop = loop


