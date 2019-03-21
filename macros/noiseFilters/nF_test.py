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
	branchList = tree.GetListOfBranches()
	passedDict = {}
	passedDict["None"] = 1
	passedDict["All"] = -1
	for branch in branchList:
		if "Filter" in branch.GetName():
			passedDict[branch.GetName()] = -3
			print(branch.GetName())
		else:
			tree.SetBranchStatus(branch.GetName(),0)

	tree.SetBranchStatus("JetsAK8",1)
	tree.SetBranchStatus("MT_AK8",1)
	tree.SetBranchStatus("NElectrons",1)
	tree.SetBranchStatus("NMuons",1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("TriggerPass",1)

	if (("Jets" in self.fileID) or ("QCD" in self.fileID) or self.fileID == "G16"): # only need to do this for MC bkg
		tree.SetBranchStatus("Weight",1)
		if "16" in self.fileID:
			lumi = 35900.
			print("2016 Lumi")
		elif "17" in self.fileID:
			lumi = 41500.
			print("2017 Lumi")
		else:
			print("Dont know what total lumi to use. default to 38.7 fb-1 (average of 16,17)")
			lumi = 38700.
	else:
		lumi = 1

	histList_MT = []
	histList_nEvents = []
	histList_Jet0Pt = []
	histList_Jet1Pt = []

	MTbins, MTstart, MTend = 50, 1000,6000
	nEbins, nEstart, nEend = 2, 0, 2
	Pt0bins, Pt0start, Pt0end = 50, 150, 3000
	Pt1bins, Pt1start, Pt1end = 50, 150, 2000

	for branchName in passedDict.keys():
		histList_MT.append(self.makeTH1F("hist_MT_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";MT;count/au",MTbins,MTstart,MTend))
		histList_nEvents.append(self.makeTH1F("hist_nEvents_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";pp;count/au",nEbins, nEstart, nEend))
		histList_Jet0Pt.append(self.makeTH1F("hist_Jet0Pt_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";Jet0Pt;count/au",Pt0bins, Pt0start, Pt0end))
		histList_Jet1Pt.append(self.makeTH1F("hist_Jet1Pt_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";Jet1Pt;count/au",Pt1bins, Pt1start, Pt1end))
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		#centeral
		
		passedDict["None"] = 1
		passedDict["BadChargedCandidateFilter"] = -3
		passedDict["BadPFMuonFilter"] = -3
		passedDict["CSCTightHaloFilter"] = 1
		passedDict["ecalBadCalibFilter"] = 1
		passedDict["EcalDeadCellTriggerPrimitiveFilter"] = -3
		passedDict["eeBadScFilter"] = -3
		passedDict["globalSuperTightHalo2016Filter"] = -3
		passedDict["globalTightHalo2016Filter"] = 1
		passedDict["HBHEIsoNoiseFilter"] = -3
		passedDict["HBHENoiseFilter"] = -3
		passedDict["PrimaryVertexFilter"] = -3
		#custom
		passedDict["METRatioFilter"] = -3
		passedDict["MuonJetFilter"] = -3
		passedDict["EcalNoiseJetFilter"] = -3
		passedDict["HTRatioFilter"] = -3
		passedDict["HTRatioTightFilter"] = -3
		passedDict["HTRatioDPhiFilter"] = -3
		passedDict["HTRatioDPhiTightFilter"] = -3
		passedDict["LowNeutralJetFilter"] = -3
		passedDict["LowNeutralJetTightFilter"] = -3
		passedDict["HEMVetoFilter"] = 1

		passedDict["BadChargedCandidateFilter"] = tree.BadChargedCandidateFilter
		passedDict["BadPFMuonFilter"] = tree.BadPFMuonFilter
		#passedDict["CSCTightHaloFilter"] = tree.CSCTightHaloFilter
		#passedDict["ecalBadCalibFilter"] = tree.ecalBadCalibFilter
		passedDict["EcalDeadCellTriggerPrimitiveFilter"] = tree.EcalDeadCellTriggerPrimitiveFilter
		passedDict["eeBadScFilter"] = tree.eeBadScFilter
		passedDict["globalSuperTightHalo2016Filter"] = tree.globalSuperTightHalo2016Filter
		#passedDict["globalTightHalo2016Filter"] = tree.globalTightHalo2016Filter
		passedDict["HBHEIsoNoiseFilter"] = tree.HBHEIsoNoiseFilter
		passedDict["HBHENoiseFilter"] = tree.HBHENoiseFilter
		passedDict["PrimaryVertexFilter"] = tree.PrimaryVertexFilter
		passedDict["METRatioFilter"] = tree.METRatioFilter
		passedDict["MuonJetFilter"] = tree.MuonJetFilter
		passedDict["EcalNoiseJetFilter"] = tree.EcalNoiseJetFilter
		passedDict["HTRatioFilter"] = tree.HTRatioFilter
		passedDict["HTRatioTightFilter"] = tree.HTRatioTightFilter
		passedDict["HTRatioDPhiFilter"] = tree.HTRatioDPhiFilter
		passedDict["HTRatioDPhiTightFilter"] = tree.HTRatioDPhiTightFilter
		passedDict["LowNeutralJetFilter"] = tree.LowNeutralJetFilter
		passedDict["LowNeutralJetTightFilter"] = tree.LowNeutralJetTightFilter
		#passedDict["HEMVetoFilter"] = tree.HEMVetoFilter
		
		nFilters = 0
		nPassed = 0
		for key, value in passedDict.items():
			if ("All" in key) or ("None" in key):
				continue
			else:
				nFilters +=1
				nPassed += value
		if nFilters == nPassed:
			passedDict["All"] = 1
		else:
			passedDict["All"] = 0 

		if "Data" in self.fileID: # need to apply Triggers to Data
			weight = 1.
			if "16" in self.fileID:
				if not ((tree.TriggerPass[10] == 1) or (tree.TriggerPass[13] == 1) or (tree.TriggerPass[103] == 1) or (tree.TriggerPass[105] == 1) or (tree.TriggerPass[106] == 1)):
					continue
			elif "17" in self.fileID:
				if not ((tree.TriggerPass[11] == 1) or (tree.TriggerPass[13] == 1) or (tree.TriggerPass[67] == 1) or (tree.TriggerPass[107] == 1)):
					continue
			else:
				exit("Data ID doesnt have year.")
		elif not (("Data" in self.fileID) or ("Jets" in self.fileID) or ("QCD" in self.fileID)): # need to apply preselection to SignalMC
			if not self.passedPreselection(tree.JetsAK8, tree.MET, tree.MT_AK8, tree.NElectrons, tree.NMuons):
				continue
			weight = 1.
		else: # no need to apply additional cuts to BackgroundMC, but do need to use weight to scale it to lumi properly
			weight = tree.Weight

		for hist in histList_MT:
			if passedDict[hist.GetName().split("_")[3]] == 1:
				hist.Fill(tree.MT_AK8, weight*lumi)
		for hist in histList_Jet0Pt:
			if passedDict[hist.GetName().split("_")[3]] == 1:
				hist.Fill(tree.JetsAK8[0].Pt(), weight*lumi)
		for hist in histList_Jet1Pt:
			if passedDict[hist.GetName().split("_")[3]] == 1:
				hist.Fill(tree.JetsAK8[1].Pt(), weight*lumi)
		for hist in histList_nEvents:
			hist.Fill(passedDict[hist.GetName().split("_")[3]], weight*lumi)

def addLoop():
	baseClass.loop = loop


