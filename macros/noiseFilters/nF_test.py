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
	passedDict["All"] = 1
	for branch in branchList:
		if "Filter" in branch.GetName():
			passedDict[branch.GetName()] = -3
		else:
			tree.SetBranchStatus(branch.GetName(),0)

	tree.SetBranchStatus("JetsAK8",1)
	tree.SetBranchStatus("MT_AK8",1)
	tree.SetBranchStatus("NElectrons",1)
	tree.SetBranchStatus("NMuons",1)
	tree.SetBranchStatus("MET",1)

	if not ("Data" in self.fileID):
		tree.SetBranchStatus("Weight",1)
		if "16" in self.fileID:
			lumi = 35900.
		elif "17" in self.fileID:
			lumi = 41500.
		else:
			print("Dont know what total lumi to use. default to 38.7 fb-1 (average of 16,17)")
			lumi = 38700.
	else:
		lumi = 1

	# values of interest
		# Percent and/or total number of events effected
		# MT resolution before and after filters

	histList_MT = []
	histList_nEvents = []

	MTbins, MTstart, MTend = 50, 1000,4000
	nEbins, nEstart, nEend = 2, 0, 2

	for branchName in passedDict.keys():
		histList_MT.append(self.makeTH1F("hist_MT_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";MT;count/au",MTbins,MTstart,MTend))
		histList_nEvents.append(self.makeTH1F("hist_nEvents_"+self.fileID+"_"+branchName,branchName+" "+self.fileID+";pp;count/au",nEbins, nEstart, nEend))
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		#centeral
		
		passedDict["All"] = 1
		passedDict["BadChargedCandidateFilter"] = -3
		passedDict["BadPFMuonFilter"] = -3
		#passedDict["CSCTightHaloFilter"] = -3
		#passedDict["ecalBadCalibFilter"] = -3
		passedDict["EcalDeadCellTriggerPrimitiveFilter"] = -3
		passedDict["eeBadScFilter"] = -3
		passedDict["globalSuperTightHalo2016Filter"] = -3
		#passedDict["globalTightHalo2016Filter"] = -3
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
		#passedDict["HEMVetoFilter"] = -3

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

		if "Data" in self.fileID:
			weight = 1.
		else:
			weight = tree.Weight

		for hist in histList_MT:
			if passedDict[hist.GetName().split("_")[3]] == 1:
				hist.Fill(tree.MT_AK8, weight*lumi)
		for hist in histList_nEvents:
			hist.Fill(passedDict[hist.GetName().split("_")[3]], weight*lumi)

def addLoop():
	baseClass.loop = loop


