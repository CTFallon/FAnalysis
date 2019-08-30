from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

def loop(self):
	# set up trees
	tree = self.getChain(self.treeNameList[0])
	# added friend tree
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))

	if not ("Data" in self.fileID): # only need to do this for MC
		tree.SetBranchStatus("Weight",1)
		if "16" in self.fileID:
			lumi = 35900.
			print("2016 Lumi")
		else: # signal samples have 2017 lumi
			lumi = 41500.
			print("2017 Lumi")
	else:
		lumi = 1
	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()

	# looking to opt. minDeltaPhi(MET,j1/2), and MET/MT

	hist_MT = self.makeTH1F("hist_MT","MT;MT;Count", 100, 0, 5000)
	nFail = 0
	nPass = 0

	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		# filter cuts:
		if not ((tree.globalSuperTightHalo2016Filter==1)and(tree.HBHENoiseFilter==1)and(tree.HBHEIsoNoiseFilter==1)and(tree.eeBadScFilter==1)and(tree.EcalDeadCellTriggerPrimitiveFilter==1)and(tree.BadChargedCandidateFilter==1)and(tree.BadPFMuonFilter==1)and(tree.NVtx > 0)):
			continue
		if not ((tree.METRatioFilter==1)and(tree.MuonJetFilter==1)and(tree.HTRatioDPhiTightFilter==1)and(tree.LowNeutralJetFilter==1)):
			continue
		if not ("Data" in self.fileID):
			weight = tree.Weight
		else:
			weight = 1
		if (tree.DeltaPhiMin_AK8 < 0.3) and (float(tree.MET)/float(tree.MT_AK8) > 0.2):
			nPass += 1
			hist_MT.Fill(float(tree.MT_AK8), weight*lumi)
		else:
			nFail += 1

	print(nPass, nFail)

def addLoop():
	baseClass.loop = loop


