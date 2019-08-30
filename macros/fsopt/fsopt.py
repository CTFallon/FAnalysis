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

	dPbins, dPmin, dPmax = 50, 0, rt.TMath.Pi()
	mRbins, mRmin, mRmax = 50, 0.15, 0.65
	hist_dPhi = self.makeTH1F("hist_dPhi","minDeltaPhi;deltaPhi;Count",dPbins, dPmin, dPmax)
	hist_metR = self.makeTH1F("hist_metR","MET/MT;MET/MT;Count",mRbins, mRmin, mRmax)

	hist2d_dPhi_metR = self.makeTH2F("hist2d_dPhi_metR","2D Hist;dPhi;metR",dPbins, dPmin, dPmax,mRbins, mRmin, mRmax)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		# TT Stitching
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName(), tree.madHT, tree.GenMET, len(tree.GenElectrons)+len(tree.GenMuons)+len(tree.GenTaus)): continue
		# filter cuts:
		if not ((tree.globalSuperTightHalo2016Filter==1)and(tree.HBHENoiseFilter==1)and(tree.HBHEIsoNoiseFilter==1)and(tree.eeBadScFilter==1)and(tree.EcalDeadCellTriggerPrimitiveFilter==1)and(tree.BadChargedCandidateFilter==1)and(tree.BadPFMuonFilter==1)and(tree.NVtx > 0)):
			continue
		if not ((tree.METRatioFilter==1)and(tree.MuonJetFilter==1)and(tree.HTRatioDPhiTightFilter==1)and(tree.LowNeutralJetFilter==1)):
			continue
		if not ("Data" in self.fileID):
			weight = tree.Weight
		else:
			weight = 1
		if self.fileID[0] == 'z':
			lowMTlim = int(self.fileID[1:]) * 100.0 * 2./3.
			highMTlim = int(self.fileID[1:]) * 100.0 * 4./3.
		else:
			lowMTlim = 2000.0
			highMTlim = 4000.0
		if lowMTlim < tree.MT_AK8  and  tree.MT_AK8 < highMTlim:
			hist_dPhi.Fill(tree.DeltaPhiMin_AK8, weight*lumi)
			hist_metR.Fill(float(tree.MET)/float(tree.MT_AK8), weight*lumi)
			hist2d_dPhi_metR.Fill(tree.DeltaPhiMin_AK8,float(tree.MET)/float(tree.MT_AK8), weight*lumi)
	

def addLoop():
	baseClass.loop = loop


