from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

tdrstyle.setTDRStyle()

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	# adding friend tree
	ff = rt.TFile.Open(self.extraDir+"outFriend.root")
	friendTree = ff.Get("friend")
	tree.AddFriend(friendTree)
	# added friend tree
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))
	
	# Turn off all branches, selective turn on branches
	tree.SetBranchStatus("*", 0)
	tree.SetBranchStatus("RunNum",1)
	tree.SetBranchStatus("EvtNum",1)
	tree.SetBranchStatus("*AK8*", 1)
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	#tree.SetBranchStatus("fracPtFromHVQuarks",1)
	tree.SetBranchStatus("numHVPartsInJet",1)
	tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)
	#tree.SetBranchStatus("zPrimept",1)
	#tree.SetBranchStatus("zPrimephi",1)
	tree.SetBranchStatus("pGJ_visible",1)
	tree.SetBranchStatus("pGJ_invis",1)
	tree.SetBranchStatus("pGJ_every",1)
	tree.SetBranchStatus("fracVisPTfromVisHVQ",1)
	tree.SetBranchStatus("fracInvPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromAllHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVisHVQ",1)
	tree.SetBranchStatus("fracTotPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVis",1)
	tree.SetBranchStatus("fracVisHVQtoInvHVQ",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	# Compare distributions for:
	hist_MinDeltaR_MDP_others = self.makeTH1F(
		"hist_MinDeltaR_MDP_others",
		"Minimum Delta R Between Jet(MaxDPhi) and all other Jets;DeltaR;",
		100,0,6)
	hist_MinDeltaR_MDP_others_2FSR = self.makeTH1F(
		"hist_MinDeltaR_MDP_others_2FSR",
		"MinDeltaR(Jet(MaxDPhi), All Other Jets) 2 FSR Jets;DeltaR;",
		100,0,6)
	hist_MinDeltaR_MDP_others_3FSR = self.makeTH1F(
		"hist_MinDeltaR_MDP_others_3FSR",
		"MinDeltaR(Jet(MaxDPhi), All Other Jets) 3 FSR Jets;DeltaR;",
		100,0,6)
	hist_MinDeltaR_MDP_others_4FSR = self.makeTH1F(
		"hist_MinDeltaR_MDP_others_4FSR",
		"MinDeltaR(Jet(MaxDPhi), All Other Jets) 4+ FSR Jets;DeltaR;",
		100,0,6)

	hist_MinDeltaPhi_MDP_others = self.makeTH1F(
		"hist_MinDeltaPhi_MDP_others",
		"Minimum Delta Phi Between Jet(MaxDPhi) and all other Jets;DeltaPhi;",
		100,0,3.2)
	hist_MinDeltaPhi_MDP_others_2FSR = self.makeTH1F(
		"hist_MinDeltaPhi_MDP_others_2FSR",
		"MinDeltaPhi(Jet(MaxDPhi), All Other Jets) 2 FSR Jets;DeltaPhi;",
		100,0,3.2)
	hist_MinDeltaPhi_MDP_others_3FSR = self.makeTH1F(
		"hist_MinDeltaPhi_MDP_others_3FSR",
		"MinDeltaPhi(Jet(MaxDPhi), All Other Jets) 3 FSR Jets;DeltaPhi;",
		100,0,3.2)
	hist_MinDeltaPhi_MDP_others_4FSR = self.makeTH1F(
		"hist_MinDeltaPhi_MDP_others_4FSR",
		"MinDeltaPhi(Jet(MaxDPhi), All Other Jets) 4+ FSR Jets;DeltaPhi;",
		100,0,3.2)

	hist_minSDV_MDP_others = self.makeTH1F(
		"hist_SDV_MDP_others",
		"SDV(Jet(MaxDPhi), All Other Jets);SDV;",
		100,0,0.34)
	hist_minSDV_MDP_others_2FSR = self.makeTH1F(
		"hist_SDV_MDP_others_2FSR",
		"SDV(Jet(MaxDPhi), All Other Jets) 2 FSR Jets;SDV;",
		100,0,0.34)
	hist_minSDV_MDP_others_3FSR = self.makeTH1F(
		"hist_SDV_MDP_others_3FSR",
		"SDV(Jet(MaxDPhi), All Other Jets) 3 FSR Jets;SDV;",
		100,0,0.34)
	hist_minSDV_MDP_others_4FSR = self.makeTH1F(
		"hist_SDV_MDP_others_4FSR",
		"SDV(Jet(MaxDPhi), All Other Jets) 4+ FSR Jets;SDV;",
		100,0,0.34)



	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		
		if nJets == 2:
			continue
		#	deltaPhiMetList = [tree.DeltaPhi1, tree.DeltaPhi2]
		else:
			deltaPhiMetList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
		jetMDPIndex = deltaPhiMetList.index(max(deltaPhiMetList))
		jetPtMDP = tree.JetsAK8[jetMDPIndex].Pt()
		jetMDP = tree.JetsAK8[jetMDPIndex]

		jetDeltaRother = []
		jetDeltaPhiother = []
		sdvOther = []
		ptList = []
		nFSR = 0
		for i in range(nJets):
			ptList.append(tree.JetsAK8[i].Pt())
			if not tree.JetsAK8_isISR[i]:
				nFSR += 1
			if i != jetMDPIndex:
				jetDeltaRother.append(jetMDP.DeltaR(tree.JetsAK8[i]))
				jetDeltaPhiother.append(abs(jetMDP.DeltaPhi(tree.JetsAK8[i])))
				sdvOther.append(SDV([jetPtMDP,tree.JetsAK8[i].Pt()]))

		
		hist_MinDeltaR_MDP_others.Fill(min(jetDeltaRother))
		hist_MinDeltaPhi_MDP_others.Fill(min(jetDeltaPhiother))
		hist_minSDV_MDP_others.Fill(SDV(ptList))
		if nFSR == 2:
			hist_MinDeltaR_MDP_others_2FSR.Fill(min(jetDeltaRother))
			hist_MinDeltaPhi_MDP_others_2FSR.Fill(min(jetDeltaPhiother))
			hist_minSDV_MDP_others_2FSR.Fill(SDV(ptList))
		if nFSR == 3:
			hist_MinDeltaR_MDP_others_3FSR.Fill(min(jetDeltaRother))
			hist_MinDeltaPhi_MDP_others_3FSR.Fill(min(jetDeltaPhiother))
			hist_minSDV_MDP_others_3FSR.Fill(SDV(ptList))
		if nFSR >= 4:
			hist_MinDeltaR_MDP_others_4FSR.Fill(min(jetDeltaRother))
			hist_MinDeltaPhi_MDP_others_4FSR.Fill(min(jetDeltaPhiother))
			hist_minSDV_MDP_others_4FSR.Fill(SDV(ptList))

	makePlots([
		hist_MinDeltaR_MDP_others,
		hist_MinDeltaR_MDP_others_2FSR,
		hist_MinDeltaR_MDP_others_3FSR,
		hist_MinDeltaR_MDP_others_4FSR],
		self.extraDir,
		"DeltaR_MDPJet_Others_nFSRJets",
		log = True)
	makePlots([
		hist_MinDeltaPhi_MDP_others,
		hist_MinDeltaPhi_MDP_others_2FSR,
		hist_MinDeltaPhi_MDP_others_3FSR,
		hist_MinDeltaPhi_MDP_others_4FSR],
		self.extraDir,
		"DeltaPhi_MDPJet_Others_nFSRJets",
		log = True)
	makePlots([
		hist_minSDV_MDP_others,
		hist_minSDV_MDP_others_2FSR,
		hist_minSDV_MDP_others_3FSR,
		hist_minSDV_MDP_others_4FSR],
		self.extraDir,
		"minSDV_MDPJet_Others_nFSRJets",
		log = True)

def addLoop():
	baseClass.loop = loop

def makePlots(hList,direct, name, log = False, lC = [0.3,0.2,0.7,0.3]):
	c1 = rt.TCanvas("c1","c1",900,600)
	if log:
		c1.SetLogy()
	for iH in range(len(hList)):
		if iH == 0:
			hList[iH].SetLineColor(1)
			hList[iH].Draw()
		else:
			hList[iH].SetLineColor(iH+1)
			hList[iH].Draw('same')
	c1.BuildLegend(lC[0],lC[1],lC[2],lC[3])
	c1.SaveAs(direct+name+".png")

def SDV(jetPtList):
	minPt = min(jetPtList)
	return (minPt)/(sum(jetPtList))


