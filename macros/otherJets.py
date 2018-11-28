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
	# look at delta jet varialbes for jets that arn't maxDPhi
	# delta R, delta Phi, delta eta, delta Pt
	
	hist_deltaR = self.makeTH1F("hist_deltaR",
		"DeltaR Between Non-MaxDPhiMet Jets",
		100,0.5,6)
	hist_deltaR_2FSR = self.makeTH1F("hist_deltaR_2FSR",
		"DeltaR Between Non-MaxDPhiMet Jets, 2FSR",
		100,0.5,6)
	hist_deltaR_3FSR = self.makeTH1F("hist_deltaR_3FSR",
		"DeltaR Between Non-MaxDPhiMet Jets, 3FSR",
		100,0.5,6)
	hist_deltaR_4FSR = self.makeTH1F("hist_deltaR_4FSR",
		"DeltaR Between Non-MaxDPhiMet Jets, 4FSR",
		100,0.5,6)

	hist_deltaPhi = self.makeTH1F("hist_deltaPhi",
		"DeltaPhi Between Non-MaxDPhiMet Jets",
		100,0,rt.TMath.Pi())
	hist_deltaPhi_2FSR = self.makeTH1F("hist_deltaPhi_2FSR",
		"DeltaPhi Between Non-MaxDPhiMet Jets, 2FSR",
		100,0,rt.TMath.Pi())
	hist_deltaPhi_3FSR = self.makeTH1F("hist_deltaPhi_3FSR",
		"DeltaPhi Between Non-MaxDPhiMet Jets, 3FSR",
		100,0,rt.TMath.Pi())
	hist_deltaPhi_4FSR = self.makeTH1F("hist_deltaPhi_4FSR",
		"DeltaPhi Between Non-MaxDPhiMet Jets, 4FSR",
		100,0,rt.TMath.Pi())

	hist_deltaEta = self.makeTH1F("hist_deltaEta",
		"DeltaEta Between Non-MaxDPhiMet Jets",
		100,0,5)
	hist_deltaEta_2FSR = self.makeTH1F("hist_deltaEta_2FSR",
		"DeltaEta Between Non-MaxDPhiMet Jets, 2FSR",
		100,0,5)
	hist_deltaEta_3FSR = self.makeTH1F("hist_deltaEta_3FSR",
		"DeltaEta Between Non-MaxDPhiMet Jets, 3FSR",
		100,0,5)
	hist_deltaEta_4FSR = self.makeTH1F("hist_deltaEta_4FSR",
		"DeltaEta Between Non-MaxDPhiMet Jets, 4FSR",
		100,0,5)

	hist_deltaPt = self.makeTH1F("hist_deltaPt",
		"DeltaPt Between Non-MaxDPhiMet Jets",
		100,0,2000)
	hist_deltaPt_2FSR = self.makeTH1F("hist_deltaPt_2FSR",
		"DeltaPt Between Non-MaxDPhiMet Jets, 2FSR",
		100,0,2000)
	hist_deltaPt_3FSR = self.makeTH1F("hist_deltaPt_3FSR",
		"DeltaPt Between Non-MaxDPhiMet Jets, 3FSR",
		100,0,2000)
	hist_deltaPt_4FSR = self.makeTH1F("hist_deltaPt_4FSR",
		"DeltaPt Between Non-MaxDPhiMet Jets, 4FSR",
		100,0,2000)

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		jets = tree.JetsAK8
		nJets = len(jets)
		#if nJets == 2:
		#	continue

		nFSR = 0
		for iJet in range(nJets):
			if not tree.JetsAK8_isISR[iJet]:
				nFSR += 1

		iJetMaxDPhi = tree.iJetMaxDeltaPhi
		jI = [x for x in range(nJets)]
		jI.remove(iJetMaxDPhi)
		

		for i in jI:
			for j in jI:
				if i <= j:
					continue
				dR = (jets[j].DeltaR(jets[i]))
				dP = abs(jets[j].DeltaPhi(jets[i]))
				dE = abs(jets[j].Eta() - jets[i].Eta())
				dPt = (jets[j].Pt() - jets[i].Pt())
				hist_deltaR.Fill(dR)
				hist_deltaPhi.Fill(dP)
				hist_deltaEta.Fill(dE)
				hist_deltaPt.Fill(dPt)
				if nFSR == 2:
					hist_deltaR_2FSR.Fill(dR)
					hist_deltaPhi_2FSR.Fill(dP)
					hist_deltaEta_2FSR.Fill(dE)
					hist_deltaPt_2FSR.Fill(dPt)
				elif nFSR == 3:
					hist_deltaR_3FSR.Fill(dR)
					hist_deltaPhi_3FSR.Fill(dP)
					hist_deltaEta_3FSR.Fill(dE)
					hist_deltaPt_3FSR.Fill(dPt)
				else:
					hist_deltaR_4FSR.Fill(dR)
					hist_deltaPhi_4FSR.Fill(dP)
					hist_deltaEta_4FSR.Fill(dE)
					hist_deltaPt_4FSR.Fill(dPt)
		
	makePlots([hist_deltaR,
			hist_deltaR_2FSR,
			hist_deltaR_3FSR,
			hist_deltaR_4FSR],self.extraDir,"deltaR_otherJets",log=True)
	makePlots([hist_deltaPhi,
			hist_deltaPhi_2FSR,
			hist_deltaPhi_3FSR,
			hist_deltaPhi_4FSR],self.extraDir,"deltaPhi_otherJets",log=True,
			lC = [0.0,0.7,0.3,1.0])
	makePlots([hist_deltaEta,
			hist_deltaEta_2FSR,
			hist_deltaEta_3FSR,
			hist_deltaEta_4FSR],self.extraDir,"deltaEta_otherJets",log=True)
	makePlots([hist_deltaPt,
			hist_deltaPt_2FSR,
			hist_deltaPt_3FSR,
			hist_deltaPt_4FSR],self.extraDir,"deltaPt_otherJets",log=True)
	

def addLoop():
	baseClass.loop = loop

def makePlots(hList,direct, name, log = False, lC = [0.7,0.7,1.0,1.0]):
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

