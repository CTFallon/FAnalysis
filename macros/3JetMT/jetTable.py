from analysisBase import baseClass
import ROOT as rt
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

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
	hist_ISR_nHVParts = self.makeTH1F("hist_ISR_nHVParts","nHVParts of ISR Jets;;",50,0,50)
	hist_FSR_nHVParts = self.makeTH1F("hist_FSR_nHVParts","nHVParts of FSR Jets;;",50,0,50)

	hist_nJets_HVFrac  = self.makeTH1F("hist_nJets_HVFrac","Jets per Event with HV Frac;;",10,0,10)
	hist_nJets_FSR = self.makeTH1F("hist_nJets_FSR","FSR Jets per Event;;",10,0,10)
	hist_nJets_HVParts = self.makeTH1F("hist_nJets_HVParts","Jets per Event with HV parts;;",10,0,10)

	nEventsPP = 0
	nEventsWith3JetsThatHaveHVConstiuents = 0
	nEventsWith3JetsThatAreTaggedAsFSR = 0
	nEventsWith3JetsThatHaveHVParts = 0
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nEventsPP += 1
		nJets = len(tree.JetsAK8)
		nJetsWithHV = 0
		nJetsFSR = 0
		nJetsWithHVParts = 0
		for iJet in range(nJets):
			if tree.fracTotPTfromAllHVQ[iJet] > 0.0:
				nJetsWithHV += 1
			if not tree.JetsAK8_isISR[iJet]:
				nJetsFSR += 1
			if tree.numHVPartsInJet[iJet] > 0:
				nJetsWithHVParts += 1
			if tree.JetsAK8_isISR[iJet]:
				hist_ISR_nHVParts.Fill(tree.numHVPartsInJet[iJet])	
			elif not tree.JetsAK8_isISR[iJet]:
				hist_FSR_nHVParts.Fill(tree.numHVPartsInJet[iJet])
		hist_nJets_HVFrac.Fill(nJetsWithHV)
		hist_nJets_FSR.Fill(nJetsFSR)
		hist_nJets_HVParts.Fill(nJetsWithHVParts)
		if nJetsWithHV >= 3:
			nEventsWith3JetsThatHaveHVConstiuents += 1
		if nJetsFSR >=3 :
			nEventsWith3JetsThatAreTaggedAsFSR += 1
		if nJetsWithHVParts >= 3:
			nEventsWith3JetsThatHaveHVParts += 1

	print("NEvents pass PreSelection {}".format(snEventsPP))
	print("NEvents with 3 or more jets that have HV Frac Pt {}".format(nEventsWith3JetsThatHaveHVConstiuents))
	print("NEvents with 3 or more jets that are tagged as FSR {}".format(nEventsWith3JetsThatAreTaggedAsFSR))
	print("NEvents with 3 more more jets that have HV parts {}".format(nEventsWith3JetsThatHaveHVParts))

	c1 = rt.TCanvas("c1","c1",900,600)
	c1.SetLogy()
	hist_ISR_nHVParts.SetLineColor(2)
	hist_FSR_nHVParts.Draw()
	hist_ISR_nHVParts.Draw('same')
	c1.BuildLegend()
	c1.SaveAs(self.extraDir+"Jet_nHVParts.png")

	hist_nJets_HVFrac.Draw()
	hist_nJets_FSR.SetLineColor(2)
	hist_nJets_FSR.Draw('same')
	hist_nJets_HVParts.SetLineColor(3)
	hist_nJets_HVParts.Draw('same')
	c1.BuildLegend()
	c1.SaveAs(self.extraDir+"nJets.png")
		

	

def addLoop():
	baseClass.loop = loop


