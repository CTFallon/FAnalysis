from analysisBase import baseClass
import ROOT as rt
from array import array

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
	#tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	#tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	#tree.SetBranchStatus("fracPtFromHVQuarks",1)
	#tree.SetBranchStatus("numHVPartsInJet",1)
	#tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	#tree.SetBranchStatus("dPhiMaxDeltaPhi",1)
	#tree.SetBranchStatus("zPrimept",1)
	#tree.SetBranchStatus("zPrimephi",1)
	#tree.SetBranchStatus("pGJ_visible",1)
	#tree.SetBranchStatus("pGJ_invis",1)
	#tree.SetBranchStatus("pGJ_every",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	# create histgrams that show the deltaR distribution between all 3 pairs of jets
	# based on if the jetPtMaxdPhi cut selects the event for 2 or 3 jet MT

	hist_deltaR12_3jetMT = self.makeTH1F("hist_deltaR12_3jetMT","3 Jet MT (12);deltaR; count/a.u.", 100, 0, 6)
	hist_deltaR13_3jetMT = self.makeTH1F("hist_deltaR13_3jetMT","3 Jet MT (13);deltaR; count/a.u.", 100, 0, 6)
	hist_deltaR23_3jetMT = self.makeTH1F("hist_deltaR23_3jetMT","3 Jet MT (23);deltaR; count/a.u.", 100, 0, 6)

	hist_deltaR12_2jetMT = self.makeTH1F("hist_deltaR12_2jetMT","2 Jet MT (12);deltaR; count/a.u.", 100, 0, 6)
	hist_deltaR13_2jetMT = self.makeTH1F("hist_deltaR13_2jetMT","2 Jet MT (13);deltaR; count/a.u.", 100, 0, 6)
	hist_deltaR23_2jetMT = self.makeTH1F("hist_deltaR23_2jetMT","2 Jet MT (23);deltaR; count/a.u.", 100, 0, 6)

	hist_sumdR_3jetMT = self.makeTH1F("hist_sumdR_3jetMT","3 Jet MT (sumDR);\Sigma \Delta R;count/a.u.",100,2,13)
	hist_sumdR_2jetMT = self.makeTH1F("hist_sumdR_2jetMT","2 Jet MT (sumDR);\Sigma \Delta R;count/a.u.",100,2,13)

	hist_3jetArea_3jetMT = self.makeTH1F("hist_3jetArea_3jetMT","3JetArea(3jetMT);Area;count/a.u.",100,0,7)
	hist_3jetArea_2jetMT = self.makeTH1F("hist_3jetArea_2jetMT","3JetArea(2jetMT);Area;count/a.u.",100,0,7)

	hist_deltadeltR_12_13_3jetMT = self.makeTH1F("hist_deltadeltR_12_13_3jetMT","ugh;ddR;count",100,-6,6)
	hist_deltadeltR_12_23_3jetMT = self.makeTH1F("hist_deltadeltR_12_23_3jetMT","ugh;ddR;count",100,-6,6)
	hist_deltadeltR_13_23_3jetMT = self.makeTH1F("hist_deltadeltR_13_23_3jetMT","ugh;ddR;count",100,-6,6)

	hist_deltadeltR_12_13_2jetMT = self.makeTH1F("hist_deltadeltR_12_13_2jetMT","ugh;ddR;count",100,-6,6)
	hist_deltadeltR_12_23_2jetMT = self.makeTH1F("hist_deltadeltR_12_23_2jetMT","ugh;ddR;count",100,-6,6)
	hist_deltadeltR_13_23_2jetMT = self.makeTH1F("hist_deltadeltR_13_23_2jetMT","ugh;ddR;count",100,-6,6)

	hist_jet1Pt_3jetMT = self.makeTH1F("hist_jet1Pt_3jetMT","jetPt;pT;count", 100, 0, 2500)
	hist_jet2Pt_3jetMT = self.makeTH1F("hist_jet2Pt_3jetMT","jetPt;pT;count", 100, 0, 2500)
	hist_jet3Pt_3jetMT = self.makeTH1F("hist_jet3Pt_3jetMT","jetPt;pT;count", 100, 0, 2500)

	hist_jet1Pt_2jetMT = self.makeTH1F("hist_jet1Pt_2jetMT","jetPt;pT;count", 100, 0, 2500)
	hist_jet2Pt_2jetMT = self.makeTH1F("hist_jet2Pt_2jetMT","jetPt;pT;count", 100, 0, 2500)
	hist_jet3Pt_2jetMT = self.makeTH1F("hist_jet3Pt_2jetMT","jetPt;pT;count", 100, 0, 2500)

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		if nJets >= 3:
			dr12 = tree.JetsAK8[0].DeltaR(tree.JetsAK8[1])
			dr13 = tree.JetsAK8[0].DeltaR(tree.JetsAK8[2])
			dr23 = tree.JetsAK8[1].DeltaR(tree.JetsAK8[2])
			if tree.JetsAK8[tree.iJetMaxDeltaPhi].Pt() < 860.0: # soft jet, use 3
				hist_deltaR12_3jetMT.Fill(dr12)
				hist_deltaR13_3jetMT.Fill(dr13)
				hist_deltaR23_3jetMT.Fill(dr23)
				hist_sumdR_3jetMT.Fill(dr12+dr13+dr23)
				hist_3jetArea_3jetMT.Fill(HeronsFormula(dr12, dr13, dr23))
				hist_deltadeltR_12_13_3jetMT.Fill(dr12-dr13)
				hist_deltadeltR_12_23_3jetMT.Fill(dr12-dr23)
				hist_deltadeltR_13_23_3jetMT.Fill(dr13-dr23)
				hist_jet1Pt_3jetMT.Fill(tree.JetsAK8[0].Pt())
				hist_jet2Pt_3jetMT.Fill(tree.JetsAK8[1].Pt())
				hist_jet3Pt_3jetMT.Fill(tree.JetsAK8[2].Pt())

			else: # hard jet, use 2 
				hist_deltaR12_2jetMT.Fill(dr12)
				hist_deltaR13_2jetMT.Fill(dr13)
				hist_deltaR23_2jetMT.Fill(dr23)
				hist_sumdR_2jetMT.Fill(dr12+dr13+dr23)
				hist_3jetArea_2jetMT.Fill(HeronsFormula(dr12, dr13, dr23))
				hist_deltadeltR_12_13_2jetMT.Fill(dr12-dr13)
				hist_deltadeltR_12_23_2jetMT.Fill(dr12-dr23)
				hist_deltadeltR_13_23_2jetMT.Fill(dr13-dr23)
				hist_jet1Pt_2jetMT.Fill(tree.JetsAK8[0].Pt())
				hist_jet2Pt_2jetMT.Fill(tree.JetsAK8[1].Pt())
				hist_jet3Pt_2jetMT.Fill(tree.JetsAK8[2].Pt())
						


	

def addLoop():
	baseClass.loop = loop

def HeronsFormula(a,b,c):
	# formula to calculate area of a triangle using lengths of the three sides
	p = (a + b + c)/2
	return rt.TMath.Sqrt(p*(p-a)*(p-b)*(p-c))


