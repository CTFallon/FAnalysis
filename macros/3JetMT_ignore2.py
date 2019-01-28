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

	nDiT = 0
	nTriT = 0
	nDiCutA = 0
	nTriCutA = 0
	nDiCutB = 0
	nTriCutB = 0

	nDiR = 0
	nTriR = 0

	hist_iJetMaxDeltaPhi_dijet = self.makeTH1F("hist_iJetMaxDeltaPhi_dijet","Jet Index - Dijet MT;Jet Index;Count", 5, 0, 5)
	hist_iJetMaxDeltaPhi_trijet= self.makeTH1F("hist_iJetMaxDeltaPhi_trijet","Jet Index - Trijet MT;Jet Index;Count", 5, 0, 5)

	hist_deltaPhi12_dijet = self.makeTH1F("hist_deltaPhi12_dijet","Delta Phi 12 - Dijet MT;dPhi12;Count", 100, 0, rt.TMath.Pi())
	hist_deltaPhi12_trijet= self.makeTH1F("hist_deltaPhi12_trijet","Delta Phi 12 - Trijet MT;dPhi12;Count", 100, 0, rt.TMath.Pi())


	testVarName = "deltaPhi23"
	# failed test variables for cut C: deltaPhi23, deltaPhi13, deltaPt2MDP, pt2/(pt2_mdppt), ptmdp, pt3, pt2, jet eta (123)
	# list cont: 
	# pt1 might be useful, but probably just seperating events based on which shower split
	testVarLow = 0
	testVarHig = rt.TMath.Pi()

	hist_test_dijet = self.makeTH1F("hist_test_dijet",testVarName+" - Dijet MT;"+testVarName+";Count", 50, testVarLow, testVarHig)
	hist_test_trijet = self.makeTH1F("hist_test_trijet",testVarName+" - Trijet MT;"+testVarName+";Count", 50, testVarLow, testVarHig)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		if len(tree.JetsAK8) <= 2:
			continue

		if bool(tree.JetsAK8_isHV[2]) == True:
			isTri = True
			nTriT +=1
		else:
			nDiT += 1
			isTri = False

		if isTri:
			hist_iJetMaxDeltaPhi_trijet.Fill(tree.iJetMaxDeltaPhi)
		else:
			hist_iJetMaxDeltaPhi_dijet.Fill(tree.iJetMaxDeltaPhi)
		if tree.iJetMaxDeltaPhi == 1:
			if isTri:
				nTriCutA += 1
			else:
				nDiCutA += 1
			continue
		else:
			deltaPhi12 = deltaPhi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi())
			if isTri:
				hist_deltaPhi12_trijet.Fill(deltaPhi12)
			else:
				hist_deltaPhi12_dijet.Fill(deltaPhi12)
			
		if deltaPhi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()) < 2.65:
			if isTri:
				nTriCutB += 1
			else:
				nDiCutB += 1
			continue

		
		testVar = deltaPhi(tree.JetsAK8[1].Phi(),tree.JetsAK8[2].Phi())
		if isTri:
			nTriR += 1
			hist_test_trijet.Fill(testVar)
		else:
			nDiR += 1
			hist_test_dijet.Fill(testVar)
			
	print("iJetMaxDeltaPhi")
	for iBin in range(1,5):
		print("{} {} {}".format(iBin,hist_iJetMaxDeltaPhi_dijet.GetBinContent(iBin),hist_iJetMaxDeltaPhi_trijet.GetBinContent(iBin)))
	print("deltaPhi Jets 1 and 2")
	for iBin in range(1,100):
		print("{} {} {}".format(iBin,hist_deltaPhi12_dijet.Integral(0,iBin),hist_deltaPhi12_trijet.Integral(0,iBin)))
	print(" Cut | nDi | nTri | Total")
	print("  T  | {} | {} | {}".format(nDiT, nTriT, nDiT+nTriT))
	print("  A  | {} | {} | {}".format(nDiCutA, nTriCutA, nDiCutA+nTriCutA))
	print("  B  | {} | {} | {}".format(nDiCutB, nTriCutB, nDiCutB+nTriCutB))
	print("  R  | {} | {} | {}".format(nDiR, nTriR, nDiR+nTriR))

	self.makePng([hist_iJetMaxDeltaPhi_dijet,hist_iJetMaxDeltaPhi_trijet],"iJetMaxDeltaPhi")
	self.makePng([hist_deltaPhi12_dijet,hist_deltaPhi12_trijet],"deltaPhi12")
	self.makePng([hist_test_dijet,hist_test_trijet],testVarName)
	
	

def addLoop():
	baseClass.loop = loop

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def deltaPhi(phi1, phi2):
	delta = phi1 - phi2
	while delta > rt.TMath.Pi():
		delta -= 2*rt.TMath.Pi()
	while delta < -rt.TMath.Pi():
		delta += 2*rt.TMath.Pi()
	return abs(delta)

