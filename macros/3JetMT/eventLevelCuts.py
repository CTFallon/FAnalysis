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
	tree.SetBranchStatus("DijetEvent",1)
	tree.SetBranchStatus("TrijetEvent",1)
	tree.SetBranchStatus("OtherjetEvent",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()

	hist_SDvar12 = self.makeTH1F("hist_SDvar12","SDvar12;;",100,0.5,1)
	hist_SDvar12_di = self.makeTH1F("hist_SDvar12_di","SDvar12_di;;",100,0.5,1)
	hist_SDvar12_tr = self.makeTH1F("hist_SDvar12_tr","SDvar12_tr;;",100,0.5,1)

	hist_SDvar13 = self.makeTH1F("hist_SDvar13","SDvar13;;",100,0.5,1)
	hist_SDvar13_di = self.makeTH1F("hist_SDvar13_di","SDvar13_di;;",100,0.5,1)
	hist_SDvar13_tr = self.makeTH1F("hist_SDvar13_tr","SDvar13_tr;;",100,0.5,1)

	hist_SDvar23 = self.makeTH1F("hist_SDvar23","SDvar23;;",100,0.5,1)
	hist_SDvar23_di = self.makeTH1F("hist_SDvar23_di","SDvar23_di;;",100,0.5,1)
	hist_SDvar23_tr = self.makeTH1F("hist_SDvar23_tr","SDvar23_tr;;",100,0.5,1)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		SDvar12 = tree.JetsAK8[0].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt())
		if len(tree.JetsAK8) > 2:
			SDvar13 = tree.JetsAK8[0].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt())
			SDvar23 = tree.JetsAK8[1].Pt()/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())
			hist_SDvar13.Fill(SDvar13)
			hist_SDvar23.Fill(SDvar23)

		hist_SDvar12.Fill(SDvar12)
		
		if tree.DijetEvent:
			hist_SDvar12_di.Fill(SDvar12)
			if len(tree.JetsAK8) > 2:
				hist_SDvar13_di.Fill(SDvar13)
				hist_SDvar23_di.Fill(SDvar23)

		if tree.TrijetEvent:
			hist_SDvar12_tr.Fill(SDvar12)
			if len(tree.JetsAK8) > 2:
				hist_SDvar13_tr.Fill(SDvar13)
				hist_SDvar23_tr.Fill(SDvar23)

	
	makePlots(hist_SDvar12, hist_SDvar13, hist_SDvar23 ,self.extraDir)
	makePlots(hist_SDvar12_di, hist_SDvar13_di, hist_SDvar23_di ,self.extraDir)
	makePlots(hist_SDvar12_tr, hist_SDvar13_tr, hist_SDvar23_tr ,self.extraDir)
	makePlots(hist_SDvar12, hist_SDvar12_di, hist_SDvar12_tr ,self.extraDir,extraName='2')
	makePlots(hist_SDvar13, hist_SDvar13_di, hist_SDvar13_tr ,self.extraDir)
	makePlots(hist_SDvar23, hist_SDvar23_di, hist_SDvar23_tr ,self.extraDir)

def addLoop():
	baseClass.loop = loop

def makePlots(All, ISR, FSR, direc, logY = True,extraName=''):
	c = rt.TCanvas("canvas","canvas",900,600)
	if logY:
		c.SetLogy()
	histStack = rt.THStack()
	histStack.SetTitle(All.GetTitle())
	All.SetLineColor(1)
	ISR.SetLineColor(2)
	FSR.SetLineColor(3)
	histStack.Add(All)
	histStack.Add(ISR)
	histStack.Add(FSR)
	histStack.Draw("NOSTACK")
	histStack.GetXaxis().SetTitle(All.GetXaxis().GetTitle())
	c.BuildLegend(0.75,0.75,0.9,0.9,"")
	c.SaveAs(direc+All.GetTitle()+extraName+".png")


