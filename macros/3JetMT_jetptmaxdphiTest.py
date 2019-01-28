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

	#histograms of njets ==2, maxdphi index and maxdphiPt
	hist_2Jets_mDPhiIndex = self.makeTH1F("hist_2Jets_mDPhiIndex",
		"Index of MaxDPhi Jet, events with only 2 jets;Index;Count",2,0,2)
	hist_2Jets_ptMaxDPhi = self.makeTH1F("hist_2Jets_ptMaxDPhi",
		"Pt of MaxDPhi Jet, events with only 2 jets;pT;count",100,0,2000)

	hist_2Jets_jetPt = self.makeTH1F("hist_2Jets_jetPt",
		"Jet Pt of all HV Jets, events with only 2 jets; pT; count",
		100,0,2000)
		
	hist_3Jets_jetPt_HV = self.makeTH1F("hist_3Jets_jetPt_HV",
		"Jet Pt of all HV jets in events with 3 or more jets;pt;count",
		100,0,2000)

	hist_3Jets_jetPt_notHV = self.makeTH1F("hist_3Jets_jetPt_notHV",
		"Jet Pt of non-HV jets in events with 3 or more jets;pt;count",
		100,0,2000)

	hist_2Jets_deltaPhi = self.makeTH1F("hist_2Jets_deltaPhi",
		"DeltaPhi of all HV Jets, events with only 2 jets; deltaPhi; count",
		100,0,rt.TMath.Pi())
	hist_3Jets_deltaPhi_HV = self.makeTH1F("hist_3Jets_deltaPhi_HV",
		"deltaPhi of all HV jets in events with 3;deltaPhi;count",
		100,0,rt.TMath.Pi())
	hist_3Jets_deltaPhi_notHV = self.makeTH1F("hist_3Jets_deltaPhi_notHV",
		"deltaPhi of non-HV jets in events with 3;deltaPhi;count",
		100,0,rt.TMath.Pi())
	
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
			hist_2Jets_mDPhiIndex.Fill(tree.iJetMaxDeltaPhi)
			hist_2Jets_ptMaxDPhi.Fill(tree.pTMaxDeltaPhi)
			for iJet in range(2):
				if bool(tree.JetsAK8_isHV[iJet]):
					hist_2Jets_jetPt.Fill(tree.JetsAK8[iJet].Pt())
					if iJet == 0:
						hist_2Jets_deltaPhi.Fill(tree.DeltaPhi1)
					if iJet == 1:
						hist_2Jets_deltaPhi.Fill(tree.DeltaPhi2)
		else:
			for iJet in range(nJets):
				if bool(tree.JetsAK8_isHV[iJet]):
					hist_3Jets_jetPt_HV.Fill(tree.JetsAK8[iJet].Pt())
					if iJet == 0:
						hist_3Jets_deltaPhi_HV.Fill(tree.DeltaPhi1)
					if iJet == 1:
						hist_3Jets_deltaPhi_HV.Fill(tree.DeltaPhi2)
					if iJet == 2:
						hist_3Jets_deltaPhi_HV.Fill(tree.DeltaPhi3)
				else:
					hist_3Jets_jetPt_notHV.Fill(tree.JetsAK8[iJet].Pt())
					if iJet == 0:
						hist_3Jets_deltaPhi_notHV.Fill(tree.DeltaPhi1)
					if iJet == 1:
						hist_3Jets_deltaPhi_notHV.Fill(tree.DeltaPhi2)
					if iJet == 2:
						hist_3Jets_deltaPhi_notHV.Fill(tree.DeltaPhi3)

	makePlot(
		[hist_2Jets_jetPt,hist_3Jets_jetPt_HV,hist_3Jets_jetPt_notHV],
		self.extraDir,
		"jetPts",
		"Pt of Jets;pT;count",logY=True
		)
	makePlot(
		[hist_2Jets_deltaPhi,hist_3Jets_deltaPhi_HV,hist_3Jets_deltaPhi_notHV],
		self.extraDir,
		"deltaPhis",
		"Deltaphi of Jets;deltaphi;count",logY=True
		)
	

def addLoop():
	baseClass.loop = loop

def makePlot(ListOfHistos, location, name, stackTitle, logY = False, gs = False):
	c1 = rt.TCanvas("c1","c1",900,600)
	c1.SetLineWidth(2)
	if gs:
		c1.SetGrayscale()
	if logY:
		c1.SetLogy()
	stack = rt.THStack("stack",stackTitle)
	stackTitle = stackTitle.split(";")
	#stack.SetMinimum(10) # use these two lines to set Y axis range
	#stack.SetMaximum(1100)
	i = 0
	for histo in ListOfHistos:
		histo.SetLineColor(i+1)
		histo.SetLineWidth(2)
		stack.Add(histo)
		i += 1
	if len(ListOfHistos) > 1:
		stack.Draw("nostack")
		c1.BuildLegend()
	else:
		stack.Draw("hist TEXT0")
	stack.SetTitle(stackTitle[0])
	stack.GetXaxis().SetTitle(stackTitle[1])
	stack.GetYaxis().SetTitle(stackTitle[2])
	rt.gPad.Modified()
	#stack.GetXaxis().SetLimits(0,4000) # use these two lines to set X axis range
	#stack.Draw("nostack")
	#t1 = rt.TText(0,80,"CMS Simulation") # use axis numbers for this
	#t1.Draw("same")
	c1.SaveAs(location+name+".png")
