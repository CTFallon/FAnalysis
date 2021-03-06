from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()
# macro to look at how using the ratio of pT from a jet can be used to 
# discriminate against FSRvsISR jets

# ideas:
	# numerator:
		# jet pT (max dPhi)
		# jet0 pT
	# denominators:
		# scalar sum of pT's of leading 3 jets
		# scalar sum of all jets pT's

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

	# var4 - jetMaxDPhi/3Jets


	hist_2_var4 = self.makeTH1F("hist_2_var4",
		"True Dijet PtFraction MaxDPhi Jet / Leading 3 Jets",
		100,-0.01,0.85)
	hist_3_var4 = self.makeTH1F("hist_3_var4",
		"True Trijet PtFraction MaxDPhi Jet / Leading 3 Jets",
		100,-0.01,0.85)
	hist_o_var4 = self.makeTH1F("hist_o_var4",
		"True Other PtFraction MaxDPhi Jet / Leading 3 Jets",
		100,-0.01,0.85)
	hist_eventMTcat = self.makeTH1F("hist_eventMTcat",
		"MT Category of Events",
		5,0,5)


	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			hist_eventMTcat.Fill(0)
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		if len(tree.JetsAK8) < 3:
			hist_eventMTcat.Fill(1)
			continue

		allJetsPt = 0
		lead3JetsPt = 0
		for iJet in range(len(tree.JetsAK8)):
			allJetsPt += tree.JetsAK8[iJet].Pt()
			if iJet <= 2:
				lead3JetsPt += tree.JetsAK8[iJet].Pt()
		var4 = tree.JetsAK8[tree.iJetMaxDeltaPhi].Pt()/lead3JetsPt

		if not tree.JetsAK8_isHV[2]:
			hist_eventMTcat.Fill(2)
			hist_2_var4.Fill(var4)
		elif tree.JetsAK8_isHV[2]:
			hist_eventMTcat.Fill(3)
			hist_3_var4.Fill(var4)


	makePlot([hist_2_var4, hist_3_var4],self.extraDir,"mdphiOver3","PtFraction (Jet Max DPhi / Leading 3 Jets);Fraction;Count")


	hist_2_var4.Scale(1/hist_2_var4.Integral())
	hist_3_var4.Scale(1/hist_3_var4.Integral())
	

	makePlot([hist_2_var4, hist_3_var4],self.extraDir,"mdphiOver3_NORM","PtFraction (Jet Max DPhi / Leading 3 Jets);Fraction;a.u.")
	hist_eventMTcat.GetXaxis().SetBinLabel(1, "failed")
	hist_eventMTcat.GetXaxis().SetBinLabel(2, "only 2 jets")
	hist_eventMTcat.GetXaxis().SetBinLabel(3, "dijet MT")
	hist_eventMTcat.GetXaxis().SetBinLabel(4, "trijet MT")
	hist_eventMTcat.GetXaxis().SetBinLabel(5, "other MT")
	makePlot([hist_eventMTcat],self.extraDir,"eventMTcat","Event MT Category;Category;Count", logY= True)
	

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
	stack.GetXaxis().SetTitle(stackTitle[1])
	stack.GetYaxis().SetTitle(stackTitle[2])
	rt.gPad.Modified()
	#stack.GetXaxis().SetLimits(0,4000) # use these two lines to set X axis range
	#stack.Draw("nostack")
	#t1 = rt.TText(0,80,"CMS Simulation") # use axis numbers for this
	#t1.Draw("same")
	c1.SaveAs(location+name+".png")
