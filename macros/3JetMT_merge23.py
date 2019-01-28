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
	
	hist_Pt_jet2 = self.makeTH1F("hist_Pt_jet2","Jet 2 MT;MT;count",100,0,4000)
	hist_Pt_jet3 = self.makeTH1F("hist_Pt_jet3","Jet 3 MT;MT;count",100,0,4000)
	hist_Pt_jetB = self.makeTH1F("hist_Pt_jetB","Jet B MT;MT;count",100,0,4000)
	hist_Pt_merged = self.makeTH1F("hist_Pt_merged","Jet 2+3 MT;MT;count",100,0,4000)
	hist_Pt_jet1 = self.makeTH1F("hist_Pt_jet1","Jet 1 MT;MT;count",100,0,4000)
	hist_Pt_HV = self.makeTH1F("hist_Pt_HV","HV MT;MT;count",100,0,4000)

	hist_deltaR_jet3HV = self.makeTH1F("hist_deltaR_jet3HV","DeltaR(j3HV);deltaR;count",
		100,0,5)
	hist_deltaR_jet3notHV = self.makeTH1F("hist_deltaR_jet3notHV","DeltaR(j3NotHV);deltaR;count",
		100,0,5)

	hist2d_jetHVPercent = self.makeTH2F("hist2d_jetHVPercent","Percent of Jets that are HV;Jet Number; HV T/F",
		7,0,7,
		2,0,2)

	nIdJets = [0.,0.,0.,0.,0.,0.,0.]
	nHVJets = [0.,0.,0.,0.,0.,0.,0.]
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15


		nJets = len(tree.JetsAK8)
		try:
			merged23 = tree.JetsAK8[1]+tree.JetsAK8[2]
		except IndexError:
			merged23 = tree.JetsAK8[1]
		hist_Pt_jet1.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		hist_Pt_jet2.Fill(trans_mass_Njet([merged23], tree.MET, tree.METPhi))
		try:
			hist_Pt_jet3.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[2]], tree.MET, tree.METPhi))
		except IndexError:
			hist_Pt_jet3.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		hist_Pt_merged.Fill(trans_mass_Njet([tree.JetsAK8[0],merged23], tree.MET, tree.METPhi))
		try:
			deltaR = tree.JetsAK8[1].DeltaR(tree.JetsAK8[2])
			if bool(tree.JetsAK8_isHV[2]):
				hist_deltaR_jet3HV.Fill(deltaR)
			else:
				hist_deltaR_jet3notHV.Fill(deltaR)
			if deltaR > 1.0:
				hist_Pt_jetB.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
			else:
				hist_Pt_jetB.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
		except IndexError:
			hist_Pt_jetB.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		HVJets = []
		for iJet in range(min(7,nJets)):
			nIdJets[iJet] += 1
			if bool(tree.JetsAK8_isHV[iJet]):
				nHVJets[iJet] += 1
				HVJets.append(tree.JetsAK8[iJet])
		hist_Pt_HV.Fill(trans_mass_Njet(HVJets, tree.MET, tree.METPhi))
	c1 = rt.TCanvas("c1","c1",1200,900)
	hist2d_jetHVPercent.SetAxisRange(0.0,1.0,"Z")
	for iJet in range(7):
		if nIdJets[iJet] > 0:
			hist2d_jetHVPercent.SetBinContent(iJet+1,1,1.-(nHVJets[iJet]/nIdJets[iJet]))
			hist2d_jetHVPercent.SetBinContent(iJet+1,2,(nHVJets[iJet]/nIdJets[iJet]))
	hist2d_jetHVPercent.Draw("colz")
	c1.SaveAs(self.extraDir+"2d_jetNumberHVPercent.png")



	c1.SetLogy()
	hist_Pt_jet2.SetLineColor(2)
	hist_Pt_jet3.SetLineColor(3)
	hist_Pt_merged.SetLineColor(4)
	hist_Pt_jetB.SetLineColor(5)
	hist_Pt_HV.SetLineStyle(2)
	
	#hist_Pt_jet2.Draw()
	hist_Pt_HV.Draw()
	hist_Pt_merged.Draw("same")
	hist_Pt_jet3.Draw("same")
	hist_Pt_jet1.Draw("same")
	hist_Pt_jetB.Draw("same")
	c1.BuildLegend(0.16,0.13,0.5,0.3)
	c1.SaveAs(self.extraDir+"meged23_MT.png")

	hist_deltaR_jet3HV.Draw()
	hist_deltaR_jet3notHV.SetLineColor(2)
	hist_deltaR_jet3notHV.Draw("same")
	c1.SaveAs(self.extraDir+"deltaRJet23.png")

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

