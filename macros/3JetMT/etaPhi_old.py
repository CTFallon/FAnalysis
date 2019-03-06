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
	#tree.SetBranchStatus("*CA11*", 1)
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
	tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
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
	nPass = 0
	hist_njets = self.makeTH1F("hist_njets","Njets;nJets;Count",4,2,6)
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		hist_njets.Fill(len(tree.JetsAK8))
		if int(iEvent) != 53:
			continue

		nPass += 1
		partsList = []
		partsListPdgId = []
		partsListisFromHVQuark = []
		for iPart in range(len(tree.GenParticles)):
			if (
				( (tree.numberOfDaughtersAParticleHas[iPart] == 0)  and (abs(tree.GenParticles_PdgId[iPart]) != 4900101)) or
				( (abs(tree.GenParticles_PdgId[iPart]) == 4900101) and (tree.GenParticles_Status[iPart] == 71) )
				):
				partsList.append(tree.GenParticles[iPart])
				partsListPdgId.append(tree.GenParticles_PdgId[iPart])
				partsListisFromHVQuark.append(tree.genParticleIsFromHVQuark[iPart])
		drawEventEtaPhiPlot(
			tree.JetsAK8,
			partsList,
			partsListPdgId,
			tree.METPhi,
			partsListisFromHVQuark,
			self.extraDir+"etaPhi", str(iEvent)+"_prunedGenParticles")
		drawEventEtaPhiPlot(
			tree.JetsAK8,
			tree.GenParticles,
			tree.GenParticles_PdgId,
			tree.METPhi,
			tree.genParticleIsFromHVQuark,
			self.extraDir+"etaPhi", str(iEvent)+"_allGenParticles")
	

	c1 = rt.TCanvas()
	hist_njets.GetYaxis().SetTitleOffset(2)
	hist_njets.Draw()
	hist_njets.GetYaxis().SetTitleOffset(2)
	hist_njets.Draw()
	c1.SaveAs(self.extraDir+"nJets.png")
def addLoop():
	baseClass.loop = loop

def drawEventEtaPhiPlot(jetCollectionAK8, partCol, particlePDGID,METPhi, isFromHVQuark, direc, plotNumber):
	vertPix = 800
	etaLim = 3.
	phiLim = rt.TMath.Pi()
	tM = 0.05
	bM = 0.13
	lM = 0.26
	rM = 0.02
	canv = rt.TCanvas("canv","canv",int((etaLim/phiLim)/(1-lM-rM)*(1-tM-bM)*vertPix),vertPix)
	canv.SetTopMargin(tM)
	canv.SetBottomMargin(bM)
	canv.SetLeftMargin(lM)
	canv.SetRightMargin(rM)
	histAxis = rt.TH2F("axisHsito", "",10,-etaLim,etaLim,10,-phiLim,phiLim)
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionAK8)):
		objectList.append(rt.TEllipse(jetCollectionAK8[iJet].Eta(), jetCollectionAK8[iJet].Phi(), 0.8, 0.8))
		objectList[-1].SetLineColor(iJet+1) # jet color matches PT ordering
		objectList[-1].SetLineWidth(2)
	for iPart in range(len(partCol)):
		if (abs(partCol[iPart].Eta()) > etaLim) or (abs(partCol[iPart].Phi()) > phiLim):
			continue
		objectList.append(rt.TMarker(partCol[iPart].Eta(),partCol[iPart].Phi(),2))
		objectList[-1].SetMarkerSize(2)
		if isFromHVQuark[iPart]: # apply only to particles that are decedants of HV-quarks
			objectList[-1].SetMarkerStyle(5) # make them X's
		if abs(particlePDGID[iPart]) == 4900101: # apply only to the HV-quarks
			objectList[-1].SetMarkerColor(3) # turn them green
			objectList[-1].SetMarkerStyle(22) # make them triangles
			objectList[-1].SetMarkerSize(3) # make them big
		elif abs(particlePDGID[iPart]) > 4900000: # apply only to invisible partilces that arn't the HVquarks
			objectList[-1].SetMarkerColor(2) # make them red
			objectList[-1].SetMarkerStyle(5)
	objectList.append(rt.TLine(-etaLim,METPhi,etaLim,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	leg = rt.TLegend(0.01,0.3,.21,0.80)
	leg1Jet = rt.TH1F("leg1Jet","1st",2,-100,-99)
	leg1Jet.SetLineWidth(2)
	leg1Jet.Draw("same")
	leg2Jet = rt.TH1F("leg2Jet","2nd",2,-100,-99)
	leg2Jet.SetLineColor(2)
	leg2Jet.SetLineWidth(2)
	leg2Jet.Draw("same")
	leg3Jet = rt.TH1F("leg3Jet","3rd",2,-100,-99)
	leg3Jet.SetLineWidth(2)
	leg3Jet.SetLineColor(3)
	leg3Jet.Draw("same")
	leg4Jet = rt.TH1F("leg4Jet","4th",2,-100,-99)
	leg4Jet.SetLineWidth(2)
	leg4Jet.SetLineColor(4)
	legMET = rt.TH1F("legMET","MET \phi",2,-100,-99)
	legMET.SetLineWidth(2)
	legMET.SetLineStyle(2)
	legMET.SetLineColor(4)
	legMET.Draw("same")
	legHVQ = rt.TH1F("legHVQ","\chi_{1,2}",2,-100,-99)
	legHVQ.SetMarkerColor(3)
	legHVQ.SetMarkerStyle(22)
	legHVQ.SetMarkerSize(2)
	legHVQ.Draw("same")
	legInvFSR = rt.TH1F("legInvFSR","Invisible, from Z\'",2,-100,-99)
	legInvFSR.SetMarkerStyle(5)
	legInvFSR.SetMarkerColor(2)
	legInvFSR.SetMarkerSize(2)
	legInvFSR.Draw("same")
	legVisFSR = rt.TH1F("legVisFSR","Visible, from Z\'",2,-100,-99)
	legVisFSR.SetMarkerStyle(5)
	legVisFSR.SetMarkerSize(2)
	legVisFSR.Draw("same")
	legISR = rt.TH1F("legISR","Not From Z\'",2,-100,-99)
	legISR.SetMarkerStyle(2)
	legISR.SetMarkerSize(2)
	legISR.Draw("same")

	# One Column
	leg.AddEntry(rt.TObject(), "Jets, p_{T} Rank","")
	leg.AddEntry("leg1Jet","","l")
	leg.AddEntry("leg2Jet","","l")
	if len(jetCollectionAK8) >= 3:
		leg.AddEntry("leg3Jet","","l")
	if len(jetCollectionAK8) >= 4:
		leg.AddEntry("leg4Jet","","l")
	leg.AddEntry(rt.TObject()," ","")
	leg.AddEntry("legMET","","l")
	leg.AddEntry(rt.TObject()," ","")
	leg.AddEntry(rt.TObject(),"Particles","")
	leg.AddEntry("legHVQ","","p")
	leg.AddEntry("legInvFSR","","p")
	leg.AddEntry("legVisFSR","","p")
	leg.AddEntry("legISR","","p")
	leg.Draw()
	# axis labels
	t = rt.TLatex()
	t.DrawLatex(-1.2*etaLim,0.9*phiLim,"\phi")
	t.DrawLatex(0.9*etaLim,-1.2*phiLim,"\eta")

	canv.SaveAs(direc+"/etaPhi_"+plotNumber+".png")
