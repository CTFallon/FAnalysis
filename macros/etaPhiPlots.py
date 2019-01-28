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
	tree.SetBranchStatus("Jets",1)
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
	for iEvent in range(nEvents):
		if nPass > 50:
			continue
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		nJetsHV = 0
		for i in range(len(tree.JetsAK8)):	
			nJetsHV += int(bool(tree.JetsAK8_isHV[i]))
		#if nJetsHV > 1:
		#	continue
		if not (len(tree.JetsAK8) > 2 and bool(tree.JetsAK8_isHV[2]) == True):
			continue
		nPass += 1
		partsList = []
		partsListPdgId = []
		partsListisFromHVQuark = []
		HVQparts = []
		for iPart in range(len(tree.GenParticles)):
			if (  
				( (tree.numberOfDaughtersAParticleHas[iPart] == 0) and (abs(tree.GenParticles_PdgId[iPart]) != 4900101)) or
				( (abs(tree.GenParticles_PdgId[iPart]) == 4900101) and (tree.GenParticles_Status[iPart] == 71) )
				):
				partsList.append(tree.GenParticles[iPart])
				partsListPdgId.append(tree.GenParticles_PdgId[iPart])
				partsListisFromHVQuark.append(tree.genParticleIsFromHVQuark[iPart])
			if ( ( abs( tree.GenParticles_PdgId[iPart] ) == 4900101 ) and ( tree.GenParticles_Status[iPart] == 71 ) ):
				HVQparts.append(tree.GenParticles[iPart])
		drawEventEtaPhiPlot(
			tree.JetsAK8,
			tree.JetsAK8_isHV,
			partsList,
			partsListPdgId,
			tree.METPhi,
			partsListisFromHVQuark,
			self.extraDir+"/etaPhi", str(nPass)+"_prunedGen")
		#drawEventEtaPhiPlot(
		#	tree.JetsAK8,
		#	tree.JetsAK8_isHV,
		#	tree.GenParticles,
		#	tree.GenParticles_PdgId,
		#	tree.METPhi,
		#	tree.genParticleIsFromHVQuark,
		#	self.extraDir+"/etaPhi", str(nPass)+"_allGen")
		drawEventEtaPhiPlot(
			tree.JetsAK8,
			tree.JetsAK8_isHV,
			partsList,
			partsListPdgId,
			tree.METPhi,
			partsListisFromHVQuark,
			self.extraDir+"/etaPhi", str(nPass)+"_AK4", AK4Jets = tree.Jets)
		preProcess(
			tree.JetsAK8,
			tree.JetsAK8_isHV,
			tree.METPhi,
			tree.MET,
			HVQparts,
			self.extraDir+"/etaPhi", str(nPass))

def addLoop():
	baseClass.loop = loop

def drawEventEtaPhiPlot(jetCollectionAK8, FSRJets, partCol, particlePDGID, METPhi, isFromHVQuark, direc, plotNumber, AK4Jets = []):
	hozPix = 800
	etaLim = 4.
	canv = rt.TCanvas("canv","canv",int(etaLim/rt.TMath.Pi()*hozPix),hozPix)
	histAxis = rt.TH2F("axisHsito", ";\eta;\phi",100,-etaLim,etaLim,100,-rt.TMath.Pi(),rt.TMath.Pi())
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionAK8)):
		objectList.append(rt.TEllipse(jetCollectionAK8[iJet].Eta(), jetCollectionAK8[iJet].Phi(), 0.8, 0.8))
		objectList[-1].SetLineColor(iJet+1) # jet color matches PT ordering
		objectList[-1].SetLineWidth(2)
		if not FSRJets[iJet]:
			objectList[-1].SetLineStyle(2) # non-HV jets are dotted circles
	if len(AK4Jets) > 0:
		for iJet in range(len(AK4Jets)):
			objectList.append(rt.TEllipse(AK4Jets[iJet].Eta(), AK4Jets[iJet].Phi(), 0.4, 0.4))
			objectList[-1].SetLineColor(iJet+1) # jet color matches PT ordering
			objectList[-1].SetLineWidth(2)
			objectList[-1].SetLineStyle(3) # AK4 jets are finely dashed circles
	for iPart in range(len(partCol)):
		#if abs(particlePDGID[iPart]) == 4900023 or abs(particlePDGID[iPart]) == 4900021: # skip Z' and hv-gluons
		#	continue
		objectList.append(rt.TMarker(partCol[iPart].Eta(),partCol[iPart].Phi(),2))
		objectList[-1].SetMarkerSize(2)
		if isFromHVQuark[iPart]: # apply only to particles that are decedants of HV-quarks
			objectList[-1].SetMarkerStyle(5) # make them X's
		if abs(particlePDGID[iPart]) == 4900101: # apply only to the HV-quarks
			objectList[-1].SetMarkerColor(3) # turn them green
			objectList[-1].SetMarkerStyle(22) # make them triangles
			objectList[-1].SetMarkerSize(3) # make them big
		elif abs(particlePDGID[iPart]) > 4900000: # apply only to invisible partilces that arn't the HVquarks
			objectList[-1].SetMarkerColor(2) # tuen them red
	objectList.append(rt.TLine(-etaLim,METPhi,etaLim,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	leg = rt.TLegend(0.01,0.4,0.12,0.85)
	#leg.SetNColumns(2)
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
	leg4Jet.Draw("same")
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
			leg.AddEntry("legMET","","l")
		else:
			leg.AddEntry("legMET","","l")
	else:
		leg.AddEntry("legMET","","l")
	leg.AddEntry(rt.TObject()," ","")
	leg.AddEntry(rt.TObject(),"Particles","")
	leg.AddEntry("legHVQ","","p")
	leg.AddEntry("legInvFSR","","p")
	leg.AddEntry("legVisFSR","","p")
	leg.AddEntry("legISR","","p")
	leg.Draw()

	canv.SaveAs(direc+"/etaPhi_"+plotNumber+".png")

def preProcess(jetCollectionAK8, FSRJets, METPhi, MET, partCol,direc, plotNumber, AK4Jets = []):
	# function to draw a "pre-processed" version of the event
	# in x-y space (phi-pho space)
	# MET should be pointing UP (+y)
	# jets are set so that jetMaxDPhi is RIGHT (+x)
	c1 = rt.TCanvas("c1","c1",800,800)

	dPhi = -METPhi+rt.TMath.Pi()/2.
	objectList = []
	metX, metY = rotateVector(dPhi, MET*rt.TMath.Cos(METPhi), MET*rt.TMath.Sin(METPhi))
	objectList.append(rt.TArrow(0,0,metX,metY,0.05,"|>"))
	objectList[-1].SetLineColor(rt.kBlue)
	objectList[-1].SetFillColor(rt.kBlue)
	objectList[-1].SetLineWidth(2)
	objectList[-1].SetFillStyle(3003)
	objectList[-1].SetLineStyle(2)
	for iJet in range(len(jetCollectionAK8)):
		jetX, jetY = rotateVector(dPhi,jetCollectionAK8[iJet].Px(),jetCollectionAK8[iJet].Py())
		if iJet == 0:
			signJetX = jetX/abs(jetX)
		objectList.append(rt.TArrow(0,0,signJetX*jetX,jetY))
		objectList[-1].SetLineColor(iJet+1)
		objectList[-1].SetLineWidth(2)
		if not FSRJets[iJet]:
			objectList[-1].SetLineStyle(2)
	for particle in partCol:
		partX, partY = rotateVector(dPhi,particle.Px(), particle.Py())
		objectList.append(rt.TArrow(0,0,partX,partY,0.025,"|>"))
		objectList[-1].SetLineColor(rt.kGreen)
		objectList[-1].SetFillColor(rt.kGreen)

	linVals = []
	for thing in objectList:
		linVals.append(abs(thing.GetX2()))
		linVals.append(abs(thing.GetY2()))
	Lim = max(linVals)
	Lim = Lim*1.1
	histAxis = rt.TH2F("axisHsito", "Transverse Momentum;P(x');P(y')",100,-Lim,Lim,100,-Lim,Lim)
	histAxis.SetStats(False)
	histAxis.Draw()
	for thing in objectList:
		thing.Draw()
	c1.SaveAs(direc+"/etaPhi"+plotNumber+"_polar.png")

def rotateVector(phi, oldX, oldY):
	newX = oldX*rt.TMath.Cos(phi) - oldY*rt.TMath.Sin(phi)
	newY = oldX*rt.TMath.Sin(phi) + oldY*rt.TMath.Cos(phi)
	return newX, newY

	
