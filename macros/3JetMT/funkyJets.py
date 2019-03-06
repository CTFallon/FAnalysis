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
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
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
	hasPrintedParticles = False
	funkyEventDict = {}
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		listOfFunkyJets = []
		eventHasFunkyJets = False
		for iJet in range(len(tree.JetsAK8)):
			if tree.JetsAK8_isISR[iJet] and tree.JetsAK8_HVfracPt[iJet] > 0.9:
				eventHasFunkyJets = True
				listOfFunkyJets.append(iJet)
	
		if eventHasFunkyJets:
			funkyEventDict[int(tree.EvtNum)] = listOfFunkyJets
			partsList = []
			partsListPdgId = []
			partsListisFromHVQuark = []
			for iPart in range(len(tree.GenParticles)):
				if (tree.GenParticles_Status[iPart] == 1):
					partsList.append(tree.GenParticles[iPart])
					partsListPdgId.append(tree.GenParticles_PdgId[iPart])
					partsListisFromHVQuark.append(tree.genParticleIsFromHVQuark[iPart])
			drawEventEtaPhiPlot(tree.JetsAK8, tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, tree.genParticleIsFromHVQuark, tree.EvtNum+1, self.extraDir)
			drawEventEtaPhiPlot(tree.JetsAK8, partsList, partsListPdgId, tree.METPhi, partsListisFromHVQuark, tree.EvtNum, self.extraDir)
			

			if hasPrintedParticles:
				continue
			hasPrintedParticles = True
			print(str(tree.EvtNum) + " GenParticlePdgID, MompdgID, Status, Jet#")
			for iPart in range(len(tree.GenParticles)):
				inJet = -1
				for iJet in range(len(tree.JetsAK8)):
					if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8:
						inJet = iJet
				iMom = tree.GenParticles_ParentIdx[iPart]
				print(iPart,tree.GenParticles_PdgId[iPart], tree.GenParticles_PdgId[iMom], tree.GenParticles_Status[iPart], inJet)

	print(funkyEventDict)
	

def addLoop():
	baseClass.loop = loop

def drawEventEtaPhiPlot(jetCollectionAK8, partCol, particlePDGID,METPhi, isFromHVQuark, plotNumber, edir):
	canv = rt.TCanvas("canv","canv",1600,800)
	histAxis = rt.TH2F("axisHsito", ";\eta;\phi",100,-6,6,100,-rt.TMath.Pi(),rt.TMath.Pi())
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionAK8)):
		objectList.append(rt.TEllipse(jetCollectionAK8[iJet].Eta(), jetCollectionAK8[iJet].Phi(), 0.8, 0.8))
		objectList[-1].SetLineColor(iJet+1)
		objectList[-1].SetLineWidth(2)
	for iPart in range(len(partCol)):
		objectList.append(rt.TMarker(partCol[iPart].Eta(),partCol[iPart].Phi(),2))
		objectList[-1].SetMarkerSize(2)
		if abs(particlePDGID[iPart]) == 4900101: # HV Quarks are
			print('Green, triangle, large')
			objectList[-1].SetMarkerColor(3) # Green
			objectList[-1].SetMarkerStyle(20) # Solid Triangles
			objectList[-1].SetMarkerSize(3) # Large
		elif abs(particlePDGID[iPart] == 4900023): # Z' is
			objectList[-1].SetMarkerStyle(4) #  thick cross
			objectList[-1].SetMarkerSize(3) # large
		elif abs(particlePDGID[iPart]) > 4900000: # HV Particles are
			objectList[-1].SetMarkerColor(2) # Red
		if isFromHVQuark[iPart]: # decendants of HV quarks are 
			objectList[-1].SetMarkerStyle(5) # X's
	objectList.append(rt.TLine(-6,METPhi,6,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	canv.SaveAs(edir+"funkyJets_etaPhi_"+str(plotNumber)+".png")


