from analysisBase import baseClass
import ROOT as rt
from array import array
# this macro makes the friend tree to our analysis ntuple tree
def loop(self):
	print("Starting Loop")
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))
	
	# Turn off all branches, selective turn on branches
	tree.SetBranchStatus("*", 0)
	tree.SetBranchStatus("*AK8*", 1)
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticle*",1)
	tree.SetBranchStatus("Electrons",1)
	tree.SetBranchStatus("Muons",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	friend = rt.TTree("friend","friend")
	self.objects.append(friend)
	# create branch variables
	passedPreSelection = array('i', [0])
	numGenParts = array('i',[0])
	numJets = array('i',[0])
	genParticleInAK8Jet = array('i' ,[-1 for x in range(400)])
	genParticleIsFromHVQuark = array('f' ,[0. for x in range(400)])
	numberOfDaughtersAParticleHas = array('f', [0. for x in range(400)])
	numHVPartsInJet = array("i",[-1 for x in range(10)])
	numSMPartsInJet = array("i",[-1 for x in range(10)])
	pGJ_visible = array('f',[-1. for x in range(10)])
	pGJ_invis = array('f',[-1. for x in range(10)])
	pGJ_every = array('f',[-1. for x in range(10)])
	iJetMaxDeltaPhi = array('i',[-1])
	pTMaxDeltaPhi = array('f',[0.])
	dPhiMaxDeltaPhi = array('f',[0.])
	MTFromParticles = array('f',[0.])
	zPrimept = array('f',[0.])
	zPrimephi = array('f',[0.])
	numPartsStat1Daugher0 = array('i',[0])
	numPartsStat1Daughternot0 = array('i',[0])
	numPartsStatnot1Daughter0 = array('i',[0])
	numPartsStatnot1Daughternot0 = array('i',[0])

	fracVisPTfromVisHVQ = array("f",[0. for x in range(10)])
	fracInvPTfromInvHVQ = array("f",[0. for x in range(10)])
	fracTotPTfromAllHVQ = array("f",[0. for x in range(10)])
	fracTotPTfromVisHVQ = array("f",[0. for x in range(10)])
	fracTotPTfromInvHVQ = array("f",[0. for x in range(10)])
	fracTotPTfromVis = array("f",[0. for x in range(10)])
	fracVisHVQtoInvHVQ = array("f",[0. for x in range(10)])

	pt_VisHVQ = array("f",[0. for x in range(10)])
	pt_Vis = array("f",[0. for x in range(10)])

	dijetEvent = array("i",[0])
	trijetEvent = array("i",[0])
	otherEvent = array("i",[0])
	
	

	friend.Branch("passedPreSelection",passedPreSelection, 'passedPreSelection/I')
	friend.Branch("numGenParts",numGenParts, 'numGenParts/I')
	friend.Branch("genParticleInAK8Jet",genParticleInAK8Jet,'genParticleInAK8Jet[400]/I')
	friend.Branch("genParticleIsFromHVQuark",genParticleIsFromHVQuark,'genParticleIsFromHVQuark[400]/F')
	friend.Branch("numberOfDaughtersAParticleHas",numberOfDaughtersAParticleHas,'numberOfDaughtersAParticleHas[400]/F')
	friend.Branch("numHVPartsInJet",numHVPartsInJet,'numHVPartsInJet[10]/I')
	friend.Branch("numSMPartsInJet",numSMPartsInJet,'numSMPartsInJet[10]/I')
	friend.Branch("iJetMaxDeltaPhi",iJetMaxDeltaPhi,'iJetMaxDeltaPhi/I')
	friend.Branch("pTMaxDeltaPhi",pTMaxDeltaPhi,'pTMaxDeltaPhi/F')
	friend.Branch("dPhiMaxDeltaPhi",dPhiMaxDeltaPhi,'dPhiMaxDeltaPhi/F')
	friend.Branch("MTFromParticles",MTFromParticles,'MTFromParticles/F')
	friend.Branch("zPrimept",zPrimept, 'zPrimept/F')
	friend.Branch("zPrimephi",zPrimephi, 'zPrimephi/F')
	friend.Branch("pGJ_visible",pGJ_visible,"pGJ_visible[10]/F")
	friend.Branch("pGJ_invis",pGJ_invis,"pGJ_invis[10]/F")
	friend.Branch("pGJ_every",pGJ_every,"pGJ_every[10]/F")
	friend.Branch("numPartsStat1Daugher0",numPartsStat1Daugher0,"numPartsStat1Daugher0/I")
	friend.Branch("numPartsStat1Daughternot0",numPartsStat1Daughternot0,"numPartsStat1Daughternot0/I")
	friend.Branch("numPartsStatnot1Daughter0",numPartsStatnot1Daughter0,"numPartsStatnot1Daughter0/I")
	friend.Branch("numPartsStatnot1Daughternot0",numPartsStatnot1Daughternot0,"numPartsStatnot1Daughternot0/I")

	friend.Branch("fracVisPTfromVisHVQ",fracVisPTfromVisHVQ,"fracVisPTfromVisHVQ[10]/F")
	friend.Branch("fracInvPTfromInvHVQ",fracInvPTfromInvHVQ,"fracInvPTfromInvHVQ[10]/F")
	friend.Branch("fracTotPTfromAllHVQ",fracTotPTfromAllHVQ,"fracTotPTfromAllHVQ[10]/F")
	friend.Branch("fracTotPTfromVisHVQ",fracTotPTfromVisHVQ,"fracTotPTfromVisHVQ[10]/F")
	friend.Branch("fracTotPTfromInvHVQ",fracTotPTfromInvHVQ,"fracTotPTfromInvHVQ[10]/F")
	friend.Branch("fracTotPTfromVis",fracTotPTfromVis,"fracTotPTfromVis[10]/F")
	friend.Branch("fracVisHVQtoInvHVQ",fracVisHVQtoInvHVQ,"fracVisHVQtoInvHVQ[10]/F")
	
	friend.Branch("pt_VisHVQ",pt_VisHVQ,"pt_VisHVQ[10]/F")
	friend.Branch("pt_Vis",pt_Vis,"pt_Vis[10]/F")
	
	friend.Branch("DijetEvent", dijetEvent,"DijetEvent/I")
	friend.Branch("TrijetEvent", trijetEvent,"TrijetEvent/I")
	friend.Branch("OtherjetEvent", otherEvent,"otherEvent/I")

	tree.AddFriend(friend)
	maxNofParticle = 0
	maxNofJets = 0 
	for iEvent in range(nEvents):
		#if iEvent !=  39854:
		#	continue
		if iEvent%1000==0:
			print("Event "+ str(iEvent)+"/"+str(nEvents))
		tree.GetEvent(iEvent)
		if len(tree.GenParticles) > 400:
			print(len(tree.GenParticles))
		passedPreSelection[0] = 0
		numGenParts[0] = len(tree.GenParticles)
		numJets[0] = len(tree.JetsAK8)
		iJetMaxDeltaPhi[0] = -1
		pTMaxDeltaPhi[0] = 0.
		dPhiMaxDeltaPhi[0] = 0.
		zPrimept[0] = 0.
		zPrimephi[0] = 0.
		numPartsStat1Daugher0[0] = 0
		numPartsStat1Daughternot0[0] = 0
		numPartsStatnot1Daughter0[0] = 0
		numPartsStatnot1Daughternot0[0] = 0
		dijetEvent[0] = 0
		trijetEvent[0] = 0
		otherEvent[0] = 0
		if maxNofParticle < numGenParts[0]:
			maxNofParticle = numGenParts[0]
		if maxNofJets < numJets[0]:
			maxNofJets = numJets[0]
		for i in range(400):
			genParticleInAK8Jet[i] = -1
			genParticleIsFromHVQuark[i] = 0
			numberOfDaughtersAParticleHas[i] = 0.
		for i in range(10):
			fracVisPTfromVisHVQ[i] = -1.
			fracInvPTfromInvHVQ[i] = -1.
			fracTotPTfromAllHVQ[i] = -1.
			fracTotPTfromVisHVQ[i] = -1.
			fracTotPTfromInvHVQ[i] = -1.
			fracTotPTfromVis[i] = -1.
			fracVisHVQtoInvHVQ[i] = -1.
			pt_VisHVQ[i] = -1.
			pt_Vis[i] = -1.
			numHVPartsInJet[i] = -1
			numSMPartsInJet[i] = -1
		for i in range(numJets[0]):
			fracVisPTfromVisHVQ[i] = 0.
			fracInvPTfromInvHVQ[i] = 0.
			fracTotPTfromAllHVQ[i] = 0.
			fracTotPTfromVisHVQ[i] = 0.
			fracTotPTfromInvHVQ[i] = 0.
			fracTotPTfromVis[i] = 0.
			fracVisHVQtoInvHVQ[i] = 0.
			pt_VisHVQ[i] = 0.
			pt_Vis[i] = 0.
			numHVPartsInJet[i] = 0
			numSMPartsInJet[i] = 0
		for i in range(10):
			pGJ_visible[i] = -1.
			pGJ_invis[i] = -1.
			pGJ_every[i] = -1.
		pass1 = 0 # atleast 2 AK8 jets
		pass2 = 0 # lead 2 jets each have pt > 200 GeV
		pass3 = 0 # lead 2 jets each have eta < 2.4
		pass4 = 0 # deltaEta(jet1, jet2) < 1.5
		pass5 = 0 # MET/MT > 0.15
		pass6 = 1 # MT > 1500 GeV, TEMP TURN OFF
		pass7 = 0 # lepton veto
		# At least 2 jets in the event, temp 3 for 3JetMT purposes
		if (len(tree.JetsAK8)>=2): 
			pass1 = 1
		else: # special case because next check can cause Index Error
			passedPreSelection[0] = 0
			friend.Fill()
			continue
		# Both leading jets must have pt > 170
			# updated to 200 GeV on Feb 13
		if ((tree.JetsAK8[0].Pt() > 200.0) and (tree.JetsAK8[1].Pt() > 200.0)): 
			pass2 = 1
		# both leading jets must have |eta| < 2.4, Feb 13 update
		if ((abs(tree.JetsAK8[0].Eta()) < 2.4) and (abs(tree.JetsAK8[1].Eta()) < 2.4)): 
			pass3 = 1
		# leading jets must be within 1.5 delta eta of eachother, feb 13
		if abs(tree.JetsAK8[0].Eta()-tree.JetsAK8[1].Eta()) < 1.5:
			pass4 = 1
		# MET/MT ratio must be greater than 0.15
		if (tree.MET/tree.MT_AK8 > 0.15):
			pass5 = 1
		# MT must be greater than 1500 GeV, feb 13
		if tree.MT_AK8 > 1500:
			pass6 = 1
		# must not have any leptons
		if ((len(tree.Electrons) + len(tree.Muons)) == 0):
			pass7 = 1
		if pass1 and pass2 and pass3 and pass4 and pass5 and pass6 and pass7:
			passedPreSelection[0] = 1
		else:
			passedPreSelection[0] = 0
			friend.Fill()
			continue
		
		# for each PARTICLE:
		# record:
		#   if it is decendant of an HV quark
		
		for iPart in range(2,len(tree.GenParticles)):
			iParent = tree.GenParticles_ParentIdx[iPart]
			if iParent != -1:
				numberOfDaughtersAParticleHas[iParent] += 1.
			if (abs(tree.GenParticles_PdgId[iParent]) == 4900023) or genParticleIsFromHVQuark[iParent]:
				genParticleIsFromHVQuark[iPart] = 1
			if tree.GenParticles_PdgId[iPart] == 4900023:
				zPrimept[0] = tree.GenParticles[iPart].Pt()
				zPrimephi[0] = tree.GenParticles[iPart].Phi()
		
		# calculate the 'final-state' MT from all visible, daughterless particles that are decendants of HV particls
		# also determin what jet particles belong to
		particlesForMT = []
		for iPart in range(2,len(tree.GenParticles)):
			## ==== # Results from this block:
				# If Status == 1, nDaughters == 0. 
				# But Status != 1, !=> nDaughters != 0 (i.e. there exists a particle that is not stable (stat == 1) but has 0 daugters.)
			# Check if particle status == 1 correlates with numberOfDaughtersAParticleHas
			if tree.GenParticles_Status[iPart] == 1 and numberOfDaughtersAParticleHas[iPart] == 0:
				numPartsStat1Daugher0[0] += 1
			if tree.GenParticles_Status[iPart] == 1 and numberOfDaughtersAParticleHas[iPart] != 0:
				numPartsStat1Daughternot0[0] += 1
			if tree.GenParticles_Status[iPart] != 1 and numberOfDaughtersAParticleHas[iPart] == 0:
				numPartsStatnot1Daughter0[0] += 1
			if tree.GenParticles_Status[iPart] != 1 and numberOfDaughtersAParticleHas[iPart] != 0:
				numPartsStatnot1Daughternot0[0] += 1
			## ====
			if genParticleIsFromHVQuark[iPart] == 1 and tree.GenParticles_Status[iPart] == 1:
				particlesForMT.append(tree.GenParticles[iPart])
			for iJet in range(len(tree.JetsAK8)-1,-1,-1):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8:
					genParticleInAK8Jet[iPart] = iJet
		MTFromParticles[0] = trans_mass_Njet(particlesForMT, tree.MET, tree.METPhi)
		# for each JET record:
		# how much of the total PT from daughter-less particles is from HV decendants
		# number of particles in the jet, total, visible, and invisible
		# also record which jet is furthest from the METPhi
		for iJet in range(len(tree.JetsAK8)):
			pGJ_vis = rt.TLorentzVector(0.,0.,0.,0.)
			pGJ_inv = rt.TLorentzVector(0.,0.,0.,0.)
			pGJ_all = rt.TLorentzVector(0.,0.,0.,0.)
			vectorPt = rt.TLorentzVector(0.,0.,0.,0.)
			HVvecPt = rt.TLorentzVector(0.,0.,0.,0.)
			deltaPhi = deltaPhi_calc(tree.JetsAK8[iJet].Phi(),tree.METPhi)
			if deltaPhi > dPhiMaxDeltaPhi[0]:
				iJetMaxDeltaPhi[0] = iJet
				pTMaxDeltaPhi[0] = tree.JetsAK8[iJet].Pt()
				dPhiMaxDeltaPhi[0] = deltaPhi
			visPt = 0.
			invPt = 0.
			totPt = 0.
			HVQallPt = 0.
			HVQvisPt = 0.
			HVQinvPt = 0.
			for iPart in range(2,len(tree.GenParticles)):
				# only want final-state particles (stable)
				# only want particles that belong to the jet
				# only want particles that descend from HVQuarks
				# ignore particles that a) anr't stable or b) arn't in the jet
				# ====  old way
				#if (numberOfDaughtersAParticleHas[iPart] != 0):
				#	continue
				# ===
				if (tree.GenParticles_Status[iPart] != 1):
					continue
				if (genParticleInAK8Jet[iPart] != iJet):
					continue
				# right now, we have all particles that are stable and in the jet
				# so we want to fill our pseudoGenJets:
				pGJ_all = pGJ_all + tree.GenParticles[iPart]
				totPt += tree.GenParticles[iPart].Pt()
				if genParticleIsFromHVQuark[iPart] == 1:
					HVQallPt += tree.GenParticles[iPart].Pt()
				if abs(tree.GenParticles_PdgId[iPart]) < 4900000: #visible particles
					pGJ_vis = pGJ_vis + tree.GenParticles[iPart]
					numSMPartsInJet[iJet] += 1
					visPt += tree.GenParticles[iPart].Pt()
					if genParticleIsFromHVQuark[iPart] == 1: # isolate particles from HVQuarks
						HVQvisPt += tree.GenParticles[iPart].Pt()
				else:# invisible
					numHVPartsInJet[iJet] += 1
					pGJ_inv = pGJ_inv + tree.GenParticles[iPart]
					invPt += tree.GenParticles[iPart].Pt()
					if genParticleIsFromHVQuark[iPart] == 1: # isolate particles from HVQuarks
						HVQinvPt += tree.GenParticles[iPart].Pt()
			pt_VisHVQ[iJet] = HVQvisPt
			pt_Vis[iJet] = visPt
			try:
				fracTotPTfromAllHVQ[iJet] = float(HVQallPt)/float(totPt)
				fracTotPTfromVisHVQ[iJet] = float(HVQvisPt)/float(totPt)
				fracTotPTfromInvHVQ[iJet] = float(HVQinvPt)/float(totPt)
				fracTotPTfromVis[iJet] = float(visPt)/float(totPt)
			except ZeroDivisionError:
				fracTotPTfromAllHVQ[iJet] = -1.
				fracTotPTfromVisHVQ[iJet] = -1.
				fracTotPTfromInvHVQ[iJet] = -1.
				fracTotPTfromVis[iJet] = -1.
			try:
				fracVisPTfromVisHVQ[iJet] = float(HVQvisPt)/float(visPt)
			except ZeroDivisionError:
				fracVisPTfromVisHVQ[iJet] = -1.
			try:
				fracInvPTfromInvHVQ[iJet] = float(HVQinvPt)/float(invPt)
			except ZeroDivisionError:
				fracInvPTfromInvHVQ[iJet] = -1.
			try:
				fracVisHVQtoInvHVQ[iJet] = float(HVQvisPt)/float(HVQinvPt)
			except ZeroDivisionError:
				fracVisHVQtoInvHVQ[iJet] = -1.
				
			pGJ_visible[iJet] = pGJ_vis.Pt()
			pGJ_invis[iJet] = pGJ_inv.Pt()
			pGJ_every[iJet] = pGJ_all.Pt()

	
		jetCode = [0,0,0,0,0]
		for iJet in range(len(tree.JetsAK8)):
			if tree.JetsAK8_isHV[iJet]:
				if iJet <= 2:
					jetCode[-iJet-1] = 1
		jetValue = jetCode[-1]*1+jetCode[-2]*2+jetCode[-3]*4+jetCode[-4]*8+jetCode[-5]*16
		if jetValue == 3:
			dijetEvent[0] = 1
		elif jetValue == 7:
			trijetEvent[0] = 1
		else:
			otherEvent[0] = 1

		#if iEvent == 39854:
		#	jets = tree.JetsAK8
		#	partsList = []
		#	partsListPdgId = []
		#	partsListisFromHVQuark = []
		#	for iPart in range(len(tree.GenParticles)):
		#		if tree.GenParticles_Status[iPart] == 1:
		#			partsList.append(tree.GenParticles[iPart])
		#			partsListPdgId.append(tree.GenParticles_PdgId[iPart])
		#			partsListisFromHVQuark.append(genParticleIsFromHVQuark[iPart])
		#	drawEventEtaPhiPlot(tree.JetsAK8, tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, genParticleIsFromHVQuark, iEvent+1)
		#	drawEventEtaPhiPlot(tree.JetsAK8, partsList, partsListPdgId, tree.METPhi, partsListisFromHVQuark, iEvent)
		#	print("index, PDG ID, parent index, iJet, Status, HVQD?, ")
		#	for iPart in range(len(tree.GenParticles)):
		#		print(iPart, tree.GenParticles_PdgId[iPart], tree.GenParticles_ParentIdx[iPart], genParticleInAK8Jet[iPart], tree.GenParticles_Status[iPart],genParticleIsFromHVQuark[iPart])
		friend.Fill()

	print(maxNofParticle)
	print(maxNofJets)
	print("Loop Ending")


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

def drawEventEtaPhiPlot(jetCollectionAK8, partCol, particlePDGID,METPhi, isFromHVQuark, plotNumber):
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
		if abs(particlePDGID[iPart]) == 4900101:
			objectList[-1].SetMarkerColor(3)
			objectList[-1].SetMarkerStyle(22)
			objectList[-1].SetMarkerSize(3)
		elif abs(particlePDGID[iPart] == 4900023):
			objectList[-1].SetMarkerStyle(43)
			objectList[-1].SetMarkerSize(3)
		elif abs(particlePDGID[iPart]) > 4900000:
			objectList[-1].SetMarkerColor(2)
		if isFromHVQuark[iPart]:
			objectList[-1].SetMarkerStyle(5)
	objectList.append(rt.TLine(-6,METPhi,6,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	canv.SaveAs("etaPhi_"+str(plotNumber)+".png")

def deltaPhi_calc(phi1, phi2):
	delta = phi1 - phi2
	while delta >= rt.TMath.Pi():
		delta -= 2*rt.TMath.Pi()
	while delta < -rt.TMath.Pi():
		delta += 2*rt.TMath.Pi()
	return abs(delta)


