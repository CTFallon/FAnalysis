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
	genParticleInAK8Jet = array('f' ,[-1. for x in range(155)])
	genParticleIsFromHVQuark = array('f' ,[0. for x in range(155)])
	numberOfDaughtersAParticleHas = array('f', [0. for x in range(155)])
	fracPtFromHVQuarks = array("f",[0. for x in range(10)])
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
	

	friend.Branch("passedPreSelection",passedPreSelection, 'passedPreSelection/I')
	friend.Branch("numGenParts",numGenParts, 'numGenParts/I')
	friend.Branch("genParticleInAK8Jet",genParticleInAK8Jet,'genParticleInAK8Jet[155]/F')
	friend.Branch("genParticleIsFromHVQuark",genParticleIsFromHVQuark,'genParticleIsFromHVQuark[155]/F')
	friend.Branch("numberOfDaughtersAParticleHas",numberOfDaughtersAParticleHas,'numberOfDaughtersAParticleHas[155]/F')
	friend.Branch("fracPtFromHVQuarks",fracPtFromHVQuarks,'fracPtFromHVQuarks[10]/F')
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
	
	tree.AddFriend(friend)
	maxNofParticle = 0
	maxNofJets = 0 
	for iEvent in range(nEvents):
		if iEvent%1000==0:
			print("Event "+ str(iEvent)+"/"+str(nEvents))
		tree.GetEvent(iEvent)
		passedPreSelection[0] = 0
		numGenParts[0] = len(tree.GenParticles)
		numJets[0] = len(tree.JetsAK8)
		iJetMaxDeltaPhi[0] = -1
		pTMaxDeltaPhi[0] = 0.
		dPhiMaxDeltaPhi[0] = 0.
		zPrimept[0] = 0.
		zPrimephi[0] = 0.
		if maxNofParticle < numGenParts[0]:
			maxNofParticle = numGenParts[0]
		if maxNofJets < numJets[0]:
			maxNofJets = numJets[0]
		for i in range(155):
			genParticleInAK8Jet[i] = -1.
			genParticleIsFromHVQuark[i] = 0.
			numberOfDaughtersAParticleHas[i] = 0.
		for i in range(10):
			fracPtFromHVQuarks[i] = -1.
			numHVPartsInJet[i] = -1
			numSMPartsInJet[i] = -1
		for i in range(numJets[0]):
			fracPtFromHVQuarks[i] = 0.
			numHVPartsInJet[i] = 0
			numSMPartsInJet[i] = 0
		for i in range(10):
			pGJ_visible[i] = -1.
			pGJ_invis[i] = -1.
			pGJ_every[i] = -1.
		pass1 = 0
		pass2 = 0
		pass3 = 0
		pass4 = 0
		# At least 2 jets in the event, temp 3 for 3JetMT purposes
		if (len(tree.JetsAK8)>=2): 
			pass1 = 1
		else: # special case because next chcek can cause Index Error
			passedPreSelection[0] = 0
			friend.Fill()
			continue
		# Both leading jets must have pt > 170
		if ((tree.JetsAK8[0].Pt() > 170.0) and (tree.JetsAK8[1].Pt() > 170.0)): 
			pass2 = 1
		# MET/MT ratio must be greater than 0.15
		if (tree.MET/tree.MT_AK8 > 0.15):
			pass3 = 1
		# must not have any leptons
		if ((len(tree.Electrons) + len(tree.Muons)) == 0):
			pass4 = 1
		if pass1 and pass2 and pass3 and pass4:
			passedPreSelection[0] = 1
		else:
			passedPreSelection[0] = 0
			friend.Fill()
			continue
		
		# for each PARTICLE:
		# record:
	 	#   if it is inside which AK8 jet
		#   if it is decendant of an HV quark
		for iPart in range(2,len(tree.GenParticles)):
			iParent = tree.GenParticles_ParentIdx[iPart]
			if iParent != -1:
				numberOfDaughtersAParticleHas[iParent] += 1.
			if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or genParticleIsFromHVQuark[iParent]:
				genParticleIsFromHVQuark[iPart] = float(1)
			if tree.GenParticles_PdgId[iPart] == 4900023:
				zPrimept[0] = tree.GenParticles[iPart].Pt()
				zPrimephi[0] = tree.GenParticles[iPart].Phi()
		
		# calculate the 'final-state' MT from all visible, daughterless particles that are decendants of HV particls
		# also determin what jets particles belong to
		particlesForMT = []
		for iPart in range(2,len(tree.GenParticles)):
			if genParticleIsFromHVQuark[iPart] and numberOfDaughtersAParticleHas[iPart] == 0:
				particlesForMT.append(tree.GenParticles[iPart])
			for iJet in range(len(tree.JetsAK8)-1,-1,-1):
				#print(iJet)
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and numberOfDaughtersAParticleHas[iPart] == 0:
					genParticleInAK8Jet[iPart] = float(iJet)
		
		MTFromParticles[0] = trans_mass_Njet(particlesForMT, tree.MET, tree.METPhi)
		# for each JET record:
		# how much of the total PT from daughter-less particles is from HV decendants
		# number of particles in the jet, total, visible, and invisible
		# also record which jet is furthest from the METPhi
		for iJet in range(len(tree.JetsAK8)):
			pGJ_vis = rt.TLorentzVector(0.,0.,0.,0.)
			pGJ_inv = rt.TLorentzVector(0.,0.,0.,0.)
			pGJ_all = rt.TLorentzVector(0.,0.,0.,0.)
			deltaPhi = abs(tree.JetsAK8[iJet].Phi()-tree.METPhi)%rt.TMath.Pi()
			if deltaPhi > dPhiMaxDeltaPhi[0]:
				iJetMaxDeltaPhi[0] = iJet
				pTMaxDeltaPhi[0] = tree.JetsAK8[iJet].Pt()
				dPhiMaxDeltaPhi[0] = deltaPhi
			totalPt = 0.
			HVPt = 0.
			for iPart in range(2,len(tree.GenParticles)):
				# only want final-state particles
				# only want particles that belong to the jet
				# only want particles that descend from HVQuarks
				# ignore particles that a) have daughters or b) arn't in the jet
				if (numberOfDaughtersAParticleHas[iPart] != 0):
					continue
				if (genParticleInAK8Jet[iPart] != iJet):
					continue
				# right now, we have all particles that are final state and in the jet
				# so we want to fill our pseudoGenJets:
				pGJ_all += tree.GenParticles[iPart]
				if abs(tree.GenParticles_PdgId[iPart] < 4900000): # visible
					pGJ_vis += tree.GenParticles[iPart]
				else: # invisible
					pGJ_inv = pGJ_inv + tree.GenParticles[iPart]
				# now, we want to ignore any particle that isn't from an HVQuark:
				if (genParticleIsFromHVQuark[iPart] != 1):
					continue
				if abs(tree.GenParticles_PdgId[iPart]) < 4900000: # finally, ignore all invisible particles
					numSMPartsInJet[iJet] += 1
					totalPt += tree.GenParticles[iPart].Pt()
					if genParticleIsFromHVQuark[iPart] == 1.:
						HVPt += tree.GenParticles[iPart].Pt()
				else:
					numHVPartsInJet[iJet] += 1
			try:
				fracPtFromHVQuarks[iJet] = HVPt/totalPt
			except ZeroDivisionError:
				fracPtFromHVQuarks[iJet] = -1
			pGJ_visible[iJet] = pGJ_vis.Pt()
			pGJ_invis[iJet] = pGJ_inv.Pt()
			pGJ_every[iJet] = pGJ_all.Pt()
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


