from analysisBase import baseClass
import ROOT as rt
from array import array
# this macro makes the friend tree to our analysis ntuple tree
def loop(self):
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
	iJetMaxDeltaPhi = array('i',[-1])
	pTMaxDeltaPhi = array('f',[0.])
	dPhiMaxDeltaPhi = array('f',[0.])
	

	friend.Branch("passedPreSelection",passedPreSelection, 'passedPreSelection/I')
	friend.Branch("numGenParts",numGenParts, 'numGenParts/I')
	friend.Branch("genParticleInAK8Jet",genParticleInAK8Jet,'genParticleInAK8Jet[155]/F')
	friend.Branch("genParticleIsFromHVQuark",genParticleIsFromHVQuark,'genParticleIsFromHVQuark[155]/F')
	friend.Branch("numberOfDaughtersAParticleHas",numberOfDaughtersAParticleHas,'numberOfDaughtersAParticleHas[155]/F')
	friend.Branch("fracPtFromHVQuarks",fracPtFromHVQuarks,'fracPtFromHVQuarks[7]/F')
	friend.Branch("numHVPartsInJet",numHVPartsInJet,'numHVPartsInJet[7]/I')
	friend.Branch("numSMPartsInJet",numSMPartsInJet,'numSMPartsInJet[7]/I')
	friend.Branch("iJetMaxDeltaPhi",iJetMaxDeltaPhi,'iJetMaxDeltaPhi/I')
	friend.Branch("pTMaxDeltaPhi",pTMaxDeltaPhi,'pTMaxDeltaPhi/F')
	friend.Branch("dPhiMaxDeltaPhi",dPhiMaxDeltaPhi,'dPhiMaxDeltaPhi/F')
	
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
		if maxNofParticle < numGenParts[0]:
			maxNofParticle = numGenParts[0]
		if maxNofJets < numJets[0]:
			maxNofJets = numJets[0]
		for i in range(155):
			genParticleInAK8Jet[i] = -1.
			genParticleIsFromHVQuark[i] = 0.
			numberOfDaughtersAParticleHas[i] = 0.
		for i in range(7):
			fracPtFromHVQuarks[i] = 0.
			numHVPartsInJet[i] = -1
			numSMPartsInJet[i] = -1
		for i in range(numJets[0]):
			numHVPartsInJet[i] = 0
			numSMPartsInJet[i] = 0
		# PreSelection Cuts
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
		for iPart in range(2,len(tree.GenParticles)):
			for iJet in range(len(tree.JetsAK8)-1,-1,-1):
				#print(iJet)
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and numberOfDaughtersAParticleHas[iPart] == 0:
					genParticleInAK8Jet[iPart] = float(iJet)

		# for each JET record:
		# how much of the total PT from daughter-less particles is from HV decendants
		# number of particles in the jet, total, visible, and invisible
		# also record which jet is furthest from the METPhi
		for iJet in range(len(tree.JetsAK8)):
			deltaPhi = abs(tree.JetsAK8[iJet].Phi()-tree.METPhi)%rt.TMath.Pi()
			if deltaPhi > dPhiMaxDeltaPhi[0]:
				iJetMaxDeltaPhi[0] = iJet
				pTMaxDeltaPhi[0] = tree.JetsAK8[iJet].Pt()
				dPhiMaxDeltaPhi[0] = deltaPhi
			totalPt = 0.
			HVPt = 0.
			for iPart in range(2,len(tree.GenParticles)):
				# only want final-state particles
				if numberOfDaughtersAParticleHas[iPart] >= 1:
					continue
				# only want particles that belong to the jet
				if genParticleInAK8Jet[iPart] != iJet:
					continue
				if abs(tree.GenParticles_PdgId[iPart]) < 4900000:
					numSMPartsInJet[iJet] += 1
					totalPt += tree.GenParticles[iPart].Pt()
					if genParticleIsFromHVQuark[iPart] == 1:
						HVPt += tree.GenParticles[iPart].Pt()
				else:
					numHVPartsInJet[iJet] += 1
			try:
				fracPtFromHVQuarks[iJet] = HVPt/totalPt
			except ZeroDivisionError:
				fracPtFromHVQuarks[iJet] = -1
		
					
					
		friend.Fill()
	print(maxNofParticle)
	print(maxNofJets)


def addLoop():
	baseClass.loop = loop


