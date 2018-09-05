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
	passednJets = array('i', [0])
	passedHighPt = array('i', [0])
	passedMETMTRatio = array('i', [0])
	passedLeptonVeto = array('i', [0])
	passedPreSelection = array('i', [0])
	
	genParticleInAK8Jet = array('i' ,[-10 for x in range(100)])
	genParticleIsFromHVQuark = array('i' ,[-10 for x in range(100)])
	
	friend.Branch("passednJets", passednJets, 'passednJets/I')
	friend.Branch("passedHighPt", passedHighPt, 'passedHighPt/I')
	friend.Branch("passedMETMTRatio", passedMETMTRatio, 'passedMETMTRatio/I')
	friend.Branch("passedLeptonVeto", passedLeptonVeto, 'passedLeptonVeto/I')
	friend.Branch("passedPreSelection", passedPreSelection, 'passedPreSelection/I')
	friend.Branch("ganParticleInAK8Jet", genParticleInAK8Jet,'genParticleInAK8Jet[iPart]/I')
	friend.Branch("genParticleIsFromHVQuark", genParticleIsFromHVQuark,'genParticleIsFromHVQuark[iPart]/I')
	
	tree.AddFriend(friend)
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		passednJets[0] = 0
		passedHighPt[0] = 0
		passedMETMTRatio[0] = 0
		passedLeptonVeto[0] = 0
		passedPreSelection[0] = 0
		for i in range(100):
			genParticleInAK8Jet[i] = =10
			genParticleIsFromHVQuark[i] = -10
		# PreSelection Cuts
		# At least 2 jets in the event
		if not (len(tree.JetsAK8)>=2): 
			continue
		passednJets[0] = 1
		# Both leading jets must have pt > 170
		if not ((tree.JetsAK8[0].Pt() > 170.0) and (tree.JetsAK8[1].Pt() > 170.0)): 
			continue
		passedHighPt[0] = 1
		# MET/MT ratio must be greater than 0.15
		if not (tree.MET/tree.MT_AK8 > 0.15):
			continue
		passedMETMTRatio[0] = 1
		# must not have any leptons
		if not ((len(tree.Electrons) + len(tree.Muons)) == 0):
			continue
		passedLeptonVeto[0] = 1
		passedPreSelection[0] = 1
		
		# for each PARTICLE:
		# record:
	 	#   if it is inside which AK8 jet
		#   if it is decendant of an HV quark
		numberOfDaughtersAParticleHas = [0 for part in tree.GenParticles]
		for iPart in range(2,len(tree.GenParticles)):
			iParent = tree.GenParticles_ParentIdx[iPart]
			if iParent != -1:
				numberOfDaughtersAParticleHas[iParent] += 1
			if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or genParticleIsFromHVQuark[iParent]:
				genParticleIsFromHVQuark[iPart] = 1
		for iPart in range(len(tree.GenParticles)):
			for iJet in range(len(tree.JetsAK8)-1,-1,-1):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) > 0.08 and numberOfDaughtersAParticleHas[iPart] == 0:
					genParticleInAK8Jet[iPart] = iJet
		friend.Fill()

	

def addLoop():
	baseClass.loop = loop


