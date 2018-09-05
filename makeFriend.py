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
	tree.SetBranchStatus("GenParticles*",1)
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
	
	#genParticleInAK8Jet = [0 for partile in len(tree.GenParticles)]
	
	friend.Branch("passednJets", passednJets, 'passednJets/I')
	friend.Branch("passedHighPt", passedHighPt, 'passedHighPt/I')
	friend.Branch("passedMETMTRatio", passedMETMTRatio, 'passedMETMTRatio/I')
	friend.Branch("passedLeptonVeto", passedLeptonVeto, 'passedLeptonVeto/I')
	friend.Branch("passedPreSelection", passedPreSelection, 'passedPreSelection/I')
	
	tree.AddFriend(friend)
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		passednJets[0] = 0
		passedHighPt[0] = 0
		passedMETMTRatio[0] = 0
		passedLeptonVeto[0] = 0
		passedPreSelection[0] = 0
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
		# copied from old stuff, probably lots of errors in terms of varible names
		# UPDATE
		for iPart in range(2, len(tree.GenParticles)):
			iParent = tree.GenParticles_ParentIdx[iPart]
			if iParent != -1: 
				numberOfDaughtersAParticleHas[iParent] += 1
		AK8jetsWithHVDecendants = ["0","0","0","0","0"]
		AK8_nHVParts = [0,0,0,0,0]
		AK8_nParts = [0,0,0,0,0]
		AK8_ptHVParts = [0,0,0,0,0]
		AK8_ptParts = [0,0,0,0,0]
		#make vector of length nGenParts that is 1 if the particle came from a HV quark
		isFromHVQuark = [0 for x in range(len(tree.GenParticles))]
		listOfHVQuarks = []
		for iPart in range(2,len(tree.GenParticles)):
			iParent = tree.GenParticles_ParentIdx[iPart]
			if abs(tree.GenParticles_PdgId[iPart]) == 4900101:
				listOfHVQuarks.append(tree.GenParticles[iPart])
			if iParent >= iPart:
				print("Ut-oh, the parent has a higher index than the child...")
			if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or (isFromHVQuark[iParent]):
				isFromHVQuark[iPart] = 1
			#finding what Jet a particle is in, only if it decends from a HVQuark:
			for iJet in range(min(len(tree.JetsAK8),5)):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
					AK8_ptParts[-iJet-1] += tree.GenParticles[iPart].Pt()
					AK8_nParts[-iJet-1] += 1
			if isFromHVQuark[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
				for iJet in range(min(len(tree.JetsAK8),5)):
					if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart]:
						AK8jetsWithHVDecendants[-iJet-1] = "1"
						AK8_ptHVParts[-iJet-1] += tree.GenParticles[iPart].Pt()
						AK8_nHVParts[-iJet-1] += 1
			#UPDATE
		friend.Fill()

	

def addLoop():
	baseClass.loop = loop


