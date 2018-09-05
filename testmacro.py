from analysisBase import baseClass
import ROOT as rt
from array import array

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

	# initalize histograms to be made, or create Friend tree to be filled
	friend = rt.TTree("friend","friend")
	self.objects.append(friend)
	# create branch variables
	passednJets = array(f, [0.])
	passedHighPt = array(f, [0.])
	passedMETMTRatio = array(f, [0.])
	passedLeptonVeto = array(f, [0.])
	passedPreSelection = array(f, [0.])
	
	#genParticleInAK8Jet = [0 for partile in len(tree.GenParticles)]
	
	friend.Branch("passednJets", passednJets, 'passednJets/I')
	friend.Branch("passedHighPt", passedHighPt, 'passedHighPt/I')
	friend.Branch("passedMETMTRatio", passedMETMTRatio, 'passedMETMTRatio/I')
	friend.Branch("passedLeptonVeto", passedLeptonVeto, 'passedLeptonVeto/I')
	friend.Branch("passedPreSelection", passedPreSelection, 'passedPreSelection/I')
	
	tree.AddFriend(friend)
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		passednJets[0] = 0.
		passedHighPt[0] = 0.
		passedMETMTRatio[0] = 0.
		passedLeptonVeto[0] = 0.
		passedPreSelection[0] = 0.
		# PreSelection Cuts
		# At least 2 jets in the event
		if not (len(jetCollection)>=2): 
			continue
		passednJets[0] = 1.
		# Both leading jets must have pt > 170
		if not ((jetCollection[0].Pt() > 170.0) and (jetCollection[1].Pt() > 170.0)): 
			continue
		passedHighPt[0] = 1.
		# MET/MT ratio must be greater than 0.15
		if not (tree.MET/tree.MT_AK8 > 0.15):
			continue
		passedMETMTRatio[0] = 1.
		# must not have any leptons
		if not ((len(tree.Electrons) + len(tree.Muons)) == 0):
			continue
		passedLeptonVeto[0] = 1.
		passedPreSelection[0] = 1.
		
		friend.Fill()

	

def addLoop():
	baseClass.loop = loop


