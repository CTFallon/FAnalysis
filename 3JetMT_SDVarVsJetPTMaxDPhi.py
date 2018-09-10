from analysisBase import baseClass
import ROOT as rt
from array import array

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	ff = rt.TFile.Open("outFriend.root")
	friendTree = ff.Get("friend")
	tree.AddFriend(friendTree)
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
	tree.SetBranchStatus("Electrons",1)
	tree.SetBranchStatus("Muons",1)
	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	#tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
	#tree.SetBranchStatus("numHVPartsInJet",1)
	#tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	# check overlap between jetPt(maxdPhi) and SDVar
	hist_overlap = self.makeTH2F("hist_overlap","Overlap;123 from jetPt(maxdPhi);123 from SDVar_13",2,0,2,2,0,2)
	
	# MT distributions
	hist_MT_12_all = self.makeTH1F("hist_MT_12_all","Dijet MT - all;MT;count/a.u.",100,0,4000)
	hist_MT_123_all = self.makeTH1F("hist_MT_123_all","Trijet MT - all;MT;count/a.u.",100,0,4000)
	
	hist_MT_12_jetPtMaxDPhi = self.makeTH1F("hist_MT_12_jetPtMaxDPhi","Dijet MT - pT(dPhi);MT;count/a.u.",100,0,4000)
	hist_MT_123_jetPtMaxDPhi = self.makeTH1F("hist_MT_123_jetPtMaxDPhi","Trijet MT - pT(dPhi);MT;count/a.u.",100,0,4000)

	hist_MT_12_SDVar13 = self.makeTH1F("hist_MT_12_SDVar13","Dijet MT - SDVar13;MT;count/a.u.",100,0,4000)
	hist_MT_123_SDVar13 = self.makeTH1F("hist_MT_123_SDVar13","Trijet MT - SDVar13;MT;count/a.u.",100,0,4000)
	
	
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		nJets = len(tree.JetsAK8)
		if nJets == 2:
			# fill all histograms for just 2 jet MT
			hist_MT_12_all.Fill()
		else:
			# process algorithm for 12 vs 123 jet MT


	

def addLoop():
	baseClass.loop = loop


