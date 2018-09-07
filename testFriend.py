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
	tree.SetBranchStatus("passedPreSelection",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	passed = self.makeTH1F("jet1Pt_passed", 100, 0, 4000)
	failed = self.makeTH1F("jet1Pt_failed", 100, 0, 4000)
	total = self.makeTH1F("jet1PT_total", 100, 0, 4000)
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		if len(tree.JetsAK8) >= 1:
			total.Fill(tree.JetsAK8[0].Pt())
			if tree.passedPreSelection == 1:
				passed.Fill(tree.JetsAK8[0].Pt())
			else:
				failed.Fill(tree.JetsAK8[0].Pt())


	

def addLoop():
	baseClass.loop = loop


