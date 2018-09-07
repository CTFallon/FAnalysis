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
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	#Step 1, make iJet vs FracPt for events with nJets
	histList_2d_iJetvsFracPt = [0,0]
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_2Jets", "Events with 2 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 2, 0, 2, 100, -.1, 1.1))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_3Jets", "Events with 3 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 3, 0, 3, 100, -.1, 1.1))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_4Jets", "Events with 4 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 4, 0, 4, 100, -.1, 1.1))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_5Jets", "Events with 5 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 5, 0, 5, 100, -.1, 1.1))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_6Jets", "Events with 6 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 6, 0, 6, 100, -.1, 1.1))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_7Jets", "Events with 7 Jets; Jet Number; Fraction of Pt from Visible HV Decendants", 7, 0, 7, 100, -.1, 1.1))
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		if tree.passedPreSelection == 1:
			nJets = len(tree.JetsAK8)
			for iJet in range(nJets):
				histList_2d_iJetvsFracPt[nJets].Fill(iJet+0.5, tree.fracPtFromHVQuarks[iJet])


	

def addLoop():
	baseClass.loop = loop


