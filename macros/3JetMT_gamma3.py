from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

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
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	#tree.SetBranchStatus("fracPtFromHVQuarks",1)
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

	# make plots of gamma3 and iJetMaxDeltaPhi
	# distinguish between Trijet and Dijet Events
	
	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	hist_g3_All = self.makeTH1F("hist_g3_All","All;#gamma_{3};count",100,0,15)
	hist_g3_Dij = self.makeTH1F("hist_g3_Dij","Dijet;#gamma_{3};count",100,0,15)
	hist_g3_Tri = self.makeTH1F("hist_g3_Tri","Trijet;#gamma_{3};count",100,0,15)

	hist_ijmdp_All = self.makeTH1F("hist_ijmdp_All","All;iJetMaxDeltaPhi;count",5,0,5)
	hist_ijmdp_Dij = self.makeTH1F("hist_ijmdp_Dij","Dijet;iJetMaxDeltaPhi;count",5,0,5)
	hist_ijmdp_Tri = self.makeTH1F("hist_ijmdp_Tri","Trijet;iJetMaxDeltaPhi;count",5,0,5)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		if len(tree.JetsAK8) < 3:
			continue
		hist_g3_All.Fill(tree.JetsAK8[2].Gamma())
		hist_ijmdp_All.Fill(tree.iJetMaxDeltaPhi)
		if bool(tree.JetsAK8_isHV[2]) == False:
			hist_g3_Dij.Fill(tree.JetsAK8[2].Gamma())
			hist_ijmdp_Dij.Fill(tree.iJetMaxDeltaPhi)
		elif bool(tree.JetsAK8_isHV[2]) == True:
			hist_g3_Tri.Fill(tree.JetsAK8[2].Gamma())
			hist_ijmdp_Tri.Fill(tree.iJetMaxDeltaPhi)
		else:
			print("Event is neighter dijet nor trijet!")

	self.makePng([hist_g3_All, hist_g3_Dij, hist_g3_Tri],"3JetMT_gamma3_gamma3")
	self.makePng([hist_ijmdp_All, hist_ijmdp_Dij, hist_ijmdp_Tri],"3JetMT_gamma3_iJetMaxDeltaPhi")

def addLoop():
	baseClass.loop = loop


