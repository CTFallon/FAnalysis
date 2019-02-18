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

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()

	# What this macro should do:
		# create two histograms for various variables, one for trijet, one for dijet
		# compute the optimal cut to seperate the two
		# report in such a way as to easily be exported to Excel

	# variabels to consider
		# gamma3
		# jetIndex (special case, two sided cut)
		# jetPtMaxDetlaPhi
		# jetPt3
		# deltaEta13

	hist_gamma3_Tri = self.makeTH1F("hist_gamma3_Tri","Trijet;gamma3;count",100,0,17)
	hist_gamma3_Dij = self.makeTH1F("hist_gamma3_Dij","Dijet;gamma3;count",100,0,17)

	hist_jetIdx_Tri = self.makeTH1F("hist_jetIdx_Tri","Trijet;jetIdx;count",5,0,5)
	hist_jetIdx_Dij = self.makeTH1F("hist_jetIdx_Dij","Dijet;jetIdx;count",5,0,5)
	
	hist_jptmdp_Tri = self.makeTH1F("hist_jptmdp_Tri","Trijet;jptmdp;count",100,0,2000)
	hist_jptmdp_Dij = self.makeTH1F("hist_jptmdp_Dij","Dijet;jptmdp;count",100,0,2000)

	hist_jetPt3_Tri = self.makeTH1F("hist_jetPt3_Tri","Trijet;jetPt3;count",100,0,1200)
	hist_jetPt3_Dij = self.makeTH1F("hist_jetPt3_Dij","Dijet;jetPt3;count",100,0,1200)
	
	hist_dEta23_Tri = self.makeTH1F("hist_dEta23_Tri","Trijet;dEta23;count",100,-0.01,5)
	hist_dEta23_Dij = self.makeTH1F("hist_dEta23_Dij","Dijet;dEta23;count",100,-0.01,5)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		if len(tree.JetsAK8) <= 2: # skip events where there arn't atleast 3 jets
			continue
		if tree.iJetMaxDeltaPhi == 1: # skip events taken out by the idx cut
			continue

		if bool(tree.JetsAK8_isHV[2]) == True:
			hist_gamma3_Tri.Fill(tree.JetsAK8[2].Gamma())
			hist_jetIdx_Tri.Fill(tree.iJetMaxDeltaPhi)
			hist_jptmdp_Tri.Fill(tree.pTMaxDeltaPhi)
			hist_jetPt3_Tri.Fill(tree.JetsAK8[2].Pt())
			hist_dEta23_Tri.Fill(abs(tree.JetsAK8[1].Eta()-tree.JetsAK8[2].Eta()))
		elif bool(tree.JetsAK8_isHV[2]) == False:
			hist_gamma3_Dij.Fill(tree.JetsAK8[2].Gamma())
			hist_jetIdx_Dij.Fill(tree.iJetMaxDeltaPhi)
			hist_jptmdp_Dij.Fill(tree.pTMaxDeltaPhi)
			hist_jetPt3_Dij.Fill(tree.JetsAK8[2].Pt())
			hist_dEta23_Dij.Fill(abs(tree.JetsAK8[1].Eta()-tree.JetsAK8[2].Eta()))
		else:
			print("Jet 3 is neither hv nor not Hv")

	self.makePng([hist_gamma3_Tri,hist_gamma3_Dij],"3JetMT_Final_vars_gamma3", doCum = True)
	#self.makePng([hist_jetIdx_Tri,hist_jetIdx_Dij],"3JetMT_Final_vars_jetIdx")
	self.makePng([hist_jptmdp_Tri,hist_jptmdp_Dij],"3JetMT_Final_vars_jptmdp", doCum = True)
	self.makePng([hist_jetPt3_Tri,hist_jetPt3_Dij],"3JetMT_Final_vars_jetPt3", doCum = True)
	self.makePng([hist_dEta23_Tri,hist_dEta23_Dij],"3JetMT_Final_vars_dEta23", doCum = True)

	hist_csd1 = hist_gamma3_Tri.GetCumulative().Clone()
	hist_csd1.Add(-1*hist_gamma3_Dij.GetCumulative())
	print("signal var min max minCut maxCut nDij nTri")
	print(self.fileID + " gamma3 {} {} {} {} {} {}".format(
			hist_csd1.GetMinimum(),hist_csd1.GetMaximum(),
			hist_csd1.GetMinimumBin(),hist_csd1.GetMaximumBin(),
			hist_gamma3_Dij.GetEntries(),hist_gamma3_Tri.GetEntries()))

	hist_csd2 = hist_jptmdp_Tri.GetCumulative().Clone()
	hist_csd2.Add(-1*hist_jptmdp_Dij.GetCumulative())
	print(self.fileID + " jptmdp {} {} {} {} {} {}".format(
			hist_csd2.GetMinimum(),hist_csd2.GetMaximum(),
			hist_csd2.GetMinimumBin(),hist_csd2.GetMaximumBin(),
			hist_jptmdp_Dij.GetEntries(),hist_jptmdp_Tri.GetEntries()))

	hist_csd3 = hist_jetPt3_Tri.GetCumulative().Clone()
	hist_csd3.Add(-1*hist_jetPt3_Dij.GetCumulative())
	print(self.fileID + " jetPt3 {} {} {} {} {} {}".format(
			hist_csd3.GetMinimum(),hist_csd3.GetMaximum(),
			hist_csd3.GetMinimumBin(),hist_csd3.GetMaximumBin(),
			hist_jetPt3_Dij.GetEntries(),hist_jetPt3_Tri.GetEntries()))

	hist_csd4 = hist_dEta23_Tri.GetCumulative().Clone()
	hist_csd4.Add(-1*hist_dEta23_Dij.GetCumulative())
	print(self.fileID + " dEta23 {} {} {} {} {} {}".format(
			hist_csd4.GetMinimum(),hist_csd4.GetMaximum(),
			hist_csd4.GetMinimumBin(),hist_csd4.GetMaximumBin(),
			hist_dEta23_Dij.GetEntries(),hist_dEta23_Tri.GetEntries()))
	self.makePng([hist_csd1],"3JetMT_Final_vars_cumDiff_gamma3")
	self.makePng([hist_csd2],"3JetMT_Final_vars_cumDiff_jptmdp")
	self.makePng([hist_csd3],"3JetMT_Final_vars_cumDiff_jetPt3")
	self.makePng([hist_csd4],"3JetMT_Final_vars_cumDiff_dEta23")
	

def addLoop():
	baseClass.loop = loop
