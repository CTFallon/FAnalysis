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
	tree.SetBranchStatus("DijetEvent",1)
	tree.SetBranchStatus("TrijetEvent",1)
	tree.SetBranchStatus("OtherjetEvent",1)

	# make plots of gamma3 and iJetMaxDeltaPhi
	# distinguish between Trijet and Dijet Events
	
	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	hist_g3_All = self.makeTH1F("hist_g3_All","All;#gamma_{3};count",200,1,17)
	hist_g3_Dij = self.makeTH1F("hist_g3_Dij","Dijet;#gamma_{3};count",200,1,17)
	hist_g3_Tri = self.makeTH1F("hist_g3_Tri","Trijet;#gamma_{3};count",200,1,17)

	hist_ijmdp_All = self.makeTH1F("hist_ijmdp_All","All;iJetMaxDeltaPhi;count",5,0,5)
	hist_ijmdp_Dij = self.makeTH1F("hist_ijmdp_Dij","Dijet;iJetMaxDeltaPhi;count",5,0,5)
	hist_ijmdp_Tri = self.makeTH1F("hist_ijmdp_Tri","Trijet;iJetMaxDeltaPhi;count",5,0,5)

	hist_pt3frac_All = self.makeTH1F("hist_pt3frac_All","All;pt3frac;count",200,-0.01,1.01)
	hist_pt3frac_Dij = self.makeTH1F("hist_pt3frac_Dij","Dijet;pt3frac;count",200,-0.01,1.01)
	hist_pt3frac_Tri = self.makeTH1F("hist_pt3frac_Tri","Trijet;pt3frac;count",200,-0.01,1.01)

	hist_sosqrtb_g3 = self.makeTH1F("hist_sosqrtb_g3","S/sqrt(B) for Gamma3 cut;Gamma3;SoSqrtB",200,1,17)

	hist_2d_g3_vs_ijmdp_Dij = self.makeTH2F("hist_2d_g3_vs_ijmdp_Dij", ";jetPt3Frac;isHV3", 200,0,1,2,0,2)
	hist_2d_g3_vs_ijmdp_Tri = self.makeTH2F("hist_2d_g3_vs_ijmdp_Tri", ";jetPt3Frac;isHV3", 200,0,1,2,0,2)
	
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
		hist_pt3frac_All.Fill(tree.fracTotPTfromAllHVQ[2])
		if bool(tree.JetsAK8_isHV[2]) == False:
			hist_g3_Dij.Fill(tree.JetsAK8[2].Gamma())
			hist_ijmdp_Dij.Fill(tree.iJetMaxDeltaPhi)
			hist_2d_g3_vs_ijmdp_Dij.Fill(tree.fracTotPTfromAllHVQ[2],int(bool(tree.JetsAK8_isHV[2])))
			hist_pt3frac_Dij.Fill(tree.fracTotPTfromAllHVQ[2])
		elif bool(tree.JetsAK8_isHV[2]) == True:
			hist_g3_Tri.Fill(tree.JetsAK8[2].Gamma())
			hist_ijmdp_Tri.Fill(tree.iJetMaxDeltaPhi)
			hist_2d_g3_vs_ijmdp_Tri.Fill(tree.fracTotPTfromAllHVQ[2],int(bool(tree.JetsAK8_isHV[2])))
			hist_pt3frac_Tri.Fill(tree.fracTotPTfromAllHVQ[2])
		else:
			print("Event is neighter dijet nor trijet")

	for iBin in range(1,hist_ijmdp_All.GetNbinsX()+1):
		try:
			ratio = float(hist_ijmdp_Dij.GetBinContent(iBin))/float(hist_ijmdp_Tri.GetBinContent(iBin))
		except ZeroDivisionError:
			ratio = -1
		if ratio != 0:
			print("When jet #{} is maxDeltaPhi, it has a dijet/trijet (inverse) ratio of {:.2f} ({:.2f})".format(iBin, ratio, 1./ratio))

	for iBin in range(1,hist_g3_All.GetNbinsX()):
		s = float(hist_g3_Tri.Integral(0,iBin))
		b = float(hist_g3_Dij.Integral(0,iBin))
		sqrtB = rt.TMath.Sqrt(b)
		try:
			hist_sosqrtb_g3.SetBinContent(iBin,s/sqrtB)
		except ZeroDivisionError:
			continue

	self.makePng([hist_g3_Dij, hist_g3_Tri],"3JetMT_gamma3_gamma3", doCum = True)
	self.makePng([hist_ijmdp_Dij, hist_ijmdp_Tri],"3JetMT_gamma3_iJetMaxDeltaPhi")
	self.makePng([hist_pt3frac_Dij, hist_pt3frac_Tri],"3JetMT_gamma3_pt3frac", log = True)

	hist_g3_normCumDiff = hist_g3_Tri.GetCumulative().Clone()
	hist_g3_normCumDiff.SetTitle("Difference of Cumulative Distributions of Normalized Histograms (tri-dij)")
	hist_g3_normCumDiff.Add(-1*hist_g3_Dij.GetCumulative())
	ncdMax = hist_g3_normCumDiff.GetMaximum()
	ncdMaxBin = hist_g3_normCumDiff.GetMaximumBin()
	self.makePng([hist_g3_normCumDiff,hist_sosqrtb_g3],"3JetMT_gamma3_FOMs")
	print("Optimal normcumdiff is {} at bin {}".format(ncdMax, ncdMaxBin))

	hist_2d_g3_vs_ijmdp_Dij.SetMarkerStyle(29)
	hist_2d_g3_vs_ijmdp_Dij.SetMarkerColorAlpha(rt.kRed, .5)
	hist_2d_g3_vs_ijmdp_Tri.SetMarkerStyle(29)
	hist_2d_g3_vs_ijmdp_Tri.SetMarkerColorAlpha(rt.kGreen, .5)
	c1 = rt.TCanvas()
	hist_2d_g3_vs_ijmdp_Tri.Draw()
	hist_2d_g3_vs_ijmdp_Dij.Draw("same")
	c1.BuildLegend()
	c1.SaveAs(self.extraDir+"_aransrequest.png")

def addLoop():
	baseClass.loop = loop


