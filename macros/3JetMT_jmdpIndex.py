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


	hist_2d_jpmdp_vs_sum3jetpt = self.makeTH2F(
		"hist_2d_jpmdp_vs_sum3jetpt",
		"Jet Pt Max DPhi vs Sum Pt 3 Leading Jets;MaxDPhi Pt; Sum 3 Pt",
		100,0,2000,100,0,4000)
	hist_2d_jpmdp_vs_sum3jetpt_dijet = self.makeTH2F(
		"hist_2d_jpmdp_vs_sum3jetpt_dijet",
		"Dijet Jet Pt Max DPhi vs Sum Pt 3 Leading Jets;MaxDPhi Pt; Sum 3 Pt",
		100,0,2000,100,0,4000)
	hist_2d_jpmdp_vs_sum3jetpt_trijet = self.makeTH2F(
		"hist_2d_jpmdp_vs_sum3jetpt_trijet",
		"Trijet Jet Pt Max DPhi vs Sum Pt 3 Leading Jets;MaxDPhi Pt; Sum 3 Pt",
		100,0,2000,100,0,4000)

	hist_2d_jpmdp_vs_MT = self.makeTH2F(
		"hist_2d_jpmdp_vs_MT",
		"Jet Pt Max DPhi vs MT;MaxDPhi Pt; MT",
		100,0,2000,100,0,4000)
	hist_2d_jpmdp_vs_MT_dijet = self.makeTH2F(
		"hist_2d_jpmdp_vs_MT_dijet",
		"Dijet Jet Pt Max DPhi vs MT;MaxDPhi Pt; MT",
		100,0,2000,100,0,4000)
	hist_2d_jpmdp_vs_MT_trijet = self.makeTH2F(
		"hist_2d_jpmdp_vs_MT_trijet",
		"Trijet Jet Pt Max DPhi vs MT;MaxDPhi Pt; MT",
		100,0,2000,100,0,4000)

	hist_2d_sum3jetpt_vs_MT = self.makeTH2F(
		"hist_2d_sum3jetpt_vs_MT",
		"sum3jetpt vs MT;sum3jetpt; MT",
		100,0,4000,100,0,4000)
	hist_2d_sum3jetpt_vs_MT_dijet = self.makeTH2F(
		"hist_2d_sum3jetpt_vs_MT_dijet",
		"Dijet sum3jetpt vs MT;sum3jetpt; MT",
		100,0,4000,100,0,4000)
	hist_2d_sum3jetpt_vs_MT_trijet = self.makeTH2F(
		"hist_2d_sum3jetpt_vs_MT_trijet",
		"Trijet sum3jetpt vs MT;sum3jetpt; MT",
		100,0,4000,100,0,4000)

	hist_ptRatio_dijet = self.makeTH1F("hist_ptRatio_dijet",
		"Jet Pt Max DPhi Over Sum Pt 3 Leading Jets, Dijet Events;ratio;count",
		100,0,1)
	hist_ptRatio_trijet = self.makeTH1F("hist_ptRatio_trijet",
		"Jet Pt Max DPhi Over Sum Pt 3 Leading Jets, Trijet Events;ratio;count",
		100,0,1)


	var1Delta = 5
	var1Low = 0
	var1Up = 5
	var2Delta = 5
	var2Low = 0
	var2Up = 5
	var3Delta = 100
	var3Low = -5
	var3Up = 5
	var4Delta = 100
	var4Low = -rt.TMath.Pi()
	var4Up = rt.TMath.Pi()

	# for jet 1 is maxdPhi
	hist_2d_comboA_jet1mdphi = self.makeTH2F(
		"hist_2d_comboA_jet1mdphi",
		"Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_dijet_jet1mdphi = self.makeTH2F(
		"hist_2d_comboA_dijet_jet1mdphi",
		"Dijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_trijet_jet1mdphi = self.makeTH2F(
		"hist_2d_comboA_trijet_jet1mdphi",
		"Trijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	# for jet 2 is maxdPhi
	hist_2d_comboA_jet2mdphi = self.makeTH2F(
		"hist_2d_comboA_jet2mdphi",
		"Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_dijet_jet2mdphi = self.makeTH2F(
		"hist_2d_comboA_dijet_jet2mdphi",
		"Dijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_trijet_jet2mdphi = self.makeTH2F(
		"hist_2d_comboA_trijet_jet2mdphi",
		"Trijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	# for jet 3 is maxdPhi
	hist_2d_comboA_jet3mdphi = self.makeTH2F(
		"hist_2d_comboA_jet3mdphi",
		"Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_dijet_jet3mdphi = self.makeTH2F(
		"hist_2d_comboA_dijet_jet3mdphi",
		"Dijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)
	hist_2d_comboA_trijet_jet3mdphi = self.makeTH2F(
		"hist_2d_comboA_trijet_jet3mdphi",
		"Trijet Jet Pt Max DPhi vs sum3jetpt;v1; v2",
		var1Delta,var1Low,var1Up,var2Delta,var2Low,var2Up)

	# for jet 1 is maxdPhi
	hist_2d_comboB_jet1mdphi = self.makeTH2F(
		"hist_2d_comboB_jet1mdphi",
		"Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_dijet_jet1mdphi = self.makeTH2F(
		"hist_2d_comboB_dijet_jet1mdphi",
		"Dijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_trijet_jet1mdphi = self.makeTH2F(
		"hist_2d_comboB_trijet_jet1mdphi",
		"Trijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	# for jet 2 is maxdPhi
	hist_2d_comboB_jet2mdphi = self.makeTH2F(
		"hist_2d_comboB_jet2mdphi",
		"Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_dijet_jet2mdphi = self.makeTH2F(
		"hist_2d_comboB_dijet_jet2mdphi",
		"Dijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_trijet_jet2mdphi = self.makeTH2F(
		"hist_2d_comboB_trijet_jet2mdphi",
		"Trijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	# for jet 3 is maxdPhi
	hist_2d_comboB_jet3mdphi = self.makeTH2F(
		"hist_2d_comboB_jet3mdphi",
		"Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_dijet_jet3mdphi = self.makeTH2F(
		"hist_2d_comboB_dijet_jet3mdphi",
		"Dijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)
	hist_2d_comboB_trijet_jet3mdphi = self.makeTH2F(
		"hist_2d_comboB_trijet_jet3mdphi",
		"Trijet Jet Pt Max DPhi vs MT;v3;v4",
		var3Delta,var3Low,var3Up,var4Delta,var4Low,var4Up)

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		if nJets < 3:
			continue
		nHVJets = 0
		for iJet in range(nJets):
			if bool(tree.JetsAK8_isHV[iJet]):
				nHVJets += 1
		
		
		sum3jetPt = tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt()
		jpmdpIdx = tree.iJetMaxDeltaPhi
		mt = trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi)
		hist_2d_jpmdp_vs_sum3jetpt.Fill(tree.pTMaxDeltaPhi,sum3jetPt)
		hist_2d_sum3jetpt_vs_MT.Fill(sum3jetPt,mt)
		hist_2d_jpmdp_vs_MT.Fill(tree.pTMaxDeltaPhi,mt)
		
		var1 = nJets
		var2 = nHVJets
		var3 = tree.JetsAK8[1].Eta()
		var4 = tree.JetsAK8[1].Phi()
		print(var1)
		
		if jpmdpIdx == 0:
			hist_2d_comboB_jet1mdphi.Fill(var1,var2)
			hist_2d_comboA_jet1mdphi.Fill(var3,var4)
		elif jpmdpIdx == 1:
			hist_2d_comboB_jet2mdphi.Fill(var1,var2)
			hist_2d_comboA_jet2mdphi.Fill(var3,var4)
		elif jpmdpIdx == 2:
			hist_2d_comboB_jet3mdphi.Fill(var1,var2)
			hist_2d_comboA_jet3mdphi.Fill(var3,var4)
		if bool(tree.JetsAK8_isHV[2]) == False:
			hist_2d_jpmdp_vs_sum3jetpt_dijet.Fill(tree.pTMaxDeltaPhi,sum3jetPt)
			hist_2d_jpmdp_vs_MT_dijet.Fill(tree.pTMaxDeltaPhi,mt)
			hist_2d_sum3jetpt_vs_MT_dijet.Fill(sum3jetPt,mt)
			hist_ptRatio_dijet.Fill(tree.pTMaxDeltaPhi/sum3jetPt)
			if jpmdpIdx == 0:
				hist_2d_comboB_dijet_jet1mdphi.Fill(var1,var2)
				hist_2d_comboA_dijet_jet1mdphi.Fill(var3,var4)
			elif jpmdpIdx == 1:
				hist_2d_comboB_dijet_jet2mdphi.Fill(var1,var2)
				hist_2d_comboA_dijet_jet2mdphi.Fill(var3,var4)
			elif jpmdpIdx == 2:
				hist_2d_comboB_dijet_jet3mdphi.Fill(var1,var2)
				hist_2d_comboA_dijet_jet3mdphi.Fill(var3,var4)
		if bool(tree.JetsAK8_isHV[2]) == True:
			hist_2d_jpmdp_vs_sum3jetpt_trijet.Fill(tree.pTMaxDeltaPhi,sum3jetPt)
			hist_2d_jpmdp_vs_MT_trijet.Fill(tree.pTMaxDeltaPhi,mt)
			hist_2d_sum3jetpt_vs_MT_trijet.Fill(sum3jetPt,mt)
			hist_ptRatio_trijet.Fill(tree.pTMaxDeltaPhi/sum3jetPt)
			if jpmdpIdx == 0:
				hist_2d_comboB_trijet_jet1mdphi.Fill(var1,var2)
				hist_2d_comboA_trijet_jet1mdphi.Fill(var3,var4)
			elif jpmdpIdx == 1:
				hist_2d_comboB_trijet_jet2mdphi.Fill(var1,var2)
				hist_2d_comboA_trijet_jet2mdphi.Fill(var3,var4)
			elif jpmdpIdx == 2:
				hist_2d_comboB_trijet_jet3mdphi.Fill(var1,var2)
				hist_2d_comboA_trijet_jet3mdphi.Fill(var3,var4)

	c1 = rt.TCanvas("c1","c1",1200,900)
	c1.Clear()
	c1.Divide(3,2)
	c1.cd(1)
	leg = rt.TLegend(0.2,0.2,0.8,0.8)
	dijetLeg = rt.TH1F("dijetLeg","Dijet Events",2,-100,-99)
	trijetLeg = rt.TH1F("trijetLeg","Trijet Events",2,-100,-99)
	dijetLeg.SetLineColor(rt.kBlue)
	trijetLeg.SetLineColor(rt.kRed)
	leg.AddEntry(dijetLeg)
	leg.AddEntry(trijetLeg)
	leg.Draw()
	c1.cd(4)
	hist_dijet_xProj = hist_2d_jpmdp_vs_sum3jetpt_dijet.ProjectionX()
	hist_dijet_xProj.SetLineColor(rt.kBlue)
	hist_trijet_xProj = hist_2d_jpmdp_vs_sum3jetpt_trijet.ProjectionX()
	hist_trijet_xProj.SetLineColor(rt.kRed)
	hist_trijet_xProj.Draw()
	hist_dijet_xProj.Draw("same")
	c1.cd(6)
	hist_dijet_yProj = hist_2d_jpmdp_vs_sum3jetpt_dijet.ProjectionY()
	hist_dijet_yProj.SetLineColor(rt.kBlue)
	hist_trijet_yProj = hist_2d_jpmdp_vs_sum3jetpt_trijet.ProjectionY()
	hist_trijet_yProj.SetLineColor(rt.kRed)
	hist_trijet_yProj.Draw()
	hist_dijet_yProj.Draw("same")
	c1.cd(5)
	hist_ptRatio_dijet.SetLineColor(rt.kBlue)
	hist_ptRatio_trijet.SetLineColor(rt.kRed)
	hist_ptRatio_trijet.Draw()
	hist_ptRatio_dijet.Draw("same")
	c1.cd(2)
	hist_2d_jpmdp_vs_sum3jetpt_dijet.Draw("colz")
	c1.cd(3)
	hist_2d_jpmdp_vs_sum3jetpt_trijet.Draw("colz")
	c1.SaveAs(self.extraDir+"multiplot_jpmdpVssum3jetpt.png")

	# sum3jets and MT
	c1.Clear()
	c1.Divide(3,2)
	c1.cd(1)
	leg = rt.TLegend(0.2,0.2,0.8,0.8)
	dijetLeg = rt.TH1F("dijetLeg","Dijet Events",2,-100,-99)
	trijetLeg = rt.TH1F("trijetLeg","Trijet Events",2,-100,-99)
	dijetLeg.SetLineColor(rt.kBlue)
	trijetLeg.SetLineColor(rt.kRed)
	leg.AddEntry(dijetLeg)
	leg.AddEntry(trijetLeg)
	leg.Draw()
	c1.cd(4)
	hist_dijet_xProj = hist_2d_sum3jetpt_vs_MT_dijet.ProjectionX()
	hist_dijet_xProj.SetLineColor(rt.kBlue)
	hist_trijet_xProj = hist_2d_sum3jetpt_vs_MT_trijet.ProjectionX()
	hist_trijet_xProj.SetLineColor(rt.kRed)
	hist_trijet_xProj.Draw()
	hist_dijet_xProj.Draw("same")
	c1.cd(6)
	hist_dijet_yProj = hist_2d_sum3jetpt_vs_MT_dijet.ProjectionY()
	hist_dijet_yProj.SetLineColor(rt.kBlue)
	hist_trijet_yProj = hist_2d_sum3jetpt_vs_MT_trijet.ProjectionY()
	hist_trijet_yProj.SetLineColor(rt.kRed)
	hist_trijet_yProj.Draw()
	hist_dijet_yProj.Draw("same")
	#c1.cd(5)
	#hist_ptRatio_dijet.SetLineColor(rt.kBlue)
	#hist_ptRatio_trijet.SetLineColor(rt.kRed)
	#hist_ptRatio_trijet.Draw()
	#hist_ptRatio_dijet.Draw("same")
	c1.cd(2)
	hist_2d_sum3jetpt_vs_MT_dijet.Draw("colz")
	c1.cd(3)
	hist_2d_sum3jetpt_vs_MT_trijet.Draw("colz")
	c1.SaveAs(self.extraDir+"multiplot_sum3jetptvsMT.png")

	# jmdp and MT
	c1.Clear()
	c1.Divide(3,2)
	c1.cd(1)
	leg = rt.TLegend(0.2,0.2,0.8,0.8)
	dijetLeg = rt.TH1F("dijetLeg","Dijet Events",2,-100,-99)
	trijetLeg = rt.TH1F("trijetLeg","Trijet Events",2,-100,-99)
	dijetLeg.SetLineColor(rt.kBlue)
	trijetLeg.SetLineColor(rt.kRed)
	leg.AddEntry(dijetLeg)
	leg.AddEntry(trijetLeg)
	leg.Draw()
	c1.cd(4)
	hist_dijet_xProj = hist_2d_jpmdp_vs_MT_dijet.ProjectionX()
	hist_dijet_xProj.SetLineColor(rt.kBlue)
	hist_trijet_xProj = hist_2d_jpmdp_vs_MT_trijet.ProjectionX()
	hist_trijet_xProj.SetLineColor(rt.kRed)
	hist_trijet_xProj.Draw()
	hist_dijet_xProj.Draw("same")
	c1.cd(6)
	hist_dijet_yProj = hist_2d_jpmdp_vs_MT_dijet.ProjectionY()
	hist_dijet_yProj.SetLineColor(rt.kBlue)
	hist_trijet_yProj = hist_2d_jpmdp_vs_MT_trijet.ProjectionY()
	hist_trijet_yProj.SetLineColor(rt.kRed)
	hist_trijet_yProj.Draw()
	hist_dijet_yProj.Draw("same")
	#c1.cd(5)
	#hist_ptRatio_dijet.SetLineColor(rt.kBlue)
	#hist_ptRatio_trijet.SetLineColor(rt.kRed)
	#hist_ptRatio_trijet.Draw()
	#hist_ptRatio_dijet.Draw("same")
	c1.cd(2)
	hist_2d_jpmdp_vs_MT_dijet.Draw("colz")
	c1.cd(3)
	hist_2d_jpmdp_vs_MT_trijet.Draw("colz")
	c1.SaveAs(self.extraDir+"multiplot_jpmdpvsMT.png")
	
	c1.Clear()
	c1.Divide(3,3)
	c1.cd(1)
	hist_2d_comboB_jet1mdphi.Draw("colz")
	c1.cd(2)
	hist_2d_comboB_dijet_jet1mdphi.Draw("colz")
	c1.cd(3)
	hist_2d_comboB_trijet_jet1mdphi.Draw("colz")
	c1.cd(4)
	hist_2d_comboB_jet2mdphi.Draw("colz")
	c1.cd(5)
	hist_2d_comboB_dijet_jet2mdphi.Draw("colz")
	c1.cd(6)
	hist_2d_comboB_trijet_jet2mdphi.Draw("colz")
	c1.cd(7)
	hist_2d_comboB_jet3mdphi.Draw("colz")
	c1.cd(8)
	hist_2d_comboB_dijet_jet3mdphi.Draw("colz")
	c1.cd(9)
	hist_2d_comboB_trijet_jet3mdphi.Draw("colz")
	c1.SaveAs(self.extraDir+"multiplot_SplitonJetIndexMaxDPhi_comboB.png")
	c1.Clear()
	c1.Divide(3,3)
	c1.cd(1)
	hist_2d_comboA_jet1mdphi.Draw("colz")
	c1.cd(2)
	hist_2d_comboA_dijet_jet1mdphi.Draw("colz")
	c1.cd(3)
	hist_2d_comboA_trijet_jet1mdphi.Draw("colz")
	c1.cd(4)
	hist_2d_comboA_jet2mdphi.Draw("colz")
	c1.cd(5)
	hist_2d_comboA_dijet_jet2mdphi.Draw("colz")
	c1.cd(6)
	hist_2d_comboA_trijet_jet2mdphi.Draw("colz")
	c1.cd(7)
	hist_2d_comboA_jet3mdphi.Draw("colz")
	c1.cd(8)
	hist_2d_comboA_dijet_jet3mdphi.Draw("colz")
	c1.cd(9)
	hist_2d_comboA_trijet_jet3mdphi.Draw("colz")
	c1.SaveAs(self.extraDir+"multiplot_SplitonJetIndexMaxDPhi_comboA.png")



def addLoop():
	baseClass.loop = loop

def makePlot(ListOfHistos, location, name, stackTitle, logY = False, gs = False):
	c1 = rt.TCanvas("c1","c1",900,600)
	c1.SetLineWidth(2)
	if gs:
		c1.SetGrayscale()
	if logY:
		c1.SetLogy()
	stack = rt.THStack("stack",stackTitle)
	stackTitle = stackTitle.split(";")
	#stack.SetMinimum(10) # use these two lines to set Y axis range
	#stack.SetMaximum(1100)
	i = 0
	for histo in ListOfHistos:
		histo.SetLineColor(i+1)
		histo.SetLineWidth(2)
		stack.Add(histo)
		i += 1
	if len(ListOfHistos) > 1:
		stack.Draw("nostack")
		c1.BuildLegend()
	else:
		stack.Draw("hist TEXT0")
	stack.SetTitle(stackTitle[0])
	stack.GetXaxis().SetTitle(stackTitle[1])
	stack.GetYaxis().SetTitle(stackTitle[2])
	rt.gPad.Modified()
	#stack.GetXaxis().SetLimits(0,4000) # use these two lines to set X axis range
	#stack.Draw("nostack")
	#t1 = rt.TText(0,80,"CMS Simulation") # use axis numbers for this
	#t1.Draw("same")
	c1.SaveAs(location+name+".png")

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))
