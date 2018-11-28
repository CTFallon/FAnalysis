from analysisBase import baseClass
import ROOT as rt
from array import array

rt.gROOT.SetBatch(True)
rt.gROOT.ForceStyle()
rt.gStyle.SetOptStat(0)

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

	histList_sdvar13 = []
	histList_jetPtMaxDPhi = []

	if self.extraDir[-3] == "p":
		zMass = 3000
	else:
		zMass = int(self.extraDir[-5])*1000
	
	nMTBins = 250
	hist_MT_FSR = self.makeTH1F("hist_MT_FSR","FSR Truth;MT;",nMTBins,0,2*zMass)
	hist_MT_Dijet = self.makeTH1F("hist_MT_Dijet","Dijet;MT;",nMTBins,0,2*zMass)
	hist_MT_Trijet = self.makeTH1F("hist_MT_Trijet","Trijet;MT;",nMTBins,0,2*zMass)
	hist_chi2_SDvar13 = self.makeTH1F("hist_chi2_SDVar13","chi2;SDvar13;Chi2",100,0,0.5)
	hist_chi2_jetPtMaxDPhi = self.makeTH1F("hist_chi2_jetPtMaxDPhi","chi2;jetPtMaxDPhi;Chi2",100,0,2000)
	hist_KS_SDvar13 = self.makeTH1F("hist_KS_SDVar13","KS;SDvar13;KS",100,0,0.5)
	hist_KS_jetPtMaxDPhi = self.makeTH1F("hist_KS_jetPtMaxDPhi","KS;jetPtMaxDPhi;KS",100,0,2000)
	hist_Reso_SDvar13 = self.makeTH1F("hist_Reso_SDVar13","Reso;SDvar13;Reso",99,0,0.5)
	hist_Reso_jetPtMaxDPhi = self.makeTH1F("hist_Reso_jetPtMaxDPhi","Reso;jetPtMaxDPhi;Reso",99,0,2000)
	hist_DKL_SDvar13 = self.makeTH1F("hist_DKL_SDvar13","DKL;SDvarCut;DKL",100,0,0.5)
	hist_DKL_jetPtMaxDPhi = self.makeTH1F("hist_DKL_jetPtMaxDPhi","DKL",100,0,2000)
	hist_pTP_SDvar13 = self.makeTH1F("hist_pTP_SDvar13","Percent True Positive;SDvar13;",100,0,0.5)
	hist_pTN_SDvar13 = self.makeTH1F("hist_pTN_SDvar13","Percent True Negative;SDvar13;",100,0,0.5)
	hist_pTP_jetPtMaxDPhi = self.makeTH1F("hist_pTP_jetPtMaxDPhi","Percent True Positive;jetPtMaxDPhi;",100,0,2000)
	hist_pTN_jetPtMaxDPhi = self.makeTH1F("hist_pTN_jetPtMaxDPhi","Percent True Negative;jetPtMaxDPhi;",100,0,2000)

	hist_SDvar13 = self.makeTH1F("hist_SDvar13","SDvar13;Sdvar13;",100,0,1.0)


	nEvents2Jets = 0
	nEvents_SDvar13_Dijet = [0.]*100
	nEvents_SDvar13_Trijet = [0.]*100
	nEvents_JetPt_Dijet = [0.]*100
	nEvents_JetPt_Trijet = [0.]*100

	nTP_SDVar13 = [0.]*100 # true positive, when we should use 123 and we do
	nFP_SDVar13 = [0.]*100 # false positive, when we should use 12 but we use 123
	nTN_SDVar13 = [0.]*100 # true negative, when we should use 12 and we do
	nFN_SDVar13 = [0.]*100 # false negative, when we should use 123 but we use 12

	nTP_JPMDP = [0.]*100 # true positive, when we should use 123 and we do
	nFP_JPMDP = [0.]*100 # false positive, when we should use 12 but we use 123
	nTN_JPMDP = [0.]*100 # true negative, when we should use 12 and we do
	nFN_JPMDP = [0.]*100 # false negative, when we should use 123 but we use 12


	nEventsFSR12 = 0
	nEventsFSR123 = 0
	nEventsFSROther = 0

	hist2d_SDVar13_vs_JetPtMaxDPhi = self.makeTH2F(
		"hist2d_SDVar13_vs_JetPtMaxDPhi",
		"SDVar13 vs JetPtMaxDPhi",
		100,0,0.5,100,0,2000)
	hist2d_deltaMT_vs_JetPtMaxDPhi = self.makeTH2F(
		"hist2d_deltaMT_vs_JetPtMaxDPhi",
		"deltaMT vs JetPtMaxDPhi",
		100,0,2000,100,0,2000)

	#hist2d_stdDev_vs_deltaMu_jetPt = self.makeTH2F(
	#"hist2d_stdDev_vs_deltaMu_jetPt",
	#"MT StdDev vs dMu(FSR, jetPt);stdDev;deltaMu",
	#100,400,700,
	#100,-200,200
	#)
	#hist2d_stdDev_vs_deltaMu_sdvar13 = self.makeTH2F(
	#"hist2d_stdDev_vs_deltaMu_sdvar13",
	#"MT StdDev vs dMu(FSR, sdvar13);stdDev;deltaMu",
	#100,400,700,
	#100,-200,200
	#)
	
	hist_maxDPhiIndex = self.makeTH1F("hist_maxDPhiIndex","maxDPhiIndex;index;count",
		6,0,6)
	hist_maxDPhiIndex_2jets = self.makeTH1F("hist_maxDPhiIndex_2jets","maxDPhiIndex_2jets;index;count",
		6,0,6)
	for iBin in range(100):
		binEdge = hist_chi2_SDvar13.GetXaxis().GetBinUpEdge(iBin)
		histList_sdvar13.append(rt.TH1F("hist_sdVar13_cut_"+str(binEdge),"sdVar("+str(binEdge)+");MT;",nMTBins,0,2*zMass))
		binEdge = hist_chi2_jetPtMaxDPhi.GetXaxis().GetBinUpEdge(iBin)
		histList_jetPtMaxDPhi.append(rt.TH1F("hist_jetPtMaxDPhi_cut_"+str(binEdge),"jpmdp("+str(binEdge)+");MT;",nMTBins,0,2*zMass))
			
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		#Alright, so FSR is 'truth'
		# SDVar13 is based on the softdrop variable between jets 1 and 3 (pt(1)/(pt(1)+pt(3)))
		# jetPtMaxDPhi is based on the pt of the jet which is furthest from the MET
		# jet3DPhi is based on whether or not jet 3 is closest to MET

		# task1: make FOM curves for MT Resolution and Chi2Test w.r.t. FSR MT

		nJets = len(tree.JetsAK8)
		MT12 = trans_mass_Njet(tree.JetsAK8[0:2], tree.MET, tree.METPhi)

		FSRJets = []
		FSRIndex = []
		for iJet in range(nJets):
			if tree.JetsAK8_isHV[iJet]:
				FSRJets.append(tree.JetsAK8[iJet])
				FSRIndex.append(iJet)
		MTFSR = trans_mass_Njet(FSRJets, tree.MET, tree.METPhi)
		#if nJets > 2:
		hist_MT_FSR.Fill(MTFSR)
		hist_MT_Dijet.Fill(MT12)

		if nJets > 2:
			if (0 in FSRIndex) and (1 in FSRIndex) and not (2 in FSRIndex):
				nEventsFSR12 += 1
			elif (0 in FSRIndex) and (1 in FSRIndex) and (2 in FSRIndex):
				nEventsFSR123 += 1
			else:
				nEventsFSROther += 1

		if nJets == 2:
			nEvents2Jets += 1
			for histo in histList_sdvar13:
				histo.Fill(MT12)
			for histo in histList_jetPtMaxDPhi:
				histo.Fill(MT12)
			hist_MT_Trijet.Fill(MT12)
			hist_maxDPhiIndex_2jets.Fill([tree.DeltaPhi1, tree.DeltaPhi2].index(max(tree.DeltaPhi1, tree.DeltaPhi2)))
		if nJets > 2:
			MT123 = trans_mass_Njet(tree.JetsAK8[0:3], tree.MET, tree.METPhi)
			hist_MT_Trijet.Fill(MT123)
			sdvar13 = tree.JetsAK8[2].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt())
			mindPhiIndex = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3].index(min(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3))
			jetPtMaxDPhi = tree.JetsAK8[[tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3].index(max(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3))].Pt()
			hist2d_SDVar13_vs_JetPtMaxDPhi.Fill(sdvar13,jetPtMaxDPhi)
			hist_maxDPhiIndex.Fill([tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3].index(max(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3)))
			for iHisto in range(len(histList_sdvar13)):
				cutValue = hist_chi2_SDvar13.GetXaxis().GetBinUpEdge(iHisto)
				if sdvar13 < cutValue: # Low SDvar13 means that pt3 is low, optimized (read: mean when fit to gaus) is .321
					histList_sdvar13[iHisto].Fill(MT12)
					nEvents_SDvar13_Dijet[iHisto] += 1
					if (0 in FSRIndex) and (1 in FSRIndex) and not (2 in FSRIndex):
						nTN_SDVar13[iHisto] += 1
					elif (0 in FSRIndex) and (1 in FSRIndex) and (2 in FSRIndex):
						nFN_SDVar13[iHisto] += 1
				else:
					histList_sdvar13[iHisto].Fill(MT123)
					nEvents_SDvar13_Trijet[iHisto] += 1
					if (0 in FSRIndex) and (1 in FSRIndex) and not (2 in FSRIndex):
						nFP_SDVar13[iHisto] += 1
					elif (0 in FSRIndex) and (1 in FSRIndex) and (2 in FSRIndex):
						nTP_SDVar13[iHisto] += 1
			for iHisto in range(len(histList_jetPtMaxDPhi)):
				cutValue = hist_chi2_jetPtMaxDPhi.GetXaxis().GetBinUpEdge(iHisto)
				if jetPtMaxDPhi > cutValue: # High Pt of the jet furthest from MET means there is a higher probability of the jets closer to MET to have been a single SVJ that split
							# optimized (read: mean when fit to gaus) is 861
					histList_jetPtMaxDPhi[iHisto].Fill(MT12)
					nEvents_JetPt_Dijet[iHisto] += 1
					if (0 in FSRIndex) and (1 in FSRIndex) and not (2 in FSRIndex):
						nTN_JPMDP[iHisto] += 1
					elif (0 in FSRIndex) and (1 in FSRIndex) and (2 in FSRIndex):
						nFN_JPMDP[iHisto] += 1
				else:
					histList_jetPtMaxDPhi[iHisto].Fill(MT123)
					nEvents_JetPt_Trijet[iHisto] += 1
					if (0 in FSRIndex) and (1 in FSRIndex) and not (2 in FSRIndex):
						nFP_JPMDP[iHisto] += 1
					elif (0 in FSRIndex) and (1 in FSRIndex) and (2 in FSRIndex):
						nTP_JPMDP[iHisto] += 1
			hist2d_deltaMT_vs_JetPtMaxDPhi.Fill(MT123-MT12, jetPtMaxDPhi)

	for iHisto in range(len(histList_sdvar13)):
		hist_chi2_SDvar13.SetBinContent(iHisto,hist_MT_FSR.Chi2Test(histList_sdvar13[iHisto]))
		hist_KS_SDvar13.SetBinContent(iHisto,hist_MT_FSR.KolmogorovTest(histList_sdvar13[iHisto]))		
		hist_Reso_SDvar13.SetBinContent(iHisto,histList_sdvar13[iHisto].GetRMS()/histList_sdvar13[iHisto].GetMean())
		hist_DKL_SDvar13.SetBinContent(iHisto,Dkl(hist_MT_FSR, histList_sdvar13[iHisto]))
		#hist2d_stdDev_vs_deltaMu_sdvar13.Fill(histList_sdvar13[iHisto].GetRMS(), (histList_sdvar13[iHisto].GetMean()-hist_MT_FSR.GetMean()))
		if (nTP_SDVar13[iHisto]+nFP_SDVar13[iHisto]) > 0:
			print(nTP_SDVar13[iHisto]/(nTP_SDVar13[iHisto]+nFP_SDVar13[iHisto]))
			hist_pTP_SDvar13.SetBinContent(iHisto,nTP_SDVar13[iHisto]/(nTP_SDVar13[iHisto]+nFP_SDVar13[iHisto]))
		if (nTN_SDVar13[iHisto]+nFN_SDVar13[iHisto]) > 0:
			hist_pTN_SDvar13.SetBinContent(iHisto,nTN_SDVar13[iHisto]/(nTN_SDVar13[iHisto]+nFN_SDVar13[iHisto]))
	for iHisto in range(len(histList_jetPtMaxDPhi)):
		hist_chi2_jetPtMaxDPhi.SetBinContent(iHisto,hist_MT_FSR.Chi2Test(histList_jetPtMaxDPhi[iHisto]))
		hist_KS_jetPtMaxDPhi.SetBinContent(iHisto,hist_MT_FSR.KolmogorovTest(histList_jetPtMaxDPhi[iHisto]))
		hist_Reso_jetPtMaxDPhi.SetBinContent(iHisto,histList_jetPtMaxDPhi[iHisto].GetRMS()/histList_jetPtMaxDPhi[iHisto].GetMean())
		hist_DKL_jetPtMaxDPhi.SetBinContent(iHisto, Dkl(hist_MT_FSR, histList_jetPtMaxDPhi[iHisto]))
		#hist2d_stdDev_vs_deltaMu_jetPt.Fill(histList_jetPtMaxDPhi[iHisto].GetRMS(), (histList_jetPtMaxDPhi[iHisto].GetMean()-hist_MT_FSR.GetMean()))
		if (nTP_JPMDP[iHisto]+nFP_JPMDP[iHisto]) > 0:
			hist_pTP_jetPtMaxDPhi.SetBinContent(iHisto,nTP_JPMDP[iHisto]/(nTP_JPMDP[iHisto]+nFP_JPMDP[iHisto]))
		if (nTN_JPMDP[iHisto]+nFN_JPMDP[iHisto]) > 0:
			hist_pTN_jetPtMaxDPhi.SetBinContent(iHisto,nTN_JPMDP[iHisto]/(nTN_JPMDP[iHisto]+nFN_JPMDP[iHisto]))




	sdVar13_maxChi2Bin = hist_chi2_SDvar13.GetMaximumBin()
	jetPtMaxDPhi_maxChi2Bin = hist_chi2_jetPtMaxDPhi.GetMaximumBin()

	sdVar13_maxKSBin = hist_KS_SDvar13.GetMaximumBin()
	jetPtMaxDPhi_maxKSBin = hist_KS_jetPtMaxDPhi.GetMaximumBin()

	sdVar13_minResoBin = hist_Reso_SDvar13.GetMinimumBin()
	jetPtMaxDPhi_minResoBin = hist_Reso_jetPtMaxDPhi.GetMinimumBin()
	
	print("Best Chi2 SDvar13 Bin = " + str(sdVar13_maxChi2Bin))
	print("Best Chi2 JetPtMaxDPhi Bin = " + str(jetPtMaxDPhi_maxChi2Bin))
	print("Chi2 SDvar13 = " +str(hist_MT_FSR.Chi2Test(histList_sdvar13[sdVar13_maxChi2Bin])))
	print("Chi2 JetPtMaxDPhi = " +str(hist_MT_FSR.Chi2Test(histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin])))

	print("Best KS SDvar13 Bin = " + str(sdVar13_maxKSBin))
	print("Best Ks JetPtMaxDPhi Bin = " + str(jetPtMaxDPhi_maxKSBin))
	print("KS SDvar13 = " +str(hist_MT_FSR.KolmogorovTest(histList_sdvar13[sdVar13_maxKSBin])))
	print("KS JetPtMaxDPhi = " +str(hist_MT_FSR.KolmogorovTest(histList_jetPtMaxDPhi[jetPtMaxDPhi_maxKSBin])))

	print("Best Reso SDvar13 Bin = " + str(sdVar13_minResoBin))
	print("Best Reso JetPtMaxDPhi Bin = " + str(jetPtMaxDPhi_minResoBin))
	
	# Jet Code & KS & Chi2 & Reso & KS & Chi2 & Reso
	print("jet Pt & {} & {} & {} & {} & {} & {} & {}\\".format(
		jetPtMaxDPhi_maxKSBin*20,
		jetPtMaxDPhi_maxChi2Bin*20,
		jetPtMaxDPhi_minResoBin*20,
		hist_KS_jetPtMaxDPhi.GetBinContent(jetPtMaxDPhi_maxKSBin),
		hist_chi2_jetPtMaxDPhi.GetBinContent(jetPtMaxDPhi_maxChi2Bin),
		hist_Reso_jetPtMaxDPhi.GetBinContent(jetPtMaxDPhi_minResoBin),
		hist_MT_FSR.GetRMS()/hist_MT_FSR.GetMean()
		))
	print("SDVar13 & {} & {} & {} & {} & {} & {} & {}\\".format(
		sdVar13_maxKSBin*0.005,
		sdVar13_maxChi2Bin*0.005,
		sdVar13_minResoBin*0.005,
		hist_KS_SDvar13.GetBinContent(sdVar13_maxKSBin),
		hist_chi2_SDvar13.GetBinContent(sdVar13_maxChi2Bin),
		hist_Reso_SDvar13.GetBinContent(sdVar13_minResoBin),
		hist_MT_FSR.GetRMS()/hist_MT_FSR.GetMean()
		))


	#print("iBin Cuts TP FP TN FN")
	#for i in range(100):
	#	print("{} SDVar13 {} {} {} {}".format(i,nTP_SDVar13[i],nFP_SDVar13[i],nTN_SDVar13[i],nFN_SDVar13[i]))
	#	print("{} JPMDP {} {} {} {}".format(i,  nTP_JPMDP[i]  ,nFP_JPMDP[i]  ,nTN_JPMDP[i]  ,nFN_JPMDP[i]))

	c1 = rt.TCanvas("c1","c1",900,600)

	hist_DKL_jetPtMaxDPhi.Draw()
	c1.SaveAs(self.extraDir+"DKL_pt.png")
	hist_DKL_SDvar13.Draw()
	c1.SaveAs(self.extraDir+"DKL_sdvar.png")


	hist_KS_SDvar13.SetLineColor(2)
	tempMax = hist_KS_SDvar13.GetMaximum()
	ksSDMaxBin = hist_KS_SDvar13.GetMaximumBin()
	print("Maximum KS SDvar is {}".format(tempMax))
	if tempMax != 0:
		for i in range(hist_KS_SDvar13.GetNbinsX()):
			tempBin = hist_KS_SDvar13.GetBinContent(i)
			hist_KS_SDvar13.SetBinContent(i,tempBin/tempMax)
	hist_KS_SDvar13.Draw()
	tempMax = hist_chi2_SDvar13.GetMaximum()
	chi2SDMaxBin = hist_chi2_SDvar13.GetMaximumBin()
	print("Maximum Chi2 SDvar is {}".format(tempMax))
	if tempMax != 0:
		for i in range(hist_chi2_SDvar13.GetNbinsX()):
			tempBin = hist_chi2_SDvar13.GetBinContent(i)
			hist_chi2_SDvar13.SetBinContent(i,tempBin/tempMax)
	hist_chi2_SDvar13.Draw("same")
	hist_Reso_SDvar13.SetLineColor(3)
	tempMax = hist_Reso_SDvar13.GetMaximum()
	resoSDMaxBin = hist_Reso_SDvar13.GetMaximumBin()
	print("Maximum Reso SDvar is {}".format(tempMax))
	print("minimum Reso SDvar is {}".format(hist_Reso_SDvar13.GetMinimum()))
	if tempMax != 0:
		for i in range(hist_Reso_SDvar13.GetNbinsX()):
			tempBin = hist_Reso_SDvar13.GetBinContent(i)
			hist_Reso_SDvar13.SetBinContent(i,tempBin/tempMax)
	hist_Reso_SDvar13.Draw("same")

	hist_DKL_SDvar13.SetLineColor(4)
	tempMax = hist_DKL_SDvar13.GetMaximum()
	dklSDMaxBin = hist_DKL_SDvar13.GetMaximumBin()
	print("Maximum Reso SDvar is {}".format(tempMax))
	print("minimum Reso SDvar is {}".format(hist_DKL_SDvar13.GetMinimum()))
	if tempMax != 0:
		for i in range(hist_DKL_SDvar13.GetNbinsX()):
			tempBin = hist_DKL_SDvar13.GetBinContent(i)
			hist_DKL_SDvar13.SetBinContent(i,tempBin/tempMax)
	hist_DKL_SDvar13.Draw("same")


	hist_pTP_SDvar13.SetLineColor(5)
	hist_pTP_SDvar13.Draw("same")
	#hist_pTN_SDvar13.SetLineColor(5)
	#hist_pTN_SDvar13.Draw("same")
	c1.BuildLegend(0.1,0.1,0.4,0.4)
	c1.SaveAs(self.extraDir+"FOM_sdvar13.png")


	hist_KS_jetPtMaxDPhi.SetLineColor(2)
	tempMax = hist_KS_jetPtMaxDPhi.GetMaximum()
	ksPtMaxBin = hist_KS_jetPtMaxDPhi.GetMaximumBin()
	print("Maximum KS JetPtMAxDPhi is {}".format(tempMax))
	if tempMax != 0:
		for i in range(hist_KS_jetPtMaxDPhi.GetNbinsX()):
			tempBin = hist_KS_jetPtMaxDPhi.GetBinContent(i)
			hist_KS_jetPtMaxDPhi.SetBinContent(i,tempBin/tempMax)
	hist_KS_jetPtMaxDPhi.Draw()
	tempMax = hist_chi2_jetPtMaxDPhi.GetMaximum()
	chi2PtMaxBin = hist_chi2_jetPtMaxDPhi.GetMaximumBin()
	print("Maximum Chi2 JetPtMAxDPhi is {}".format(tempMax))
	if tempMax != 0:
		for i in range(hist_chi2_jetPtMaxDPhi.GetNbinsX()):
			tempBin = hist_chi2_jetPtMaxDPhi.GetBinContent(i)
			hist_chi2_jetPtMaxDPhi.SetBinContent(i,tempBin/tempMax)
	hist_chi2_jetPtMaxDPhi.Draw("same")
	hist_Reso_jetPtMaxDPhi.SetLineColor(3)
	tempMax = hist_Reso_jetPtMaxDPhi.GetMaximum()
	resoPtMaxBin = hist_Reso_jetPtMaxDPhi.GetMinimumBin()
	print("Maximum Reso JetPtMAxDPhi is {}".format(tempMax))
	print("Minimum Reso JetPtMAxDPhi is {}".format(hist_Reso_jetPtMaxDPhi.GetMinimum()))
	if tempMax != 0:
		for i in range(hist_Reso_jetPtMaxDPhi.GetNbinsX()):
			tempBin = hist_Reso_jetPtMaxDPhi.GetBinContent(i)
			hist_Reso_jetPtMaxDPhi.SetBinContent(i,tempBin/tempMax)
	hist_Reso_jetPtMaxDPhi.Draw("same")

	hist_DKL_jetPtMaxDPhi.SetLineColor(4)
	tempMax = hist_DKL_jetPtMaxDPhi.GetMaximum()
	dkljetPtMaxDPhiMaxBin = hist_DKL_jetPtMaxDPhi.GetMaximumBin()
	print("Maximum Reso jetPtMaxDPhi is {}".format(tempMax))
	print("minimum Reso jetPtMaxDPhi is {}".format(hist_DKL_jetPtMaxDPhi.GetMinimum()))
	if tempMax != 0:
		for i in range(hist_DKL_jetPtMaxDPhi.GetNbinsX()):
			tempBin = hist_DKL_jetPtMaxDPhi.GetBinContent(i)
			hist_DKL_jetPtMaxDPhi.SetBinContent(i,tempBin/tempMax)
	hist_DKL_jetPtMaxDPhi.Draw("same")

	hist_pTP_jetPtMaxDPhi.SetLineColor(5)
	hist_pTP_jetPtMaxDPhi.Draw("same")
	#hist_pTN_jetPtMaxDPhi.SetLineColor(5)
	#hist_pTN_jetPtMaxDPhi.Draw("same")
	c1.BuildLegend(0.1,0.1,0.4,0.4)
	c1.SaveAs(self.extraDir+"FOM_jetPtMaxDPhi.png")

	#c1.SetLogy()
	histList_jetPtMaxDPhi[ksPtMaxBin].SetLineColor(2)
	histList_jetPtMaxDPhi[ksPtMaxBin].Draw()
	histList_jetPtMaxDPhi[chi2PtMaxBin].SetLineColor(3)
	histList_jetPtMaxDPhi[chi2PtMaxBin].Draw("same")
	histList_jetPtMaxDPhi[resoPtMaxBin].SetLineColor(4)
	histList_jetPtMaxDPhi[resoPtMaxBin].Draw("same")
	hist_MT_Dijet.SetLineColor(1)
	hist_MT_Dijet.SetLineStyle(2)
	hist_MT_Dijet.Draw('same')
	hist_MT_Trijet.SetLineColor(1)
	hist_MT_Trijet.SetLineStyle(3)
	hist_MT_Trijet.Draw('same')
	hist_MT_FSR.SetLineColor(1)
	hist_MT_FSR.SetLineStyle(4)
	hist_MT_FSR.Draw("same")
	c1.BuildLegend(0.0,0.7,0.3,0.9)
	c1.SaveAs(self.extraDir+"MTDistro_jetPtOptimal.png")

	histList_sdvar13[ksSDMaxBin].SetLineColor(2)
	histList_sdvar13[ksSDMaxBin].Draw()
	histList_sdvar13[chi2SDMaxBin].SetLineColor(3)
	histList_sdvar13[chi2SDMaxBin].Draw("same")
	histList_sdvar13[resoSDMaxBin].SetLineColor(4)
	histList_sdvar13[resoSDMaxBin].Draw("same")
	hist_MT_Dijet.SetLineColor(1)
	hist_MT_Dijet.SetLineStyle(2)
	hist_MT_Dijet.Draw('same')
	hist_MT_Trijet.SetLineColor(1)
	hist_MT_Trijet.SetLineStyle(3)
	hist_MT_Trijet.Draw('same')
	hist_MT_FSR.SetLineColor(1)
	hist_MT_FSR.SetLineStyle(4)
	hist_MT_FSR.Draw("same")
	c1.BuildLegend(0.0,0.7,0.3,0.9)
	c1.SaveAs(self.extraDir+"MTDistro_sdVar13Optimal.png")

	hist2d_SDVar13_vs_JetPtMaxDPhi.Draw("colz")
	c1.SaveAs(self.extraDir+"2d_sdvar13_vs_jetptmaxdphi.png")
	hist2d_deltaMT_vs_JetPtMaxDPhi.Draw("colz")
	c1.SaveAs(self.extraDir+"2d_deltaMT_vs_jetptmaxdphi.png")


	hist_maxDPhiIndex.Draw()
	hist_maxDPhiIndex_2jets.SetLineColor(2)
	hist_maxDPhiIndex_2jets.Draw("same")
	c1.BuildLegend()
	c1.SaveAs(self.extraDir+"maxDPhiIndex.png")

	#hist2d_stdDev_vs_deltaMu_sdvar13.Draw("colz")
	#c1.SaveAs(self.extraDir+"2d_stddevdMu_sdvar.png")
	#hist2d_stdDev_vs_deltaMu_jetPt.Draw("colz")
	#c1.SaveAs(self.extraDir+"2d_stddevdMu_jetPt.png")
	"""
	#hist_MT_FSR.GetXaxis().SetLabelSize(3)
	hist_MT_FSR.Draw()
	hist_MT_Dijet.SetLineColor(6)
	#hist_MT_Dijet.Draw("same")
	hist_MT_Trijet.SetLineColor(7)
	#hist_MT_Trijet.Draw("same")
	histList_sdvar13[sdVar13_maxChi2Bin].SetLineColor(3)
	#histList_sdvar13[sdVar13_maxChi2Bin].SetLineWidth(4)
	histList_sdvar13[sdVar13_maxChi2Bin].Draw("same")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].SetLineColor(4)
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].Draw("same")
	c1.BuildLegend(0.7,0.6,1.0,0.9)
	hist_MT_FSR.SetTitle("Optimized MT Distribtions")
	c1.SaveAs(self.extraDir+"VariablePlots.png")

	
	hist_MT_FSR.Draw()
	hist_MT_Dijet.Draw("same")
	hist_MT_Trijet.Draw("same")
	c1.BuildLegend(0.7,0.6,1.0,0.9)
	hist_MT_FSR.SetTitle("Simple MT Distributions")
	c1.SaveAs(self.extraDir+"DefaultPlots.png")

	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].SetLineColor(1)
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxKSBin].SetLineColor(2)
	histList_sdvar13[jetPtMaxDPhi_maxChi2Bin].SetLineColor(3)
	histList_sdvar13[jetPtMaxDPhi_maxKSBin].SetLineColor(4)
	histList_jetPtMaxDPhi[jetPtMaxDPhi_minResoBin].SetLineColor(5)
	histList_sdvar13[jetPtMaxDPhi_minResoBin].SetLineColor(6)

	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].SetTitle("JetPtMaxDPhi, Chi2")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxKSBin].SetTitle("JetPtMaxDPhi, KS")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_minResoBin].SetTitle("JetPtMaxDPhi, Reso")
	histList_sdvar13[jetPtMaxDPhi_maxChi2Bin].SetTitle("SDVar13, Chi2")
	histList_sdvar13[jetPtMaxDPhi_maxKSBin].SetTitle("SDVar13, KS")
	histList_sdvar13[jetPtMaxDPhi_minResoBin].SetTitle("SDVar13, Reso")

	hist_MT_FSR.SetLineStyle(2)
	hist_MT_FSR.Draw()
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].Draw("same")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxKSBin].Draw("same")
	histList_sdvar13[jetPtMaxDPhi_maxChi2Bin].Draw("same")
	histList_sdvar13[jetPtMaxDPhi_maxKSBin].Draw("same")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_minResoBin].Draw("same")
	histList_sdvar13[jetPtMaxDPhi_minResoBin].Draw("same")

	c1.BuildLegend(0.7,0.6,1.0,0.9)
	c1.SaveAs(self.extraDir+"BestChi2andKSPlots.png")
	
	hist_MT_FSR.Draw()
	hist_MT_Dijet.Draw("same")
	hist_MT_Trijet.Draw("same")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].Draw("same")
	c1.BuildLegend(0.7,0.6,1.0,0.9)
	hist_MT_FSR.SetTitle("Simple and Best MT Distributions")
	c1.SaveAs(self.extraDir+"DefaultAndBestPlots.png")


	c1.SetLogy()
	hist_MT_FSR.Draw()
	histList_sdvar13[sdVar13_maxChi2Bin].Draw("same")
	histList_jetPtMaxDPhi[jetPtMaxDPhi_maxChi2Bin].Draw("same")
	c1.BuildLegend(0.7,0.6,1.0,0.9)
	hist_MT_FSR.SetTitle("Optimized MT Distributions, Log Scale")
	c1.SaveAs(self.extraDir+"VariablePlotsLOG.png")

	print("")
	print("nEvents2Jets " + str(nEvents2Jets))
	print("nEvents_SDvar13_Dijet " + str(nEvents_SDvar13_Dijet))
	print("nEvents_SDvar13_Trijet " + str(nEvents_SDvar13_Trijet))
	print("nEvents_JetPt_Dijet " + str(nEvents_JetPt_Dijet))
	print("nEvents_JetPt_Trijet " + str(nEvents_JetPt_Trijet))
	print("nEventsFSR12 " + str(nEventsFSR12))
	print("nEventsFSR123 " + str(nEventsFSR123))
	print("nEventsFSROther " + str(nEventsFSROther))
	
	print("Table of KS, Chi2, and MT Resolution")

	for iBin in range(100):
		binEdge = hist_chi2_SDvar13.GetXaxis().GetBinUpEdge(iBin)
		KS = hist_KS_SDvar13.GetBinContent(iBin)
		chi2 = hist_chi2_SDvar13.GetBinContent(iBin)
		reso = hist_Reso_SDvar13.GetBinContent(iBin)
		print("{} {} {} {}".format(KS, chi2, reso, binEdge))

	for iBin in range(100):
		binEdge = hist_chi2_jetPtMaxDPhi.GetXaxis().GetBinUpEdge(iBin)
		KS = hist_KS_jetPtMaxDPhi.GetBinContent(iBin)
		chi2 = hist_chi2_jetPtMaxDPhi.GetBinContent(iBin)
		reso = hist_Reso_jetPtMaxDPhi.GetBinContent(iBin)
		print("{} {} {} {}".format(KS, chi2, reso, binEdge))
	"""
	


	

def addLoop():
	baseClass.loop = loop

def deltaPhi(phi1, phi2):
	x = phi1 - phi2
	while x >= rt.TMath.Pi():
		x = x - 2*rt.TMath.Pi()
	while x < -rt.TMath.Pi():
		x = x + 2*rt.TMath.Pi()
	return x

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def Dkl(A, B):
	# takes two histograms with same binning,
	# and computes the Kullback-Liebler Divergence
	# defined as sum over bins (i) of A[i]*log(A[i]/B[i])
	res = 0.
	for i in range(A.GetNbinsX()):
		Ai = A.GetBinContent(i)
		Bi = B.GetBinContent(i)
		if Bi > 0. and Ai > 0.:
			res += Ai*rt.TMath.Log10(Ai/Bi)
	return res








































