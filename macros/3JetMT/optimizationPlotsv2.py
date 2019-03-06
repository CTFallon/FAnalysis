from analysisBase import baseClass
import ROOT as rt
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

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


	print(self.fileID[:2])
	print(self.fileID[-1:])
	print(self.fileID[-2:])

	if self.fileID[:2] == "mz":
		if self.fileID[-1] == "5":
			zMass = int(self.fileID[-2:])*100
		else:
			zMass = int(self.fileID[-1:])*1000
	else:	
		zMass = 3000
	print("zMass is {}".format(zMass))

	mtBins = 100
	hist_MT_FSR = self.makeTH1F("hist_MT_FSR","FSR;MT;",mtBins,0,2*zMass)
	hist_MT_Dijet = self.makeTH1F("hist_MT_Dijet","Dijet;MT;",mtBins,0,2*zMass)
	hist_MT_Trijet = self.makeTH1F("hist_MT_Trijet","Trijet;MT;",mtBins,0,2*zMass)
	hist_MT_SDVar13 = self.makeTH1F("hist_MT_SDVar13","SDVar13;MT;",mtBins,0,2*zMass)
	hist_MT_JetPtMaxDPhi = self.makeTH1F("hist_MT_JetPtMaxDPhi","JetPtMaxDPhi;MT;",mtBins,0,2*zMass)
	hist_MT_Jet3DPhi = self.makeTH1F("hist_MT_Jet3DPhi","Jet3DPhi;MT;",mtBins,0,2*zMass)

	nFSRDi = 0
	nFSRTr = 0
	nFSRot = 0
	nSDDi = 0
	nSDTr = 0
	nJPDi = 0
	nJPTr = 0
	nJ3Di = 0
	nJ3Tr = 0
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		
		# optimal cuts for sdvar13 and jetptmaxdphi are 32 and 860
	
		# This macro is to:
		# 3rd jet selection: 
		#	remake true/false positive/negative tables
		#	2D plots
		#	chi^2 and resolution vs. parameters
		
		#make FSR histogram
		FSRJets = []
		FSRIndex = []
		for iJet in range(nJets):
			if tree.JetsAK8_isHV[iJet]:
				FSRJets.append(tree.JetsAK8[iJet])
				FSRIndex.append(iJet)
		hist_MT_FSR.Fill(trans_mass_Njet(FSRJets,tree.MET,tree.METPhi))


		if nJets > 2:
			if ((0 in FSRIndex)
				and (1 in FSRIndex)
				and not (2 in FSRIndex)
				and not (3 in FSRIndex)
				and not (4 in FSRIndex)
				and not (5 in FSRIndex)):
				nFSRDi += 1
			elif ((0 in FSRIndex)
				and (1 in FSRIndex)
				and (2 in FSRIndex)
				and not (3 in FSRIndex)
				and not (4 in FSRIndex)
				and not (5 in FSRIndex)):
				nFSRTr += 1
			else:
				nFSRot += 1
		

		#Step 1, deal with events that only have 2 jets
		MT12 = trans_mass_Njet(tree.JetsAK8[0:2],tree.MET,tree.METPhi)
		hist_MT_Dijet.Fill(MT12)
		if nJets == 2: # fill MT histograms
			hist_MT_SDVar13.Fill(MT12)
			hist_MT_Trijet.Fill(MT12)
			hist_MT_JetPtMaxDPhi.Fill(MT12)
			hist_MT_Jet3DPhi.Fill(MT12)
		else:
			MT123 = trans_mass_Njet(tree.JetsAK8[0:3],tree.MET,tree.METPhi)
			#make variables
			hist_MT_Trijet.Fill(MT123)
			SDVar13 = tree.JetsAK8[2].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt())
			JetPtMaxDPhi = tree.JetsAK8[[tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3].index(max(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3))].Pt()
			MinDPhiIndex = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3].index(min(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3))

			if SDVar13 < .321:
				hist_MT_SDVar13.Fill(MT12)
				nSDDi += 1
			else:
				hist_MT_SDVar13.Fill(MT123)
				nSDTr += 1

			if JetPtMaxDPhi > 840:
				hist_MT_JetPtMaxDPhi.Fill(MT12)
				nJPDi += 1
			else:
				hist_MT_JetPtMaxDPhi.Fill(MT123)
				nJPTr += 1

			if MinDPhiIndex == 2:
				hist_MT_Jet3DPhi.Fill(MT123)
				nJ3Tr += 1
			else:
				hist_MT_Jet3DPhi.Fill(MT12)
				nJ3Di += 1
	
	rt.gStyle.SetOptTitle(0)

	hist_MT_FSR.SetLineColor(1)
	hist_MT_Dijet.SetLineColor(2)
	hist_MT_Trijet.SetLineColor(3)
	hist_MT_JetPtMaxDPhi.SetLineColor(4)

	hist_MT_FSR.SetTitle("M_{T}(#chi_{1},#chi_{2})")
	hist_MT_Dijet.SetTitle("M_{T}(j_{1},j_{2})")
	hist_MT_Trijet.SetTitle("M_{T}(j_{1},j_{2},[j_{3}])")
	hist_MT_JetPtMaxDPhi.SetTitle("p_{T,#Delta#phi_{max}(Jet, MET)} Algorithm")

	hist_MT_FSR.SetLineWidth(2)
	hist_MT_Dijet.SetLineWidth(2)
	hist_MT_Trijet.SetLineWidth(2)
	hist_MT_JetPtMaxDPhi.SetLineWidth(2)

	hist_MT_FSR.SetLineStyle(1)
	hist_MT_Dijet.SetLineStyle(2)
	hist_MT_Trijet.SetLineStyle(2)
	hist_MT_JetPtMaxDPhi.SetLineStyle(1)
	
	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_solid")
	makePlot([
		hist_MT_FSR,
		hist_MT_JetPtMaxDPhi,
		hist_MT_Dijet,
		hist_MT_Trijet], self.extraDir, "DOE_prop_log_solid",logY=True)

	#hist_MT_FSR.SetLineColor(rt.kBlack)
	#hist_MT_Dijet.SetLineColor(rt.kGray)
	#hist_MT_Trijet.SetLineColor(rt.kGray+1)
	#hist_MT_JetPtMaxDPhi.SetLineColor(rt.kGray+2)

	
	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_gs_solid")
	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_gs_log_solid",logY=True)

	#hist_MT_FSR.SetLineStyle(2)
	#hist_MT_Dijet.SetLineStyle(2)
	#hist_MT_Trijet.SetLineStyle(3)
	#hist_MT_JetPtMaxDPhi.SetLineStyle(3)

	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_gs_dashed")
	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_gs_log_dashed",logY=True)

	#hist_MT_FSR.SetLineColor(1)
	#hist_MT_Dijet.SetLineColor(2)
	#hist_MT_Trijet.SetLineColor(3)
	#hist_MT_JetPtMaxDPhi.SetLineColor(4)

	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_dashed")
	#makePlot([
	#	hist_MT_FSR,
	#	hist_MT_Dijet,
	#	hist_MT_Trijet,
	#	hist_MT_JetPtMaxDPhi], self.extraDir, "DOE_prop_log_dashed",logY=True)

	rt.gStyle.SetOptTitle(1)
	#KS tests
	KS_FSR = hist_MT_FSR.KolmogorovTest(hist_MT_FSR)
	KS_Dijet = hist_MT_FSR.KolmogorovTest(hist_MT_Dijet)
	KS_SDVar13 = hist_MT_FSR.KolmogorovTest(hist_MT_SDVar13)
	KS_JetPtMaxDPhi = hist_MT_FSR.KolmogorovTest(hist_MT_JetPtMaxDPhi)
	KS_Jet3DPhi = hist_MT_FSR.KolmogorovTest(hist_MT_Jet3DPhi)
	#Chi2 tests
	Chi2_FSR = hist_MT_FSR.Chi2Test(hist_MT_FSR)
	Chi2_Dijet = hist_MT_FSR.Chi2Test(hist_MT_Dijet)
	Chi2_SDVar13 = hist_MT_FSR.Chi2Test(hist_MT_SDVar13)
	Chi2_JetPtMaxDPhi = hist_MT_FSR.Chi2Test(hist_MT_JetPtMaxDPhi)
	Chi2_Jet3DPhi = hist_MT_FSR.Chi2Test(hist_MT_Jet3DPhi)
	#MT Resolution
	Reso_FSR = hist_MT_FSR.GetRMS()/hist_MT_FSR.GetMean()
	Reso_Dijet = hist_MT_Dijet.GetRMS()/hist_MT_Dijet.GetMean()
	Reso_SDVar13 = hist_MT_SDVar13.GetRMS()/hist_MT_SDVar13.GetMean()
	Reso_JetPtMaxDPhi = hist_MT_JetPtMaxDPhi.GetRMS()/hist_MT_JetPtMaxDPhi.GetMean()
	Reso_Jet3DPhi = hist_MT_Jet3DPhi.GetRMS()/hist_MT_Jet3DPhi.GetMean()

	makePlot([hist_MT_FSR,hist_MT_Dijet,hist_MT_SDVar13,hist_MT_JetPtMaxDPhi,hist_MT_Jet3DPhi], self.extraDir, "optv2")
	makePlot([hist_MT_FSR,hist_MT_Dijet,hist_MT_SDVar13,hist_MT_JetPtMaxDPhi,hist_MT_Jet3DPhi], self.extraDir, "optv2LOG", True)

	fitToBeta(hist_MT_FSR, self.extraDir)
	fitToCB(hist_MT_FSR, self.extraDir)
	fitToCB(hist_MT_Dijet, self.extraDir)
	#fitToCBBetaBW(hist_MT_SDVar13, self.extraDir)
	#fitToCBBetaBW(hist_MT_JetPtMaxDPhi, self.extraDir)
	#fitToCBBetaBW(hist_MT_Jet3DPhi, self.extraDir)

	print("Cut nDijet nTrijet KS Chi2 Reso")
	print("FSR {} {} {} {} {}".format(nFSRDi, nFSRTr, KS_FSR, Chi2_FSR, Reso_FSR))
	print("Dijet {} {} {} {} {}".format(0, 0, KS_Dijet, Chi2_Dijet, Reso_Dijet))
	print("SDVar13 {} {} {} {} {}".format(nSDDi, nSDTr, KS_SDVar13,Chi2_SDVar13,Reso_SDVar13))
	print("JetPtMaxDPhi {} {} {} {} {}".format(nJPDi, nJPTr, KS_JetPtMaxDPhi,Chi2_JetPtMaxDPhi,Reso_JetPtMaxDPhi))
	print("Jet3DPhi {} {} {} {} {}".format(nJ3Di, nJ3Tr, KS_Jet3DPhi,Chi2_Jet3DPhi,Reso_Jet3DPhi))
	


def addLoop():
	baseClass.loop = loop

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def makePlot(ListOfHistos, location, name, logY = False, gs = False):
	c1 = rt.TCanvas("c1","c1",900,600)
	c1.SetLineWidth(2)
	if gs:
		c1.SetGrayscale()
	if logY:
		c1.SetLogy()
	stack = rt.THStack("stack","Transverse Mass for Various Jet Inclusion Methods;M_{T} [GeV];Count")
	stack.SetMinimum(10)
	stack.SetMaximum(1100)	
	#stack.GetXaxis().SetTitle("MT [GeV]")
	#stack.GetYaxis().SetTitle("Count")
	for histo in ListOfHistos:
		#histo.SetLineColor(ListOfHistos.index(histo)+1)
		#histo.SetAxisRange(100,1000,"Y")
		#histo.SetAxisRange(0,4000,"X")
		stack.Add(histo)
	stack.Draw("nostack")
	stack.GetXaxis().SetLimits(0,4000)
	stack.Draw("nostack")
	t1 = rt.TText(50,1200,"CMS Simulation")
	t1.Draw("same")
	#stack.SetAxisRange(10,1000,"Y")
	#c1.Update()
	c1.BuildLegend(0.2,0.7,0.4,0.9)
	c1.SaveAs(location+name+".pdf")

def fitToCB(histo, direc, extraName=""):

	histo.Scale(1/histo.Integral())
	
	cbFunc = rt.TF1("cbFunc","crystalball",0,histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX()))
	cbFunc.SetParameters(
		0.5, # arbitrary, half of bounds below
		histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX())/2, # midpoint of histogram
		histo.GetXaxis().GetBinUpEdge(histo.GetNbinsX())/10, # one-tenth of histogram
		2, # arbitrary
		100 # arbitrary
	)
	
	cbFunc.SetParLimits(0,0,1) # constant
	cbFunc.SetParLimits(1,0,6000) # mean
	cbFunc.SetParLimits(2,10,1000) # sigma
	cbFunc.SetParLimits(3,0,3) # alpha
	#cbFunc.SetParLimits(4,0,100) # N

	cbFunc.SetLineColor(2)

	histo.Fit('cbFunc')
	
	print("Cutoff at {}".format(cbFunc.GetParameter(1)+cbFunc.GetParameter(2)*cbFunc.GetParameter(3)))
	c2 = rt.TCanvas("c2","c2",900,600)
	histo.SetLineColor(1)
	histo.Draw()
	cbFunc.Draw("same")
	c2.SaveAs(direc+"cbFunc_"+histo.GetName()+extraName+".png")

def betaFunc(x, p):
	# p is a list of parameters
	# p[0] = xLo
	# p[1] = xHi
	# p[2] = alpha
	# p[3] = beta
	# xLo and xHi should not be fit to
	# alpha and beta should be fit
	numer = rt.TMath.Gamma(p[2]+p[3])
	denom = rt.TMath.Gamma(p[2])*rt.TMath.Gamma(p[3])
	factor = rt.TMath.Power(x[0]-p[0],p[2]-1)*rt.TMath.Power(p[1]-x[0],p[3]-1)
	return factor*(numer/denom)

def fitToBeta(histo, direc, extraName=""):
	#histo.Scale(1/histo.Integral())
	xLow = histo.GetBinLowEdge(1)
	xHigh = histo.GetBinLowEdge(histo.GetNbinsX()+1)
	btFunc = rt.TF1("btFunc",betaFunc, xLow, xHigh, 4)
	btFunc.SetParNames("xLow","xHigh","alpha","beta")
	btFunc.FixParameter(0, xLow)
	btFunc.FixParameter(1, xHigh)
	#btFunc.FixParameter(2,5)
	#btFunc.FixParameter(3,2)
	btFunc.SetParLimits(2,0,100)
	btFunc.SetParLimits(3,0,100)
	
	histo.Fit('btFunc')
	
	c3 = rt.TCanvas("c3","c3",900,600)
	histo.SetLineColor(1)
	histo.Draw()
	btFunc.Draw("same")
	c3.SaveAs(direc+"btFunc_"+histo.GetName()+extraName+".png") 
	










































