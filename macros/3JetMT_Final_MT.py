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
	tree.SetBranchStatus("Weight",1)
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

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	newTree = rt.TTree("newTree","newTree")
	self.objects.append(newTree)
	dijMT = array('f',[0])
	triMT = array('f',[0])
	algMT = array('f',[0])
	passFinal = array('i',[0])
	weight = array('f',[0])
	ak8MT = array("f",[0])
	useTriJetMT = array('i',[0])

	newTree.Branch("dijet_MT",dijMT,"dijetMT/F")
	newTree.Branch("trijet_MT",triMT,"trijetMT/F")
	newTree.Branch("algo_MT",algMT,"algoMT/F")
	newTree.Branch("ak8_MT",ak8MT,"ak8MT/F")
	newTree.Branch("passFinal",passFinal,"passFinal/I")
	newTree.Branch("weight",weight,"weight/F")
	newTree.Branch("useTri",useTriJetMT,"useTri/I")
	# make MT histograms for three cases:
		# all Dijet
		# All Trijet
		# algorithmic, JetIdx and gamma3 cuts

	hist_MT_dij = self.makeTH1F("hist_MT_dij","Dijet;MT;count",100,0,5000)
	hist_MT_tri = self.makeTH1F("hist_MT_tri","Trijet;MT;count",100,0,5000)
	hist_MT_alg = self.makeTH1F("hist_MT_alg","Algorithm;MT;count",100,0,5000)
	hist_MT_prf = self.makeTH1F("hist_MT_prf","Perfect;MT;count",100,0,5000)
	
	hist_algoCuts_nEv = self.makeTH1F("hist_algoCuts_nEv","Undecided;Cut;count",6,0,6)
	hist_algoCuts_dij = self.makeTH1F("hist_algoCuts_dij","Dijet;Cut;count",6,0,6)
	hist_algoCuts_dijNew = self.makeTH1F("hist_algoCuts_dijNew","New Dijet;Cut;count",6,0,6)
	hist_algoCuts_tri = self.makeTH1F("hist_algoCuts_tri","Trijet;Cut;count",6,0,6)
	hist_algoCuts_triNew = self.makeTH1F("hist_algoCuts_triNew","New Trijet;Cut;count",6,0,6)
	hist_algoCuts_dijTruthD = self.makeTH1F("hist_algoCuts_dijTruthD","DijetTruthP;Cut;count",6,0,6)
	hist_algoCuts_dijTruthT = self.makeTH1F("hist_algoCuts_dijTruthT","DijetTruthN;Cut;count",6,0,6)
	hist_algoCuts_dijTruthO = self.makeTH1F("hist_algoCuts_dijTruthO","DijetTruthO;Cut;count",6,0,6)
	hist_algoCuts_triTruthT = self.makeTH1F("hist_algoCuts_triTruthT","TrijetTruthP;Cut;count",6,0,6)
	hist_algoCuts_triTruthD = self.makeTH1F("hist_algoCuts_triTruthD","TrijetTruthN;Cut;count",6,0,6)
	hist_algoCuts_triTruthO = self.makeTH1F("hist_algoCuts_triTruthO","TrijetTruthO;Cut;count",6,0,6)
	# 0 if passedPreselection
	# 1 is has 3 or more AK8 Jets
	# 2 is has Jet 2 be maxDeltaPhi
	# 3 is gamma3 < 5.82
	# 4 is total
	# 5 is Truth
	# each time we fill "NEW", we need to also fill in all of the bins to the left in the not-"NEW"
	# at each cut, when we don't assign an MT, we need to fill in nEv

	nDijTruthDij_2Jet = 0
	nDijTruthOth_2Jet = 0
	for iEvent in range(nEvents):
		dijMT[0] = 0.
		triMT[0] = 0.
		algMT[0] = 0.
		passFinal[0] = 0
		weight[0] = 0.
		ak8MT[0] = 0.
		useTriJetMT[0] = -1
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		useTri = -1
		W = float(tree.Weight)
		weight[0] = W
		ak8MT[0] = tree.MT_AK8
		hist_algoCuts_nEv.Fill(0,W)
		if useTri == -1: # always check to make sure we still need to decide for this event
			if len(tree.JetsAK8) <= 2: # First check, make sure we have 3 or more jets
				useTri = 0
				hist_algoCuts_dijNew.Fill(1,W)
				hist_algoCuts_dij.Fill(2,W)
				hist_algoCuts_dij.Fill(3,W)
				hist_algoCuts_dij.Fill(4,W)
				if tree.DijetEvent == 1:
					nDijTruthDij_2Jet += 1
				elif tree.OtherjetEvent == 1:
					nDijTruthOth_2Jet += 1
			else:
				hist_algoCuts_nEv.Fill(1,W)
			
		if useTri == -1: #  # always check to make sure we still need to decide for this event
			if tree.iJetMaxDeltaPhi == 1: # Second check, if jet 2 (index 1) is maxDelta phi, we have Trijet
				useTri = 1
				hist_algoCuts_triNew.Fill(2,W)
				hist_algoCuts_tri.Fill(3,W)
				hist_algoCuts_tri.Fill(4,W)
			else:
				hist_algoCuts_nEv.Fill(2,W)
			
		if useTri == -1: #  # always check to make sure we still need to decide for this event
			if tree.JetsAK8[2].Gamma() < 5.82: # Third check, if gamma3 < 5.82, use Tri
				useTri = 1
				hist_algoCuts_triNew.Fill(3,W)
				hist_algoCuts_tri.Fill(4,W)
			else:
				useTri = 0
				hist_algoCuts_dijNew.Fill(3,W)
				hist_algoCuts_dij.Fill(4,W)

		if useTri == 0:
			if tree.DijetEvent == 1:
				hist_algoCuts_dijTruthD.Fill(5,W)
			elif tree.TrijetEvent == 1:
				hist_algoCuts_dijTruthT.Fill(5,W)
			elif tree.OtherjetEvent == 1:
				hist_algoCuts_dijTruthO.Fill(5,W)
		elif useTri == 1:
			if tree.DijetEvent == 1:
				hist_algoCuts_triTruthD.Fill(5,W)
			elif tree.TrijetEvent == 1:
				hist_algoCuts_triTruthT.Fill(5,W)
			elif tree.OtherjetEvent == 1:
				hist_algoCuts_triTruthO.Fill(5,W)
		else:
			print("Dun goofed!")

		MT2 = trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi)
		dijMT[0] = MT2
		hist_MT_dij.Fill(MT2,W)
		if len(tree.JetsAK8)  <= 2:
			useTriJetMT[0] = 0
			triMT[0] = MT2
			algMT[0] = MT2
			hist_MT_tri.Fill(MT2,W)
			hist_MT_alg.Fill(MT2,W)
			hist_MT_prf.Fill(MT2,W)
		else:
			MT3 = trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)
			hist_MT_tri.Fill(MT3,W)
			triMT[0] = MT3
			if useTri == 1:
				useTriJetMT[0] = 1
				algMT[0] = MT3
				hist_MT_alg.Fill(MT3,W)
			elif useTri == 0:
				useTriJetMT[0] = 0
				algMT[0] = MT2
				hist_MT_alg.Fill(MT2,W)
			else:
				print("Dun goofed, pt 2 electric boogaloo!")
			if tree.TrijetEvent == 1:
				hist_MT_prf.Fill(MT3,W)
			else:
				hist_MT_prf.Fill(MT2,W)
		if tree.MET/tree.MT_AK8 > 0.25 and tree.DeltaPhiMin_AK8 < 0.75 and tree.MT_AK8 > 1500:
			passFinal[0] = 1
		else:
			passFinal[0] = 0
		newTree.Fill()
			
	print("bkg n2jets nDijet nTrijet")
	print(self.fileID + " {} {} {}".format(hist_algoCuts_dijNew.GetBinContent(2),hist_algoCuts_dijNew.GetBinContent(4),hist_algoCuts_tri.GetBinContent(5)))


	# configure, draw, and save the plot for "cutflow" MT
	c1 = rt.TCanvas("c1","c1",1200, 800)
	cutStack = rt.THStack()
	hist_algoCuts_dij.SetFillColor(rt.kGreen+2)
	hist_algoCuts_dijNew.SetFillColor(rt.kGreen)
	hist_algoCuts_tri.SetFillColor(rt.kRed+2)
	hist_algoCuts_triNew.SetFillColor(rt.kRed)
	hist_algoCuts_nEv.SetFillColor(rt.kBlack)
	hist_algoCuts_nEv.SetFillStyle(3005)
	hist_algoCuts_dijTruthD.SetFillColor(rt.kGreen+2)
	hist_algoCuts_dijTruthT.SetFillColor(rt.kGreen)
	hist_algoCuts_dijTruthO.SetFillStyle(3002)
	hist_algoCuts_dijTruthO.SetFillColor(rt.kGreen)
	hist_algoCuts_triTruthO.SetFillColor(rt.kRed)
	hist_algoCuts_triTruthD.SetFillColor(rt.kRed)
	hist_algoCuts_triTruthO.SetFillStyle(3002)
	hist_algoCuts_triTruthT.SetFillColor(rt.kRed+2)
	cutStack.Add(hist_algoCuts_dij)
	cutStack.Add(hist_algoCuts_dijNew)
	cutStack.Add(hist_algoCuts_tri)
	cutStack.Add(hist_algoCuts_triNew)
	cutStack.Add(hist_algoCuts_nEv)
	cutStack.Add(hist_algoCuts_dijTruthD)
	cutStack.Add(hist_algoCuts_triTruthD)
	cutStack.Add(hist_algoCuts_dijTruthO)
	cutStack.Add(hist_algoCuts_triTruthO)
	cutStack.Add(hist_algoCuts_dijTruthT)
	cutStack.Add(hist_algoCuts_triTruthT)
	hist_algoCuts_dij.GetXaxis().SetBinLabel(1,"PreSelection")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(2,"Num. of AK8 Jets")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(3,"#Delta#phi^{max} Index")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(4,"#gamma_{3}")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(5,"Final")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(6,"Truth")
	cutStack.Draw()
	c1.SaveAs(self.extraDir+"algorithmFlow.png")


	resoDij = hist_MT_dij.GetRMS()/hist_MT_dij.GetMean()
	resoTri = hist_MT_tri.GetRMS()/hist_MT_tri.GetMean()
	resoAlg = hist_MT_alg.GetRMS()/hist_MT_alg.GetMean()
	resoPrf = hist_MT_prf.GetRMS()/hist_MT_prf.GetMean()

	ksDij = hist_MT_dij.KolmogorovTest(hist_MT_prf)
	ksTri = hist_MT_tri.KolmogorovTest(hist_MT_prf)
	ksAlg = hist_MT_alg.KolmogorovTest(hist_MT_prf)

	self.makePng([hist_MT_dij,hist_MT_tri,hist_MT_alg,hist_MT_prf], "3JetMT_Final_4MTplots", doCum = True)

	nTriTruthDij = hist_algoCuts_triTruthD.GetBinContent(6)
	nTriTruthTri = hist_algoCuts_triTruthT.GetBinContent(6)
	nTriTruthOth = hist_algoCuts_triTruthO.GetBinContent(6)
	nDijTruthDij = hist_algoCuts_dijTruthD.GetBinContent(6)
	nDijTruthTri = hist_algoCuts_dijTruthT.GetBinContent(6)
	nDijTruthOth = hist_algoCuts_dijTruthO.GetBinContent(6)

	


	print("file nTriTruthDij nTriTruthTri nTriTruthOth nDijTruthDij nDijTruthTri nDijTruthOth nDijTruthDij2Jet nDijTruthOth2Jet resoDij resoTri resoAlg resoPrf ksDij ksTri ksAlg")
	print(self.fileID + " {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
			nTriTruthDij, nTriTruthTri, nTriTruthOth,
			nDijTruthDij, nDijTruthTri, nDijTruthOth,
			nDijTruthDij_2Jet, nDijTruthOth_2Jet,
			resoDij, resoTri, resoAlg, resoPrf,
			ksDij, ksTri, ksAlg))

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

