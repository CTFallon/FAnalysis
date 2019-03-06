from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

def loop(self):
	

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	newTree = rt.TTree("newTree","newTree")
	self.objects.append(newTree)
	dijMT = array('f',[0])
	triMT = array('f',[0])
	algMT = array('f',[0])
	passFinal = array('i',[0])
	weight = array('f',[0])
	useTriJetMT = array('i',[0])

	newTree.Branch("dijet_MT",dijMT,"dijetMT/F")
	newTree.Branch("trijet_MT",triMT,"trijetMT/F")
	newTree.Branch("algo_MT",algMT,"algoMT/F")
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
	
	hist_algoCuts_nEv = self.makeTH1F("hist_algoCuts_nEv","Undecided;Cut;count",6,0,6)
	hist_algoCuts_dij = self.makeTH1F("hist_algoCuts_dij","Dijet;Cut;count",6,0,6)
	hist_algoCuts_dijNew = self.makeTH1F("hist_algoCuts_dijNew","New Dijet;Cut;count",6,0,6)
	hist_algoCuts_tri = self.makeTH1F("hist_algoCuts_tri","Trijet;Cut;count",6,0,6)
	hist_algoCuts_triNew = self.makeTH1F("hist_algoCuts_triNew","New Trijet;Cut;count",6,0,6)
	hist_algoCuts_dijTruthD = self.makeTH1F("hist_algoCuts_dijTruthD","only Two Jets;Cut;count",6,0,6)
	hist_algoCuts_dijTruthT = self.makeTH1F("hist_algoCuts_dijTruthT","3+ jets, Trijet;Cut;count",6,0,6)
	hist_algoCuts_dijTruthO = self.makeTH1F("hist_algoCuts_dijTruthO","3+ Jets, Dijet;Cut;count",6,0,6)
	# 0 if passedPreselection
	# 1 is has 3 or more AK8 Jets
	# 2 is has Jet 2 be maxDeltaPhi
	# 3 is gamma3 < 5.82
	# 4 is total
	# 5 is Truth
	# each time we fill "NEW", we need to also fill in all of the bins to the left in the not-"NEW"
	# at each cut, when we don't assign an MT, we need to fill in nEv

	nDij_2Jet = 0


	# set up trees
	nEvents = 0
	for rootID in self.inputFileList:
		f = rt.TFile.Open(rootID)
		tree = f.Get("tree")
		oldnEvents = nEvents
		nEvents += tree.GetEntries()
		print(rootID.split("_")[2:])
		for iEvent in range(tree.GetEntries()):
			dijMT[0] = 0.
			triMT[0] = 0.
			algMT[0] = 0.
			passFinal[0] = 0
			weight[0] = 0.
			useTriJetMT[0] = -1
			if iEvent%1000 == 0:
				print("Event: " + str(iEvent+oldnEvents) + "/" + str(nEvents))
			tree.GetEvent(iEvent)
			W = tree.Weight
			weight[0] = W
			# sepcial TTJets Stitching
			#TTJets                                 madHT < 600 and GenElectrons->size()==0 and GenMuons->size()==0 and GenTaus->size()==0
			#TTJets	DiLept                          madHT < 600 and GenMET < 150                                                       
			#TTJets	SingleLeptFromT                 madHT < 600 and GenMET < 150                                                         
			#TTJets	SingleLeptFromTbar              madHT < 600 and GenMET < 150                                                          
			#TTJets	DiLept				genMET150   madHT < 600 and GenMET >= 150                                                        
			#TTJets	SingleLeptFromT		genMET150   madHT < 600 and GenMET >= 150                                                         
			#TTJets	SingleLeptFromTbar	genMET150   madHT < 600 and GenMET >= 150                                                        
			#TTJets	HT600to800                      madHT >= 600                                                                         
			#TTJets	HT800to1200                     madHT >= 600                                                                         
			#TTJets	HT1200to2500                    madHT >= 600                                                                         
			#TTJets	HT2500toInf                     madHT >= 600  

			if rootID.split("_")[2] == "TTJets":
				disc = rootID.split("_")
				if disc[3][:2] == "MC":
					if not (tree.madHT < 600 and tree.GenElectrons.size()==0 and tree.GenMuons.size()==0 and tree.GenTaus.size()==0):
						continue
				elif disc[3][:2] == "HT":
					if not (tree.madHT >= 600):
						continue
				elif disc[4][:2] == "MC":
					if not (tree.madHT < 600 and tree.GenMET < 150):
						continue
				elif disc[4][:2] == "ge":
					if not (tree.madHT < 600 and tree.GenMET >= 150):
						continue
				else:
					print("TTBar stitiching error!")
					print(rootID.split("_")[2:])

			useTri = -1
			hist_algoCuts_nEv.Fill(0,W)
			if useTri == -1: # always check to make sure we still need to decide for this event
				if len(tree.JetsAK8) <= 2: # First check, make sure we have 3 or more jets
					useTri = 0
					hist_algoCuts_dijNew.Fill(1,W)
					hist_algoCuts_dij.Fill(2,W)
					hist_algoCuts_dij.Fill(3,W)
					hist_algoCuts_dij.Fill(4,W)
				else:
					hist_algoCuts_nEv.Fill(1,W)
			if useTri == -1: #  # always check to make sure we still need to decide for this event
				#find iJetMaxDeltaPhi
				deltaPhis = [deltaPhi(tree.JetsAK8[0].Phi(),tree.METPhi),deltaPhi(tree.JetsAK8[1].Phi(),tree.METPhi),deltaPhi(tree.JetsAK8[2].Phi(),tree.METPhi)]
				jetIdx = deltaPhis.index(max(deltaPhis))
				if jetIdx == 1: # Second check, if jet 2 (index 1) is maxDelta phi, we have Trijet
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
			if useTri == -1:
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
				nDij_2Jet += 1
				hist_algoCuts_dijTruthD.Fill(5,W)
			else:
				MT3 = trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)
				triMT[0] = MT3
				hist_MT_tri.Fill(MT3,W)
				if useTri == 1:
					useTriJetMT[0] = 1
					hist_MT_alg.Fill(MT3,W)
					algMT[0] = MT3
					hist_algoCuts_dijTruthT.Fill(5,W)
				elif useTri == 0:
					useTriJetMT[0] = 0
					hist_MT_alg.Fill(MT2,W)
					algMT[0] = MT2
					hist_algoCuts_dijTruthO.Fill(5,W)
				else:
					print("Dun goofed, pt 2 electric boogaloo!")
			if tree.MET/tree.MT_AK8 > 0.25 and tree.DeltaPhiMin_AK8 < 0.75:
				passFinal[0] = 1
			else:
				passFinal[0] = 0
			newTree.Fill()
			
		f.Close()
	self.outRootFile.cd()
	print("bkg n2jets nDijet nTrijet")
	print(self.fileID + " {} {} {}".format(hist_algoCuts_dijTruthD.GetBinContent(6),hist_algoCuts_dijTruthO.GetBinContent(6),hist_algoCuts_dijTruthT.GetBinContent(6)))


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
	hist_algoCuts_dijTruthT.SetFillColor(rt.kRed+2)
	hist_algoCuts_dijTruthO.SetFillColor(rt.kGreen)
	cutStack.Add(hist_algoCuts_dij)
	cutStack.Add(hist_algoCuts_dijNew)
	cutStack.Add(hist_algoCuts_tri)
	cutStack.Add(hist_algoCuts_triNew)
	cutStack.Add(hist_algoCuts_nEv)
	cutStack.Add(hist_algoCuts_dijTruthD)
	cutStack.Add(hist_algoCuts_dijTruthO)
	cutStack.Add(hist_algoCuts_dijTruthT)
	hist_algoCuts_dij.GetXaxis().SetBinLabel(1,"PreSelection")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(2,"Num. of AK8 Jets")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(3,"#Delta#phi^{max} Index")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(4,"#gamma_{3}")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(5,"Final")
	hist_algoCuts_dij.GetXaxis().SetBinLabel(6,"Truth")
	cutStack.Draw("hist")
	c1.SaveAs(self.extraDir+"algorithmFlow_bkg.png")


	resoDij = hist_MT_dij.GetRMS()/hist_MT_dij.GetMean()
	resoTri = hist_MT_tri.GetRMS()/hist_MT_tri.GetMean()
	resoAlg = hist_MT_alg.GetRMS()/hist_MT_alg.GetMean()

	self.makeRatio(hist_MT_dij, hist_MT_alg ,"3JetMT_Final_MT_bkg_ratioalgdij", doLeg = True, log = False)
	self.makeRatio(hist_MT_dij, hist_MT_tri ,"3JetMT_Final_MT_bkg_ratiotridij", doLeg = True, log = False)
	self.makePng([hist_MT_dij,hist_MT_tri,hist_MT_alg], "3JetMT_Final_bkg_4MTplots", doCum = True)

	#nTriTruthDij = hist_algoCuts_triTruthD.GetBinContent(6)
	#nTriTruthTri = hist_algoCuts_triTruthT.GetBinContent(6)
	#nTriTruthOth = hist_algoCuts_triTruthO.GetBinContent(6)

	

	#print("file nTriTruthDij nTriTruthTri nTriTruthOth nDijTruthDij nDijTruthTri nDijTruthOth nDijTruthDij2Jet nDijTruthOth2Jet resoDij resoTri resoAlg resoPrf ksDij ksTri ksAlg")
	#print(self.fileID + " {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(
	#		nTriTruthDij, nTriTruthTri, nTriTruthOth,
	#		nDijTruthDij, nDijTruthTri, nDijTruthOth,
	#		nDijTruthDij_2Jet, nDijTruthOth_2Jet,
	#		resoDij, resoTri, resoAlg, resoPrf,
	#		ksDij, ksTri, ksAlg))

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

def deltaPhi(phi1, phi2):
	dphi = phi1 - phi2
	while dphi >= rt.TMath.Pi():
		dphi -= 2*rt.TMath.Pi()
	while dphi < -rt.TMath.Pi():
		dphi += 2*rt.TMath.Pi()
	return abs(dphi)

