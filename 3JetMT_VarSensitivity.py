from analysisBase import baseClass
import ROOT as rt
from array import array

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	ff = rt.TFile.Open(self.extraDir+"outFriend.root")
	friendTree = ff.Get("friend")
	tree.AddFriend(friendTree)
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))
	
	# Turn off all branches, selective turn on branches
	tree.SetBranchStatus("*", 0)
	tree.SetBranchStatus("RunNum",1)
	tree.SetBranchStatus("EvtNum",1)
	tree.SetBranchStatus("*AK8*", 1)
	#tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	#tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)
	tree.SetBranchStatus("passedPreSelection",1)
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	#Step 1, make iJet vs FracPt for events with nJets
	histList_2d_iJetvsFracPt = [0,0]
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_2Jets", "Events with 2 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 2, 0, 2, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_3Jets", "Events with 3 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 3, 0, 3, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_4Jets", "Events with 4 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 4, 0, 4, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_5Jets", "Events with 5 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 5, 0, 5, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_6Jets", "Events with 6 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 6, 0, 6, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_7Jets", "Events with 7 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 7, 0, 7, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_8Jets", "Events with 8 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 8, 0, 8, 100, -.01, 1.01))
	histList_2d_iJetvsFracPt.append(self.makeTH2F("hist_iJetvsFracPt_9Jets", "Events with 9 Jets;Jet Number;Fraction of Pt from Visible HV Decendants", 9, 0, 9, 100, -.01, 1.01))

	#coarse grading for optimal MT resolution for fracPt cut
	cutFractions = [x*0.01 for x in range(0,100)]
	#first, only vary the cut on one jet at a time
	hist_MTLead2 = self.makeTH1F("hist_MTlead2Jets","Base MT;MT;count/a.u.",100,0,4000)
	histList_MTcut = []
	hist_MTcut = self.makeTH1F("hist_MTcut","'true' MT;MT;count/a.u.",100,0,4000)
	for cutVal in cutFractions:
		histList_MTcut.append(self.makeTH1F("hist_MT2jet_"+str(cutVal),"2jet MT;MT;count/a.u.",100,0,4000))

	# step 3, plot of 'variables of interest'
	# SDVar for various jets
	# jet pt of jet that has highest delta phi
	# others...
	hist_SDVar12_2jet = self.makeTH1F("SDVar_12_2jet","SDvar_12_2jet;sdVar;count/a.u.",100,0,0.5)
	hist_SDVar13_2jet = self.makeTH1F("SDVar_13_2jet","SDvar_13_2jet;sdVar;count/a.u.",100,0,0.5)
	hist_SDVar23_2jet = self.makeTH1F("SDVar_23_2jet","SDvar_23_2jet;sdVar;count/a.u.",100,0,0.5)
	hist_jetPtMaxdPhi_2jet = self.makeTH1F("pTMaxdPhi_2jet","Pt of MaxdPhiJet_2jet;pT;count/a.u.",100,0,3000)

	hist_SDVar12_3jet = self.makeTH1F("SDVar_12_3jet","SDvar_12_3jet;sdVar;count/a.u.",100,0,0.5)
	hist_SDVar13_3jet = self.makeTH1F("SDVar_13_3jet","SDvar_13_3jet;sdVar;count/a.u.",100,0,0.5)
	hist_SDVar23_3jet = self.makeTH1F("SDVar_23_3jet","SDvar_23_3jet;sdVar;count/a.u.",100,0,0.5)
	hist_jetPtMaxdPhi_3jet = self.makeTH1F("pTMaxdPhi_3jet","Pt of MaxdPhiJet_3jet;pT;count/a.u.",100,0,3000)

	# plot ptFrac[2] with 'sensitive variables'
	hist_ptFracJet3_vs_SDVar12 = self.makeTH2F("ist_ptFracJet3_vs_SDVar12","pTFrac of Jet 3 vs SDVar12;ptFrac;SDVar",100,-.01,1.01,100,0,0.5)
	hist_ptFracJet3_vs_SDVar13 = self.makeTH2F("ist_ptFracJet3_vs_SDVar13","pTFrac of Jet 3 vs SDVar13;ptFrac;SDVar",100,-.01,1.01,100,0,0.5)
	hist_ptFracJet3_vs_SDVar23 = self.makeTH2F("ist_ptFracJet3_vs_SDVar23","pTFrac of Jet 3 vs SDVar23;ptFrac;SDVar",100,-.01,1.01,100,0,0.5)
	hist_ptFracJet3_vs_jetPtMaxdPhi = self.makeTH2F("ist_ptFracJet3_vs_jetPtMaxdPhi","pTFrac of Jet 3 vs jetPt(maxdPhi);ptFrac;pT",100,-.01,1.1,100,0,2000)

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print(str(iEvent)+"/"+str(nEvents))
		tree.GetEvent(iEvent)
		met = tree.MET
		metPhi = tree.METPhi
		if tree.passedPreSelection == 1:
			nJets = len(tree.JetsAK8)
			if iEvent > 50000:
				print(iEvent, nJets)
			for iJet in range(nJets):
				if iEvent > 50000:
					print(iEvent, iJet, len(histList_2d_iJetvsFracPT), len(tree.fracPtFromHVQuarks))
				histList_2d_iJetvsFracPt[nJets].Fill(iJet+0.5, tree.fracPtFromHVQuarks[iJet])
			hist_MTLead2.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], met, metPhi))
			if nJets == 2:
				jetsToUse = 2
			elif tree.fracPtFromHVQuarks[2] > 0.2:
				jetsToUse = 3
			else:
				jetsToUse = 2
			if jetsToUse == 3:
				hist_MTcut.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], met, metPhi))
				hist_SDVar12_3jet.Fill(tree.JetsAK8[1].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()))
				hist_SDVar13_3jet.Fill(tree.JetsAK8[2].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()))
				hist_SDVar23_3jet.Fill(tree.JetsAK8[2].Pt()/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt()))
				hist_jetPtMaxdPhi_3jet.Fill(tree.JetsAK8[tree.iJetMaxDeltaPhi].Pt())
			else:
				hist_MTcut.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], met, metPhi))
				hist_SDVar12_2jet.Fill(tree.JetsAK8[1].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()))
				if nJets > 2:
					hist_SDVar13_2jet.Fill(tree.JetsAK8[2].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()))
					hist_SDVar23_2jet.Fill(tree.JetsAK8[2].Pt()/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt()))
				hist_jetPtMaxdPhi_2jet.Fill(tree.JetsAK8[tree.iJetMaxDeltaPhi].Pt())

			if len(tree.JetsAK8) >= 3:
				hist_ptFracJet3_vs_SDVar12.Fill(tree.fracPtFromHVQuarks[2],tree.JetsAK8[1].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()))
				hist_ptFracJet3_vs_SDVar13.Fill(tree.fracPtFromHVQuarks[2],tree.JetsAK8[2].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()))
				hist_ptFracJet3_vs_SDVar23.Fill(tree.fracPtFromHVQuarks[2],tree.JetsAK8[2].Pt()/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt()))
				hist_ptFracJet3_vs_jetPtMaxdPhi.Fill(tree.fracPtFromHVQuarks[2],tree.JetsAK8[tree.iJetMaxDeltaPhi].Pt())
				
			for iCut in range(len(cutFractions)):
				jetsForMt = []
				cutVal = cutFractions[iCut]
				if nJets == 2:
					jetsForMt.append(tree.JetsAK8[0])
					jetsForMt.append(tree.JetsAK8[1])
				else:
					if tree.fracPtFromHVQuarks[0] > 0.0:
						jetsForMt.append(tree.JetsAK8[0])
					if tree.fracPtFromHVQuarks[1] > 0.0:
						jetsForMt.append(tree.JetsAK8[1])
					if tree.fracPtFromHVQuarks[2] > cutVal:
						jetsForMt.append(tree.JetsAK8[2])
				histList_MTcut[iCut].Fill(trans_mass_Njet(jetsForMt, met, metPhi))
	
	print("No cut has Resolution " + str(hist_MTLead2.GetRMS()/hist_MTLead2.GetMean()))
	for histo in histList_MTcut:
		try:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
		except ZeroDivisionError:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is NULL")
	# optimal MT resolution is found at a cutVal of .2. So, for all events with 3 or more jets, if jet 3 has ptFrac > .2, include it. THIS IS FOR BASELINE
					


	

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


