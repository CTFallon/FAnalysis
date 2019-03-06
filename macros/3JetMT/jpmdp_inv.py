from analysisBase import baseClass
import ROOT as rt
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

# macro to investigate the jetPtMaxDPhi method to select
# events to use 3 jet MT on

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
	
	#Look at deltaPhi  and deltaR between all FSR jets

	hist_deltaPhi_FSR_all = self.makeTH1F(
		"hist_deltaPhi_FSR_all",
		"AllEvents;deltaPhi;count",
		100,0,rt.TMath.Pi())

	hist_deltaPhi_FSR_2jet = self.makeTH1F(
		"hist_deltaPhi_FSR_2jet",
		"TwoJetEvents;deltaPhi;count",
		100,0,rt.TMath.Pi())

	hist_deltaPhi_FSR_3jet = self.makeTH1F(
		"hist_deltaPhi_FSR_3jet",
		"ThreeJetsEvents;deltaPhi;count",
		100,0,rt.TMath.Pi())

	hist_deltaPhi_FSR_min = self.makeTH1F(
		"hist_deltaPhi_FSR_min",
		"MinDeltaPhi;deltaPhi;count",
		100,0,rt.TMath.Pi())

	hist_deltaR_FSR_all = self.makeTH1F(
		"hist_deltaR_FSR_all",
		"AllEvents;deltaR;count",
		100,0.8,6)

	hist_deltaR_FSR_2jet = self.makeTH1F(
		"hist_deltaR_FSR_2jet",
		"TwoJetEvents;deltaR;count",
		100,0.8,6)

	hist_deltaR_FSR_3jet = self.makeTH1F(
		"hist_deltaR_FSR_3jet",
		"ThreeJetsEvents;deltaR;count",
		100,0.8,6)

	hist_deltaR_FSR_min = self.makeTH1F(
		"hist_deltaR_FSR_min",
		"MinDeltaR;deltaR;count",
		100,0.8,6)

	# compare deltaR_FSRMin and deltaPhi_FSRMin to the results of
	# the FSR Tag
	# three categories: FSR2, FSR3, other
	# FSR2 has only 2 FSR jets
	# FSR3 has only 3 FSR jets
	# other has more than 3 FSR jets


	hist_deltaPhi_FSRmin_FSR2 = self.makeTH1F(
		"hist_deltaPhi_FSRmin_FSR2",
		"deltPhi_FSRmin_FSR2",100,0,rt.TMath.Pi())
	hist_deltaPhi_FSRmin_FSR3= self.makeTH1F(
		"hist_deltaPhi_FSRmin_FSR3",
		"deltPhi_FSRmin_FSR3",100,0,rt.TMath.Pi())
	hist_deltaPhi_FSRmin_other= self.makeTH1F(
		"hist_deltaPhi_FSRmin_other",
		"deltPhi_FSRmin_other",100,0,rt.TMath.Pi())


	hist_deltaR_FSRmin_FSR2 = self.makeTH1F(
		"hist_deltaR_FSRmin_FSR2",
		"deltR_FSRmin_FSR2",100,0.8,6)
	hist_deltaR_FSRmin_FSR3 = self.makeTH1F(
		"hist_deltaR_FSRmin_FSR3",
		"deltR_FSRmin_FSR3",100,0.8,6)
	hist_deltaR_FSRmin_other = self.makeTH1F(
		"hist_deltaR_FSRmin_other",
		"deltR_FSRmin_other",100,0.8,6)

	hist_jpmdp_all = self.makeTH1F(
		"hist_jpmdp_all",
		"jpmdp_all",100,0,2000)
	hist_jpmdp_FSR2 = self.makeTH1F(
		"hist_jpmdp_FSR2",
		"jpmdp_FSR2",100,0,2000)
	hist_jpmdp_FSR3 = self.makeTH1F(
		"hist_jpmdp_FSR3",
		"jpmdp_FSR3",100,0,2000)
	hist_jpmdp_other = self.makeTH1F(
		"hist_jpmdp_other",
		"jpmdp_other",100,0,2000)

	# Plot distributions of deltaPhi between jet(MaxDPhi) and other jets.
	# plot distributions of deltaPhi between jet(maxDPhi) and closest jet
	# for cases where both jets are FSR, and cases where only one jet is FSR

	histList_deltaPhi_jetMaxDPhi = []
	histList_deltaPhi_jetMaxDPhi.append(
		self.makeTH1F("hist_deltaPhi_jetMaxDPhi_jet1",
						"maxDPhi_deltaPhi_jet1;deltaPhi;",
						100, 0 , rt.TMath.Pi()
					))
	histList_deltaPhi_jetMaxDPhi.append(
		self.makeTH1F("hist_deltaPhi_jetMaxDPhi_jet2",
						"maxDPhi_deltaPhi_jet2;deltaPhi;",
						100, 0 , rt.TMath.Pi()
					))
	histList_deltaPhi_jetMaxDPhi.append(
		self.makeTH1F("hist_deltaPhi_jetMaxDPhi_jet3",
						"maxDPhi_deltaPhi_jet3;deltaPhi;",
						100, 0 , rt.TMath.Pi()
					))

	hist_deltaPhi_jetMaxDPhi_closestJet = self.makeTH1F(
		"hist_deltaPhi_jetMaxDPhi_closestJet",
		"jetMaxDPhi_closestJet;deltaPhi;",
		100,0,rt.TMath.Pi()
		)
	hist_deltaPhi_jetMaxDPhi_closestJet_FSR = self.makeTH1F(
		"hist_deltaPhi_jetMaxDPhi_closestJet_FSR",
		"jetMaxDPhi_closestJet_FSR;deltaPhi;",
		100,0,rt.TMath.Pi()
		)
	hist_deltaPhi_jetMaxDPhi_closestJet_ISR = self.makeTH1F(
		"hist_deltaPhi_jetMaxDPhi_closestJet_ISR",
		"jetMaxDPhi_closestJet_ISR;deltaPhi;",
		100,0,rt.TMath.Pi()
		)
	hist_deltaPhi_jetMaxDPhi_ISR_closestJet = self.makeTH1F(
		"hist_deltaPhi_jetMaxDPhi_ISR_closestJet",
		"jetMaxDPhi_ISR_closestJet;deltaPhi;",
		100,0,rt.TMath.Pi()
		)
	hist_deltaPhi_jetMaxDPhi_closestJet_2ISR = self.makeTH1F(
		"hist_deltaPhi_jetMaxDPhi_closestJet_2ISR",
		"jetMaxDPhi_closestJet_2ISR;deltaPhi;",
		100,0,rt.TMath.Pi()
		)

	hist_jetMaxDPhiIndex = self.makeTH1F(
		"hist_jetMaxDPhiIndex",
		"jetMaxDPhiIndex;index;",
		6,0,6
		)
	hist2d_jetMaxDPhiIndex_nFSRJets = self.makeTH2F(
		"hist2d_jetMaxDPhiIndex_nFSRJets",
		"JetMaxDPhiIndex_vs_nFSRJets;Jet Index;num FSR jets",
		6,0,6,
		6,0,6
		)

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		if nJets != 3:
			continue
		listOfDeltaPhi = []
		listOfDeltaR = []
		jetsFSR = []
		nFSRJets = 0
		nJetsToCheck = nJets #nJets, min(nJets,3)
		for iJet in range(nJetsToCheck):
			if not tree.JetsAK8_isISR[iJet]:
				jetsFSR.append(iJet)
				nFSRJets += 1
			if not tree.JetsAK8_isISR[iJet] and iJet+1 != min(nJets, 3):
			#if iJet+1 <= nJetsToCheck:
				for jJet in range(iJet+1,nJetsToCheck):
					if not tree.JetsAK8_isISR[jJet]:
					#if True:
						deltaPhi = abs(tree.JetsAK8[iJet].DeltaPhi(
										tree.JetsAK8[jJet]
										))
						deltaR = abs(tree.JetsAK8[iJet].DeltaR(
										tree.JetsAK8[jJet]
										))
						hist_deltaPhi_FSR_all.Fill(deltaPhi)
						hist_deltaR_FSR_all.Fill(deltaR)
						if nJets == 2:
							hist_deltaPhi_FSR_2jet.Fill(deltaPhi)
							hist_deltaR_FSR_2jet.Fill(deltaR)
						else:
							hist_deltaPhi_FSR_3jet.Fill(deltaPhi)
							hist_deltaR_FSR_3jet.Fill(deltaR)
						listOfDeltaPhi.append(deltaPhi)
						listOfDeltaR.append(deltaR)
		#if nFSRJets < 3:
		#	continue
		if nJets == 2:
		#	continue
			deltaPhiMetList = [tree.DeltaPhi1, tree.DeltaPhi2]
		else:
			deltaPhiMetList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
		jetMaxDPhiIndex = deltaPhiMetList.index(max(deltaPhiMetList))
		jetPtMaxDPhi = tree.JetsAK8[jetMaxDPhiIndex].Pt()

		# deltaPhi between maxDPhiJet and otherJets
		# deltaPhi between maxDPhiJet and closestJet, FSR/ISR

		closestJetIndex = -1
		closestJetIsFSR = False
		closestJetDPhi = 4
		for iJet in range(min(nJets,3)):
			if iJet == jetMaxDPhiIndex: # dont need to compare jet to itself
				continue
			iJetDPhi = abs(tree.JetsAK8[jetMaxDPhiIndex].DeltaPhi(
													tree.JetsAK8[iJet]))
			if iJetDPhi < closestJetDPhi:
				closestJetIndex = iJet
				closestJetIsISR = tree.JetsAK8_isISR[iJet]
				closestJetDPhi = iJetDPhi

			histList_deltaPhi_jetMaxDPhi[iJet].Fill(iJetDPhi)
		
		hist_deltaPhi_jetMaxDPhi_closestJet.Fill(closestJetDPhi)
		if ((not tree.JetsAK8_isISR[jetMaxDPhiIndex]) and (not closestJetIsISR)):
			hist_deltaPhi_jetMaxDPhi_closestJet_FSR.Fill(closestJetDPhi)
		elif ((not tree.JetsAK8_isISR[jetMaxDPhiIndex]) and closestJetIsISR):
			hist_deltaPhi_jetMaxDPhi_closestJet_ISR.Fill(closestJetDPhi)
		elif (tree.JetsAK8_isISR[jetMaxDPhiIndex] and (not closestJetIsISR)):
			hist_deltaPhi_jetMaxDPhi_ISR_closestJet.Fill(closestJetDPhi)
		else:
			hist_deltaPhi_jetMaxDPhi_closestJet_2ISR.Fill(closestJetDPhi)
		hist_jetMaxDPhiIndex.Fill(jetMaxDPhiIndex)
		hist2d_jetMaxDPhiIndex_nFSRJets.Fill(jetMaxDPhiIndex,nFSRJets)
		# determine the 'true' jets for MT
		if (0 in jetsFSR) and (1 in jetsFSR) and (len(jetsFSR) == 2):
			FSRJetCat = 'a'
		elif (
				(0 in jetsFSR)
				and (1 in jetsFSR)
				and (2 in jetsFSR)
				and (len(jetsFSR) == 3)
			):
			FSRJetCat = 'b'
		else:
			FSRJetCat = 'c'

		hist_jpmdp_all.Fill(jetPtMaxDPhi)
		if FSRJetCat == 'a':
			hist_jpmdp_FSR2.Fill(jetPtMaxDPhi)
		elif FSRJetCat == 'b':
			hist_jpmdp_FSR3.Fill(jetPtMaxDPhi)
		elif FSRJetCat == 'c':
			hist_jpmdp_other.Fill(jetPtMaxDPhi)
		try:
			hist_deltaPhi_FSR_min.Fill(min(listOfDeltaPhi))
			hist_deltaR_FSR_min.Fill(min(listOfDeltaR))
			if FSRJetCat == 'a':
				hist_deltaPhi_FSRmin_FSR2.Fill(min(listOfDeltaPhi))
				hist_deltaR_FSRmin_FSR2.Fill(min(listOfDeltaR))
			elif FSRJetCat == 'b':
				hist_deltaPhi_FSRmin_FSR3.Fill(min(listOfDeltaPhi))
				hist_deltaR_FSRmin_FSR3.Fill(min(listOfDeltaR))
			elif FSRJetCat == 'c':
				hist_deltaPhi_FSRmin_other.Fill(min(listOfDeltaPhi))
				hist_deltaR_FSRmin_other.Fill(min(listOfDeltaR))
			else:
				print("No FSR Jet Category!")
		except ValueError:
			continue


	makePlots([hist_deltaPhi_FSR_all,
				hist_deltaPhi_FSR_2jet,
				hist_deltaPhi_FSR_3jet,
				hist_deltaPhi_FSR_min],
				self.extraDir, "deltaPhi_FSR",True)

	makePlots([hist_deltaR_FSR_all,
				hist_deltaR_FSR_2jet,
				hist_deltaR_FSR_3jet,
				hist_deltaR_FSR_min],
				self.extraDir, "deltaR_FSR",True)

	makePlots([hist_deltaPhi_FSR_min,
				hist_deltaPhi_FSRmin_FSR2,
				hist_deltaPhi_FSRmin_FSR3,
				hist_deltaPhi_FSRmin_other],
				self.extraDir, "deltaPhi_FSRmin",True)
	makePlots([hist_deltaR_FSR_min,
				hist_deltaR_FSRmin_FSR2,
				hist_deltaR_FSRmin_FSR3,
				hist_deltaR_FSRmin_other],
				self.extraDir, "deltaR_FSRmin",True)
	makePlots([hist_jpmdp_all,
				hist_jpmdp_FSR2,
				hist_jpmdp_FSR3,
				hist_jpmdp_other],
				self.extraDir, "jpmdp_3omJets", True)

	makePlots([histList_deltaPhi_jetMaxDPhi[1],
				histList_deltaPhi_jetMaxDPhi[0],
				histList_deltaPhi_jetMaxDPhi[2]],
				self.extraDir,"deltaPhi_jetMaxDPhi_otherJets",True)

	makePlots([
				hist_deltaPhi_jetMaxDPhi_closestJet_FSR,
				hist_deltaPhi_jetMaxDPhi_closestJet_ISR,
				hist_deltaPhi_jetMaxDPhi_ISR_closestJet,
				hist_deltaPhi_jetMaxDPhi_closestJet_2ISR],
				self.extraDir,"deltaPhi_jetMaxDPhi_closestJet",True)
	
	makePlots([hist_jetMaxDPhiIndex],
				self.extraDir,
				"jetMaxDPhi_Index",
				log = True,
				lC = [0.7,0.7,0.9,0.9])
	make2DPlot(hist2d_jetMaxDPhiIndex_nFSRJets, self.extraDir, "2d_indexVsnFSR")

def addLoop():
	baseClass.loop = loop

def makePlots(hList,direct, name, log = False, lC = [0.3,0.2,0.7,0.3]):
	c1 = rt.TCanvas("c1","c1",900,600)
	if log:
		c1.SetLogy()
	for iH in range(len(hList)):
		if iH == 0:
			hList[iH].SetLineColor(1)
			hList[iH].Draw()
		else:
			hList[iH].SetLineColor(iH+1)
			hList[iH].Draw('same')
	c1.BuildLegend(lC[0],lC[1],lC[2],lC[3])
	c1.SaveAs(direct+name+".png")

def make2DPlot(hist, direct, name, logX = False, logY = False, logZ = False):
	c1 = rt.TCanvas("c1","c1",900,600)
	if logX:
		c1.SetLogx()
	if logY:
		c1.SetLogy()
	if logZ:
		c1.SetLogz()
	hist.Draw("colz")
	c1.SaveAs(direct+name+'.png')


