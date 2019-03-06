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
	tree.SetBranchStatus("HT",1)
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

	histDict_nJets = {}
	histDict_nISRJets = {}
	histDict_nFSRJets = {}

	for x in ["all","2","3","4","5","6","7","3+"]:
		histDict_nJets[x] = self.makeTH1F("hist_nJets_"+x+"_all",
							"nJets_"+x+";nJets;",
							7,0,7)
		histDict_nISRJets[x] = self.makeTH1F("hist_nJets_"+x+"_ISR",
							"nISRJets_"+x+";nISRJets;",
							7,0,7)
		histDict_nFSRJets[x] = self.makeTH1F("hist_nJets_"+x+"_FSR",
							"nFSRJets_"+x+";nISRJets;",
							7,0,7)
	
	
	hist2d_nJets_nFSRJets = self.makeTH2F("hist2d_nJets_nFSRJets",
							"TotalvsFSR;TotalJets;FSRJets",
							7,0,7,7,0,7)
	hist2d_nJets_nISRJets = self.makeTH2F("hist2d_nJets_nISRJets",
							"TotalvsISR;TotalJets;ISRJets",
							7,0,7,7,0,7)
	hist2d_nFSRJets_nISRJets = self.makeTH2F("hist2d_nFSRJets_nISRJets",
							"FSRvsISR;FSRJets;ISRJets",
							7,0,7,7,0,7)


	hist_jetPtMaxDPhi = self.makeTH1F("hist_jetPtMaxDPhi",
			"jetPtMaxDPhi;jetPt;",100,0,3000)
	hist_jetPtMaxDPhi_2FSR = self.makeTH1F("hist_jetPtMaxDPhi_2FSR",
			"jetPtMaxDPhi_2FSR;jetPt;",100,0,3000)
	hist_jetPtMaxDPhi_3FSR = self.makeTH1F("hist_jetPtMaxDPhi_3FSR",
			"jetPtMaxDPhi_3FSR;jetPt;",100,0,3000)
	hist_jetPtMaxDPhi_1FSR = self.makeTH1F("hist_jetPtMaxDPhi_1FSR",
			"jetPtMaxDPhi_1FSR;jetPt;",100,0,3000)

	hist_MET = self.makeTH1F("hist_MET",
			"MET;MET;",100,0,1500)
	hist_MET_2FSR = self.makeTH1F("hist_MET_2FSR",
			"MET_2FSR;MET;",100,0,1500)
	hist_MET_3FSR = self.makeTH1F("hist_MET_3FSR",
			"MET_3FSR;MET;",100,0,1500)
	hist_MET_1FSR = self.makeTH1F("hist_MET_1FSR",
			"MET_1FSR;MET;",100,0,1500)

	hist_HT = self.makeTH1F("hist_HT",
			"HT;HT;",100,0,5000)
	hist_HT_2FSR = self.makeTH1F("hist_HT_2FSR",
			"hT_2FSR;HT;",100,0,5000)
	hist_HT_3FSR = self.makeTH1F("hist_HT_3FSR",
			"hT_3FSR;HT;",100,0,5000)
	hist_HT_1FSR = self.makeTH1F("hist_HT_1FSR",
			"HT_1FSR;HT;",100,0,5000)

	hist2d_2params = self.makeTH2F("hist_2params",";MET;HT",
		100,0,1500,100,0,5000)
	hist2d_2params_1FSR = self.makeTH2F("hist_2params_1FSR",";MET;HT",
		100,0,1500,100,0,5000)
	hist2d_2params_2FSR = self.makeTH2F("hist_2params_2FSR",";MET;HT",
		100,0,1500,100,0,5000)
	hist2d_2params_3FSR = self.makeTH2F("hist_2params_3FSR",";MET;HT",
		100,0,1500,100,0,5000)

	hist_pt_jetmaxdphi = self.makeTH1F("hist_pt_jetmaxdphi",
		"JetMaxDPhi;Pt;",
		100,0,2000)
	hist_pt_jetmaxdphi_1FSR = self.makeTH1F("hist_pt_jetmaxdphi_1FSR",
		"JetMaxDPhi_1FSR;Pt;",
		100,0,2000)
	hist_pt_jetmaxdphi_2FSR = self.makeTH1F("hist_pt_jetmaxdphi_2FSR",
		"JetMaxDPhi_2FSR;Pt;",
		100,0,2000)
	hist_pt_jetmaxdphi_3FSR = self.makeTH1F("hist_pt_jetmaxdphi_3FSR",
		"JetMaxDPhi_3FSR;Pt;",
		100,0,2000)
	hist_pt_jetmindphi = self.makeTH1F("hist_pt_jetmindphi",
		"JetMinDPhi;Pt;",
		100,0,2000)
	hist_pt_jetmindphi_1FSR = self.makeTH1F("hist_pt_jetmindphi_1FSR",
		"JetMinDPhi_1FSR;Pt;",
		100,0,2000)
	hist_pt_jetmindphi_2FSR = self.makeTH1F("hist_pt_jetmindphi_2FSR",
		"JetMinDPhi_2FSR;Pt;",
		100,0,2000)
	hist_pt_jetmindphi_3FSR = self.makeTH1F("hist_pt_jetmindphi_3FSR",
		"JetMinDPhi_3FSR;Pt;",
		100,0,2000)
	hist_pt_jetmiddphi = self.makeTH1F("hist_pt_jetmiddphi",
		"JetMidDPhi;Pt;",
		100,0,2000)
	hist_pt_jetmiddphi_1FSR = self.makeTH1F("hist_pt_jetmiddphi_1FSR",
		"JetMidDPhi_1FSR;Pt;",
		100,0,2000)
	hist_pt_jetmiddphi_2FSR = self.makeTH1F("hist_pt_jetmiddphi_2FSR",
		"JetMidDPhi_2FSR;Pt;",
		100,0,2000)
	hist_pt_jetmiddphi_3FSR = self.makeTH1F("hist_pt_jetmiddphi_3FSR",
		"JetMidDPhi_3FSR;Pt;",
		100,0,2000)


	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		nJets = len(tree.JetsAK8)

		nFSRJets = 0
		nISRJets = 0
		for value in tree.JetsAK8_isISR:
			if not value:
				nFSRJets += 1
			else:
				nISRJets += 1

		histDict_nJets["all"].Fill(nJets)
		histDict_nISRJets["all"].Fill(nISRJets)
		histDict_nFSRJets["all"].Fill(nFSRJets)

		
		histDict_nJets[str(nJets)].Fill(nJets)
		histDict_nISRJets[str(nJets)].Fill(nISRJets)
		histDict_nFSRJets[str(nJets)].Fill(nFSRJets)

		if nJets >= 3:
			histDict_nJets["3+"].Fill(nJets)
			histDict_nISRJets["3+"].Fill(nISRJets)
			histDict_nFSRJets["3+"].Fill(nFSRJets)

		hist2d_nJets_nFSRJets.Fill(nJets,nFSRJets)
		hist2d_nJets_nISRJets.Fill(nJets,nISRJets)
		hist2d_nFSRJets_nISRJets.Fill(nFSRJets,nISRJets)

		if nJets == 3:
			deltaPhiMetList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
			jetMaxDPhiIndex = deltaPhiMetList.index(max(deltaPhiMetList))
			jetPtMaxDPhi = tree.JetsAK8[jetMaxDPhiIndex].Pt()
			jetMinDPhiIndex = deltaPhiMetList.index(min(deltaPhiMetList))
			jetPtMinDPhi = tree.JetsAK8[jetMinDPhiIndex].Pt()
			jetMidDPhiIndex = -(jetMaxDPhiIndex+jetMinDPhiIndex)+3
			jetPtMidDPhi = tree.JetsAK8[jetMidDPhiIndex].Pt()
			hist_jetPtMaxDPhi.Fill(jetPtMaxDPhi)
			hist_MET.Fill(tree.MET)
			hist_HT.Fill(tree.HT)
			hist2d_2params.Fill(tree.MET,tree.HT)
			hist_pt_jetmaxdphi.Fill(jetPtMaxDPhi)
			hist_pt_jetmindphi.Fill(jetPtMinDPhi)
			hist_pt_jetmiddphi.Fill(jetPtMidDPhi)
			#other possible variables to use?
			# event-level variables:
				# MET, HT, 
			if nFSRJets == 2:
				hist_jetPtMaxDPhi_2FSR.Fill(jetPtMaxDPhi)
				hist_MET_2FSR.Fill(tree.MET)
				hist_HT_2FSR.Fill(tree.HT)
				hist2d_2params_2FSR.Fill(tree.MET,tree.HT)
				hist_pt_jetmaxdphi_2FSR.Fill(jetPtMaxDPhi)
				hist_pt_jetmindphi_2FSR.Fill(jetPtMinDPhi)
				hist_pt_jetmiddphi_2FSR.Fill(jetPtMidDPhi)
			elif nFSRJets == 3:
				hist_jetPtMaxDPhi_3FSR.Fill(jetPtMaxDPhi)
				hist_MET_3FSR.Fill(tree.MET)
				hist_HT_3FSR.Fill(tree.HT)
				hist2d_2params_3FSR.Fill(tree.MET,tree.HT)
				hist_pt_jetmaxdphi_3FSR.Fill(jetPtMaxDPhi)
				hist_pt_jetmindphi_3FSR.Fill(jetPtMinDPhi)
				hist_pt_jetmiddphi_3FSR.Fill(jetPtMidDPhi)
			elif nFSRJets == 1:
				hist_jetPtMaxDPhi_1FSR.Fill(jetPtMaxDPhi)
				hist_MET_1FSR.Fill(tree.MET)
				hist_HT_1FSR.Fill(tree.HT)
				hist2d_2params_1FSR.Fill(tree.MET,tree.HT)
				hist_pt_jetmaxdphi_1FSR.Fill(jetPtMaxDPhi)
				hist_pt_jetmindphi_1FSR.Fill(jetPtMinDPhi)
				hist_pt_jetmiddphi_1FSR.Fill(jetPtMidDPhi)

	makePlots([hist_jetPtMaxDPhi,
				hist_jetPtMaxDPhi_1FSR,
				hist_jetPtMaxDPhi_2FSR,
				hist_jetPtMaxDPhi_3FSR],
				self.extraDir,"jetPtMaxDPhi_3JetsTotal",log = True)
	makePlots([hist_MET,
				hist_MET_1FSR,
				hist_MET_2FSR,
				hist_MET_3FSR],
				self.extraDir,"MET_3JetsTotal",log = True)
	makePlots([hist_HT,
				hist_HT_1FSR,
				hist_HT_2FSR,
				hist_HT_3FSR],
				self.extraDir,"HT_3JetsTotal",log = True)

	makePlots([hist_pt_jetmaxdphi,
				hist_pt_jetmaxdphi_1FSR,
				hist_pt_jetmaxdphi_2FSR,
				hist_pt_jetmaxdphi_3FSR],
			self.extraDir,"jetptmaxdphi_only3jets",log = True)
	makePlots([hist_pt_jetmindphi,
				hist_pt_jetmindphi_1FSR,
				hist_pt_jetmindphi_2FSR,
				hist_pt_jetmindphi_3FSR],
			self.extraDir,"jetptmindphi_only3jets",log = True)
	makePlots([hist_pt_jetmiddphi,
				hist_pt_jetmiddphi_1FSR,
				hist_pt_jetmiddphi_2FSR,
				hist_pt_jetmiddphi_3FSR],
			self.extraDir,"jetptmiddphi_only3jets",log = True)

	make2DPlot(hist2d_2params,self.extraDir,"METvsHT")
	make2DPlot(hist2d_2params_1FSR,self.extraDir,"METvsHT_1FSR")
	make2DPlot(hist2d_2params_2FSR,self.extraDir,"METvsHT_2FSR")
	make2DPlot(hist2d_2params_3FSR,self.extraDir,"METvsHT_3FSR")


	print("| Key | nJets | nFSRJets | nISRJets |")
	for key in list(histDict_nJets):

		scale = histDict_nJets[key].GetEntries()
		if scale == 0:
			continue
		histDict_nJets[key].Scale(1./scale)
		histDict_nFSRJets[key].Scale(1./scale)
		histDict_nISRJets[key].Scale(1./scale)
	
		makePlots([
			histDict_nJets[key],
			histDict_nFSRJets[key],
			histDict_nISRJets[key]
			], self.extraDir,"nJets_"+key, log = False, lC = [0.7,0.7,0.9,0.9])

		avg_nJets = histDict_nJets[key].GetMean()
		avg_nFSRJets = histDict_nFSRJets[key].GetMean()
		avg_nISRJets = histDict_nISRJets[key].GetMean()

		print("| {} | {} | {} | {} |".format(key, avg_nJets, avg_nFSRJets, avg_nISRJets))

	makePlots([
		histDict_nFSRJets["all"]
		], self.extraDir,"nFSRJets", log = True, lC = [0.7,0.7,0.9,0.9])
	makePlots([
		histDict_nISRJets["all"]
		], self.extraDir,"nISRJets", log = True, lC = [0.7,0.7,0.9,0.9])
		
	make2DPlot(hist2d_nJets_nFSRJets,self.extraDir,"jetsVSfsr",logZ=True)
	make2DPlot(hist2d_nJets_nISRJets,self.extraDir,"jetsVSisr",logZ=True)
	make2DPlot(hist2d_nFSRJets_nISRJets,self.extraDir,"fsrVSisr",logZ=True)


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


