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

	nEventsTot = 0
	nEventsPassPre = 0
	nEventsPP3Jets = 0

	nDiT = 0
	nTriT = 0

	nDiA = 0
	nTriA = 0
	nDiB = 0
	nTriB = 0
	nDiC = 0
	nTriC = 0

	nDiR = 0
	nTriR = 0


	hist_MT_AllEvents_dijet  = self.makeTH1F("hist_MT_AllEvents_dijet" ,"MT(12);MT; a.u."  ,100,1500,4000)
	hist_MT_AllEvents_trijet = self.makeTH1F("hist_MT_AllEvents_trijet","MT(12[3]);MT; a.u.",100,1500,4000)
	hist_MT_AllEvents_algo   = self.makeTH1F("hist_MT_AllEvents_algo"  ,"MT(algo);MT; a.u." ,100,1500,4000)
	hist_MT_AllEvents_perf   = self.makeTH1F("hist_MT_AllEvents_perf"  ,"MT(perf);MT; a.u." ,100,1500,4000)

	hist_A_dijet =  self.makeTH1F("hist_A_dijet","gamma3 - Dijet;gamma3;Count/a.u.",200,0,20)
	hist_A_trijet =  self.makeTH1F("hist_A_trijet","gamma3 - Trijet;gamma3;Count/a.u.",200,0,20)

	hist_B_dijet =  self.makeTH1F("hist_B_dijet","iJetMaxDeltaPhi - Dijet;iJetMaxDeltaPhi;Count/a.u.",5,0,5)
	hist_B_trijet =  self.makeTH1F("hist_B_trijet","iJetMaxDeltaPhi - Trijet;iJetMaxDeltaPhi;Count/a.u.",5,0,5)
	
	testVarName = "d_gamma3"
	testVarLow = 0
	testVarHig = 20
	nBins = 500

	hist_test_dijet = self.makeTH1F("hist_test_dijet",testVarName+" - Dijet MT;"+testVarName+";Count/a.u.", nBins, testVarLow, testVarHig)
	hist_test_trijet = self.makeTH1F("hist_test_trijet",testVarName+" - Trijet MT;"+testVarName+";Count/a.u.", nBins, testVarLow, testVarHig)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		nEventsTot += 1
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		nEventsPassPre += 1
		hist_MT_AllEvents_dijet.Fill(tree.MT_AK8)
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		if len(tree.JetsAK8) <= 2:
			hist_MT_AllEvents_trijet.Fill(tree.MT_AK8)
			hist_MT_AllEvents_algo.Fill(tree.MT_AK8)
			hist_MT_AllEvents_perf.Fill(tree.MT_AK8)
			continue
		nEventsPP3Jets += 1
		MT3 = trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)
		hist_MT_AllEvents_trijet.Fill(MT3)

		if bool(tree.JetsAK8_isHV[2]) == True:
			isTri = True
			nTriT +=1
			hist_MT_AllEvents_perf.Fill(MT3)
		else:
			nDiT += 1
			isTri = False
			hist_MT_AllEvents_perf.Fill(tree.MT_AK8)


		if isTri:
			hist_A_trijet.Fill(tree.JetsAK8[2].Gamma())
		else:
			hist_A_dijet.Fill(tree.JetsAK8[2].Gamma())
		#if tree.JetsAK8[2].Gamma() < 5.2:
		#	useTri = True
		#	if isTri:
		#		nTriA += 1
		#		hist_B_trijet.Fill(tree.iJetMaxDeltaPhi)
		#	else:
		#		nDiA += 1
		#		hist_B_dijet.Fill(tree.iJetMaxDeltaPhi)
		#elif tree.iJetMaxDeltaPhi == 1:
		#	useTri = True
		#	if isTri:
		#		nTriB += 1
		#	else:
		#		nDiB += 1
		#else:
		useTri = False
		#Test new cuts here
		testVar = tree.JetsAK8[2].Gamma()
		if isTri:
			nTriR += 1
			hist_test_trijet.Fill(testVar)
		else:
			nDiR += 1
			hist_test_dijet.Fill(testVar)

		if useTri:		
			hist_MT_AllEvents_algo.Fill(MT3)
		else:
			hist_MT_AllEvents_algo.Fill(tree.MT_AK8)

			

	print("-------------------------")
	print("Number of Events In MC: {}".format(nEventsTot))
	print("Number of Events Passing PreSelection: {}".format(nEventsPassPre))
	print("Number of Events Passing PreSelection with 3 or more AK8 Jets: {}".format(nEventsPP3Jets))
	print("Number of Events PP3+Jets with only 2 HV jets: {}".format(nDiT))
	print("Number of Events PP3+Jets with 3 HV jets: {}".format(nTriT))
	print("-------------------------")
	print(" Cut | nDi | nTri | Total")
	print(" Total | {} | {} | {}".format(nDiT, nTriT, nDiT+nTriT))
	print(" gamma3 | {} | {} | {}".format(nDiA, nTriA, nDiA+nTriA))
	print(" Jet Index | {} | {} | {}".format(nDiB, nTriB, nDiB+nTriB))
	print(" Remaining | {} | {} | {}".format(nDiR, nTriR, nDiR+nTriR))

	# add the MT Resolution to each MT distribution's title
	hist_MT_AllEvents_dijet.SetTitle(hist_MT_AllEvents_dijet.GetTitle()+" {:.4f}".format(hist_MT_AllEvents_dijet.GetRMS()/hist_MT_AllEvents_dijet.GetMean()))
	hist_MT_AllEvents_trijet.SetTitle(hist_MT_AllEvents_trijet.GetTitle()+" {:.4f}".format(hist_MT_AllEvents_trijet.GetRMS()/hist_MT_AllEvents_trijet.GetMean()))
	hist_MT_AllEvents_algo.SetTitle(hist_MT_AllEvents_algo.GetTitle()+" {:.4f}".format(hist_MT_AllEvents_algo.GetRMS()/hist_MT_AllEvents_algo.GetMean()))
	hist_MT_AllEvents_perf.SetTitle(hist_MT_AllEvents_perf.GetTitle()+" {:.4f}".format(hist_MT_AllEvents_perf.GetRMS()/hist_MT_AllEvents_perf.GetMean()))

	self.makeRatio(hist_MT_AllEvents_perf, hist_MT_AllEvents_dijet,"MT_dijet", doLeg = True, log = False)
	self.makeRatio(hist_MT_AllEvents_perf, hist_MT_AllEvents_trijet,"MT_trijet", doLeg = True, log = False)
	self.makeRatio(hist_MT_AllEvents_perf, hist_MT_AllEvents_algo,"MT_algo", doLeg = True, log = False)


	self.makePng([hist_MT_AllEvents_dijet,hist_MT_AllEvents_trijet,hist_MT_AllEvents_algo,hist_MT_AllEvents_perf],"MT_dtap")
	self.makePng([hist_A_dijet,hist_A_trijet],"gamma3", log= False, doCum = True)
	self.makePng([hist_A_dijet.GetCumulative()-hist_A_trijet.GetCumulative()],"gamma3_diff", log= False, doLeg = False)

	self.makePng([hist_B_dijet,hist_B_trijet],"iJetMaxDeltaPhi", log= False, doCum = True)
	self.makePng([hist_B_dijet.GetCumulative()-hist_B_trijet.GetCumulative()],"iJetMaxDeltaPhi_diff", log= False, doLeg = False)

	
	self.makePng([hist_test_dijet,hist_test_trijet],testVarName, log= False, doCum = True)
	hist_CumDiff = hist_test_dijet.GetCumulative()-hist_test_trijet.GetCumulative()
	maximum = hist_CumDiff.GetMaximum()
	maxBin = hist_CumDiff.GetMaximumBin()
	minimum = hist_CumDiff.GetMinimum()
	minBin = hist_CumDiff.GetMinimumBin()
	print("Maximum Value is {} at bin number {}, corresponding to cut value of {}".format(maximum, maxBin, testVarLow+maxBin*(testVarHig-testVarLow)/float(nBins)))
	print("Minimum Value is {} at bin number {}, corresponding to cut value of {}".format(minimum, minBin, testVarLow+minBin*(testVarHig-testVarLow)/float(nBins)))
	self.makePng([hist_CumDiff],testVarName+"_diff", log= False, doLeg = False)
	
	
	print("\\begin{frame}{Impact of the Problem, Sample Code "+self.fileID+"}")
	print("\\begin{center}")

	print("\\begin{tabular}{c|c|c|c}")
	print("Event Selection & nEvents  & \\% of PreSelec & \\% of 3+ AK8\\\\ \\hline")
	print("All MC & {} &  - & -\\\\".format(nEventsTot))
	print("PreSelection & {} & {:.2f}\\% & -\\\\ ".format(nEventsPassPre,100*float(nEventsPassPre)/nEventsPassPre))
	print("3+ AK8 Jets & {} & {:.2f}\\% & {:.2f}\\% \\\\ \\hline".format(nEventsPP3Jets,100*float(nEventsPP3Jets)/nEventsPassPre,100*float(nEventsPP3Jets)/nEventsPP3Jets))
	print("3+ AK8 Jets \\& Dijet & {} & {:.2f}\\% & {:.2f}\\%\\\\ ".format(nDiT,100*float(nDiT)/nEventsPassPre,100*float(nDiT)/nEventsPP3Jets))
	print("3+ AK8 Jets \\& Trijet & {} & {:.2f}\\% & {:.2f}\\%\\\\".format(nTriT,100*float(nTriT)/nEventsPassPre,100*float(nTriT)/nEventsPP3Jets))
	print("\end{tabular}")
	print("\end{center}")
	print("\end{frame}")
	
	print("ECPN: min minCut max maxCut nMC nPS n3+ nDi nTri nDiG3 nTriG3")
	print("ECPN: {} {} {} {} {} {} {} {} {} {} {}".format(minimum,testVarLow+minBin*(testVarHig-testVarLow)/float(nBins),
										maximum,testVarLow+maxBin*(testVarHig-testVarLow)/float(nBins),
										nEventsTot, nEventsPassPre,nEventsPP3Jets,nDiT,nTriT, nDiR, nTriR))


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
	delta = phi1 - phi2
	while delta > rt.TMath.Pi():
		delta -= 2*rt.TMath.Pi()
	while delta < -rt.TMath.Pi():
		delta += 2*rt.TMath.Pi()
	return abs(delta)

	

