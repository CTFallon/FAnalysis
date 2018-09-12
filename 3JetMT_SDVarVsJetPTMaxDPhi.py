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
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	tree.SetBranchStatus("Electrons",1)
	tree.SetBranchStatus("Muons",1)
	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	#tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
	#tree.SetBranchStatus("numHVPartsInJet",1)
	#tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	# check overlap between jetPt(maxdPhi) and SDVar
	hist_overlap = self.makeTH2F("hist_overlap","Overlap;123 from jetPt(maxdPhi);123 from SDVar_13",2,0,2,2,0,2)
	
	# MT distributions

	cutPt = [x*10. for x in range(0,300)]# current optimal is 860 -> 861
	pTOptCut = 860
	#first, only vary the cut on one jet at a time
	histList_jetPtMaxDPhicut = []
	for cutVal in cutPt:
		histList_jetPtMaxDPhicut.append(self.makeTH1F("hist_MT_jetPtMaxdPhi_"+str(cutVal),"jetPtMaxdPhiCut;MT;count/a.u.",100,0,4000))

	cutSD = [x*0.01 for x in range(0,50)] # current optimal is .22 -> 0.242
	sdOptCut = 0.22
	#first, only vary the cut on one jet at a time
	histList_SDcut = []
	for cutVal in cutSD:
		histList_SDcut.append(self.makeTH1F("hist_MT_SD_"+str(cutVal),"SD;MT;count/a.u.",100,0,4000))

	hist_MT_optimalSDVar = self.makeTH1F("hist_MT_optimalSDVar","SDVar;MT;count/a.u.", 100, 0, 4000)
	hist_MT_optimalJetPT = self.makeTH1F("hist_MT_optimalJetPT","JetPT;MT;count/a.u.", 100, 0, 4000)
	hist_MT_optimalPTFrac = self.makeTH1F("hist_MT_optimalPTFrac","PTfrac;MT;count/a.u.", 100, 0, 4000)
	hist_MT_base = self.makeTH1F("hist_MT_base","base;MT;count/a.u.", 100, 0, 4000)
	hist_MT_just2jets = self.makeTH1F("hist_MT_just2jets","just2jets;MT;count/a.u.", 100, 0, 4000)
	
	count_shouldwas_sdVar = 0
	count_shouldwasnt_sdVar = 0
	count_shouldntwas_sdVar = 0
	count_shouldntwasnt_sdVar = 0

	count_shouldwas_jetPt = 0
	count_shouldwasnt_jetPt = 0
	count_shouldntwas_jetPt = 0
	count_shouldntwasnt_jetPt = 0
	
	count_should_truth = 0
	count_shouldnt_truth = 0
	
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print(str(iEvent) + "/"+ str(nEvents))
		tree.GetEvent(iEvent)
		nJets = len(tree.JetsAK8)
		jets = tree.JetsAK8
		met = tree.MET
		metPhi = tree.METPhi
		if not tree.passedPreSelection:
			continue
		#Optimize cuts for SDVar13 and jetPT(maxDPhi), only using RECO level information

		for iCut in range(len(cutPt)):
			jetsForMt = []
			cutVal = cutPt[iCut]
			jetsForMt.append(jets[0])
			jetsForMt.append(jets[1])
			if nJets != 2 and jets[tree.iJetMaxDeltaPhi].Pt() < cutVal:
				jetsForMt.append(tree.JetsAK8[2])
			histList_jetPtMaxDPhicut[iCut].Fill(trans_mass_Njet(jetsForMt, met, metPhi))


		for iCut in range(len(cutSD)):
			jetsForMt = []
			cutVal = cutSD[iCut]
			jetsForMt.append(jets[0])
			jetsForMt.append(jets[1])
			if nJets != 2 and jets[2].Pt()/(jets[0].Pt()+jets[2].Pt()) > cutVal:
				jetsForMt.append(tree.JetsAK8[2])
			histList_SDcut[iCut].Fill(trans_mass_Njet(jetsForMt, met, metPhi))

		# check optimal cuts for comparison to 'truth', and overlap
		if nJets != 2:
			wasSDVar = 0
			wasjetPt = 0
			if tree.fracPtFromHVQuarks[2] > 0.2:
				count_should_truth += 1
				if jets[2].Pt()/(jets[0].Pt()+jets[2].Pt()) > sdOptCut:
					count_shouldwas_sdVar += 1
					wasSDVar = 1
				else:
					count_shouldwasnt_sdVar += 1
				if jets[tree.iJetMaxDeltaPhi].Pt() < pTOptCut:
					count_shouldwas_jetPt += 1
					wasjetPt = 1
				else:
					count_shouldwasnt_jetPt += 1
			else:
				count_shouldnt_truth += 1
				if jets[2].Pt()/(jets[0].Pt()+jets[2].Pt()) > sdOptCut:
					count_shouldntwas_sdVar += 1
					wasSDVar = 1
				else:
					count_shouldntwasnt_sdVar += 1
				if jets[tree.iJetMaxDeltaPhi].Pt() < pTOptCut:
					count_shouldntwas_jetPt += 1
					wasjetPt = 1
				else:
					count_shouldntwasnt_jetPt += 1
			hist_overlap.Fill(wasjetPt, wasSDVar)

		if nJets != 2:
			if jets[2].Pt()/(jets[0].Pt()+jets[2].Pt()) > sdOptCut:
				hist_MT_optimalSDVar.Fill(trans_mass_Njet(jets[0:3], met, metPhi))
			else:
				hist_MT_optimalSDVar.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
			if jets[tree.iJetMaxDeltaPhi].Pt() < pTOptCut:
				hist_MT_optimalJetPT.Fill(trans_mass_Njet(jets[0:3], met, metPhi))
			else:
				hist_MT_optimalJetPT.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
			if tree.fracPtFromHVQuarks[2] < 0.2:
				hist_MT_optimalPTFrac.Fill(trans_mass_Njet(jets[0:3], met, metPhi))
			else:
				hist_MT_optimalPTFrac.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
		else:
			hist_MT_just2jets.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
			hist_MT_optimalSDVar.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
			hist_MT_optimalJetPT.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
			hist_MT_optimalPTFrac.Fill(trans_mass_Njet(jets[0:2], met, metPhi))
		hist_MT_base.Fill(trans_mass_Njet(jets[0:2], met, metPhi))

	print("Truth: should = " + str(count_should_truth))
	print("Truth: shouldnt = " + str(count_shouldnt_truth))

	print("sdVar: should, was = " + str(count_shouldwas_sdVar))
	print("sdVar: should, wasnt = " + str(count_shouldwasnt_sdVar))
	print("sdVar: shouldnt, was = " + str(count_shouldntwas_sdVar))
	print("sdVar: shouldnt, wasnt = " + str(count_shouldntwasnt_sdVar))

	print("jetPt: should, was = " + str(count_shouldwas_jetPt))
	print("jetPt: should, wasnt = " + str(count_shouldwasnt_jetPt))
	print("jetPt: shouldnt, was = " + str(count_shouldntwas_jetPt))
	print("jetPt: shouldnt, wasnt = " + str(count_shouldntwasnt_jetPt))
	
	for histo in histList_jetPtMaxDPhicut:
		try:
			print("Cut at " + histo.GetName() + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
		except ZeroDivisionError:
			print("Cut at " + histo.GetName() + " Resolution is NULL")
	for histo in histList_SDcut:
		try:
			print("Cut at " + histo.GetName() + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
		except ZeroDivisionError:
			print("Cut at " + histo.GetName() + " Resolution is NULL")

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


