from analysisBase import baseClass
import ROOT as rt
from array import array

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	ff = rt.TFile.Open("outFriend.root")
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

	cutPt = [x*10. for x in range(0,300)]
	#first, only vary the cut on one jet at a time
	histList_jetPtMaxDPhicut = []
	for cutVal in cutPt:
		histList_jetPtMaxDPhicut.append(self.makeTH1F("hist_MT_jetPtMaxdPhi_"+str(cutVal),"jetPtMaxdPhiCut;MT;count/a.u.",100,0,4000))

	cutSD = [x*0.01 for x in range(0,50)]
	#first, only vary the cut on one jet at a time
	histList_SDcut = []
	for cutVal in cutSD:
		histList_SDcut.append(self.makeTH1F("hist_MT_SD_"+str(cutVal),"SD;MT;count/a.u.",100,0,4000))
	
	
	
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
		#Optimize cuts for SDVar13 and jetPT(maxDPhi)

		for iCut in range(len(cutPt)):
			jetsForMt = []
			cutVal = cutPt[iCut]
			if nJets == 2:
				jetsForMt.append(jets[0])
				jetsForMt.append(jets[1])
			else:
				if tree.fracPtFromHVQuarks[0] > 0.0:
					jetsForMt.append(tree.JetsAK8[0])
				if tree.fracPtFromHVQuarks[1] > 0.0:
					jetsForMt.append(tree.JetsAK8[1])
				if tree.fracPtFromHVQuarks[2] > 0.02 and jets[tree.iJetMaxDeltaPhi].Pt() > cutVal:
					jetsForMt.append(tree.JetsAK8[2])
			histList_jetPtMaxDPhicut[iCut].Fill(trans_mass_Njet(jetsForMt, met, metPhi))


		for iCut in range(len(cutSD)):
			jetsForMt = []
			cutVal = cutSD[iCut]
			if nJets == 2:
				jetsForMt.append(jets[0])
				jetsForMt.append(jets[1])
			else:
				if tree.fracPtFromHVQuarks[0] > 0.0:
					jetsForMt.append(tree.JetsAK8[0])
				if tree.fracPtFromHVQuarks[1] > 0.0:
					jetsForMt.append(tree.JetsAK8[1])
				if tree.fracPtFromHVQuarks[2] > 0.02 and jets[2].Pt()/(jets[0].Pt()+jets[2].Pt()) > cutVal:
					jetsForMt.append(tree.JetsAK8[2])
			histList_SDcut[iCut].Fill(trans_mass_Njet(jetsForMt, met, metPhi))
	
	for histo in histList_jetPtMaxDPhicut:
		try:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
		except ZeroDivisionError:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is NULL")
	for histo in histList_SDcut:
		try:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
		except ZeroDivisionError:
			print("Cut at " + histo.GetName()[-3:] + " Resolution is NULL")

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


