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

	#coarse grading for optimal MT resolution for fracPt cut
	cutFractions = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	#first, only vary the cut on one jet at a time
	hist_MTLead2 = self.makeTH1F("hist_MTlead2Jets","Base MT;MT;count/a.u.",100,0,4000)
	histList_MTcut = []
	for cutVal in cutFractions:
		histList_MTcut.append(self.makeTH1F("hist_MT2jet_"+str(cutVal),"2jet MT;MT;count/a.u.",100,0,4000))
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		met = tree.MET
		metPhi = tree.METPhi
		if tree.passedPreSelection == 1:
			nJets = len(tree.JetsAK8)
			for iJet in range(nJets):
				histList_2d_iJetvsFracPt[nJets].Fill(iJet+0.5, tree.fracPtFromHVQuarks[iJet])
			hist_MTLead2.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], met, metPhi))
			for iCut in range(len(cutFractions)):
				cutVal = cutFractions[iCut]
				if nJets == 2:
					histList_MTcut[iCut].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], met, metPhi))
				elif tree.fracPtFromHVQuarks[2] > cutVal:
					histList_MTcut[iCut].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], met, metPhi))
				else:
					histList_MTcut[iCut].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], met, metPhi))
	
	print("No cut has Resolution " + str(hist_MTLead2.GetRMS()/hist_MTLead2.GetMean()))
	for histo in histList_MTcut:
		print("Cut at " + histo.GetName[-3:] + " Resolution is " + str(histo.GetRMS()/histo.GetMean()))
					


	

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


