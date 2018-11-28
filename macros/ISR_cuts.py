from analysisBase import baseClass
import ROOT as rt
from array import array

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
	#tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	#tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	tree.SetBranchStatus("fracPtFromHVQuarks",1)
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

	#JetsAK8_* : axismajor, axisminor, doubleBDiscriminator, girth, HVfracPt, HVminDeltaR, isISR, ISRfracPt, ISRminDeltaR, momenthalf, multiplicity, nHVISRparts, nSMISRparts, NsubjettinessTau1, NsubjettinessTau2, NsubjettinessTau3, NumBhadrons, NumChadrons, overflow, prunedMass, ptD, softDropMass, kinematicVariables

	#ptD is a good first cut, greater than 0.9 is ISR

	hist_phi = self.makeTH1F("hist_phi","phi;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())
	hist_phi_ISR = self.makeTH1F("hist_phi_ISR","phi;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())
	hist_phi_FSR = self.makeTH1F("hist_phi_FSR","phi;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())

	hist_eta = self.makeTH1F("hist_eta","eta;eta",100,-5,5)
	hist_eta_ISR = self.makeTH1F("hist_eta_ISR","eta;eta",100,-5,5)
	hist_eta_FSR = self.makeTH1F("hist_eta_FSR","eta;eta",100,-5,5)

	hist_pt = self.makeTH1F("hist_pt","pt;pt",100,0,2500)
	hist_pt_ISR = self.makeTH1F("hist_pt_ISR","pt;pt",100,0,2500)
	hist_pt_FSR = self.makeTH1F("hist_pt_FSR","pt;pt",100,0,2500)

	hist_px = self.makeTH1F("hist_px","px;px",100,-2000,2000)
	hist_px_ISR = self.makeTH1F("hist_px_ISR","px;px",100,-2000,2000)
	hist_px_FSR = self.makeTH1F("hist_px_FSR","px;px",100,-2000,2000)

	hist_py = self.makeTH1F("hist_py","py;py",100,-2000,2000)
	hist_py_ISR = self.makeTH1F("hist_py_ISR","py;py",100,-2000,2000)
	hist_py_FSR = self.makeTH1F("hist_py_FSR","py;py",100,-2000,2000)

	hist_pz = self.makeTH1F("hist_pz","pz;pz",100,-6000,6000)
	hist_pz_ISR = self.makeTH1F("hist_pz_ISR","pz;pz",100,-6000,6000)
	hist_pz_FSR = self.makeTH1F("hist_pz_FSR","pz;pz",100,-6000,6000)

	hist_axismajor = self.makeTH1F("hist_axismajor","axisMajor;axismajor;",100,0,0.5)
	hist_axismajor_ISR = self.makeTH1F("hist_axismajor_ISR","axisMajor;axismajor;",100,0,0.5)
	hist_axismajor_FSR = self.makeTH1F("hist_axismajor_FSR","axisMajor;axismajor;",100,0,0.5)

	hist_axisminor = self.makeTH1F("hist_axisminor","axisMinor;axisminor;",100,0,0.3)
	hist_axisminor_ISR = self.makeTH1F("hist_axisminor_ISR","axisMinor;axisminor;",100,0,0.3)
	hist_axisminor_FSR = self.makeTH1F("hist_axisminor_FSR","axisMinor;axisminor;",100,0,0.3)

	hist_doubleBDiscriminator = self.makeTH1F("hist_doubleBDiscriminator","doubleBDiscriminator;doubleBDiscriminator;",100,-1.1,1.1)
	hist_doubleBDiscriminator_ISR = self.makeTH1F("hist_doubleBDiscriminator_ISR","doubleBDiscriminator;doubleBDiscriminator;",100,-1.1,1.1)
	hist_doubleBDiscriminator_FSR = self.makeTH1F("hist_doubleBDiscriminator_FSR","doubleBDiscriminator;doubleBDiscriminator;",100,-1.1,1.1)

	hist_girth = self.makeTH1F("hist_girth","girth;girth;",100,0,0.6)
	hist_girth_ISR = self.makeTH1F("hist_girth_ISR","girth;girth;",100,0,0.6)
	hist_girth_FSR = self.makeTH1F("hist_girth_FSR","girth;girth;",100,0,0.6)

	hist_HVfracPt = self.makeTH1F("hist_HVfracPt","HVfracPt;HVfracPt;",100,-0.01,1.01)
	hist_HVfracPt_ISR = self.makeTH1F("hist_HVfracPt_ISR","HVfracPt;HVfracPt;",100,-0.01,1.01)
	hist_HVfracPt_FSR = self.makeTH1F("hist_HVfracPt_FSR","HVfracPt;HVfracPt;",100,-0.01,1.01)

	hist_HVminDeltaR = self.makeTH1F("hist_HVminDeltaR","HVminDeltaR;HVminDeltaR;",100,0,15.01)
	hist_HVminDeltaR_ISR = self.makeTH1F("hist_HVminDeltaR_ISR","HVminDeltaR;HVminDeltaR;",100,0,15.01)
	hist_HVminDeltaR_FSR = self.makeTH1F("hist_HVminDeltaR_FSR","HVminDeltaR;HVminDeltaR;",100,0,15.01)

	hist_ISRfracPt = self.makeTH1F("hist_ISRfracPt","ISRfracPt;ISRfracPt;",100,-0.01,1.01)
	hist_ISRfracPt_ISR = self.makeTH1F("hist_ISRfracPt_ISR","ISRfracPt;ISRfracPt;",100,-0.01,1.01)
	hist_ISRfracPt_FSR = self.makeTH1F("hist_ISRfracPt_FSR","ISRfracPt;ISRfracPt;",100,-0.01,1.01)

	hist_ISRminDeltaR = self.makeTH1F("hist_ISRminDeltaR","ISRminDeltaR;ISRminDeltaR;",100,0,15.01)
	hist_ISRminDeltaR_ISR = self.makeTH1F("hist_ISRminDeltaR_ISR","ISRminDeltaR;ISRminDeltaR;",100,0,15.01)
	hist_ISRminDeltaR_FSR = self.makeTH1F("hist_ISRminDeltaR_FSR","ISRminDeltaR;ISRminDeltaR;",100,0,15.01)

	hist_momenthalf = self.makeTH1F("hist_momenthalf","momenthalf;momenthalf;",100,0,0.8)
	hist_momenthalf_ISR = self.makeTH1F("hist_momenthalf_ISR","momenthalf;momenthalf;",100,0,0.8)
	hist_momenthalf_FSR = self.makeTH1F("hist_momenthalf_FSR","momenthalf;momenthalf;",100,0,0.8)

	hist_multiplicity = self.makeTH1F("hist_multiplicity","multiplicity;multiplicity;",400,0,400)
	hist_multiplicity_ISR = self.makeTH1F("hist_multiplicity_ISR","multiplicity;multiplicity;",400,0,400)
	hist_multiplicity_FSR = self.makeTH1F("hist_multiplicity_FSR","multiplicity;multiplicity;",400,0,400)

	hist_nHVISRparts = self.makeTH1F("hist_nHVISRparts","nHVISRparts;nHVISRparts;",15,0,15)
	hist_nHVISRparts_ISR = self.makeTH1F("hist_nHVISRparts_ISR","nHVISRparts;nHVISRparts;",15,0,15)
	hist_nHVISRparts_FSR = self.makeTH1F("hist_nHVISRparts_FSR","nHVISRparts;nHVISRparts;",15,0,15)

	hist_nSMISRparts = self.makeTH1F("hist_nSMISRparts","nSMISRparts;nSMISRparts;",35,0,35)
	hist_nSMISRparts_ISR = self.makeTH1F("hist_nSMISRparts_ISR","nSMISRparts;nSMISRparts;",35,0,35)
	hist_nSMISRparts_FSR = self.makeTH1F("hist_nSMISRparts_FSR","nSMISRparts;nSMISRparts;",35,0,35)

	hist_NsubjettinessTau1 = self.makeTH1F("hist_NsubjettinessTau1","NsubjettinessTau1;NsubjettinessTau1;",100,0,0.7)
	hist_NsubjettinessTau1_ISR = self.makeTH1F("hist_NsubjettinessTau1_ISR","NsubjettinessTau1;NsubjettinessTau1;",100,0,0.7)
	hist_NsubjettinessTau1_FSR = self.makeTH1F("hist_NsubjettinessTau1_FSR","NsubjettinessTau1;NsubjettinessTau1;",100,0,0.7)

	hist_NsubjettinessTau2 = self.makeTH1F("hist_NsubjettinessTau2","NsubjettinessTau2;NsubjettinessTau2;",100,0,0.5)
	hist_NsubjettinessTau2_ISR = self.makeTH1F("hist_NsubjettinessTau2_ISR","NsubjettinessTau2;NsubjettinessTau2;",100,0,0.5)
	hist_NsubjettinessTau2_FSR = self.makeTH1F("hist_NsubjettinessTau2_FSR","NsubjettinessTau2;NsubjettinessTau2;",100,0,0.5)

	hist_NsubjettinessTau3 = self.makeTH1F("hist_NsubjettinessTau3","NsubjettinessTau3;NsubjettinessTau3;",100,0,0.5)
	hist_NsubjettinessTau3_ISR = self.makeTH1F("hist_NsubjettinessTau3_ISR","NsubjettinessTau3;NsubjettinessTau3;",100,0,0.5)
	hist_NsubjettinessTau3_FSR = self.makeTH1F("hist_NsubjettinessTau3_FSR","NsubjettinessTau3;NsubjettinessTau3;",100,0,0.5)

	hist_Tau12 = self.makeTH1F("hist_Tau12","Tau12;Tau12;",100,0,10)
	hist_Tau12_ISR = self.makeTH1F("hist_Tau12_ISR","Tau12;Tau12;",100,0,10)
	hist_Tau12_FSR = self.makeTH1F("hist_Tau12_FSR","Tau12;Tau12;",100,0,10)

	hist_Tau23 = self.makeTH1F("hist_Tau23","Tau23;Tau23;",100,0,10)
	hist_Tau23_ISR = self.makeTH1F("hist_Tau23_ISR","Tau23;Tau23;",100,0,10)
	hist_Tau23_FSR = self.makeTH1F("hist_Tau23_FSR","Tau23;Tau23;",100,0,10)

	hist_NumBhadrons = self.makeTH1F("hist_NumBhadrons","NumBhadrons;NumBhadrons;",5,0,5)
	hist_NumBhadrons_ISR = self.makeTH1F("hist_NumBhadrons_ISR","NumBhadrons;NumBhadrons;",5,0,5)
	hist_NumBhadrons_FSR = self.makeTH1F("hist_NumBhadrons_FSR","NumBhadrons;NumBhadrons;",5,0,5)

	hist_NumChadrons = self.makeTH1F("hist_NumChadrons","NumChadrons;NumChadrons;",9,0,9)
	hist_NumChadrons_ISR = self.makeTH1F("hist_NumChadrons_ISR","NumChadrons;NumChadrons;",9,0,9)
	hist_NumChadrons_FSR = self.makeTH1F("hist_NumChadrons_FSR","NumChadrons;NumChadrons;",9,0,9)

	hist_overflow = self.makeTH1F("hist_overflow","overflow;overflow;",100,0,1)
	hist_overflow_ISR = self.makeTH1F("hist_overflow_ISR","overflow;overflow;",100,0,1)
	hist_overflow_FSR = self.makeTH1F("hist_overflow_FSR","overflow;overflow;",100,0,1)

	hist_prunedMass = self.makeTH1F("hist_prunedMass","prunedMass;prunedMass;",200,-10,800)
	hist_prunedMass_ISR = self.makeTH1F("hist_prunedMass_ISR","prunedMass;prunedMass;",200,-10,800)
	hist_prunedMass_FSR = self.makeTH1F("hist_prunedMass_FSR","prunedMass;prunedMass;",200,-10,800)

	hist_ptD = self.makeTH1F("hist_ptD","ptD;ptD;",100,0,1)
	hist_ptD_ISR = self.makeTH1F("hist_ptD_ISR","ptD;ptD;",100,0,1)
	hist_ptD_FSR = self.makeTH1F("hist_ptD_FSR","ptD;ptD;",100,0,1)

	hist_softDropMass = self.makeTH1F("hist_softDropMass","softDropMass;softDropMass;",200,-10,900)
	hist_softDropMass_ISR = self.makeTH1F("hist_softDropMass_ISR","softDropMass;softDropMass;",200,-10,900)
	hist_softDropMass_FSR = self.makeTH1F("hist_softDropMass_FSR","softDropMass;softDropMass;",200,-10,900)

	hist_jetNumber = self.makeTH1F("hist_jetNumber","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_ISR = self.makeTH1F("hist_jetNumber_ISR","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_FSR = self.makeTH1F("hist_jetNumber_FSR","jetNumber;jetNumber",10,0,10)

	hist_jetNumber_nj2 = self.makeTH1F("hist_jetNumber_nj2","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj2_ISR = self.makeTH1F("hist_jetNumber_nj2_ISR","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj2_FSR = self.makeTH1F("hist_jetNumber_nj2_FSR","jetNumber;jetNumber",10,0,10)

	hist_jetNumber_nj3 = self.makeTH1F("hist_jetNumber_nj3","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj3_ISR = self.makeTH1F("hist_jetNumber_nj3_ISR","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj3_FSR = self.makeTH1F("hist_jetNumber_nj3_FSR","jetNumber;jetNumber",10,0,10)

	hist_jetNumber_nj4 = self.makeTH1F("hist_jetNumber_nj4","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj4_ISR = self.makeTH1F("hist_jetNumber_nj4_ISR","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_nj4_FSR = self.makeTH1F("hist_jetNumber_nj4_FSR","jetNumber;jetNumber",10,0,10)
	

	



	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		nJets = len(tree.JetsAK8)
		for iJet in range(nJets):
			if(tree.JetsAK8_ptD[iJet] > 0.9):
				continue
			#Fill All Jet Histograms
			hist_phi.Fill(tree.JetsAK8[iJet].Phi())
			hist_eta.Fill(tree.JetsAK8[iJet].Eta())
			hist_pt.Fill(tree.JetsAK8[iJet].Pt())
			hist_px.Fill(tree.JetsAK8[iJet].Px())
			hist_py.Fill(tree.JetsAK8[iJet].Py())
			hist_pz.Fill(tree.JetsAK8[iJet].Pz())
			hist_axismajor.Fill(tree.JetsAK8_axismajor[iJet])
			hist_axisminor.Fill(tree.JetsAK8_axisminor[iJet])
			hist_doubleBDiscriminator.Fill(tree.JetsAK8_doubleBDiscriminator[iJet])
			hist_girth.Fill(tree.JetsAK8_girth[iJet])
			hist_HVfracPt.Fill(tree.JetsAK8_HVfracPt[iJet])
			hist_HVminDeltaR.Fill(tree.JetsAK8_HVminDeltaR[iJet])
			hist_ISRfracPt.Fill(tree.JetsAK8_ISRfracPt[iJet])
			hist_ISRminDeltaR.Fill(tree.JetsAK8_ISRminDeltaR[iJet])
			hist_momenthalf.Fill(tree.JetsAK8_momenthalf[iJet])
			hist_multiplicity.Fill(tree.JetsAK8_multiplicity[iJet])
			hist_nHVISRparts.Fill(tree.JetsAK8_nHVISRparts[iJet])
			hist_nSMISRparts.Fill(tree.JetsAK8_nSMISRparts[iJet])
			hist_NsubjettinessTau1.Fill(tree.JetsAK8_NsubjettinessTau1[iJet])
			hist_NsubjettinessTau2.Fill(tree.JetsAK8_NsubjettinessTau2[iJet])
			hist_NsubjettinessTau3.Fill(tree.JetsAK8_NsubjettinessTau3[iJet])
			hist_Tau12.Fill(tree.JetsAK8_NsubjettinessTau1[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
			hist_Tau23.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau3[iJet])
			hist_NumBhadrons.Fill(tree.JetsAK8_NumBhadrons[iJet])
			hist_NumChadrons.Fill(tree.JetsAK8_NumChadrons[iJet])
			hist_overflow.Fill(tree.JetsAK8_overflow[iJet])
			hist_prunedMass.Fill(tree.JetsAK8_prunedMass[iJet])
			hist_ptD.Fill(tree.JetsAK8_ptD[iJet])
			hist_softDropMass.Fill(tree.JetsAK8_softDropMass[iJet])
			hist_jetNumber.Fill(iJet)
			if nJets == 2:
				hist_jetNumber_nj2.Fill(iJet)
			elif nJets == 3:
				hist_jetNumber_nj3.Fill(iJet)
			elif nJets == 4:
				hist_jetNumber_nj4.Fill(iJet)
			if tree.JetsAK8_isISR[iJet]: #Fill ISR jet histograms
				hist_phi_ISR.Fill(tree.JetsAK8[iJet].Phi())
				hist_eta_ISR.Fill(tree.JetsAK8[iJet].Eta())
				hist_pt_ISR.Fill(tree.JetsAK8[iJet].Pt())
				hist_px_ISR.Fill(tree.JetsAK8[iJet].Px())
				hist_py_ISR.Fill(tree.JetsAK8[iJet].Py())
				hist_pz_ISR.Fill(tree.JetsAK8[iJet].Pz())
				hist_axismajor_ISR.Fill(tree.JetsAK8_axismajor[iJet])
				hist_axisminor_ISR.Fill(tree.JetsAK8_axisminor[iJet])
				hist_doubleBDiscriminator_ISR.Fill(tree.JetsAK8_doubleBDiscriminator[iJet])
				hist_girth_ISR.Fill(tree.JetsAK8_girth[iJet])
				hist_HVfracPt_ISR.Fill(tree.JetsAK8_HVfracPt[iJet])
				hist_HVminDeltaR_ISR.Fill(tree.JetsAK8_HVminDeltaR[iJet])
				hist_ISRfracPt_ISR.Fill(tree.JetsAK8_ISRfracPt[iJet])
				hist_ISRminDeltaR_ISR.Fill(tree.JetsAK8_ISRminDeltaR[iJet])
				hist_momenthalf_ISR.Fill(tree.JetsAK8_momenthalf[iJet])
				hist_multiplicity_ISR.Fill(tree.JetsAK8_multiplicity[iJet])
				hist_nHVISRparts_ISR.Fill(tree.JetsAK8_nHVISRparts[iJet])
				hist_nSMISRparts_ISR.Fill(tree.JetsAK8_nSMISRparts[iJet])
				hist_NsubjettinessTau1_ISR.Fill(tree.JetsAK8_NsubjettinessTau1[iJet])
				hist_NsubjettinessTau2_ISR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_NsubjettinessTau3_ISR.Fill(tree.JetsAK8_NsubjettinessTau3[iJet])
				hist_Tau12_ISR.Fill(tree.JetsAK8_NsubjettinessTau1[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_Tau23_ISR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau3[iJet])
				hist_NumBhadrons_ISR.Fill(tree.JetsAK8_NumBhadrons[iJet])
				hist_NumChadrons_ISR.Fill(tree.JetsAK8_NumChadrons[iJet])
				hist_overflow_ISR.Fill(tree.JetsAK8_overflow[iJet])
				hist_prunedMass_ISR.Fill(tree.JetsAK8_prunedMass[iJet])
				hist_ptD_ISR.Fill(tree.JetsAK8_ptD[iJet])
				hist_softDropMass_ISR.Fill(tree.JetsAK8_softDropMass[iJet])
				hist_jetNumber_ISR.Fill(iJet)
				if nJets == 2:
					hist_jetNumber_nj2_ISR.Fill(iJet)
				elif nJets == 3:
					hist_jetNumber_nj3_ISR.Fill(iJet)
				elif nJets == 4:
					hist_jetNumber_nj4_ISR.Fill(iJet)
			else: # fill FSR jet histograms
				hist_phi_FSR.Fill(tree.JetsAK8[iJet].Phi())
				hist_eta_FSR.Fill(tree.JetsAK8[iJet].Eta())
				hist_pt_FSR.Fill(tree.JetsAK8[iJet].Pt())
				hist_px_FSR.Fill(tree.JetsAK8[iJet].Px())
				hist_py_FSR.Fill(tree.JetsAK8[iJet].Py())
				hist_pz_FSR.Fill(tree.JetsAK8[iJet].Pz())
				hist_axismajor_FSR.Fill(tree.JetsAK8_axismajor[iJet])
				hist_axisminor_FSR.Fill(tree.JetsAK8_axisminor[iJet])
				hist_doubleBDiscriminator_FSR.Fill(tree.JetsAK8_doubleBDiscriminator[iJet])
				hist_girth_FSR.Fill(tree.JetsAK8_girth[iJet])
				hist_HVfracPt_FSR.Fill(tree.JetsAK8_HVfracPt[iJet])
				hist_HVminDeltaR_FSR.Fill(tree.JetsAK8_HVminDeltaR[iJet])
				hist_ISRfracPt_FSR.Fill(tree.JetsAK8_ISRfracPt[iJet])
				hist_ISRminDeltaR_FSR.Fill(tree.JetsAK8_ISRminDeltaR[iJet])
				hist_momenthalf_FSR.Fill(tree.JetsAK8_momenthalf[iJet])
				hist_multiplicity_FSR.Fill(tree.JetsAK8_multiplicity[iJet])
				hist_nHVISRparts_FSR.Fill(tree.JetsAK8_nHVISRparts[iJet])
				hist_nSMISRparts_FSR.Fill(tree.JetsAK8_nSMISRparts[iJet])
				hist_NsubjettinessTau1_FSR.Fill(tree.JetsAK8_NsubjettinessTau1[iJet])
				hist_NsubjettinessTau2_FSR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_NsubjettinessTau3_FSR.Fill(tree.JetsAK8_NsubjettinessTau3[iJet])
				hist_Tau12_FSR.Fill(tree.JetsAK8_NsubjettinessTau1[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_Tau23_FSR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau3[iJet])
				hist_NumBhadrons_FSR.Fill(tree.JetsAK8_NumBhadrons[iJet])
				hist_NumChadrons_FSR.Fill(tree.JetsAK8_NumChadrons[iJet])
				hist_overflow_FSR.Fill(tree.JetsAK8_overflow[iJet])
				hist_prunedMass_FSR.Fill(tree.JetsAK8_prunedMass[iJet])
				hist_ptD_FSR.Fill(tree.JetsAK8_ptD[iJet])
				hist_softDropMass_FSR.Fill(tree.JetsAK8_softDropMass[iJet])
				hist_jetNumber_FSR.Fill(iJet)
				if nJets == 2:
					hist_jetNumber_nj2_FSR.Fill(iJet)
				elif nJets == 3:
					hist_jetNumber_nj3_FSR.Fill(iJet)
				elif nJets == 4:
					hist_jetNumber_nj4_FSR.Fill(iJet)


	

def addLoop():
	baseClass.loop = loop


