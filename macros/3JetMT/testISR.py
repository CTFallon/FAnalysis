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

	#JetsAK8_* : axismajor, axisminor, doubleBDiscriminator, girth, HVfracPt, HVminDeltaR, isISR, ISRfracPt, ISRminDeltaR, momenthalf, multiplicity, nHVISRparts, nSMISRparts, NsubjettinessTau1, NsubjettinessTau2, NsubjettinessTau3, NumBhadrons, NumChadrons, overflow, prunedMass, ptD, softDropMass, kinematicVariables

	hist_phi = self.makeTH1F("hist_phi","phi;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())
	hist_phi_ISR = self.makeTH1F("hist_phi_ISR","phi_ISR;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())
	hist_phi_FSR = self.makeTH1F("hist_phi_FSR","phi_FSR;phi",100,-rt.TMath.Pi(), rt.TMath.Pi())

	hist_eta = self.makeTH1F("hist_eta","eta;eta",100,-5,5)
	hist_eta_ISR = self.makeTH1F("hist_eta_ISR","eta_ISR;eta",100,-5,5)
	hist_eta_FSR = self.makeTH1F("hist_eta_FSR","eta_FSR;eta",100,-5,5)

	hist_mopt = self.makeTH1F("hist_mopt","mopt;mopt;",100,0,1)
	hist_mopt_ISR = self.makeTH1F("hist_mopt_ISR","mopt_ISR;mopt;",100,0,1)
	hist_mopt_FSR = self.makeTH1F("hist_mopt_FSR","mopt_FSR;mopt;",100,0,1)

	hist_pt = self.makeTH1F("hist_pt","pt;pt",100,0,2500)
	hist_pt_ISR = self.makeTH1F("hist_pt_ISR","pt_ISR;pt",100,0,2500)
	hist_pt_FSR = self.makeTH1F("hist_pt_FSR","pt_FSR;pt",100,0,2500)

	hist_px = self.makeTH1F("hist_px","px;px",100,-2000,2000)
	hist_px_ISR = self.makeTH1F("hist_px_ISR","px_ISR;px",100,-2000,2000)
	hist_px_FSR = self.makeTH1F("hist_px_FSR","px_FSR;px",100,-2000,2000)

	hist_py = self.makeTH1F("hist_py","py;py",100,-2000,2000)
	hist_py_ISR = self.makeTH1F("hist_py_ISR","py_ISR;py",100,-2000,2000)
	hist_py_FSR = self.makeTH1F("hist_py_FSR","py_FSR;py",100,-2000,2000)

	hist_pz = self.makeTH1F("hist_pz","pz;pz",100,-6000,6000)
	hist_pz_ISR = self.makeTH1F("hist_pz_ISR","pz_ISR;pz",100,-6000,6000)
	hist_pz_FSR = self.makeTH1F("hist_pz_FSR","pz_FSR;pz",100,-6000,6000)

	hist_axismajor = self.makeTH1F("hist_axismajor","axisMajor;axismajor;",100,0,0.5)
	hist_axismajor_ISR = self.makeTH1F("hist_axismajor_ISR","axisMajor_ISR;axismajor;",100,0,0.5)
	hist_axismajor_FSR = self.makeTH1F("hist_axismajor_FSR","axisMajor_FSR;axismajor;",100,0,0.5)

	hist_axisminor = self.makeTH1F("hist_axisminor","axisMinor;axisminor;",100,0,0.3)
	hist_axisminor_ISR = self.makeTH1F("hist_axisminor_ISR","axisMinor_ISR;axisminor;",100,0,0.3)
	hist_axisminor_FSR = self.makeTH1F("hist_axisminor_FSR","axisMinor_FSR;axisminor;",100,0,0.3)

	hist_doubleBDiscriminator = self.makeTH1F("hist_doubleBDiscriminator","doubleBDiscriminator;doubleBDiscriminator;",100,-1.1,1.1)
	hist_doubleBDiscriminator_ISR = self.makeTH1F("hist_doubleBDiscriminator_ISR","doubleBDiscriminator_ISR;doubleBDiscriminator;",100,-1.1,1.1)
	hist_doubleBDiscriminator_FSR = self.makeTH1F("hist_doubleBDiscriminator_FSR","doubleBDiscriminator_FSR;doubleBDiscriminator;",100,-1.1,1.1)

	hist_girth = self.makeTH1F("hist_girth","girth;girth;",100,0,0.6)
	hist_girth_ISR = self.makeTH1F("hist_girth_ISR","girth_ISR;girth;",100,0,0.6)
	hist_girth_FSR = self.makeTH1F("hist_girth_FSR","girth_FSR;girth;",100,0,0.6)

	hist_girthopt = self.makeTH1F("hist_girthopt","girthopt;girthopt;",100,0,0.004)
	hist_girthopt_ISR = self.makeTH1F("hist_girthopt_ISR","girthopt_ISR;girthopt;",100,0,0.004)
	hist_girthopt_FSR = self.makeTH1F("hist_girthopt_FSR","girthopt_FSR;girthopt;",100,0,0.004)

	hist_HVfracPt = self.makeTH1F("hist_HVfracPt","HVfracPt;HVfracPt;",100,-0.01,1.01)
	hist_HVfracPt_ISR = self.makeTH1F("hist_HVfracPt_ISR","HVfracPt_ISR;HVfracPt;",100,-0.01,1.01)
	hist_HVfracPt_FSR = self.makeTH1F("hist_HVfracPt_FSR","HVfracPt_FSR;HVfracPt;",100,-0.01,1.01)

	hist_HVminDeltaR = self.makeTH1F("hist_HVminDeltaR","HVminDeltaR;HVminDeltaR;",100,0,15.01)
	hist_HVminDeltaR_ISR = self.makeTH1F("hist_HVminDeltaR_ISR","HVminDeltaR_ISR;HVminDeltaR;",100,0,15.01)
	hist_HVminDeltaR_FSR = self.makeTH1F("hist_HVminDeltaR_FSR","HVminDeltaR_FSR;HVminDeltaR;",100,0,15.01)

	hist_ISRfracPt = self.makeTH1F("hist_ISRfracPt","ISRfracPt;ISRfracPt;",100,-0.01,1.01)
	hist_ISRfracPt_ISR = self.makeTH1F("hist_ISRfracPt_ISR","ISRfracPt_ISR;ISRfracPt;",100,-0.01,1.01)
	hist_ISRfracPt_FSR = self.makeTH1F("hist_ISRfracPt_FSR","ISRfracPt_FSR;ISRfracPt;",100,-0.01,1.01)

	hist_ISRminDeltaR = self.makeTH1F("hist_ISRminDeltaR","ISRminDeltaR;ISRminDeltaR;",100,0,15.01)
	hist_ISRminDeltaR_ISR = self.makeTH1F("hist_ISRminDeltaR_ISR","ISRminDeltaR_ISR;ISRminDeltaR;",100,0,15.01)
	hist_ISRminDeltaR_FSR = self.makeTH1F("hist_ISRminDeltaR_FSR","ISRminDeltaR_FSR;ISRminDeltaR;",100,0,15.01)

	hist_momenthalf = self.makeTH1F("hist_momenthalf","momenthalf;momenthalf;",100,0,0.8)
	hist_momenthalf_ISR = self.makeTH1F("hist_momenthalf_ISR","momenthalf_ISR;momenthalf;",100,0,0.8)
	hist_momenthalf_FSR = self.makeTH1F("hist_momenthalf_FSR","momenthalf_FSR;momenthalf;",100,0,0.8)

	hist_multiplicity = self.makeTH1F("hist_multiplicity","multiplicity;multiplicity;",400,0,400)
	hist_multiplicity_ISR = self.makeTH1F("hist_multiplicity_ISR","multiplicity_ISR;multiplicity;",400,0,400)
	hist_multiplicity_FSR = self.makeTH1F("hist_multiplicity_FSR","multiplicity_FSR;multiplicity;",400,0,400)

	hist_nHVISRparts = self.makeTH1F("hist_nHVISRparts","nHVISRparts;nHVISRparts;",15,0,15)
	hist_nHVISRparts_ISR = self.makeTH1F("hist_nHVISRparts_ISR","nHVISRparts_ISR;nHVISRparts;",15,0,15)
	hist_nHVISRparts_FSR = self.makeTH1F("hist_nHVISRparts_FSR","nHVISRparts_FSR;nHVISRparts;",15,0,15)

	hist_nSMISRparts = self.makeTH1F("hist_nSMISRparts","nSMISRparts;nSMISRparts;",35,0,35)
	hist_nSMISRparts_ISR = self.makeTH1F("hist_nSMISRparts_ISR","nSMISRparts_ISR;nSMISRparts;",35,0,35)
	hist_nSMISRparts_FSR = self.makeTH1F("hist_nSMISRparts_FSR","nSMISRparts_FSR;nSMISRparts;",35,0,35)

	hist_NsubjettinessTau1 = self.makeTH1F("hist_NsubjettinessTau1","NsubjettinessTau1;NsubjettinessTau1;",100,0,0.7)
	hist_NsubjettinessTau1_ISR = self.makeTH1F("hist_NsubjettinessTau1_ISR","NsubjettinessTau1_ISR;NsubjettinessTau1;",100,0,0.7)
	hist_NsubjettinessTau1_FSR = self.makeTH1F("hist_NsubjettinessTau1_FSR","NsubjettinessTau1_FSR;NsubjettinessTau1;",100,0,0.7)

	hist_NsubjettinessTau2 = self.makeTH1F("hist_NsubjettinessTau2","NsubjettinessTau2;NsubjettinessTau2;",100,0,0.5)
	hist_NsubjettinessTau2_ISR = self.makeTH1F("hist_NsubjettinessTau2_ISR","NsubjettinessTau2_ISR;NsubjettinessTau2;",100,0,0.5)
	hist_NsubjettinessTau2_FSR = self.makeTH1F("hist_NsubjettinessTau2_FSR","NsubjettinessTau2_FSR;NsubjettinessTau2;",100,0,0.5)

	hist_NsubjettinessTau3 = self.makeTH1F("hist_NsubjettinessTau3","NsubjettinessTau3;NsubjettinessTau3;",100,0,0.5)
	hist_NsubjettinessTau3_ISR = self.makeTH1F("hist_NsubjettinessTau3_ISR","NsubjettinessTau3_ISR;NsubjettinessTau3;",100,0,0.5)
	hist_NsubjettinessTau3_FSR = self.makeTH1F("hist_NsubjettinessTau3_FSR","NsubjettinessTau3_FSR;NsubjettinessTau3;",100,0,0.5)

	hist_Tau21 = self.makeTH1F("hist_Tau21","Tau21;Tau21;",100,0,1.1)
	hist_Tau21_ISR = self.makeTH1F("hist_Tau21_ISR","Tau21_ISR;Tau21;",100,0,1.1)
	hist_Tau21_FSR = self.makeTH1F("hist_Tau21_FSR","Tau21_FSR;Tau21;",100,0,1.1)

	hist_Tau32 = self.makeTH1F("hist_Tau32","Tau32;Tau32;",100,0,1.1)
	hist_Tau32_ISR = self.makeTH1F("hist_Tau32_ISR","Tau32_ISR;Tau32;",100,0,1.1)
	hist_Tau32_FSR = self.makeTH1F("hist_Tau32_FSR","Tau32_FSR;Tau32;",100,0,1.1)

	hist_NumBhadrons = self.makeTH1F("hist_NumBhadrons","NumBhadrons;NumBhadrons;",5,0,5)
	hist_NumBhadrons_ISR = self.makeTH1F("hist_NumBhadrons_ISR","NumBhadrons_ISR;NumBhadrons;",5,0,5)
	hist_NumBhadrons_FSR = self.makeTH1F("hist_NumBhadrons_FSR","NumBhadrons_FSR;NumBhadrons;",5,0,5)

	hist_NumChadrons = self.makeTH1F("hist_NumChadrons","NumChadrons;NumChadrons;",9,0,9)
	hist_NumChadrons_ISR = self.makeTH1F("hist_NumChadrons_ISR","NumChadrons_ISR;NumChadrons;",9,0,9)
	hist_NumChadrons_FSR = self.makeTH1F("hist_NumChadrons_FSR","NumChadrons_FSR;NumChadrons;",9,0,9)

	hist_overflow = self.makeTH1F("hist_overflow","overflow;overflow;",100,0,1)
	hist_overflow_ISR = self.makeTH1F("hist_overflow_ISR","overflow_ISR;overflow;",100,0,1)
	hist_overflow_FSR = self.makeTH1F("hist_overflow_FSR","overflow_FSR;overflow;",100,0,1)

	hist_prunedMass = self.makeTH1F("hist_prunedMass","prunedMass;prunedMass;",200,-10,800)
	hist_prunedMass_ISR = self.makeTH1F("hist_prunedMass_ISR","prunedMass_ISR;prunedMass;",200,-10,800)
	hist_prunedMass_FSR = self.makeTH1F("hist_prunedMass_FSR","prunedMass_FSR;prunedMass;",200,-10,800)

	hist_ptD = self.makeTH1F("hist_ptD","ptD;ptD;",100,0,1)
	hist_ptD_ISR = self.makeTH1F("hist_ptD_ISR","ptD_ISR;ptD;",100,0,1)
	hist_ptD_FSR = self.makeTH1F("hist_ptD_FSR","ptD_FSR;ptD;",100,0,1)

	hist_softDropMass = self.makeTH1F("hist_softDropMass","softDropMass;softDropMass;",200,-10,900)
	hist_softDropMass_ISR = self.makeTH1F("hist_softDropMass_ISR","softDropMass_ISR;softDropMass;",200,-10,900)
	hist_softDropMass_FSR = self.makeTH1F("hist_softDropMass_FSR","softDropMass_FSR;softDropMass;",200,-10,900)

	hist_jetNumber = self.makeTH1F("hist_jetNumber","jetNumber;jetNumber",10,0,10)
	hist_jetNumber_ISR = self.makeTH1F("hist_jetNumber_ISR","jetNumber_ISR;jetNumber",10,0,10)
	hist_jetNumber_FSR = self.makeTH1F("hist_jetNumber_FSR","jetNumber_FSR;jetNumber",10,0,10)

	hist_jetNumber_nj2 = self.makeTH1F("hist_jetNumber_nj2","jetNumber_nj2;jetNumber",10,0,10)
	hist_jetNumber_nj2_ISR = self.makeTH1F("hist_jetNumber_nj2_ISR","jetNumber_nj2_ISR;jetNumber",10,0,10)
	hist_jetNumber_nj2_FSR = self.makeTH1F("hist_jetNumber_nj2_FSR","jetNumber_nj2_FSR;jetNumber",10,0,10)

	hist_jetNumber_nj3 = self.makeTH1F("hist_jetNumber_nj3","jetNumber_nj3;jetNumber",10,0,10)
	hist_jetNumber_nj3_ISR = self.makeTH1F("hist_jetNumber_nj3_ISR","jetNumber_ISR_nj3;jetNumber",10,0,10)
	hist_jetNumber_nj3_FSR = self.makeTH1F("hist_jetNumber_nj3_FSR","jetNumber_FSR_nj3;jetNumber",10,0,10)

	hist_jetNumber_nj4 = self.makeTH1F("hist_jetNumber_nj4","jetNumber_nj4;jetNumber",10,0,10)
	hist_jetNumber_nj4_ISR = self.makeTH1F("hist_jetNumber_nj4_ISR","jetNumber_nj4_ISR;jetNumber",10,0,10)
	hist_jetNumber_nj4_FSR = self.makeTH1F("hist_jetNumber_nj4_FSR","jetNumber_nj4_FSR;jetNumber",10,0,10)
	
	hist_VisPTfromVisHVQ = self.makeTH1F("hist_VisPTfromVisHVQ","VisPTfromVisHVQ;VisPTfromVisHVQ",100,-0.01,1.01)
	hist_VisPTfromVisHVQ_ISR = self.makeTH1F("hist_VisPTfromVisHVQ_ISR","VisPTfromVisHVQ_ISR;VisPTfromVisHVQ",100,-0.01,1.01)
	hist_VisPTfromVisHVQ_FSR = self.makeTH1F("hist_VisPTfromVisHVQ_FSR","VisPTfromVisHVQ_FSR;VisPTfromVisHVQ",100,-0.01,1.01)

	hist_InvPTfromInvHVQ = self.makeTH1F("hist_InvPTfromInvHVQ","InvPTfromInvHVQ;InvPTfromInvHVQ",100,-0.01,1.01)
	hist_InvPTfromInvHVQ_ISR = self.makeTH1F("hist_InvPTfromInvHVQ_ISR","InvPTfromInvHVQ_ISR;InvPTfromInvHVQ",100,-0.01,1.01)
	hist_InvPTfromInvHVQ_FSR = self.makeTH1F("hist_InvPTfromInvHVQ_FSR","InvPTfromInvHVQ_FSR;InvPTfromInvHVQ",100,-0.01,1.01)

	hist_TotPTfromAllHVQ = self.makeTH1F("hist_TotPTfromAllHVQ","TotPTfromAllHVQ;TotPTfromAllHVQ",100,-0.01,1.01)
	hist_TotPTfromAllHVQ_ISR = self.makeTH1F("hist_TotPTfromAllHVQ_ISR","TotPTfromAllHVQ_ISR;TotPTfromAllHVQ",100,-0.01,1.01)
	hist_TotPTfromAllHVQ_FSR = self.makeTH1F("hist_TotPTfromAllHVQ_FSR","TotPTfromAllHVQ_FSR;TotPTfromAllHVQ",100,-0.01,1.01)

	hist_TotPTfromVisHVQ = self.makeTH1F("hist_TotPTfromVisHVQ","TotPTfromVisHVQ;TotPTfromVisHVQ",100,-0.01,1.01)
	hist_TotPTfromVisHVQ_ISR = self.makeTH1F("hist_TotPTfromVisHVQ_ISR","TotPTfromVisHVQ_ISR;TotPTfromVisHVQ",100,-0.01,1.01)
	hist_TotPTfromVisHVQ_FSR = self.makeTH1F("hist_TotPTfromVisHVQ_FSR","TotPTfromVisHVQ_FSR;TotPTfromVisHVQ",100,-0.01,1.01)

	hist_TotPTfromInvHVQ = self.makeTH1F("hist_TotPTfromInvHVQ","TotPTfromInvHVQ;TotPTfromInvHVQ",100,-0.01,1.01)
	hist_TotPTfromInvHVQ_ISR = self.makeTH1F("hist_TotPTfromInvHVQ_ISR","TotPTfromInvHVQ_ISR;TotPTfromInvHVQ",100,-0.01,1.01)
	hist_TotPTfromInvHVQ_FSR = self.makeTH1F("hist_TotPTfromInvHVQ_FSR","TotPTfromInvHVQ_FSR;TotPTfromInvHVQ",100,-0.01,1.01)

	hist_TotPTfromVis = self.makeTH1F("hist_TotPTfromVis","TotPTfromVis;TotPTfromVis",100,-0.01,1.01)
	hist_TotPTfromVis_ISR = self.makeTH1F("hist_TotPTfromVis_ISR","TotPTfromVis_ISR;TotPTfromVis",100,-0.01,1.01)
	hist_TotPTfromVis_FSR = self.makeTH1F("hist_TotPTfromVis_FSR","TotPTfromVis_FSR;TotPTfromVis",100,-0.01,1.01)

	hist_MT_dijet = self.makeTH1F("hist_MT_dijet","MT_dijet;MT;",100,0,6000)
	hist_MT_trijet = self.makeTH1F("hist_MT_trujet","MT_trijet;MT;",100,0,6000)
	hist_MT_FSR = self.makeTH1F("hist_MT_FSR","MT_FSR;MT;",100,0,6000)
	hist_MT_AllJets = self.makeTH1F("hist_MT_AllJets","MT_AllJets;MT;",100,0,6000)
	hist_MT_jet3dPhi = self.makeTH1F("hist_MT_jet3dPhi","MT_jet3dPhi;MT;",100,0,6000)
	
	hist_deltaPhi = self.makeTH1F("hist_deltaPhi","deltaPhi;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi_ISR = self.makeTH1F("hist_deltaPhi_ISR","deltaPhi_ISR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi_FSR = self.makeTH1F("hist_deltaPhi_FSR","deltaPhi_FSR;dPhi;",100,0,rt.TMath.Pi())

	hist_deltaPhiopt = self.makeTH1F("hist_deltaPhiopt","deltaPhiopt;dPhi;",100,0,0.02)
	hist_deltaPhiopt_ISR = self.makeTH1F("hist_deltaPhiopt_ISR","deltaPhiopt_ISR;dPhi;",100,0,0.02)
	hist_deltaPhiopt_FSR = self.makeTH1F("hist_deltaPhiopt_FSR","deltaPhiopt_FSR;dPhi;",100,0,0.02)

	hist_deltaPhi1 = self.makeTH1F("hist_deltaPhi1","deltaPhi1;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi2 = self.makeTH1F("hist_deltaPhi2","deltaPhi2;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi3 = self.makeTH1F("hist_deltaPhi3","deltaPhi3;dPhi;",100,0,rt.TMath.Pi())

	hist_deltaPhi1_ISR = self.makeTH1F("hist_deltaPhi1_ISR","deltaPhi1_ISR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi1_FSR = self.makeTH1F("hist_deltaPhi1_FSR","deltaPhi1_FSR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi2_ISR = self.makeTH1F("hist_deltaPhi2_ISR","deltaPhi2_ISR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi2_FSR = self.makeTH1F("hist_deltaPhi2_FSR","deltaPhi2_FSR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi3_ISR = self.makeTH1F("hist_deltaPhi3_ISR","deltaPhi3_ISR;dPhi;",100,0,rt.TMath.Pi())
	hist_deltaPhi3_FSR = self.makeTH1F("hist_deltaPhi3_FSR","deltaPhi3_FSR;dPhi;",100,0,rt.TMath.Pi())
	
	hist_FSRJets = self.makeTH1F("hist_FSRJets","JetCode;Code;",32,0,32)
	hist_FSRJets_nj2 = self.makeTH1F("hist_FSRJets_nj2","JetCode_nj2;Code;",32,0,32)
	hist_FSRJets_nj3 = self.makeTH1F("hist_FSRJets_nj3","JetCode_nj3;Code;",32,0,32)
	hist_FSRJets_nj4 = self.makeTH1F("hist_FSRJets_nj4","JetCode_nj4;Code;",32,0,32)
	hist_FSRJets_nj5 = self.makeTH1F("hist_FSRJets_nj5","JetCode_nj5;Code;",32,0,32)

	hist_nHVParts = self.makeTH1F("hist_nHVParts","nHVparts;nHVParts;", 15,0,15)
	hist_nHVParts_ISR = self.makeTH1F("hist_nHVParts_ISR","nHVparts_ISR;nHVParts;", 15,0,15)
	hist_nHVParts_FSR = self.makeTH1F("hist_nHVParts_FSR","nHVparts_FSR;nHVParts;", 15,0,15)
	
	hist_nSMParts = self.makeTH1F("hist_nSMParts","nSMparts;nSMParts;", 30,0,30)
	hist_nSMParts_ISR = self.makeTH1F("hist_nSMParts_ISR","nSMparts_ISR;nSMParts;", 30,0,30)
	hist_nSMParts_FSR = self.makeTH1F("hist_nSMParts_FSR","nSMparts_FSR;nSMParts;", 30,0,30)

	hist2d_jetCodevsminDeltaPhi = self.makeTH2F("hist2d_jetCodevsminDeltaPhi",";jetCode;MinDeltaPhiIndex",32,0,32,3,0,3)

	hist_minDPhi = self.makeTH1F("hist_minDPhi","mindPhi;dPhi;",100,0,rt.TMath.Pi())
	hist_minDPhi_ISR = self.makeTH1F("hist_minDPhi_ISR","mindPhi_ISR;dPhi;",100,0,rt.TMath.Pi())
	hist_minDPhi_FSR = self.makeTH1F("hist_minDPhi_FSR","mindPhi_FSR;dPhi;",100,0,rt.TMath.Pi())
	



	dicOfFunkyJets = {}
	totJets = 0
	totFunkyJets = 0

	etaPhiTracker = 0

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		#if len(tree.JetsAK8) > 2:
		#	continue
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		#Fill MT histograms
		hist_MT_dijet.Fill(trans_mass_Njet(tree.JetsAK8[0:2], tree.MET, tree.METPhi))
		hist_MT_AllJets.Fill(trans_mass_Njet(tree.JetsAK8[0:min(4, nJets)], tree.MET, tree.METPhi))
		if nJets >=3:
			hist_MT_trijet.Fill(trans_mass_Njet(tree.JetsAK8[0:3], tree.MET, tree.METPhi))
		else:
			hist_MT_trijet.Fill(trans_mass_Njet(tree.JetsAK8[0:2], tree.MET, tree.METPhi))
		FSRJets = []
		jetCode = [0,0,0,0,0]
		eventHasFunkyJets = False

		jetClosestToDeltaPhi = -1
		minDeltaPhi = min(tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3)
		if tree.DeltaPhi1 == minDeltaPhi:
			jetClosestToDeltaPhi = 0
		elif tree.DeltaPhi2 == minDeltaPhi:
			jetClosestToDeltaPhi = 1
		elif tree.DeltaPhi3 == minDeltaPhi:
			jetClosestToDeltaPhi = 2

		for iJet in range(nJets):
			#if tree.JetsAK8[iJet].Pt() > 1000.:
			#	continue
			if not tree.JetsAK8_isISR[iJet]:
				FSRJets.append(tree.JetsAK8[iJet])
				if iJet <= 4:
					jetCode[-iJet-1] = 1
			totJets +=1
			if (tree.JetsAK8_isISR[iJet] and tree.JetsAK8_HVfracPt[iJet]>0.9):
				totFunkyJets += 1
				eventHasFunkyJets = True


			nHVParts = 0
			nSMParts = 0
			#for iPart in range(len(tree.GenParticles)):
			#	if tree.GenParticles_Status[iPart] == 1 and tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart])< 0.8:
			#		if tree.GenParticles_PdgId[iPart] >= 4900000:
			#			nHVParts += 1
			#		if tree.GenParticles_PdgId[iPart] < 4900000:
			#			nSMParts += 1
	
			
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
			hist_Tau21.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau1[iJet])
			hist_Tau32.Fill(tree.JetsAK8_NsubjettinessTau3[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
			hist_NumBhadrons.Fill(tree.JetsAK8_NumBhadrons[iJet])
			hist_NumChadrons.Fill(tree.JetsAK8_NumChadrons[iJet])
			hist_overflow.Fill(tree.JetsAK8_overflow[iJet])
			hist_prunedMass.Fill(tree.JetsAK8_prunedMass[iJet])
			hist_ptD.Fill(tree.JetsAK8_ptD[iJet])
			hist_softDropMass.Fill(tree.JetsAK8_softDropMass[iJet])
			hist_jetNumber.Fill(iJet)
			hist_VisPTfromVisHVQ.Fill(tree.fracVisPTfromVisHVQ[iJet])
			hist_InvPTfromInvHVQ.Fill(tree.fracInvPTfromInvHVQ[iJet])
			hist_TotPTfromAllHVQ.Fill(tree.fracTotPTfromAllHVQ[iJet])
			hist_TotPTfromVisHVQ.Fill(tree.fracTotPTfromVisHVQ[iJet])
			hist_TotPTfromInvHVQ.Fill(tree.fracTotPTfromInvHVQ[iJet])
			hist_TotPTfromVis.Fill(tree.fracTotPTfromVis[iJet])
			hist_deltaPhi.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi())))
			hist_mopt.Fill(tree.JetsAK8[iJet].M()/tree.JetsAK8[iJet].Pt())
			hist_nHVParts.Fill(nHVParts)
			hist_nSMParts.Fill(nSMParts)
			hist_girthopt.Fill(tree.JetsAK8_girth[iJet]/tree.JetsAK8[iJet].Pt())
			hist_deltaPhiopt.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi()))/tree.JetsAK8[iJet].Pt())
			if nJets == 2:
				hist_jetNumber_nj2.Fill(iJet)
			elif nJets == 3:
				hist_jetNumber_nj3.Fill(iJet)
			elif nJets == 4:
				hist_jetNumber_nj4.Fill(iJet)
			if iJet == 0:
				hist_deltaPhi1.Fill(tree.DeltaPhi1)
			if iJet == 1:
				hist_deltaPhi2.Fill(tree.DeltaPhi2)
			if iJet == 2:
				hist_deltaPhi3.Fill(tree.DeltaPhi3)
			if iJet == jetClosestToDeltaPhi:
				hist_minDPhi.Fill(minDeltaPhi)
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
				hist_Tau21_ISR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau1[iJet])
				hist_Tau32_ISR.Fill(tree.JetsAK8_NsubjettinessTau3[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_NumBhadrons_ISR.Fill(tree.JetsAK8_NumBhadrons[iJet])
				hist_NumChadrons_ISR.Fill(tree.JetsAK8_NumChadrons[iJet])
				hist_overflow_ISR.Fill(tree.JetsAK8_overflow[iJet])
				hist_prunedMass_ISR.Fill(tree.JetsAK8_prunedMass[iJet])
				hist_ptD_ISR.Fill(tree.JetsAK8_ptD[iJet])
				hist_softDropMass_ISR.Fill(tree.JetsAK8_softDropMass[iJet])
				hist_jetNumber_ISR.Fill(iJet)
				hist_VisPTfromVisHVQ_ISR.Fill(tree.fracVisPTfromVisHVQ[iJet])
				hist_InvPTfromInvHVQ_ISR.Fill(tree.fracInvPTfromInvHVQ[iJet])
				hist_TotPTfromAllHVQ_ISR.Fill(tree.fracTotPTfromAllHVQ[iJet])
				hist_TotPTfromVisHVQ_ISR.Fill(tree.fracTotPTfromVisHVQ[iJet])
				hist_TotPTfromInvHVQ_ISR.Fill(tree.fracTotPTfromInvHVQ[iJet])
				hist_TotPTfromVis_ISR.Fill(tree.fracTotPTfromVis[iJet])
				hist_deltaPhi_ISR.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi())))
				hist_mopt_ISR.Fill(tree.JetsAK8[iJet].M()/tree.JetsAK8[iJet].Pt())
				hist_nHVParts_ISR.Fill(nHVParts)
				hist_nSMParts_ISR.Fill(nSMParts)
				hist_girthopt_ISR.Fill(tree.JetsAK8_girth[iJet]/tree.JetsAK8[iJet].Pt())
				hist_deltaPhiopt_ISR.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi()))/tree.JetsAK8[iJet].Pt())
				if nJets == 2:
					hist_jetNumber_nj2_ISR.Fill(iJet)
				elif nJets == 3:
					hist_jetNumber_nj3_ISR.Fill(iJet)
				elif nJets == 4:
					hist_jetNumber_nj4_ISR.Fill(iJet)
				if iJet == 0:
					hist_deltaPhi1_ISR.Fill(tree.DeltaPhi1)
				if iJet == 1:
					hist_deltaPhi2_ISR.Fill(tree.DeltaPhi2)
				if iJet == 2:
					hist_deltaPhi3_ISR.Fill(tree.DeltaPhi3)
				if iJet == jetClosestToDeltaPhi:
					hist_minDPhi_ISR.Fill(minDeltaPhi)
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
				hist_Tau21_FSR.Fill(tree.JetsAK8_NsubjettinessTau2[iJet]/tree.JetsAK8_NsubjettinessTau1[iJet])
				hist_Tau32_FSR.Fill(tree.JetsAK8_NsubjettinessTau3[iJet]/tree.JetsAK8_NsubjettinessTau2[iJet])
				hist_NumBhadrons_FSR.Fill(tree.JetsAK8_NumBhadrons[iJet])
				hist_NumChadrons_FSR.Fill(tree.JetsAK8_NumChadrons[iJet])
				hist_overflow_FSR.Fill(tree.JetsAK8_overflow[iJet])
				hist_prunedMass_FSR.Fill(tree.JetsAK8_prunedMass[iJet])
				hist_ptD_FSR.Fill(tree.JetsAK8_ptD[iJet])
				hist_softDropMass_FSR.Fill(tree.JetsAK8_softDropMass[iJet])
				hist_jetNumber_FSR.Fill(iJet)
				hist_VisPTfromVisHVQ_FSR.Fill(tree.fracVisPTfromVisHVQ[iJet])
				hist_InvPTfromInvHVQ_FSR.Fill(tree.fracInvPTfromInvHVQ[iJet])
				hist_TotPTfromAllHVQ_FSR.Fill(tree.fracTotPTfromAllHVQ[iJet])
				hist_TotPTfromVisHVQ_FSR.Fill(tree.fracTotPTfromVisHVQ[iJet])
				hist_TotPTfromInvHVQ_FSR.Fill(tree.fracTotPTfromInvHVQ[iJet])
				hist_TotPTfromVis_FSR.Fill(tree.fracTotPTfromVis[iJet])
				hist_deltaPhi_FSR.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi())))
				hist_mopt_FSR.Fill(tree.JetsAK8[iJet].M()/tree.JetsAK8[iJet].Pt())
				hist_nHVParts_FSR.Fill(nHVParts)
				hist_nSMParts_FSR.Fill(nSMParts)
				hist_girthopt_FSR.Fill(tree.JetsAK8_girth[iJet]/tree.JetsAK8[iJet].Pt())
				hist_deltaPhiopt_FSR.Fill(abs(deltaPhi(tree.METPhi, tree.JetsAK8[iJet].Phi()))/tree.JetsAK8[iJet].Pt())
				if nJets == 2:
					hist_jetNumber_nj2_FSR.Fill(iJet)
				elif nJets == 3:
					hist_jetNumber_nj3_FSR.Fill(iJet)
				elif nJets == 4:
					hist_jetNumber_nj4_FSR.Fill(iJet)
				if iJet == 0:
					hist_deltaPhi1_FSR.Fill(tree.DeltaPhi1)
				if iJet == 1:
					hist_deltaPhi2_FSR.Fill(tree.DeltaPhi2)
				if iJet == 2:
					hist_deltaPhi3_FSR.Fill(tree.DeltaPhi3)
				if iJet == jetClosestToDeltaPhi:
					hist_minDPhi_FSR.Fill(minDeltaPhi)
				
		
		hist_MT_FSR.Fill(trans_mass_Njet(FSRJets, tree.MET, tree.METPhi))
		jetCodeDec = jetCode[-1]*1+jetCode[-2]*2+jetCode[-3]*4+jetCode[-4]*8+jetCode[-5]*16
		hist_FSRJets.Fill(jetCodeDec)
		if nJets == 2:
			hist_FSRJets_nj2.Fill(jetCodeDec)
		if nJets == 3:
			hist_FSRJets_nj3.Fill(jetCodeDec)
		if nJets == 4:
			hist_FSRJets_nj4.Fill(jetCodeDec)
		if nJets == 5:
			hist_FSRJets_nj5.Fill(jetCodeDec)
		
		

		if nJets == 3:
			hist2d_jetCodevsminDeltaPhi.Fill(jetCodeDec, jetClosestToDeltaPhi)

		if nJets == 2:
			hist_MT_jet3dPhi.Fill(trans_mass_Njet(tree.JetsAK8, tree.MET, tree.METPhi))
		else:
			if jetClosestToDeltaPhi == 2:
				hist_MT_jet3dPhi.Fill(trans_mass_Njet(tree.JetsAK8[0:3], tree.MET, tree.METPhi))
			else:
				hist_MT_jet3dPhi.Fill(trans_mass_Njet(tree.JetsAK8[0:2], tree.MET, tree.METPhi))
		
		#if eventHasFunkyJets:
		#	#funkyEventDict[int(tree.EvtNum)] = listOfFunkyJets
		#	partsList = []
		#	partsListPdgId = []
		#	partsListisFromHVQuark = []
		#	for iPart in range(len(tree.GenParticles)):
		#		if (tree.GenParticles_Status[iPart] == 1):
		#			partsList.append(tree.GenParticles[iPart])
		#			partsListPdgId.append(tree.GenParticles_PdgId[iPart])
		#			partsListisFromHVQuark.append(tree.genParticleIsFromHVQuark[iPart])
		#	drawEventEtaPhiPlot(tree.JetsAK8, tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, tree.genParticleIsFromHVQuark, tree.EvtNum, self.extraDir+"etaPhi/all")
		#	drawEventEtaPhiPlot(tree.JetsAK8, partsList, partsListPdgId, tree.METPhi, partsListisFromHVQuark, tree.EvtNum, self.extraDir+"etaPhi/stat")

	makePlots(hist_phi, hist_phi_ISR, hist_phi_FSR,  self.extraDir)
	makePlots(hist_eta, hist_eta_ISR, hist_eta_FSR,   self.extraDir)
	makePlots(hist_pt, hist_pt_ISR, hist_pt_FSR,   self.extraDir)
	makePlots(hist_px, hist_px_ISR, hist_px_FSR,   self.extraDir)
	makePlots(hist_py, hist_py_ISR, hist_py_FSR,   self.extraDir)
	makePlots(hist_pz, hist_pz_ISR, hist_pz_FSR,   self.extraDir)
	makePlots(hist_axismajor, hist_axismajor_ISR, hist_axismajor_FSR,   self.extraDir)
	makePlots(hist_axisminor, hist_axisminor_ISR, hist_axisminor_FSR,   self.extraDir)
	makePlots(hist_doubleBDiscriminator, hist_doubleBDiscriminator_ISR, hist_doubleBDiscriminator_FSR,   self.extraDir)
	makePlots(hist_girth, hist_girth_ISR, hist_girth_FSR,   self.extraDir)
	makePlots(hist_HVfracPt, hist_HVfracPt_ISR, hist_HVfracPt_FSR,   self.extraDir)
	makePlots(hist_HVminDeltaR, hist_HVminDeltaR_ISR, hist_HVminDeltaR_FSR,   self.extraDir)
	makePlots(hist_ISRfracPt, hist_ISRfracPt_ISR, hist_ISRfracPt_FSR,   self.extraDir)
	makePlots(hist_ISRminDeltaR, hist_ISRminDeltaR_ISR, hist_ISRminDeltaR_FSR,   self.extraDir)
	makePlots(hist_momenthalf, hist_momenthalf_ISR, hist_momenthalf_FSR,   self.extraDir)
	makePlots(hist_multiplicity, hist_multiplicity_ISR, hist_multiplicity_FSR,   self.extraDir)
	makePlots(hist_nHVISRparts, hist_nHVISRparts_ISR, hist_nHVISRparts_FSR,   self.extraDir)
	makePlots(hist_nSMISRparts, hist_nSMISRparts_ISR, hist_nSMISRparts_FSR,   self.extraDir)
	makePlots(hist_NsubjettinessTau1, hist_NsubjettinessTau1_ISR, hist_NsubjettinessTau1_FSR,   self.extraDir)
	makePlots(hist_NsubjettinessTau2, hist_NsubjettinessTau2_ISR, hist_NsubjettinessTau2_FSR,   self.extraDir)
	makePlots(hist_NsubjettinessTau3, hist_NsubjettinessTau3_ISR, hist_NsubjettinessTau3_FSR,   self.extraDir)
	makePlots(hist_Tau21, hist_Tau21_ISR, hist_Tau21_FSR,   self.extraDir)
	makePlots(hist_Tau32, hist_Tau32_ISR , hist_Tau32_FSR,   self.extraDir)
	makePlots(hist_NumBhadrons, hist_NumBhadrons_ISR, hist_NumBhadrons_FSR,   self.extraDir)
	makePlots(hist_NumChadrons, hist_NumChadrons_ISR, hist_NumChadrons_FSR,   self.extraDir)
	makePlots(hist_overflow, hist_overflow_ISR, hist_overflow_FSR,   self.extraDir)
	makePlots(hist_prunedMass, hist_prunedMass_ISR, hist_prunedMass_FSR,   self.extraDir)
	makePlots(hist_ptD, hist_ptD_ISR, hist_ptD_FSR,   self.extraDir)
	makePlots(hist_softDropMass, hist_softDropMass_ISR, hist_softDropMass_FSR,   self.extraDir)
	makePlots(hist_jetNumber, hist_jetNumber_ISR, hist_jetNumber_FSR,   self.extraDir)
	makePlots(hist_jetNumber_nj2, hist_jetNumber_nj2_ISR, hist_jetNumber_nj2_FSR,   self.extraDir)
	makePlots(hist_jetNumber_nj3, hist_jetNumber_nj3_ISR, hist_jetNumber_nj3_FSR,   self.extraDir)
	makePlots(hist_jetNumber_nj4, hist_jetNumber_nj4_ISR, hist_jetNumber_nj4_FSR,   self.extraDir)	
	makePlots(hist_VisPTfromVisHVQ, hist_VisPTfromVisHVQ_ISR, hist_VisPTfromVisHVQ_FSR,   self.extraDir)
	makePlots(hist_InvPTfromInvHVQ, hist_InvPTfromInvHVQ_ISR, hist_InvPTfromInvHVQ_FSR,   self.extraDir)
	makePlots(hist_TotPTfromAllHVQ, hist_TotPTfromAllHVQ_ISR, hist_TotPTfromAllHVQ_FSR,   self.extraDir)
	makePlots(hist_TotPTfromVisHVQ, hist_TotPTfromVisHVQ_ISR, hist_TotPTfromVisHVQ_FSR,   self.extraDir)
	makePlots(hist_TotPTfromInvHVQ, hist_TotPTfromInvHVQ_ISR, hist_TotPTfromInvHVQ_FSR,   self.extraDir)
	makePlots(hist_TotPTfromVis, hist_TotPTfromVis_ISR, hist_TotPTfromVis_FSR,   self.extraDir)
	makePlots(hist_deltaPhi, hist_deltaPhi_ISR, hist_deltaPhi_FSR, self.extraDir)
	makePlots(hist_MT_dijet, hist_MT_trijet, hist_MT_FSR, self.extraDir, extraName="_trijet_FSR")
	makePlots(hist_MT_dijet,hist_MT_trijet,  hist_MT_FSR, self.extraDir,False,extraName="_trijet_FSR_Linear")
	makePlots(hist_MT_dijet,hist_MT_trijet,  hist_MT_AllJets, self.extraDir,extraName="_trijet_AllJets")
	makePlots(hist_MT_dijet,hist_MT_FSR,  hist_MT_AllJets, self.extraDir,extraName="_FSR_AllJets")
	makePlots(hist_MT_dijet,hist_MT_trijet,  hist_MT_jet3dPhi, self.extraDir,extraName="_trijet_jet3dPhi")
	makePlots(hist_MT_dijet,hist_MT_FSR,  hist_MT_jet3dPhi, self.extraDir,extraName="_FSR_jet3dPhi")
	makePlots(hist_MT_dijet,hist_MT_trijet,  hist_MT_AllJets, self.extraDir,False,extraName="_trijet_AllJets_Linear")
	makePlots(hist_MT_dijet,hist_MT_FSR,  hist_MT_AllJets, self.extraDir,False,extraName="_FSR_AllJets_Linear")
	makePlots(hist_MT_dijet,hist_MT_trijet,  hist_MT_jet3dPhi, self.extraDir,False,extraName="_trijet_jet3dPhi_Linear")
	makePlots(hist_MT_dijet,hist_MT_FSR,  hist_MT_jet3dPhi, self.extraDir,False,extraName="_FSR_jet3dPhi_Linear")
	makePlots(hist_FSRJets, hist_FSRJets, hist_FSRJets, self.extraDir)
	makePlots(hist_FSRJets_nj2, hist_FSRJets_nj2, hist_FSRJets_nj2, self.extraDir)
	makePlots(hist_FSRJets_nj3, hist_FSRJets_nj3, hist_FSRJets_nj3, self.extraDir)
	makePlots(hist_FSRJets_nj4, hist_FSRJets_nj4, hist_FSRJets_nj4, self.extraDir)
	makePlots(hist_FSRJets_nj5, hist_FSRJets_nj5, hist_FSRJets_nj5, self.extraDir)
	makePlots(hist_FSRJets_nj2, hist_FSRJets_nj3, hist_FSRJets_nj4, self.extraDir, extraName="34")
	makePlots(hist_FSRJets_nj3, hist_FSRJets_nj4, hist_FSRJets_nj5, self.extraDir, extraName="45")
	makePlots(hist_mopt, hist_mopt_ISR,hist_mopt_FSR, self.extraDir)
	makePlots(hist_nHVParts,hist_nHVParts_ISR,hist_nHVParts_FSR,self.extraDir)
	makePlots(hist_nSMParts,hist_nSMParts_ISR,hist_nSMParts_FSR,self.extraDir)
	makePlots(hist_deltaPhi1,hist_deltaPhi2,hist_deltaPhi3,self.extraDir,extraName="23")
	makePlots(hist_deltaPhi1_ISR,hist_deltaPhi2_ISR,hist_deltaPhi3_ISR,self.extraDir,extraName="23_ISR")
	makePlots(hist_deltaPhi1_FSR,hist_deltaPhi2_FSR,hist_deltaPhi3_FSR,self.extraDir,extraName="23_FSR")
	makePlots(hist_deltaPhi1,hist_deltaPhi1_ISR,hist_deltaPhi1_FSR,self.extraDir,extraName="_ISR_FSR")
	makePlots(hist_deltaPhi2,hist_deltaPhi2_ISR,hist_deltaPhi2_FSR,self.extraDir,extraName="_ISR_FSR")
	makePlots(hist_deltaPhi3,hist_deltaPhi3_ISR,hist_deltaPhi3_FSR,self.extraDir,extraName="_ISR_FSR")
	makePlots(hist_girthopt, hist_girthopt_ISR, hist_girthopt_FSR,   self.extraDir)
	makePlots(hist_deltaPhiopt, hist_deltaPhiopt_ISR, hist_deltaPhiopt_FSR, self.extraDir)
	makePlots(hist_minDPhi,hist_minDPhi_ISR,hist_minDPhi_FSR,self.extraDir)
	fitToCrystalBall(hist_MT_dijet, self.extraDir)
	fitToCrystalBall(hist_MT_trijet, self.extraDir)
	fitToCrystalBall(hist_MT_FSR, self.extraDir)
	fitToCrystalBall(hist_MT_AllJets, self.extraDir)
	fitToCrystalBall(hist_MT_jet3dPhi, self.extraDir)
	c1 = rt.TCanvas('c1','c1',900,600)
	hist2d_jetCodevsminDeltaPhi.Draw("colz")
	c1.SetLogz()
	c1.SaveAs(self.extraDir+"jetCodevsMinDeltaPhiIndex.png")

	print("Jets | MT Resolution")
	print("DiJet | " + str(hist_MT_dijet.GetRMS()/hist_MT_dijet.GetMean()))
	print("TriJet | " + str(hist_MT_trijet.GetRMS()/hist_MT_trijet.GetMean()))
	print("FSR jets | " + str(hist_MT_FSR.GetRMS()/hist_MT_FSR.GetMean()))
	print("All Jets | " + str(hist_MT_AllJets.GetRMS()/hist_MT_AllJets.GetMean()))
	print("jet3dPhi | " + str(hist_MT_jet3dPhi.GetRMS()/hist_MT_jet3dPhi.GetMean()))
	print("Tot Jets and nJetsFunky: " + str(totJets) + " " + str(totFunkyJets))
	for x in range(15):
		print(bin(x), x)
	for e, jL in dicOfFunkyJets.items():
		print(e, jL)
		
def addLoop():
	baseClass.loop = loop

def makePlots(All, ISR, FSR, direc, logY = True,extraName = ""):
	c = rt.TCanvas("canvas","canvas",900,600)
	if logY:
		c.SetLogy()
	histStack = rt.THStack()
	histStack.SetTitle(All.GetTitle())
	All.SetLineColor(1)
	ISR.SetLineColor(2)
	FSR.SetLineColor(3)
	histStack.Add(All)
	histStack.Add(ISR)
	histStack.Add(FSR)
	histStack.Draw("NOSTACK")
	histStack.GetXaxis().SetTitle(All.GetXaxis().GetTitle())
	c.BuildLegend(0.75,0.75,0.9,0.9,"")
	c.SaveAs(direc+All.GetTitle()+extraName+".png")

def deltaPhi(phi1, phi2):
	x = phi1 - phi2
	while x >= rt.TMath.Pi():
		x = x - 2*rt.TMath.Pi()
	while x < -rt.TMath.Pi():
		x = x + 2*rt.TMath.Pi()
	return x

def drawEventEtaPhiPlot(jetCollectionAK8, partCol, particlePDGID,METPhi, isFromHVQuark, plotNumber, edir):
	canv = rt.TCanvas("canv","canv",1600,800)
	histAxis = rt.TH2F("axisHsito", ";\eta;\phi",100,-6,6,100,-rt.TMath.Pi(),rt.TMath.Pi())
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionAK8)):
		objectList.append(rt.TEllipse(jetCollectionAK8[iJet].Eta(), jetCollectionAK8[iJet].Phi(), 0.8, 0.8))
		objectList[-1].SetLineColor(iJet+1)
		objectList[-1].SetLineWidth(2)
	for iPart in range(len(partCol)):
		objectList.append(rt.TMarker(partCol[iPart].Eta(),partCol[iPart].Phi(),2))
		objectList[-1].SetMarkerSize(2)
		if abs(particlePDGID[iPart]) == 4900101:
			objectList[-1].SetMarkerColor(3)
			objectList[-1].SetMarkerStyle(22)
			objectList[-1].SetMarkerSize(3)
		elif abs(particlePDGID[iPart] == 4900023):
			objectList[-1].SetMarkerStyle(43)
			objectList[-1].SetMarkerSize(3)
		elif abs(particlePDGID[iPart]) > 4900000:
			objectList[-1].SetMarkerColor(2)
		if isFromHVQuark[iPart]:
			objectList[-1].SetMarkerStyle(5)
	objectList.append(rt.TLine(-6,METPhi,6,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	canv.SaveAs(edir+"etaPhi_"+str(plotNumber)+".png")

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def fitToCrystalBall(histo, direc, extraName=""):
	cbFunc = rt.TF1("cbFunc","crystalball",0,6000)
	cbFunc.SetParameters(800, 2500, 400, 0.8,1000000)
	cbFunc.SetLineColor(2)
	histo.Fit('cbFunc')
	print(histo.GetEntries())
	c2 = rt.TCanvas("c2","c2",900,600)
	histo.SetLineColor(1)
	histo.Draw()
	cbFunc.Draw("same")
	c2.SaveAs(direc+"cbFunc_"+histo.GetName()+extraName+".png")
	


