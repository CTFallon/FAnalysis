from analysisBase import baseClass
import ROOT as rt
from array import array
# try to find GEN level variables that can be used to discriminate which jets should be used in the MT calculation
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
	tree.SetBranchStatus("GenJets",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	#tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	tree.SetBranchStatus("numHVPartsInJet",1)
	tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)
	tree.SetBranchStatus("zPrimept",1)
	tree.SetBranchStatus("zPrimephi",1)
	tree.SetBranchStatus("pGJ_visible",1)
	tree.SetBranchStatus("pGJ_invis",1)
	tree.SetBranchStatus("pGJ_every",1)
	tree.SetBranchStatus("fracVisPTfromVisHVQ",1)
	#tree.SetBranchStatus("fracInvPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromAllHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVisHVQ",1)
	tree.SetBranchStatus("fracTotPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVis",1)
	tree.SetBranchStatus("fracVisHVQtoInvHVQ",1)
	tree.SetBranchStatus("pt_VisHVQ",1)
	tree.SetBranchStatus("pt_Vis",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	# make histograms for all events with 2 and 3 jets
	hist_MT_2jets_all = self.makeTH1F("hist_MT_2jets_all","2 Jets;MT;count/a.u.",100,0,4000)
	hist_MT_3jets_all = self.makeTH1F("hist_MT_3jets_all","3 Jets;MT;count/a.u.",100,0,4000)
	hist_MT_MTcut = self.makeTH1F("hist_MT_MTcut","MTcut;MT;count/a.u.",100,0,4000)
	hist_deltaMT_2jet3jet = self.makeTH1F("hist_deltaMT_2jet3jet","Difference of MT;\Delta MT;count",100,-100,4000)


	# fractionalPT plots
	hist_fracVisPTfromVisHVQ = self.makeTH2F("hist_fracVisPTfromVisHVQ","VisHVQ/Vis;iJet;FracPt",10,0,10,100,-0.01,1.01)
	hist_fracTotPTfromAllHVQ = self.makeTH2F("hist_fracTotPTfromAllHVQ","AllHVQ/Tot;iJet;FracPt",10,0,10,100,-0.01,1.01)
	hist_fracTotPTfromVisHVQ = self.makeTH2F("hist_fracTotPTfromVisHVQ","VisHVQ/Tot;iJet;FracPt",10,0,10,100,-0.01,1.01)
	hist_fracTotPTfromInvHVQ = self.makeTH2F("hist_fracTotPTfromInvHVQ","InvHVQ/Tot;iJet;FracPt",10,0,10,100,-0.01,1.01)
	hist_fracTotPTfromVis = self.makeTH2F("hist_fracTotPTfromVis","Vis/Tot;iJet;FracPt",10,0,10,100,-0.01,1.01)
	hist_fracVisHVQtoInvHVQ = self.makeTH2F("hist_fracVisHVQtoInvHVQ","VisHVQ/InvHVQ;iJet;FractPt",10,0,10,1000,-0.01,1.01)
	hist_VisPt = self.makeTH2F("hist_VisPt","Vis Pt;Jet;Pt",10,0,10,1210,-1,1200)
	hist_VisPtHVQ = self.makeTH2F("hist_VisPtHVQ","Vis Pt HVQ;Jet;Pt",10,0,10,6010,-1,600)
	
	hist_2d_fracPartsHV_fracPtInv = self.makeTH2F("hist_2d_fracPartsHV_fracPtHV","% Invisible Particles vs % Pt from Invisible Particles;% Inv Parts;% Inv Pt",100,-0.01,1.01,100,-0.01,1.01)

	hist_iJetvsnParts = self.makeTH2F("hist_iJetvsnParts","All Particles;iJet;number Particles",10,0,10,28,0,28)
	hist_iJetvsnSMParts = self.makeTH2F("hist_iJetvsnSMParts","SM Particles;iJet;number Particles",10,0,10,28,0,28)
	hist_iJetvsnHVParts = self.makeTH2F("hist_iJetvsnHVParts","HV Particles;iJet;number Particles",10,0,10,28,0,28)

	hist_JetEta = self.makeTH1F("hist_JetEta","Jet Eta;eta;count",100,-5.0,5.0)
	hist_JetPhi = self.makeTH1F("hist_JetPhi","Jet Phi;phi;count",100,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_JetEta_vs_fracVisParts = self.makeTH2F("hist_JetEta_vs_FracVisParts","Jet Eta vs Frac Vis Parts;eta;% of particles that are visible",100,-5.0,5.0,100,0,1)
	#coarse grading for optimal MT resolution for fracPt cut
	#cutFractions = [x*0.01 for x in range(0,100)]
	#histList_MTcut = []
	#hist_MTcut = self.makeTH1F("hist_MTcut","'true' MT;MT;count/a.u.",100,0,4000)
	#for cutVal in cutFractions:
	#	histList_MTcut.append(self.makeTH1F("hist_MT2jet_"+str(cutVal),"2jet MT;MT;count/a.u.",100,0,4000))

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print(str(iEvent)+"/"+str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1:
			continue
		# set some commonly used things
		Jets = tree.JetsAK8
		nJets = len(Jets)
		met = tree.MET
		metPhi = tree.METPhi
		mt12 = trans_mass_Njet(Jets[:2], met, metPhi)
		
		hist_MT_2jets_all.Fill(mt12)
		if nJets > 2:
			mt123 = trans_mass_Njet(Jets[:3], met, metPhi)
			hist_deltaMT_2jet3jet.Fill(mt123-mt12)
			hist_MT_3jets_all.Fill(mt123)
			if mt123 > 3000.:
				hist_MT_MTcut.Fill(mt12)
			else:
				hist_MT_MTcut.Fill(mt123)
		else:
			hist_MT_3jets_all.Fill(mt12)
			hist_MT_MTcut.Fill(mt12)
		#if nJets != 4: #for checking specific nJet cases
		#	continue
		for iJet in range(nJets):
			hist_JetEta.Fill(tree.JetsAK8[iJet].Eta())
			hist_JetPhi.Fill(tree.JetsAK8[iJet].Phi())
			numParts = tree.numHVPartsInJet[iJet]+tree.numSMPartsInJet[iJet]
			hist_iJetvsnParts.Fill(iJet, numParts)
			hist_iJetvsnSMParts.Fill(iJet, tree.numSMPartsInJet[iJet])
			hist_iJetvsnHVParts.Fill(iJet, tree.numHVPartsInJet[iJet])
			try:
				fracHVParts = float(tree.numHVPartsInJet[iJet])/numParts
				hist_JetEta_vs_fracVisParts.Fill(tree.JetsAK8[iJet].Eta(),float(tree.numSMPartsInJet[iJet])/numParts)
				hist_2d_fracPartsHV_fracPtInv.Fill(fracHVParts,tree.fracTotPTfromInvHVQ[iJet])
			except ZeroDivisionError:
				print("Jet {} in Event {} has {} visible particles and {} invisible particles.".format(iJet, iEvent, tree.numSMPartsInJet[iJet], tree.numHVPartsInJet[iJet]))
			hist_VisPt.Fill(iJet, tree.pt_Vis[iJet])
			hist_VisPtHVQ.Fill(iJet, tree.pt_VisHVQ[iJet])
			hist_fracVisPTfromVisHVQ.Fill(iJet,tree.fracVisPTfromVisHVQ[iJet])
			hist_fracTotPTfromAllHVQ.Fill(iJet,tree.fracTotPTfromAllHVQ[iJet])
			hist_fracTotPTfromVisHVQ.Fill(iJet,tree.fracTotPTfromVisHVQ[iJet])
			hist_fracTotPTfromInvHVQ.Fill(iJet,tree.fracTotPTfromInvHVQ[iJet])
			hist_fracTotPTfromVis.Fill(iJet,tree.fracTotPTfromVis[iJet])
			hist_fracVisHVQtoInvHVQ.Fill(iJet,tree.fracVisHVQtoInvHVQ[iJet])

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


