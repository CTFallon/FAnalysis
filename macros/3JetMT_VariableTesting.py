from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array
from sys import exit

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

	# Define parameters for the test variable
	varName = "gamma2overgamma1" # dont forget to change the value of varVal in the loop!
	varLo = 1
	varHi = 10
	varSplit = 300

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	# Need five sets of DKL:
		# 1 - every event
		# 2 - Events with jet 1 as MDP
		# 3 - Events with jet 2 as MDP
		# 4 - Events with jet 3 as MDP
		# 5 - Events with jet 4(+) as MDP

	#New Set here
	hist_MT_HV = self.makeTH1F("hist_MT_HV",
			"MT HV Events;MT;count",100,0,4000)
	hist_FOM_DKL = self.makeTH1F("hist_"+varName+"_DKL",
			"FOM All;"+varName+";value", varSplit, varLo, varHi)
	histList_MT = []
	for iBin in range(varSplit+1):
		histList_MT.append(rt.TH1F("histList_MT_{}".format(iBin),
			";;",100,0,4000))
	#New Set here
	hist_MT_HV_mdp1 = self.makeTH1F("hist_MT_HV_mdp1",
			"MT HV Events;MT;count",100,0,4000)
	hist_FOM_DKL_mdp1 = self.makeTH1F("hist_"+varName+"_DKL_mdp1",
			"FOM MDP1;"+varName+";value", varSplit, varLo, varHi)
	histList_MT_mdp1 = []
	for iBin in range(varSplit+1):
		histList_MT_mdp1.append(rt.TH1F("histList_MT_mdp1_{}".format(iBin),
			";;",100,0,4000))
	#New Set here
	hist_MT_HV_mdp2 = self.makeTH1F("hist_MT_HV_mdp2",
			"MT HV Events;MT;count",100,0,4000)
	hist_FOM_DKL_mdp2 = self.makeTH1F("hist_"+varName+"_DKL_mdp2",
			"FOM MDP2;"+varName+";value", varSplit, varLo, varHi)
	histList_MT_mdp2 = []
	for iBin in range(varSplit+1):
		histList_MT_mdp2.append(rt.TH1F("histList_MT_mdp2_{}".format(iBin),
			";;",100,0,4000))
	#New Set here
	hist_MT_HV_mdp3 = self.makeTH1F("hist_MT_HV_mdp3",
			"MT HV Events;MT;count",100,0,4000)
	hist_FOM_DKL_mdp3 = self.makeTH1F("hist_"+varName+"_DKL_mdp3",
			"FOM MDP3;"+varName+";value", varSplit, varLo, varHi)
	histList_MT_mdp3 = []
	for iBin in range(varSplit+1):
		histList_MT_mdp3.append(rt.TH1F("histList_MT_mdp3_{}".format(iBin),
			";;",100,0,4000))
	#New Set here
	hist_MT_HV_mdp4 = self.makeTH1F("hist_MT_HV_mdp4",
			"MT HV Events;MT;count",100,0,4000)
	hist_FOM_DKL_mdp4 = self.makeTH1F("hist_"+varName+"_DKL_mdp4",
			"FOM MDP4+;"+varName+";value", varSplit, varLo, varHi)
	histList_MT_mdp4 = []
	for iBin in range(varSplit+1):
		histList_MT_mdp4.append(rt.TH1F("histList_MT_mdp4_{}".format(iBin),
			";;",100,0,4000))

	losBins = 300
	losLo = 1
	losHi = 10
	losName = "gamma2overgamma1" # dont forget to change LOS in the loop
	hist_LOS_all = self.makeTH1F("hist_"+losName+"_all",losName + " - All Events;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_dijet = self.makeTH1F("hist_"+losName+"_dijet",losName + " - Dijet Events;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_trijet = self.makeTH1F("hist_"+losName+"_trijet",losName + " - Trijet Events;"+losName+";Count",
		losBins, losLo, losHi)

	hist_LOS_all_mdp1 = self.makeTH1F("hist_"+losName+"_all_mdp1",losName + " - All Events - MDP1;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_dijet_mdp1 = self.makeTH1F("hist_"+losName+"_dijet_mdp1",losName + " - Dijet Events - MDP1;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_trijet_mdp1 = self.makeTH1F("hist_"+losName+"_trijet_mdp1",losName + " - Trijet Events - MDP1;"+losName+";Count",
		losBins, losLo, losHi)

	hist_LOS_all_mdp2 = self.makeTH1F("hist_"+losName+"_all_mdp2",losName + " - All Events - MDP2;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_dijet_mdp2 = self.makeTH1F("hist_"+losName+"_dijet_mdp2",losName + " - Dijet Events - MDP2;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_trijet_mdp2 = self.makeTH1F("hist_"+losName+"_trijet_mdp2",losName + " - Trijet Events - MDP2;"+losName+";Count",
		losBins, losLo, losHi)

	hist_LOS_all_mdp3 = self.makeTH1F("hist_"+losName+"_all_mdp3",losName + " - All Events - MDP3;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_dijet_mdp3 = self.makeTH1F("hist_"+losName+"_dijet_mdp3",losName + " - Dijet Events - MDP3;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_trijet_mdp3 = self.makeTH1F("hist_"+losName+"_trijet_mdp3",losName + " - Trijet Events - MDP3;"+losName+";Count",
		losBins, losLo, losHi)

	hist_LOS_all_mdp4 = self.makeTH1F("hist_"+losName+"_all_mdp4",losName + " - All Events - MDP4;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_dijet_mdp4 = self.makeTH1F("hist_"+losName+"_dijet_mdp4",losName + " - Dijet Events - MDP4;"+losName+";Count",
		losBins, losLo, losHi)
	hist_LOS_trijet_mdp4 = self.makeTH1F("hist_"+losName+"_trijet_mdp4",losName + " - Trijet Events - MDP4;"+losName+";Count",
		losBins, losLo, losHi)

	#hist2d_ptRatio_vs_deltaR_all = self.makeTH2F("hist2d_ptRatio_vs_deltaR_all",
	#	"Pt Ratio vs DeltaR; ptRatio mdp; deltaR13",100,0,1,100,0,5)
	#hist2d_ptRatio_vs_deltaR_dijet = self.makeTH2F("hist2d_ptRatio_vs_deltaR_dijet",
	#	"Pt Ratio vs DeltaR; ptRatio mdp; deltaR13",100,0,1,100,0,5)
	#hist2d_ptRatio_vs_deltaR_trijet = self.makeTH2F("hist2d_ptRatio_vs_deltaR_trijet",
	#	"Pt Ratio vs DeltaR; ptRatio mdp; deltaR13",100,0,1,100,0,5)
	
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15

		# trying different cuts depending on what jet is furhtest from MET
		#if tree.iJetMaxDeltaPhi != 2:
		#	continue

		try:
			varVal = tree.JetsAK8[1].Gamma()/tree.JetsAK8[0].Gamma()
		except IndexError:
			continue

		if len(tree.JetsAK8) >2:
			#dis12 = tree.JetsAK8[0].DeltaR(tree.JetsAK8[1])
			#dis23 = tree.JetsAK8[1].DeltaR(tree.JetsAK8[2])
			dis31 = tree.JetsAK8[2].DeltaR(tree.JetsAK8[0])
			threeJets = [tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]]
			#skip events that get caught by other cuts
			if (tree.iJetMaxDeltaPhi == 1) or (deltaPhi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()) < 2.65) or (tree.JetsAK8[2].Gamma() < 5.24):
				continue
			ptRatio = tree.pTMaxDeltaPhi/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())
			#losRatio = tree.JetsAK8[0].Pt()/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())
			losRatio = tree.JetsAK8[1].Gamma()/tree.JetsAK8[0].Gamma()
			hist_LOS_all.Fill(losRatio)
			#hist2d_ptRatio_vs_deltaR_all.Fill(ptRatio,dis31)
			if bool(tree.JetsAK8_isHV[2]):
				hist_LOS_trijet.Fill(losRatio)
				#hist2d_ptRatio_vs_deltaR_trijet.Fill(ptRatio,dis31)
			else:
				hist_LOS_dijet.Fill(losRatio)
				#hist2d_ptRatio_vs_deltaR_dijet.Fill(ptRatio,dis31)
			if tree.iJetMaxDeltaPhi == 0:
				hist_LOS_all_mdp1.Fill(losRatio)
				if bool(tree.JetsAK8_isHV[2]):
					hist_LOS_trijet_mdp1.Fill(losRatio)
				else:
					hist_LOS_dijet_mdp1.Fill(losRatio)
			elif tree.iJetMaxDeltaPhi == 1:
				hist_LOS_all_mdp2.Fill(losRatio)
				if bool(tree.JetsAK8_isHV[2]):
					hist_LOS_trijet_mdp2.Fill(losRatio)
				else:
					hist_LOS_dijet_mdp2.Fill(losRatio)
			elif tree.iJetMaxDeltaPhi == 2:
				hist_LOS_all_mdp3.Fill(losRatio)
				if bool(tree.JetsAK8_isHV[2]):
					hist_LOS_trijet_mdp3.Fill(losRatio)
				else:
					hist_LOS_dijet_mdp3.Fill(losRatio)
			else:
				hist_LOS_all_mdp4.Fill(losRatio)
				if bool(tree.JetsAK8_isHV[2]):
					hist_LOS_trijet_mdp4.Fill(losRatio)
				else:
					hist_LOS_dijet_mdp4.Fill(losRatio)

		FSRJets = []
		for iJet in range(len(tree.JetsAK8)):
			if bool(tree.JetsAK8_isHV[iJet]):
				FSRJets.append(tree.JetsAK8[iJet])
		FSR_MT = trans_mass_Njet(FSRJets, tree.MET, tree.METPhi)
		hist_MT_HV.Fill(FSR_MT)
		if tree.iJetMaxDeltaPhi == 0:
			hist_MT_HV_mdp1.Fill(FSR_MT)
		if tree.iJetMaxDeltaPhi == 1:
			hist_MT_HV_mdp2.Fill(FSR_MT)
		if tree.iJetMaxDeltaPhi == 2:
			hist_MT_HV_mdp3.Fill(FSR_MT)
		if tree.iJetMaxDeltaPhi >= 3:
			hist_MT_HV_mdp4.Fill(FSR_MT)

		# next, fill the MT histograms for various levels of cutting on the varible
		varBin = int(varVal/(varHi-varLo)*varSplit)
		fillListOfMtCuts(histList_MT,tree,varBin)
		if tree.iJetMaxDeltaPhi == 0:
			fillListOfMtCuts(histList_MT_mdp1,tree,varBin)
		if tree.iJetMaxDeltaPhi == 1:
			fillListOfMtCuts(histList_MT_mdp2,tree,varBin)
		if tree.iJetMaxDeltaPhi == 2:
			fillListOfMtCuts(histList_MT_mdp3,tree,varBin)
		if tree.iJetMaxDeltaPhi >= 3:
			fillListOfMtCuts(histList_MT_mdp4,tree,varBin)

	populateDKL(histList_MT, hist_MT_HV, hist_FOM_DKL)
	populateDKL(histList_MT_mdp1, hist_MT_HV_mdp1, hist_FOM_DKL_mdp1)
	populateDKL(histList_MT_mdp2, hist_MT_HV_mdp2, hist_FOM_DKL_mdp2)
	populateDKL(histList_MT_mdp3, hist_MT_HV_mdp3, hist_FOM_DKL_mdp3)
	populateDKL(histList_MT_mdp4, hist_MT_HV_mdp4, hist_FOM_DKL_mdp4)

	#self.makePng([hist_FOM_DKL,hist_FOM_DKL_mdp1,hist_FOM_DKL_mdp2,hist_FOM_DKL_mdp3,hist_FOM_DKL_mdp4],"varTest_"+varName)
	#self.makePng([hist_LOS_all,hist_LOS_dijet,hist_LOS_trijet],"varTest_"+losName, doLeg = False)
	#self.makePng([hist_LOS_all_mdp1,hist_LOS_dijet_mdp1,hist_LOS_trijet_mdp1],"varTest_"+losName+"_mdp1", doLeg = False)
	#self.makePng([hist_LOS_all_mdp2,hist_LOS_dijet_mdp2,hist_LOS_trijet_mdp2],"varTest_"+losName+"_mdp2", doLeg = False)
	#self.makePng([hist_LOS_all_mdp3,hist_LOS_dijet_mdp3,hist_LOS_trijet_mdp3],"varTest_"+losName+"_mdp3", doLeg = False)
	#self.makePng([hist_LOS_all_mdp4,hist_LOS_dijet_mdp4,hist_LOS_trijet_mdp4],"varTest_"+losName+"_mdp4", doLeg = False)

	#c1 = rt.TCanvas("c1","c1",1200,900)
	#hist2d_ptRatio_vs_deltaR_all.Draw("colz")
	#c1.SaveAs(self.extraDir+"2d_ptRatiovsDetlaR13_all.png")
	#hist2d_ptRatio_vs_deltaR_dijet.Draw("colz")
	#c1.SaveAs(self.extraDir+"2d_ptRatiovsDetlaR13_dijet.png")
	#hist2d_ptRatio_vs_deltaR_trijet.Draw("colz")
	#c1.SaveAs(self.extraDir+"2d_ptRatiovsDetlaR13_trijet.png")
	

def addLoop():
	baseClass.loop = loop

def makePlots(LoH, name):
	c = rt.TCanvas("c1","c1",1200,900)
	stack = rt.THStack()
	for i in range(len(LoH)):
		LoH[i].SetLineColor(i+1)
		stack.Add(LoH[i])
	stack.Draw("nostack")
	c.BuildLegend(0.6,0.2,0.8,0.4)
	c.SaveAs(self.extraDir+name+".png")

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def Dkl(A, B):
	# takes two histograms with same binning,
	# and computes the Kullback-Liebler Divergence
	# defined as sum over bins (i) of A[i]*log(A[i]/B[i])
	res = 0.
	for i in range(A.GetNbinsX()):
		Ai = A.GetBinContent(i)
		Bi = B.GetBinContent(i)
		if Bi > 0. and Ai > 0.:
			res += Ai*rt.TMath.Log10(Ai/Bi)
	return res


def fillListOfMtCuts(hL, tR, vB):
	for histo in hL:
		mt = 0.
		if len(tR.JetsAK8) == 2:
			mt = tR.MT_AK8
		else:
			histCut = int(histo.GetName().split("_")[-1])
			if vB > histCut:
				mt = trans_mass_Njet([tR.JetsAK8[0],tR.JetsAK8[1],tR.JetsAK8[2]], tR.MET, tR.METPhi)
			else:
				mt = tR.MT_AK8
		histo.Fill(mt)

def populateDKL(hL, mtHist, dklHist):
	for histo in hL:
		histBin = int(histo.GetName().split("_")[-1])
		dkl = mtHist.KolmogorovTest(histo,"N")
		dklHist.SetBinContent(histBin, dkl)
	maxVal = dklHist.GetMaximum()
	for iBin in range(1,dklHist.GetNbinsX()+1):
		try:
			dklHist.SetBinContent(iBin, dklHist.GetBinContent(iBin)/maxVal)	
		except ZeroDivisionError:
			dklHist.SetBinContent(iBin, 0)
	print("Maximum Bin is #" + str(dklHist.GetMaximumBin()))

def LOS(a, b, c):
	var = (b**2+c**2-a**2)/(2*b*c)
	ang = rt.TMath.ACos(var)
	den = rt.TMath.Sin(ang)
	return a/den

def deltaPhi(phi1, phi2):
	delta = phi1 - phi2
	while delta > rt.TMath.Pi():
		delta -= 2*rt.TMath.Pi()
	while delta < -rt.TMath.Pi():
		delta += 2*rt.TMath.Pi()
	return abs(delta)
	
	
