import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)


def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutations_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

def trans_mass_Njet(jets, met, metPhi):
	visible = rt.TLorentzVector()
	for jet in jets:
		visible += jet
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))

def drawHistos(histo, plotName, logY = False, stack = False):
	canvas = rt.TCanvas("canvas","canvas",900,600)
	histStack = rt.THStack()
	histStack.SetTitle(histo[0].GetTitle())
	for hist in histo:
		hist.SetLineColor(histo.index(hist)+1)
		histStack.Add(hist)
		print("Stats: " + hist.GetName() + " " + str(hist.GetMean()) + " " + str(hist.GetRMS()))
	if logY == True:
		canvas.SetLogy()
	histStack.SetTitle(plotName)
	if not stack:
		histStack.Draw("NOSTACK")
	else:
		histStack.Draw("")
	histStack.GetXaxis().SetTitle(histo[0].GetXaxis().GetTitle())
	if len(histo) > 1:
		canvas.BuildLegend(0.75,0.75,0.9,0.9,"")
	canvas.SaveAs(createImageName(mZprime, mDark, rinv, alpha, plotName))

def pseudoAK(jetList, R):
	endJetList = []
	#print("Doing Reclustering...")
	while len(jetList) > 0:
		nJets = len(jetList)
		#print(nJets)
		dMat = [[1000 for i in range(nJets)] for j in range(nJets)]
		dMatMin, iMin, jMin = 1000., 10, 10
		for i in range(nJets):
			for j in range(nJets):
				if i == j:
					dMat[i][j] = 1/(jetList[i].Pt())**2
				elif i > j:
					dMat[i][j] = min(1/(jetList[i].Pt())**2,1/(jetList[j].Pt())**2)*(jetList[i].DeltaR(jetList[j]))/R
				if dMatMin > dMat[i][j]:
					dMatMin, iMin, jMin = dMat[i][j], i, j
		if iMin != jMin:
			newJet = jetList[iMin] + jetList[jMin]
			jetList.pop(iMin)
			jetList.pop(jMin)
			jetList.append(newJet)
		if iMin == jMin:
			endJetList.append(jetList.pop(iMin))
	return endJetList

def drawEventEtaPhiPlot(jetCollectionAK8, jetCollectionCA11, partCol, particlePDGID,METPhi, isFromHVQuark, plotNumber):
	canv = rt.TCanvas("canv","canv",1600,800)
	histAxis = rt.TH2F("axisHsito", ";\eta;\phi",100,-6,6,100,-rt.TMath.Pi(),rt.TMath.Pi())
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionCA11)):
		objectList.append(rt.TEllipse(jetCollectionCA11[iJet].Eta(), jetCollectionCA11[iJet].Phi(), 1.1, 1.1))
		objectList[-1].SetLineStyle(2)
		objectList[-1].SetLineColor(iJet+1)
		objectList[-1].SetLineWidth(2)
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
	canv.SaveAs("etaPhi_"+str(plotNumber)+".png")


#step1, know what data set we're looking at

mZprime = 3000
mDark = '20'
rinv = '0p3'
alpha = '0p2'
f_rvis = 1.0 - 0.3

inFile = rt.TFile("../../"+createInFileName(mZprime, mDark,rinv,alpha),"read")
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
nEventsWith2OrMoreJets = 0
nEventsWithLeadingTwoJetsHavingPtGreaterThan170 = 0
nEventsWithMETMTRatioGreaterThanp15 = 0
nEventsPassedPreSelection = 0


# initilize histograms

hist_AK8DiJetMt = rt.TH1F("AK8DiJetMt",";M;count/a.u.",100,0,5000)
hist_CA11DiJetMt = rt.TH1F("CA11DiJetMt",";M;count/a.u.",100,0,5000)
hist_HVQuarksMt = rt.TH1F("HVQuarksMt",";M;count/a.u.",100,0,5000)
hist_ZprimeMt = rt.TH1F("ZprimeMt",";M;count/a.u.",100,0,5000)
hist_AK8DiJetMjj = rt.TH1F("AK8DiJetMjj",";M;count/a.u.",100,0,5000)
hist_CA11DiJetMjj = rt.TH1F("CA11DiJetMjj",";M;count/a.u.",100,0,5000)

hist_AK8AllJetsMt = rt.TH1F("AK8AllJetMt",";M;count/a.u.",100,0,5000)
hist_CA11AllJetsMt = rt.TH1F("CA11AllJetMt",";M;count/a.u.",100,0,5000)

hist_JetsWithHVParts = rt.TH1F("JetsWithHVParts",";JetDict;count/a.u.", 8, 0, 8)

hist_JetsWithHVPartsvsNjets = rt.TH2F("2D_HVJetsvsnJets",";HVJets;nJets", 8, 0, 8, 10, 0, 10)

hist_nAK8_vs_nCA11 = rt.TH2F("2D_AK8vsCA11",";nAK8;nCA11",10,0,10,10,0,10)
hist_nAK8 = rt.TH1F("nAK8",";nJets;",10,0,10)
hist_nCA11 = rt.TH1F("nCA11",";nJets;",10,0,10)

hist_SDMass_withHV = rt.TH1F("SDMass_HV",";Mass_{SD};count/a.u.", 100, 0, 800)
hist_SDMass_withoutHV = rt.TH1F("SDMass_noHV",";Mass_{SD};count/a.u.", 100, 0, 800)

hist_deltaPhi1 = rt.TH1F("dPhi1", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi2 = rt.TH1F("dPhi2", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi3 = rt.TH1F("dPhi3", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_deltaPhi1_HV = rt.TH1F("dPhi1_HV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi2_HV = rt.TH1F("dPhi2_HV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi3_HV = rt.TH1F("dPhi3_HV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_deltaPhi1_noHV = rt.TH1F("dPhi1_noHV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi2_noHV = rt.TH1F("dPhi2_noHV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi3_noHV = rt.TH1F("dPhi3_noHV", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_deltaPhi1_3jetMT = rt.TH1F("dPhi1_3jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi2_3jetMT = rt.TH1F("dPhi2_3jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi3_3jetMT = rt.TH1F("dPhi3_3jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_deltaPhi1_2jetMT = rt.TH1F("dPhi1_2jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi2_2jetMT = rt.TH1F("dPhi2_2jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_deltaPhi3_2jetMT = rt.TH1F("dPhi3_2jetMT", ";dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_deltaR_12_3jetMT = rt.TH1F("dR_12_3jetMT", ";dR;count/a.u.", 100, 0, 6)
hist_deltaR_13_3jetMT = rt.TH1F("dR_13_3jetMT", ";dR;count/a.u.", 100, 0, 6)
hist_deltaR_23_3jetMT = rt.TH1F("dR_23_3jetMT", ";dR;count/a.u.", 100, 0, 6)

hist_deltaR_12_2jetMT = rt.TH1F("dR_12_2jetMT", ";dR;count/a.u.", 100, 0, 6)
hist_deltaR_13_2jetMT = rt.TH1F("dR_13_2jetMT", ";dR;count/a.u.", 100, 0, 6)
hist_deltaR_23_2jetMT = rt.TH1F("dR_13_2jetMT", ";dR;count/a.u.", 100, 0, 6)

hist2d_deltaR_1_3jetMT = rt.TH2F("dr_1X_3jetMT",";dR12;dR13", 100, 0, 6, 100, 0, 6)
hist2d_deltaR_1_2jetMT = rt.TH2F("dr_1X_2jetMT",";dR12;dR13", 100, 0, 6, 100, 0, 6)

hist2d_deltaR_2_3jetMT = rt.TH2F("dr_2X_3jetMT",";dR21;dR23", 100, 0, 6, 100, 0, 6)
hist2d_deltaR_2_2jetMT = rt.TH2F("dr_2X_2jetMT",";dR21;dR23", 100, 0, 6, 100, 0, 6)

hist2d_deltaR_3_3jetMT = rt.TH2F("dr_3X_3jetMT",";dR31;dR32", 100, 0, 6, 100, 0, 6)
hist2d_deltaR_3_2jetMT = rt.TH2F("dr_3X_2jetMT",";dR31;dR32", 100, 0, 6, 100, 0, 6)

hist2d_SDmass_12_3jetMT = rt.TH2F("mSD_12_3jetMT",";mSD1;mSD2", 100, 0, 700, 100, 0, 400)
hist2d_SDmass_12_2jetMT = rt.TH2F("mSD_12_2jetMT",";mSD1;mSD2", 100, 0, 700, 100, 0, 400)
hist2d_SDmass_13_3jetMT = rt.TH2F("mSD_13_3jetMT",";mSD1;mSD3", 100, 0, 700, 100, 0, 300)
hist2d_SDmass_13_2jetMT = rt.TH2F("mSD_13_2jetMT",";mSD1;mSD3", 100, 0, 700, 100, 0, 300)
hist2d_SDmass_23_3jetMT = rt.TH2F("mSD_23_3jetMT",";mSD2;mSD3", 100, 0, 400, 100, 0, 300)
hist2d_SDmass_23_2jetMT = rt.TH2F("mSD_23_2jetMT",";mSD2;mSD3", 100, 0, 400, 100, 0, 300)

hist_sdVar_12 = rt.TH1F("sdVar_12",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_13 = rt.TH1F("sdVar_13",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_23 = rt.TH1F("sdVar_23",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_3Jets = rt.TH1F("sdVar_3jets",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)

hist_sdVar_12_3jetMT = rt.TH1F("sdVar_12_3jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_13_3jetMT = rt.TH1F("sdVar_13_3jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_23_3jetMT = rt.TH1F("sdVar_23_3jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_3Jets_3jetMT = rt.TH1F("sdVar_3jets_3jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)

hist_sdVar_12_2jetMT = rt.TH1F("sdVar_12_2jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_13_2jetMT = rt.TH1F("sdVar_13_2jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_23_2jetMT = rt.TH1F("sdVar_23_2jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_3Jets_2jetMT = rt.TH1F("sdVar_3jets_2jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)


hist2d_sdVardR_13 = rt.TH2F("sdVardR_13",";min(p_T_i)/\Sigma(p_T_i);dR",100,0,0.5,100,0,6)
hist2d_sdVardR_13_3jetMT = rt.TH2F("sdVardR_13_3jetMT",";min(p_T_i)/\Sigma(p_T_i);dR",100,0,0.5,100,0,6)
hist2d_sdVardR_13_2jetMT = rt.TH2F("sdVardR_13_2jetMT",";min(p_T_i)/\Sigma(p_T_i);dR",100,0,0.5,100,0,6)

hist_jetPt_mindPhi = rt.TH1F("jetPt_mindPhi",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_middPhi = rt.TH1F("jetPt_middPhi",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_maxdPhi = rt.TH1F("jetPt_maxdPhi",";pT;count/a.u.", 100, 0, 2000)

hist_dPhi_mindPhi = rt.TH1F("dPhi_mindPhi",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_middPhi = rt.TH1F("dPhi_middPhi",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_maxdPhi = rt.TH1F("dPhi_maxdPhi",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())

hist_jetPt_mindPhi_3jetMT = rt.TH1F("jetPt_mindPhi_3jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_middPhi_3jetMT = rt.TH1F("jetPt_middPhi_3jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_maxdPhi_3jetMT = rt.TH1F("jetPt_maxdPhi_3jetMT",";pT;count/a.u.", 100, 0, 2000)

hist_dPhi_mindPhi_3jetMT = rt.TH1F("dPhi_mindPhi_3jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_middPhi_3jetMT = rt.TH1F("dPhi_middPhi_3jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_maxdPhi_3jetMT = rt.TH1F("dPhi_maxdPhi_3jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())

hist_jetPt_mindPhi_2jetMT = rt.TH1F("jetPt_mindPhi_2jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_middPhi_2jetMT = rt.TH1F("jetPt_middPhi_2jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_maxdPhi_2jetMT = rt.TH1F("jetPt_maxdPhi_2jetMT",";pT;count/a.u.", 100, 0, 2000)

hist_dPhi_mindPhi_2jetMT = rt.TH1F("dPhi_mindPhi_2jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_middPhi_2jetMT = rt.TH1F("dPhi_middPhi_2jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_dPhi_maxdPhi_2jetMT = rt.TH1F("dPhi_maxdPhi_2jetMT",";dPhi;count/a.u.", 100, 0, rt.TMath.Pi())

hist_deltadPhi_mindPhimiddPhi = rt.TH1F("deltadPhi_minmid",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_middPhimaxdPhi = rt.TH1F("deltadPhi_midmax",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_mindPhimaxdPhi = rt.TH1F("deltadPhi_minmax",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())

hist_deltadPhi_mindPhimiddPhi_3jetMT = rt.TH1F("deltadPhi_minmid_3jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_middPhimaxdPhi_3jetMT = rt.TH1F("deltadPhi_midmax_3jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_mindPhimaxdPhi_3jetMT = rt.TH1F("deltadPhi_minmax_3jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())

hist_deltadPhi_mindPhimiddPhi_2jetMT = rt.TH1F("deltadPhi_minmid_2jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_middPhimaxdPhi_2jetMT = rt.TH1F("deltadPhi_midmax_2jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())
hist_deltadPhi_mindPhimaxdPhi_2jetMT = rt.TH1F("deltadPhi_minmax_2jetMT",";deltadPhi;count/a.u.", 100, 0, rt.TMath.Pi())





histList = []
histList.append(hist_AK8DiJetMt)
histList.append(hist_CA11DiJetMt)
histList.append(hist_HVQuarksMt)
histList.append(hist_ZprimeMt)
histList.append(hist_AK8DiJetMjj)
histList.append(hist_CA11DiJetMjj)
#histList.append(hist_JetsWithHVParts)
histList.append(hist_AK8AllJetsMt)
histList.append(hist_CA11AllJetsMt)
histList.append(hist_nAK8)
histList.append(hist_nCA11)
histList.append(hist_SDMass_withHV)
histList.append(hist_SDMass_withoutHV)
histList.append(hist_deltaPhi1)
histList.append(hist_deltaPhi2)
histList.append(hist_deltaPhi3)
histList.append(hist_deltaPhi1_HV)
histList.append(hist_deltaPhi2_HV)
histList.append(hist_deltaPhi3_HV)
histList.append(hist_deltaPhi1_noHV)
histList.append(hist_deltaPhi2_noHV)
histList.append(hist_deltaPhi3_noHV)
histList.append(hist_jetPt_mindPhi)
histList.append(hist_jetPt_middPhi)
histList.append(hist_jetPt_maxdPhi)

nPlotsMade = 0


mat_JetsHVParts_nJets = [[0 for x in range(32)] for y in range(10)]

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	jetCollection = tree.JetsAK8
	nAK8Jets = len(jetCollection)

	#some calculations to determine which histos to fill
	# PreSelection Cuts
	if not (len(jetCollection)==3): #moreThan2Jets, temp set to ==3 for 3JetMt classification
		continue
	nEventsWith2OrMoreJets += 1
	if not ((jetCollection[0].Pt() > 170.0) and (jetCollection[1].Pt() > 170.0)): #leadingJetsPtOver170
		continue
	nEventsWithLeadingTwoJetsHavingPtGreaterThan170 += 1
	if not (tree.MET/tree.MT_AK8 > 0.15): #METoverMTgreaterThan0p15
		continue
	nEventsWithMETMTRatioGreaterThanp15 += 1
	if not ((len(tree.Electrons) + len(tree.Muons)) == 0): #leptonNotPresent
		continue
	nEventsPassedPreSelection += 1

	jetsWithHVDecendants = ["0","0","0","0","0"]
	#make vector of length nGenParts that is 1 if the particle came from a HV quark
	isFromHVQuark = [0 for x in range(len(tree.GenParticles))]
	listOfHVQuarks = []
	for iPart in range(len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if abs(tree.GenParticles_PdgId[iPart]) == 4900101:
			listOfHVQuarks.append(tree.GenParticles[iPart])
		if tree.GenParticles_PdgId[iPart] == 4900023:
			hist_ZprimeMt.Fill(trans_mass_Njet([tree.GenParticles[iPart]], 0, tree.METPhi))
		if iParent >= iPart:
			print("Ut-oh, the parent has a higher index than the child...")
		#print(str(tree.GenParticles_PdgId[iParent]) + " (" + str(iParent) + ") --> " + str(tree.GenParticles_PdgId[iPart])+ " (" + str(iPart) + ")")
		if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or (isFromHVQuark[iParent]):
			isFromHVQuark[iPart] = 1
		#finding what Jet a particle is in, only if it decends from a HVQuark:
		if isFromHVQuark[iPart]:
			for iJet in range(min(len(tree.JetsAK8),5)):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8:
					jetsWithHVDecendants[-iJet-1] = "1"
	mat_JetsHVParts_nJets[nAK8Jets][int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2)] += 1

	hist_JetsWithHVParts.Fill(int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2))
	hist_JetsWithHVPartsvsNjets.Fill(int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2), nAK8Jets)
	#print(isFromHVQuark)
	# figureout which order the jets rae in distance from METPhi
	deltaPhiList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
	IdxList = [0,1,2]
	mindPhiIdx = deltaPhiList.index(min(deltaPhiList))
	IdxList.remove(mindPhiIdx)
	maxdPhiIdx = deltaPhiList.index(max(deltaPhiList))
	IdxList.remove(maxdPhiIdx)
	middPhiIdx = IdxList[0]
	
	hist_jetPt_mindPhi.Fill(tree.JetsAK8[mindPhiIdx].Pt())
	hist_jetPt_middPhi.Fill(tree.JetsAK8[middPhiIdx].Pt())
	hist_jetPt_maxdPhi.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
	hist_dPhi_mindPhi.Fill(deltaPhiList[mindPhiIdx])
	hist_dPhi_middPhi.Fill(deltaPhiList[middPhiIdx])
	hist_dPhi_maxdPhi.Fill(deltaPhiList[maxdPhiIdx])
	hist_deltadPhi_mindPhimiddPhi.Fill(deltaPhiList[middPhiIdx]-deltaPhiList[mindPhiIdx])
	hist_deltadPhi_middPhimaxdPhi.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[middPhiIdx])
	hist_deltadPhi_mindPhimaxdPhi.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[mindPhiIdx])
	

	if nPlotsMade < 20 and len(tree.JetsAK8) >= 3:
		nPlotsMade += 1
		drawEventEtaPhiPlot(tree.JetsAK8, tree.JetsCA11, tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, isFromHVQuark, iEvent)

	for iJet in range(min(len(tree.JetsAK8), 5)):
		if int(jetsWithHVDecendants[-iJet-1]):
			hist_SDMass_withHV.Fill(tree.JetsAK8_softDropMass[iJet])
		else:
			hist_SDMass_withoutHV.Fill(tree.JetsAK8_softDropMass[iJet])
	sdVar = [min(tree.JetsAK8[0].Pt(),tree.JetsAK8[1].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()),min(tree.JetsAK8[0].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()),min(tree.JetsAK8[1].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())]
	sdVar3Jets = min(tree.JetsAK8[0].Pt(),tree.JetsAK8[1].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())
	
	hist_sdVar_12.Fill(sdVar[0])
	hist_sdVar_13.Fill(sdVar[1])
	hist_sdVar_23.Fill(sdVar[2])
	hist_sdVar_3Jets.Fill(sdVar3Jets)
	hist2d_sdVardR_13.Fill(sdVar[1],tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))

	if int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2) == 3:
		hist_deltaPhi1_2jetMT.Fill(tree.DeltaPhi1)
		hist_deltaPhi2_2jetMT.Fill(tree.DeltaPhi2)
		hist_deltaPhi3_2jetMT.Fill(tree.DeltaPhi3)
		hist_deltaR_12_2jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]))
		hist_deltaR_13_2jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
		hist_deltaR_23_2jetMT.Fill(tree.JetsAK8[1].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_1_2jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]),tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_2_2jetMT.Fill(tree.JetsAK8[1].DeltaR(tree.JetsAK8[0]),tree.JetsAK8[1].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_3_2jetMT.Fill(tree.JetsAK8[2].DeltaR(tree.JetsAK8[0]),tree.JetsAK8[2].DeltaR(tree.JetsAK8[1]))
		hist2d_SDmass_12_2jetMT.Fill(tree.JetsAK8_softDropMass[0],tree.JetsAK8_softDropMass[1])
		hist2d_SDmass_13_2jetMT.Fill(tree.JetsAK8_softDropMass[0],tree.JetsAK8_softDropMass[2])
		hist2d_SDmass_23_2jetMT.Fill(tree.JetsAK8_softDropMass[1],tree.JetsAK8_softDropMass[2])
		hist_sdVar_12_2jetMT.Fill(sdVar[0])
		hist_sdVar_13_2jetMT.Fill(sdVar[1])
		hist_sdVar_23_2jetMT.Fill(sdVar[2])
		hist_jetPt_mindPhi_2jetMT.Fill(tree.JetsAK8[mindPhiIdx].Pt())
		hist_jetPt_middPhi_2jetMT.Fill(tree.JetsAK8[middPhiIdx].Pt())
		hist_jetPt_maxdPhi_2jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
		hist_dPhi_mindPhi_2jetMT.Fill(deltaPhiList[mindPhiIdx])
		hist_dPhi_middPhi_2jetMT.Fill(deltaPhiList[middPhiIdx])
		hist_dPhi_maxdPhi_2jetMT.Fill(deltaPhiList[maxdPhiIdx])
		hist_deltadPhi_mindPhimiddPhi_2jetMT.Fill(deltaPhiList[middPhiIdx]-deltaPhiList[mindPhiIdx])
		hist_deltadPhi_middPhimaxdPhi_2jetMT.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[middPhiIdx])
		hist_deltadPhi_mindPhimaxdPhi_2jetMT.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[mindPhiIdx])
		hist_sdVar_3Jets_2jetMT.Fill(sdVar3Jets)
		hist2d_sdVardR_13_2jetMT.Fill(sdVar[1],tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
	else:
		hist_deltaPhi1_3jetMT.Fill(tree.DeltaPhi1)
		hist_deltaPhi2_3jetMT.Fill(tree.DeltaPhi2)
		hist_deltaPhi3_3jetMT.Fill(tree.DeltaPhi3)
		hist_deltaR_12_3jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]))
		hist_deltaR_13_3jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
		hist_deltaR_23_3jetMT.Fill(tree.JetsAK8[1].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_1_3jetMT.Fill(tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]),tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_2_3jetMT.Fill(tree.JetsAK8[1].DeltaR(tree.JetsAK8[0]),tree.JetsAK8[1].DeltaR(tree.JetsAK8[2]))
		hist2d_deltaR_3_3jetMT.Fill(tree.JetsAK8[2].DeltaR(tree.JetsAK8[0]),tree.JetsAK8[2].DeltaR(tree.JetsAK8[1]))
		hist2d_SDmass_12_3jetMT.Fill(tree.JetsAK8_softDropMass[0],tree.JetsAK8_softDropMass[1])
		hist2d_SDmass_13_3jetMT.Fill(tree.JetsAK8_softDropMass[0],tree.JetsAK8_softDropMass[2])
		hist2d_SDmass_23_3jetMT.Fill(tree.JetsAK8_softDropMass[1],tree.JetsAK8_softDropMass[2])
		hist_sdVar_12_3jetMT.Fill(sdVar[0])
		hist_sdVar_13_3jetMT.Fill(sdVar[1])
		hist_sdVar_23_3jetMT.Fill(sdVar[2])
		hist_jetPt_mindPhi_3jetMT.Fill(tree.JetsAK8[mindPhiIdx].Pt())
		hist_jetPt_middPhi_3jetMT.Fill(tree.JetsAK8[middPhiIdx].Pt())
		hist_jetPt_maxdPhi_3jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
		hist_dPhi_mindPhi_3jetMT.Fill(deltaPhiList[mindPhiIdx])
		hist_dPhi_middPhi_3jetMT.Fill(deltaPhiList[middPhiIdx])
		hist_dPhi_maxdPhi_3jetMT.Fill(deltaPhiList[maxdPhiIdx])
		hist_deltadPhi_mindPhimiddPhi_3jetMT.Fill(deltaPhiList[middPhiIdx]-deltaPhiList[mindPhiIdx])
		hist_deltadPhi_middPhimaxdPhi_3jetMT.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[middPhiIdx])
		hist_deltadPhi_mindPhimaxdPhi_3jetMT.Fill(deltaPhiList[maxdPhiIdx]-deltaPhiList[mindPhiIdx])
		hist_sdVar_3Jets_3jetMT.Fill(sdVar3Jets)
		hist2d_sdVardR_13_3jetMT.Fill(sdVar[1],tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]))
	
	hist_deltaPhi1.Fill(tree.DeltaPhi1)
	hist_deltaPhi2.Fill(tree.DeltaPhi2)
	hist_deltaPhi3.Fill(tree.DeltaPhi3)
	if jetsWithHVDecendants[4] == '1':
		hist_deltaPhi1_HV.Fill(tree.DeltaPhi1)
	else:
		hist_deltaPhi1_noHV.Fill(tree.DeltaPhi1)
	if jetsWithHVDecendants[3] == '1':
		hist_deltaPhi2_HV.Fill(tree.DeltaPhi2)
	else:
		hist_deltaPhi2_noHV.Fill(tree.DeltaPhi2)
	#if jetsWithHVDecendants[2] == '1':
	#	hist_deltaPhi3_HV.Fill(tree.DeltaPhi3)
	#else:
	#	hist_deltaPhi3_noHV.Fill(tree.DeltaPhi3)
		
			


	if len(listOfHVQuarks) == 2:
		hist_HVQuarksMt.Fill(trans_mass_Njet(listOfHVQuarks, tree.MET, tree.METPhi))
	hist_AK8DiJetMt.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
	hist_AK8DiJetMjj.Fill((tree.JetsAK8[0]+tree.JetsAK8[1]).M())
	#print((tree.JetsAK8[0]+tree.JetsAK8[1]).M())
	if len(tree.JetsCA11) >= 2:
		hist_CA11DiJetMt.Fill(trans_mass_Njet([tree.JetsCA11[0],tree.JetsCA11[1]], tree.MET, tree.METPhi))
		hist_CA11DiJetMjj.Fill((tree.JetsCA11[0]+tree.JetsCA11[1]).M())

	hist_nAK8_vs_nCA11.Fill(len(tree.JetsAK8), len(tree.JetsCA11))
	hist_nAK8.Fill(len(tree.JetsAK8))
	hist_nCA11.Fill(len(tree.JetsCA11))
	hist_AK8AllJetsMt.Fill(trans_mass_Njet(tree.JetsAK8, tree.MET, tree.METPhi))
	hist_CA11AllJetsMt.Fill(trans_mass_Njet(tree.JetsCA11, tree.MET, tree.METPhi))

#for histo in histList:
#	total = histo.GetEntries()
#	if total != 0:
#		histo.Scale(1/(total))

#setting bin labels for the 'binary' plot
binLabels = ["0", "1", "2","12", "3","13","23","123", "4","14","24","124","34","134","234","1234", "5","15","25","125","35","135","235","1235","45","145","245","1245","345","1345","2345","12345"]
hist_JetsWithHVParts.SetStats(0)
for x in range(1,9):
	hist_JetsWithHVParts.GetXaxis().SetBinLabel(x,binLabels[x-1])
c1 = rt.TCanvas()
hist_JetsWithHVPartsvsNjets.Draw("colz")
c1.SaveAs("2DPlot.png")
hist_nAK8_vs_nCA11.Draw("colz")
c1.SaveAs("nAK8vsnCA11.png")
hist2d_deltaR_1_3jetMT.Draw("colz")
c1.SaveAs("deltaR_1_3jetMT.png")
hist2d_deltaR_2_3jetMT.Draw("colz")
c1.SaveAs("deltaR_2_3jetMT.png")
hist2d_deltaR_3_3jetMT.Draw("colz")
c1.SaveAs("deltaR_3_3jetMT.png")
hist2d_deltaR_1_2jetMT.Draw("colz")
c1.SaveAs("deltaR_1_2jetMT.png")
hist2d_deltaR_2_2jetMT.Draw("colz")
c1.SaveAs("deltaR_2_2jetMT.png")
hist2d_deltaR_3_2jetMT.Draw("colz")
c1.SaveAs("deltaR_3_2jetMT.png")
hist2d_SDmass_12_3jetMT.Draw("colz")
c1.SaveAs("SDmass_12_3jetMT.png")
hist2d_SDmass_13_3jetMT.Draw("colz")
c1.SaveAs("SDmass_13_3jetMT.png")
hist2d_SDmass_23_3jetMT.Draw("colz")
c1.SaveAs("SDmass_23_3jetMT.png")
hist2d_SDmass_12_2jetMT.Draw("colz")
c1.SaveAs("SDmass_12_2jetMT.png")
hist2d_SDmass_13_2jetMT.Draw("colz")
c1.SaveAs("SDmass_13_2jetMT.png")
hist2d_SDmass_23_2jetMT.Draw("colz")
c1.SaveAs("SDmass_23_2jetMT.png")
hist2d_sdVardR_13.Draw("colz")
c1.SaveAs("sdVarvsdR_13.png")
hist2d_sdVardR_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarvsdR_13_3jetMT.png")
hist2d_sdVardR_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarvsdR_13_2jetMT.png")

for x in range(32):
	print(str(mat_JetsHVParts_nJets[0][x]) + " " +str(mat_JetsHVParts_nJets[1][x]) + " " +str(mat_JetsHVParts_nJets[2][x]) + " " +str(mat_JetsHVParts_nJets[3][x]) + " " +str(mat_JetsHVParts_nJets[4][x]) + " " +str(mat_JetsHVParts_nJets[5][x]) + " " +str(mat_JetsHVParts_nJets[6][x]) + " " +str(mat_JetsHVParts_nJets[7][x]) + " " +str(mat_JetsHVParts_nJets[8][x]) + " " +str(mat_JetsHVParts_nJets[9][x]))

drawHistos([hist_nAK8, hist_nCA11],"nJets")
drawHistos([hist_SDMass_withHV,hist_SDMass_withoutHV],"SoftDropMassLOG", True)
drawHistos([hist_JetsWithHVParts],"JetsWithHVParts")
drawHistos([hist_SDMass_withHV,hist_SDMass_withoutHV],"SoftDropMass")
drawHistos([hist_JetsWithHVParts],"JetsWithHVPartsLOG", True)
drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt, hist_HVQuarksMt], "DiJetMt")
drawHistos([hist_AK8DiJetMt, hist_AK8DiJetMjj,hist_AK8AllJetsMt],"AK8")
drawHistos([hist_CA11DiJetMt, hist_CA11DiJetMjj,hist_CA11AllJetsMt],"CA11")
drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt,hist_AK8AllJetsMt,hist_CA11AllJetsMt], "allJetsLOG", True)
drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt,hist_AK8AllJetsMt,hist_CA11AllJetsMt], "allJets")

drawHistos([hist_deltaPhi1,hist_deltaPhi2,hist_deltaPhi3],"DeltaPhi")
drawHistos([hist_deltaPhi1_HV,hist_deltaPhi2_HV,hist_deltaPhi3_HV],"DeltaPhi_HV")
drawHistos([hist_deltaPhi1_noHV,hist_deltaPhi2_noHV,hist_deltaPhi3_noHV],"DeltaPhi_noHV")

drawHistos([hist_deltaPhi1,hist_deltaPhi1_HV,hist_deltaPhi1_noHV],"DeltaPhi1")
drawHistos([hist_deltaPhi2,hist_deltaPhi2_HV,hist_deltaPhi2_noHV],"DeltaPhi2")
drawHistos([hist_deltaPhi3,hist_deltaPhi3_HV,hist_deltaPhi3_noHV],"DeltaPhi3")

drawHistos([hist_deltaPhi1_3jetMT,hist_deltaPhi2_3jetMT,hist_deltaPhi3_3jetMT],"DeltaPhi_3jetMT")
drawHistos([hist_deltaPhi1_2jetMT,hist_deltaPhi2_2jetMT,hist_deltaPhi3_2jetMT],"DeltaPhi_2jetMT")

drawHistos([hist_deltaPhi1_2jetMT,hist_deltaPhi1_3jetMT],"DeltaPhi1_XjetMT")
drawHistos([hist_deltaPhi2_2jetMT,hist_deltaPhi2_3jetMT],"DeltaPhi2_XjetMT")
drawHistos([hist_deltaPhi3_2jetMT,hist_deltaPhi3_3jetMT],"DeltaPhi3_XjetMT")

drawHistos([hist_deltaR_12_3jetMT,hist_deltaR_13_3jetMT,hist_deltaR_23_3jetMT], "deltaR_3jetMT")
drawHistos([hist_deltaR_12_2jetMT,hist_deltaR_13_2jetMT,hist_deltaR_23_2jetMT], "deltaR_2jetMT")

drawHistos([hist_deltaR_12_2jetMT,hist_deltaR_12_3jetMT], "deltaR_12")
drawHistos([hist_deltaR_13_2jetMT,hist_deltaR_13_3jetMT], "deltaR_13")
drawHistos([hist_deltaR_23_2jetMT,hist_deltaR_23_3jetMT], "deltaR_23")

drawHistos([hist_sdVar_12,hist_sdVar_13,hist_sdVar_23],"sdVar")
drawHistos([hist_sdVar_12_3jetMT,hist_sdVar_13_3jetMT,hist_sdVar_23_3jetMT],"sdVar_3jetMT")
drawHistos([hist_sdVar_12_2jetMT,hist_sdVar_13_2jetMT,hist_sdVar_23_2jetMT],"sdVar_2jetMT")

drawHistos([hist_sdVar_12,hist_sdVar_12_3jetMT,hist_sdVar_12_2jetMT],"sdVar_12")
drawHistos([hist_sdVar_13,hist_sdVar_13_3jetMT,hist_sdVar_13_2jetMT],"sdVar_13")
drawHistos([hist_sdVar_23,hist_sdVar_23_3jetMT,hist_sdVar_23_2jetMT],"sdVar_23")

drawHistos([hist_jetPt_mindPhi,hist_jetPt_middPhi,hist_jetPt_maxdPhi],"jetPt_dPhiOrder")
drawHistos([hist_dPhi_mindPhi,hist_dPhi_middPhi,hist_dPhi_maxdPhi],"dPhi_dPhiOrder")

drawHistos([hist_jetPt_mindPhi_3jetMT,hist_jetPt_middPhi_3jetMT,hist_jetPt_maxdPhi_3jetMT],"jetPt_dPhiOrder_3jetMT")
drawHistos([hist_dPhi_mindPhi_3jetMT,hist_dPhi_middPhi_3jetMT,hist_dPhi_maxdPhi_3jetMT],"dPhi_dPhiOrder_3jetMT")

drawHistos([hist_jetPt_mindPhi_2jetMT,hist_jetPt_middPhi_2jetMT,hist_jetPt_maxdPhi_2jetMT],"jetPt_dPhiOrder_2jetMT")
drawHistos([hist_dPhi_mindPhi_2jetMT,hist_dPhi_middPhi_2jetMT,hist_dPhi_maxdPhi_2jetMT],"dPhi_dPhiOrder_2jetMT")

drawHistos([hist_jetPt_mindPhi,hist_jetPt_mindPhi_3jetMT,hist_jetPt_mindPhi_2jetMT],"jetPt_mindPhi")
drawHistos([hist_jetPt_middPhi,hist_jetPt_middPhi_3jetMT,hist_jetPt_middPhi_2jetMT],"jetPt_middPhi")
drawHistos([hist_jetPt_maxdPhi,hist_jetPt_maxdPhi_3jetMT,hist_jetPt_maxdPhi_2jetMT],"jetPt_maxdPhi")

drawHistos([hist_dPhi_mindPhi,hist_dPhi_mindPhi_3jetMT,hist_dPhi_mindPhi_2jetMT],"dPhi_mindPhi")
drawHistos([hist_dPhi_middPhi,hist_dPhi_middPhi_3jetMT,hist_dPhi_middPhi_2jetMT],"dPhi_middPhi")
drawHistos([hist_dPhi_maxdPhi,hist_dPhi_maxdPhi_3jetMT,hist_dPhi_maxdPhi_2jetMT],"dPhi_maxdPhi")

drawHistos([hist_deltadPhi_mindPhimiddPhi,hist_deltadPhi_middPhimaxdPhi],"deltadPhi")
drawHistos([hist_deltadPhi_mindPhimiddPhi_3jetMT,hist_deltadPhi_middPhimaxdPhi_3jetMT],"deltadPhi_3jetMT")
drawHistos([hist_deltadPhi_mindPhimiddPhi_2jetMT,hist_deltadPhi_middPhimaxdPhi_2jetMT],"deltadPhi_2jetMT")

drawHistos([hist_deltadPhi_mindPhimiddPhi,hist_deltadPhi_mindPhimiddPhi_3jetMT,hist_deltadPhi_mindPhimiddPhi_2jetMT],"deltadPhi_minmid")
drawHistos([hist_deltadPhi_middPhimaxdPhi,hist_deltadPhi_middPhimaxdPhi_3jetMT,hist_deltadPhi_middPhimaxdPhi_2jetMT],"deltadPhi_midmax")
drawHistos([hist_deltadPhi_mindPhimaxdPhi,hist_deltadPhi_mindPhimaxdPhi_3jetMT,hist_deltadPhi_mindPhimaxdPhi_2jetMT],"deltadPhi_minmax")
drawHistos([hist_sdVar_3Jets,hist_sdVar_3Jets_3jetMT,hist_sdVar_3Jets_2jetMT],"sdVar_3jets")
print("Num Events with 2 or more jets = " + str(nEventsWith2OrMoreJets))
print("Num Events where leading two jet Pt are both above 170 GeV = " + str(nEventsWithLeadingTwoJetsHavingPtGreaterThan170))
print("Num Events with MET to MT ratio above .15 = " + str(nEventsWithMETMTRatioGreaterThanp15))
print("Num Events without leptons = " + str(nEventsPassedPreSelection))

inFile.Close()
