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
hist_jetPt_maxdPhi = rt.TH1F("jetPt_maxdPhi",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_maxdPhi_3jetMT = rt.TH1F("jetPt_maxdPhi_3jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_jetPt_maxdPhi_2jetMT = rt.TH1F("jetPt_maxdPhi_2jetMT",";pT;count/a.u.", 100, 0, 2000)
hist_correct_jetPt = rt.TH1F("jetPt_correct",";pT;count/a.u.", 100, 0, 2000)
hist_total_jetPt = rt.TH1F("jetPt_total",";pT;count/a.u.", 100, 0, 2000)

hist_sdVar_13 = rt.TH1F("sdVar_13",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_13_3jetMT = rt.TH1F("sdVar_13_3jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_sdVar_13_2jetMT = rt.TH1F("sdVar_13_2jetMT",";min(p_T_i)/\Sigma(p_T_i);count/a.u.",100,0,0.5)
hist_correct_sdVar = rt.TH1F("sdVar_correct",";sdVar;count/a.u.", 100,0,0.5)
hist_total_sdVar = rt.TH1F("sdVar_total",";sdVar;count/a.u.", 100,0,0.5)

hist_deltaR_12 = rt.TH1F("deltaR_12",";deltaR;count/a.u.",100,0,5.5)
hist_deltaR_12_3jetMT = rt.TH1F("deltaR_12_3jetMT",";deltaR;count/a.u.",100,0,5.5)
hist_deltaR_12_2jetMT = rt.TH1F("deltaR_12_2jetMT",";deltaR;count/a.u.",100,0,5.5)
hist_correct_deltaR_12 = rt.TH1F("deltaR_12_correct",";deltaR;count/a.u.", 100,0,5.5)
hist_total_deltaR_12 = rt.TH1F("deltaR_12_total",";deltaR;count/a.u.", 100,0,5.5)

hist_maxDphi = rt.TH1F("maxDphi",";dPhi;count/a.u.",100,0,rt.TMath.Pi())
hist_maxDphi_3jetMT = rt.TH1F("maxDphi_3jetMT",";dPhi;count/a.u.",100,0,rt.TMath.Pi())
hist_maxDphi_2jetMT = rt.TH1F("maxDphi_2jetMT",";dPhi;count/a.u.",100,0,rt.TMath.Pi())
hist_correct_maxDphi = rt.TH1F("maxDphi_correct",";dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_total_maxDphi = rt.TH1F("maxDphi_total",";dPhi;count/a.u.", 100,0,rt.TMath.Pi())



hist2d_sdVarDeltaR_12 = rt.TH2F("sdVarDeltaR_12","Jets 12;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_13 = rt.TH2F("sdVarDeltaR_13","Jets 13;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_23 = rt.TH2F("sdVarDeltaR_23","Jets 23;sdVar;DeltaR",100,0,0.5,100,0,7)

hist2d_sdVarDeltaR_12_2jetMT = rt.TH2F("sdVarDeltaR_12_2jetMT","Jets 12_2jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_13_2jetMT = rt.TH2F("sdVarDeltaR_13_2jetMT","Jets 13_2jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_23_2jetMT = rt.TH2F("sdVarDeltaR_23_2jetMT","Jets 23_2jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)

hist2d_sdVarDeltaR_12_3jetMT = rt.TH2F("sdVarDeltaR_12_3jetMT","Jets 12_3jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_13_3jetMT = rt.TH2F("sdVarDeltaR_13_3jetMT","Jets 13_3jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)
hist2d_sdVarDeltaR_23_3jetMT = rt.TH2F("sdVarDeltaR_23_3jetMT","Jets 23_3jetMT;sdVar;DeltaR",100,0,0.5,100,0,7)

hist2d_sdVarpT_12 = rt.TH2F("sdVarpT_12","Jets 12;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_13 = rt.TH2F("sdVarpT_13","Jets 13;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_23 = rt.TH2F("sdVarpT_23","Jets 23;sdVar;DpT",100,0,0.6,100,0,1600)

hist2d_sdVarpT_12_2jetMT = rt.TH2F("sdVarpT_12_2jetMT","Jets 12_2jetMT;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_13_2jetMT = rt.TH2F("sdVarpT_13_2jetMT","Jets 13_2jetMT;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_23_2jetMT = rt.TH2F("sdVarpT_23_2jetMT","Jets 23_2jetMT;sdVar;DpT",100,0,0.6,100,0,1600)

hist2d_sdVarpT_12_3jetMT = rt.TH2F("sdVarpT_12_3jetMT","Jets 12_3jetMT;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_13_3jetMT = rt.TH2F("sdVarpT_13_3jetMT","Jets 13_3jetMT;sdVar;DpT",100,0,0.6,100,0,1600)
hist2d_sdVarpT_23_3jetMT = rt.TH2F("sdVarpTR_23_3jetMT","Jets 23_3jetMT;sdVar;DpT",100,0,0.6,100,0,1600)


hist2d_sdVarpT1_12 = rt.TH2F("sdVarpT1_12","Jets 12;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_13 = rt.TH2F("sdVarpT1_13","Jets 13;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_23 = rt.TH2F("sdVarpT1_23","Jets 23;sdVar;pT1",100,0,0.6,100,0,3000)

hist2d_sdVarpT1_12_2jetMT = rt.TH2F("sdVarpT1_12_2jetMT","Jets 12_2jetMT;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_13_2jetMT = rt.TH2F("sdVarpT1_13_2jetMT","Jets 13_2jetMT;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_23_2jetMT = rt.TH2F("sdVarpT1_23_2jetMT","Jets 23_2jetMT;sdVar;pT1",100,0,0.6,100,0,3000)

hist2d_sdVarpT1_12_3jetMT = rt.TH2F("sdVarpT1_12_3jetMT","Jets 12_3jetMT;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_13_3jetMT = rt.TH2F("sdVarpT1_13_3jetMT","Jets 13_3jetMT;sdVar;pT1",100,0,0.6,100,0,3000)
hist2d_sdVarpT1_23_3jetMT = rt.TH2F("sdVarpT1_23_3jetMT","Jets 23_3jetMT;sdVar;pT1",100,0,0.6,100,0,3000)

hist2d_sdVarpT2_12 = rt.TH2F("sdVarpT2_12","Jets 12;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_13 = rt.TH2F("sdVarpT2_13","Jets 13;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_23 = rt.TH2F("sdVarpT2_23","Jets 23;sdVar;pT2",100,0,0.6,100,0,3000)

hist2d_sdVarpT2_12_2jetMT = rt.TH2F("sdVarpT2_12_2jetMT","Jets 12_2jetMT;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_13_2jetMT = rt.TH2F("sdVarpT2_13_2jetMT","Jets 13_2jetMT;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_23_2jetMT = rt.TH2F("sdVarpT2_23_2jetMT","Jets 23_2jetMT;sdVar;pT2",100,0,0.6,100,0,3000)

hist2d_sdVarpT2_12_3jetMT = rt.TH2F("sdVarpT2_12_3jetMT","Jets 12_3jetMT;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_13_3jetMT = rt.TH2F("sdVarpT2_13_3jetMT","Jets 13_3jetMT;sdVar;pT2",100,0,0.6,100,0,3000)
hist2d_sdVarpT2_23_3jetMT = rt.TH2F("sdVarpT2_23_3jetMT","Jets 23_3jetMT;sdVar;pT2",100,0,0.6,100,0,3000)

hist2d_sdVarpT3_12 = rt.TH2F("sdVarpT3_12","Jets 12;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_13 = rt.TH2F("sdVarpT3_13","Jets 13;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_23 = rt.TH2F("sdVarpT3_23","Jets 23;sdVar;pT3",100,0,0.6,100,0,3000)

hist2d_sdVarpT3_12_2jetMT = rt.TH2F("sdVarpT3_12_2jetMT","Jets 12_2jetMT;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_13_2jetMT = rt.TH2F("sdVarpT3_13_2jetMT","Jets 13_2jetMT;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_23_2jetMT = rt.TH2F("sdVarpT3_23_2jetMT","Jets 23_2jetMT;sdVar;pT3",100,0,0.6,100,0,3000)

hist2d_sdVarpT3_12_3jetMT = rt.TH2F("sdVarpT3_12_3jetMT","Jets 12_3jetMT;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_13_3jetMT = rt.TH2F("sdVarpT3_13_3jetMT","Jets 13_3jetMT;sdVar;pT3",100,0,0.6,100,0,3000)
hist2d_sdVarpT3_23_3jetMT = rt.TH2F("sdVarpT3_23_3jetMT","Jets 23_3jetMT;sdVar;pT3",100,0,0.6,100,0,3000)

hist2_jetPtMaxDPhiVSsdVar13 = rt.TH2F("jetPtMaxDPhiVSsdVar13",";pT of maxdPhi Jet;sdVar_13",100,0,3000,100,0,0.5)
hist2_jetPtMaxDPhiVSsdVar13_2jetMT = rt.TH2F("jetPtMaxDPhiVSsdVar13_2jetMT",";pT of maxdPhi Jet;sdVar_13",100,0,3000,100,0,0.5)
hist2_jetPtMaxDPhiVSsdVar13_3jetMT = rt.TH2F("jetPtMaxDPhiVSsdVar13_3jetMT",";pT of maxdPhi Jet;sdVar_13",100,0,3000,100,0,0.5)

hist_maxDphiIdx = rt.TH1F("maxdPhiIdx",";Jet Number;count/a.u.",3,0,3)
hist_maxDphiIdx_2jetMT = rt.TH1F("maxdPhiIdx_2jetMT",";Jet Number;count/a.u.",3,0,3)
hist_maxDphiIdx_3jetMT = rt.TH1F("maxdPhiIdx_3jetMT",";Jet Number;count/a.u.",3,0,3)
	

histList = []
histList.append(hist_jetPt_maxdPhi)
histList.append(hist_jetPt_maxdPhi_3jetMT)
histList.append(hist_jetPt_maxdPhi_2jetMT)
histList.append(hist_sdVar_13)
histList.append(hist_sdVar_13_3jetMT)
histList.append(hist_sdVar_13_2jetMT)
histList.append(hist_deltaR_12)
histList.append(hist_deltaR_12_3jetMT)
histList.append(hist_deltaR_12_2jetMT)
histList.append(hist_maxDphi)
histList.append(hist_maxDphi_3jetMT)
histList.append(hist_maxDphi_2jetMT)

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
	jetCode = int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2)
	mat_JetsHVParts_nJets[nAK8Jets][int(jetsWithHVDecendants[0]+jetsWithHVDecendants[1]+jetsWithHVDecendants[2]+jetsWithHVDecendants[3]+jetsWithHVDecendants[4],2)] += 1
	
	sdVar = [min(tree.JetsAK8[0].Pt(),tree.JetsAK8[1].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()),min(tree.JetsAK8[0].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()),min(tree.JetsAK8[1].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())]
	DeltaRList = [tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]),tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]),tree.JetsAK8[1].DeltaR(tree.JetsAK8[2])]
	hist_sdVar_13.Fill(sdVar[1])

	#print(isFromHVQuark)
	# figureout which order the jets rae in distance from METPhi
	deltaPhiList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
	IdxList = [0,1,2]
	mindPhiIdx = deltaPhiList.index(min(deltaPhiList))
	IdxList.remove(mindPhiIdx)
	maxdPhiIdx = deltaPhiList.index(max(deltaPhiList))
	IdxList.remove(maxdPhiIdx)
	middPhiIdx = IdxList[0]

	hist_jetPt_maxdPhi.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
	hist_deltaR_12.Fill(DeltaRList[0])
	hist_maxDphi.Fill(deltaPhiList[maxdPhiIdx])
	if jetCode == 7:
		hist_jetPt_maxdPhi_3jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
		hist2d_sdVarDeltaR_12_3jetMT.Fill(sdVar[0],DeltaRList[0])
		hist2d_sdVarDeltaR_13_3jetMT.Fill(sdVar[1],DeltaRList[1])
		hist2d_sdVarDeltaR_23_3jetMT.Fill(sdVar[2],DeltaRList[2])
		hist2d_sdVarpT_12_3jetMT.Fill(sdVar[0],tree.JetsAK8[0].Pt()-tree.JetsAK8[1].Pt())
		hist2d_sdVarpT_13_3jetMT.Fill(sdVar[1],tree.JetsAK8[0].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT_23_3jetMT.Fill(sdVar[2],tree.JetsAK8[1].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT1_12_3jetMT.Fill(sdVar[0],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_13_3jetMT.Fill(sdVar[1],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_23_3jetMT.Fill(sdVar[2],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT2_12_3jetMT.Fill(sdVar[0],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_13_3jetMT.Fill(sdVar[1],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_23_3jetMT.Fill(sdVar[2],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT3_12_3jetMT.Fill(sdVar[0],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_13_3jetMT.Fill(sdVar[1],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_23_3jetMT.Fill(sdVar[2],tree.JetsAK8[2].Pt())
		hist_deltaR_12_3jetMT.Fill(DeltaRList[0])
		hist_sdVar_13_3jetMT.Fill(sdVar[1])
		hist_maxDphi_3jetMT.Fill(deltaPhiList[maxdPhiIdx])
		hist2_jetPtMaxDPhiVSsdVar13_3jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt(),sdVar[1])
		hist_maxDphiIdx_2jetMT.Fill(maxdPhiIdx)
	elif jetCode == 3:
		hist_jetPt_maxdPhi_2jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt())
		hist2d_sdVarDeltaR_12_2jetMT.Fill(sdVar[0],DeltaRList[0])
		hist2d_sdVarDeltaR_13_2jetMT.Fill(sdVar[1],DeltaRList[1])
		hist2d_sdVarDeltaR_23_2jetMT.Fill(sdVar[2],DeltaRList[2])
		hist2d_sdVarpT_12_2jetMT.Fill(sdVar[0],tree.JetsAK8[0].Pt()-tree.JetsAK8[1].Pt())
		hist2d_sdVarpT_13_2jetMT.Fill(sdVar[1],tree.JetsAK8[0].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT_23_2jetMT.Fill(sdVar[2],tree.JetsAK8[1].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT1_12_2jetMT.Fill(sdVar[0],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_13_2jetMT.Fill(sdVar[1],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_23_2jetMT.Fill(sdVar[2],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT2_12_2jetMT.Fill(sdVar[0],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_13_2jetMT.Fill(sdVar[1],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_23_2jetMT.Fill(sdVar[2],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT3_12_2jetMT.Fill(sdVar[0],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_13_2jetMT.Fill(sdVar[1],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_23_2jetMT.Fill(sdVar[2],tree.JetsAK8[2].Pt())
		hist_deltaR_12_2jetMT.Fill(DeltaRList[0])
		hist_sdVar_13_2jetMT.Fill(sdVar[1])
		hist_maxDphi_2jetMT.Fill(deltaPhiList[maxdPhiIdx])
		hist2_jetPtMaxDPhiVSsdVar13_2jetMT.Fill(tree.JetsAK8[maxdPhiIdx].Pt(),sdVar[1])
		hist_maxDphiIdx_2jetMT.Fill(maxdPhiIdx)
	else:
		hist2d_sdVarDeltaR_12.Fill(sdVar[0],DeltaRList[0])
		hist2d_sdVarDeltaR_13.Fill(sdVar[1],DeltaRList[1])
		hist2d_sdVarDeltaR_23.Fill(sdVar[2],DeltaRList[2])
		hist2d_sdVarpT_12.Fill(sdVar[0],tree.JetsAK8[0].Pt()-tree.JetsAK8[1].Pt())
		hist2d_sdVarpT_13.Fill(sdVar[1],tree.JetsAK8[0].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT_23.Fill(sdVar[2],tree.JetsAK8[1].Pt()-tree.JetsAK8[2].Pt())
		hist2d_sdVarpT1_12.Fill(sdVar[0],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_13.Fill(sdVar[1],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT1_23.Fill(sdVar[2],tree.JetsAK8[0].Pt())
		hist2d_sdVarpT2_12.Fill(sdVar[0],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_13.Fill(sdVar[1],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT2_23.Fill(sdVar[2],tree.JetsAK8[1].Pt())
		hist2d_sdVarpT3_12.Fill(sdVar[0],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_13.Fill(sdVar[1],tree.JetsAK8[2].Pt())
		hist2d_sdVarpT3_23.Fill(sdVar[2],tree.JetsAK8[2].Pt())
		hist2_jetPtMaxDPhiVSsdVar13.Fill(tree.JetsAK8[maxdPhiIdx].Pt(),sdVar[1])
		hist_maxDphiIdx.Fill(maxdPhiIdx)
	for iBin in range(hist_jetPt_maxdPhi.GetXaxis().GetNbins()):
		ptCut = hist_jetPt_maxdPhi.GetXaxis().GetBinLowEdge(iBin)
		if (jetCode == 3 and tree.JetsAK8[maxdPhiIdx].Pt() > ptCut) or (jetCode != 3 and tree.JetsAK8[maxdPhiIdx].Pt() <= ptCut):
			hist_correct_jetPt.Fill(ptCut)
		hist_total_jetPt.Fill(ptCut)

	for iBin in range(hist_maxDphi.GetXaxis().GetNbins()):
		dPhiCut = hist_maxDphi.GetXaxis().GetBinLowEdge(iBin)
		if (jetCode == 3 and  deltaPhiList[maxdPhiIdx] > dPhiCut) or (jetCode != 3 and  deltaPhiList[maxdPhiIdx] <= dPhiCut):
			hist_correct_maxDphi.Fill(dPhiCut)
		hist_total_maxDphi.Fill(dPhiCut)
	
	for iBin in range(hist_sdVar_13.GetXaxis().GetNbins()):
		zCut = hist_sdVar_13.GetXaxis().GetBinLowEdge(iBin)
		if (jetCode == 3 and sdVar[1] < zCut) or (jetCode != 3 and sdVar[1] >= zCut):
			hist_correct_sdVar.Fill(zCut)
		hist_total_sdVar.Fill(zCut)

	for iBin in range(hist_deltaR_12.GetXaxis().GetNbins()):
		dRcut = hist_deltaR_12.GetXaxis().GetBinLowEdge(iBin)
		if (jetCode == 3 and DeltaRList[0] > dRcut) or (jetCode != 3 and DeltaRList[0] <= dRcut):
			hist_correct_deltaR_12.Fill(dRcut)
		hist_total_deltaR_12.Fill(dRcut)

for histo in histList:
	total = histo.GetEntries()
	if total != 0:
		histo.Scale(1/(total))

c1 = rt.TCanvas("c1","c1",900,600)
c1.SetRightMargin(0.1)
hist2_jetPtMaxDPhiVSsdVar13.Draw("colz")
c1.SaveAs("jetPtMaxDPhiVSsdVar13.png")
hist2_jetPtMaxDPhiVSsdVar13_2jetMT.Draw("colz")
c1.SaveAs("jetPtMaxDPhiVSsdVar13_2jetMT.png")
hist2_jetPtMaxDPhiVSsdVar13_3jetMT.Draw("colz")
c1.SaveAs("jetPtMaxDPhiVSsdVar13_3jetMT.png")
hist2d_sdVarDeltaR_12.Draw("colz")
c1.SaveAs("sdVarDeltaR_12.png")
hist2d_sdVarDeltaR_13.Draw("colz")
c1.SaveAs("sdVarDeltaR_13.png")
hist2d_sdVarDeltaR_23.Draw("colz")
c1.SaveAs("sdVarDeltaR_23.png")

hist2d_sdVarDeltaR_12_2jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_12_2jetMT.png")
hist2d_sdVarDeltaR_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_13_2jetMT.png")
hist2d_sdVarDeltaR_23_2jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_23_2jetMT.png")

hist2d_sdVarDeltaR_12_3jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_12_3jetMT.png")
hist2d_sdVarDeltaR_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_13_3jetMT.png")
hist2d_sdVarDeltaR_23_3jetMT.Draw("colz")
c1.SaveAs("sdVarDeltaR_23_3jetMT.png")

hist2d_sdVarpT_12.Draw("colz")
c1.SaveAs("sdVarpT_12.png")
hist2d_sdVarpT_13.Draw("colz")
c1.SaveAs("sdVarpT_13.png")
hist2d_sdVarpT_23.Draw("colz")
c1.SaveAs("sdVarpT_23.png")

hist2d_sdVarpT_12_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT_12_2jetMT.png")
hist2d_sdVarpT_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT_13_2jetMT.png")
hist2d_sdVarpT_23_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT_23_2jetMT.png")

hist2d_sdVarpT_12_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT_12_3jetMT.png")
hist2d_sdVarpT_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT_13_3jetMT.png")
hist2d_sdVarpT_23_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT_23_3jetMT.png")

hist2d_sdVarpT1_12.Draw("colz")
c1.SaveAs("sdVarpT1_12.png")
hist2d_sdVarpT1_13.Draw("colz")
c1.SaveAs("sdVarpT1_13.png")
hist2d_sdVarpT1_23.Draw("colz")
c1.SaveAs("sdVarpT1_23.png")

hist2d_sdVarpT1_12_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_12_2jetMT.png")
hist2d_sdVarpT1_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_13_2jetMT.png")
hist2d_sdVarpT1_23_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_23_2jetMT.png")

hist2d_sdVarpT1_12_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_12_3jetMT.png")
hist2d_sdVarpT1_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_13_3jetMT.png")
hist2d_sdVarpT1_23_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT1_23_3jetMT.png")


hist2d_sdVarpT2_12.Draw("colz")
c1.SaveAs("sdVarpT2_12.png")
hist2d_sdVarpT2_13.Draw("colz")
c1.SaveAs("sdVarpT2_13.png")
hist2d_sdVarpT2_23.Draw("colz")
c1.SaveAs("sdVarpT2_23.png")

hist2d_sdVarpT2_12_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_12_2jetMT.png")
hist2d_sdVarpT2_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_13_2jetMT.png")
hist2d_sdVarpT2_23_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_23_2jetMT.png")

hist2d_sdVarpT2_12_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_12_3jetMT.png")
hist2d_sdVarpT2_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_13_3jetMT.png")
hist2d_sdVarpT2_23_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT2_23_3jetMT.png")

hist2d_sdVarpT3_12.Draw("colz")
c1.SaveAs("sdVarpT3_12.png")
hist2d_sdVarpT3_13.Draw("colz")
c1.SaveAs("sdVarpT3_13.png")
hist2d_sdVarpT3_23.Draw("colz")
c1.SaveAs("sdVarpT3_23.png")

hist2d_sdVarpT3_12_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_12_2jetMT.png")
hist2d_sdVarpT3_13_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_13_2jetMT.png")
hist2d_sdVarpT3_23_2jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_23_2jetMT.png")

hist2d_sdVarpT3_12_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_12_3jetMT.png")
hist2d_sdVarpT3_13_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_13_3jetMT.png")
hist2d_sdVarpT3_23_3jetMT.Draw("colz")
c1.SaveAs("sdVarpT3_23_3jetMT.png")

drawHistos([hist_maxDphiIdx,hist_maxDphiIdx_3jetMT,hist_maxDphiIdx_2jetMT],"maxdPhiIdx",True)

drawHistos([hist_jetPt_maxdPhi,hist_jetPt_maxdPhi_3jetMT,hist_jetPt_maxdPhi_2jetMT],"jetPt_maxdPhi")
drawHistos([hist_correct_jetPt, hist_total_jetPt],"correctandtotal_jetPt")
hist_correct_jetPt.Divide(hist_total_jetPt)
drawHistos([hist_correct_jetPt],"FOM_jetPt")

drawHistos([hist_sdVar_13,hist_sdVar_13_3jetMT,hist_sdVar_13_2jetMT],"sdVar")
drawHistos([hist_correct_sdVar, hist_total_sdVar],"correctandtotal_sdVar")
hist_correct_sdVar.Divide(hist_total_sdVar)
drawHistos([hist_correct_sdVar],"FOM_sdVar")

drawHistos([hist_deltaR_12,hist_deltaR_12_3jetMT,hist_deltaR_12_2jetMT],"deltaR_12")
drawHistos([hist_correct_deltaR_12, hist_total_deltaR_12],"correctandtotal_deltaR_12")
hist_correct_deltaR_12.Divide(hist_total_deltaR_12)
drawHistos([hist_correct_deltaR_12],"FOM_deltaR_12")

drawHistos([hist_maxDphi,hist_maxDphi_3jetMT,hist_maxDphi_2jetMT],"maxDphi")
drawHistos([hist_correct_maxDphi, hist_total_maxDphi],"correctandtotal_maxDphi")
hist_correct_maxDphi.Divide(hist_total_maxDphi)
drawHistos([hist_correct_maxDphi],"FOM_maxDphi")


print(hist_correct_maxDphi.GetName() + " max value is " + str(hist_correct_maxDphi.GetMaximum()) + " at bin low edge " + str(hist_correct_maxDphi.GetBinLowEdge(hist_correct_maxDphi.GetMaximumBin())))
print(hist_correct_maxDphi.GetName() + " min value is " + str(hist_correct_maxDphi.GetMinimum()) + " at bin low edge " + str(hist_correct_maxDphi.GetBinLowEdge(hist_correct_maxDphi.GetMinimumBin())))

print(hist_correct_jetPt.GetName() + " max value is " + str(hist_correct_jetPt.GetMaximum()) + " at bin low edge " + str(hist_correct_jetPt.GetBinLowEdge(hist_correct_jetPt.GetMaximumBin())))
print(hist_correct_jetPt.GetName() + " min value is " + str(hist_correct_jetPt.GetMinimum()) + " at bin low edge " + str(hist_correct_jetPt.GetBinLowEdge(hist_correct_jetPt.GetMinimumBin())))

print(hist_correct_sdVar.GetName() + " max value is " + str(hist_correct_sdVar.GetMaximum()) + " at bin low edge " + str(hist_correct_sdVar.GetBinLowEdge(hist_correct_sdVar.GetMaximumBin())))
print(hist_correct_sdVar.GetName() + " min value is " + str(hist_correct_sdVar.GetMinimum()) + " at bin low edge " + str(hist_correct_sdVar.GetBinLowEdge(hist_correct_sdVar.GetMinimumBin())))


print(hist_correct_deltaR_12.GetName() + " max value is " + str(hist_correct_deltaR_12.GetMaximum()) + " at bin low edge " + str(hist_correct_deltaR_12.GetBinLowEdge(hist_correct_deltaR_12.GetMaximumBin())))
print(hist_correct_deltaR_12.GetName() + " min value is " + str(hist_correct_deltaR_12.GetMinimum()) + " at bin low edge " + str(hist_correct_deltaR_12.GetBinLowEdge(hist_correct_deltaR_12.GetMinimumBin())))


print("Num Events with 2 or more jets = " + str(nEventsWith2OrMoreJets))
print("Num Events where leading two jet Pt are both above 170 GeV = " + str(nEventsWithLeadingTwoJetsHavingPtGreaterThan170))
print("Num Events with MET to MT ratio above .15 = " + str(nEventsWithMETMTRatioGreaterThanp15))
print("Num Events without leptons = " + str(nEventsPassedPreSelection))

inFile.Close()
