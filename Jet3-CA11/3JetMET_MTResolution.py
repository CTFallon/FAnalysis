import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
#Script to compute the MT resolution at several cuts along a variable
# should also include the MT distribution, along with (mean or peak, stdDev) for the all 12, all 123, 'perfect', and cuts.

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
alpha = '1' #['0p1', '0p2', '0p5', '1']
f_rvis = 1.0 - 0.3

includeAllEvents = False

inFile = rt.TFile.Open("root://cmsxrootd.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV14/"+createInFileName(mZprime, mDark,rinv,alpha))
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
nEventsWith2OrMoreJets = 0
nEventsWithLeadingTwoJetsHavingPtGreaterThan170 = 0
nEventsWithMETMTRatioGreaterThanp15 = 0
nEventsPassedPreSelection = 0


# initilize histograms

hist_MT_12   = rt.TH1F("MT_12",";MT;count/a.u.",100,0,5000)
hist_MT_123  = rt.TH1F("MT_123",";MT;count/a.u.",100,0,5000)
hist_MT_perf = rt.TH1F("MT_perf",";MT;count/a.u.",100,0,5000)
hist_MT_ideal = rt.TH1F("MT_ideal",";MT;count/a.u.",100,0,5000)
hist_MT_ptFracCut = rt.TH1F("MT_ptFracCut",";MT;count/a.u.",100,0,5000)

hist_MT_delta =  rt.TH1F("MT_delta",";deltaMT;count/a.u.",100,0,3000)
hist_MT_delta_12 =  rt.TH1F("MT_delta_12",";deltaMT;count/a.u.",100,0,3000)
hist_MT_delta_123 =  rt.TH1F("MT_delta_123",";deltaMT;count/a.u.",100,0,3000)

hist_jetCode = rt.TH1F("jetCode",";jetCode;count",10,0,10)
hist_jetCodePassesptFracCut = rt.TH1F("jetCodePassesptFracCut",";jetCode;count",10,0,10)

hist_sdVar_13_MTReso = rt.TH1F("sdVar_13_MTReso",";sdVar_13;MTReso",100,0,0.5) # 100,0,0.5
hist_MT_sdVar13List = []
for iBin in range(hist_sdVar_13_MTReso.GetXaxis().GetNbins()):
	hist_MT_sdVar13List.append(rt.TH1F("MT_sdVarCut_"+str(iBin),";MT;count/a.u.",100,0,5000))


hist_jetPt_maxDphi_MTReso = rt.TH1F("jetPt_maxDphi_MTReso",";jetPt(maxdPhi));MTReso",100, 0, 2000)# 100, 0, 2000
hist_MT_jetPtmaxDphiList = []
for iBin in range(hist_jetPt_maxDphi_MTReso.GetXaxis().GetNbins()):
	hist_MT_jetPtmaxDphiList.append(rt.TH1F("MT_jetPTmaxDphiCut_"+str(iBin),";MT;count/a.u.",100,0,5000))

hist2D_sdVar_jetPtMaxdPhi_MTReso = rt.TH2F("sdVar_jetPtMaxdPhi_MTReso","MT Resolution;sdVarCut;jetPT(maxDPhi)",10,0.16,0.18,10,890,960)
hist2D_MT_List = [[] for iBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetNbins())]
for iBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetNbins()):
	for jBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetYaxis().GetNbins()):
		hist2D_MT_List[iBin].append(rt.TH1F("MT_"+str(iBin)+"_"+str(jBin),";MT;count/a.u.",100,0,5000))

histList = []

nPlotsMade = 0


mat_JetsHVParts_nJets = [[0 for x in range(32)] for y in range(10)]
for iEvent in range(nEvents):
	#if nEventsPassedPreSelection > 1000:
	#	continue

	tree.GetEvent(iEvent)
	if iEvent%100==0:
		print("EventNumber: "+str(iEvent)+"/"+str(nEvents)+ " ("+str(100*round(iEvent/float(nEvents),4))+"%)")

	jetCollection = tree.JetsAK8
	nAK8Jets = len(jetCollection)

	#some calculations to determine which histos to fill
	# PreSelection Cuts
	if not (len(jetCollection)>=2): #moreThan2Jets, temp set to ==3 for 3JetMt classification
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

	#if trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi) < 3000:
	#	continue
	numberOfDaughtersAParticleHas = [0 for x in range(len(tree.GenParticles))]
	for iPart in range(2, len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if iParent != -1: 
			numberOfDaughtersAParticleHas[iParent] += 1
	AK8jetsWithHVDecendants = ["0","0","0","0","0"]
	AK8_nHVParts = [0,0,0,0,0]
	AK8_nParts = [0,0,0,0,0]
	AK8_ptHVParts = [0,0,0,0,0]
	AK8_ptParts = [0,0,0,0,0]
	#make vector of length nGenParts that is 1 if the particle came from a HV quark
	isFromHVQuark = [0 for x in range(len(tree.GenParticles))]
	listOfHVQuarks = []
	for iPart in range(2,len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if abs(tree.GenParticles_PdgId[iPart]) == 4900101:
			listOfHVQuarks.append(tree.GenParticles[iPart])
		if iParent >= iPart:
			print("Ut-oh, the parent has a higher index than the child...")
		if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or (isFromHVQuark[iParent]):
			isFromHVQuark[iPart] = 1
		#finding what Jet a particle is in, only if it decends from a HVQuark:
		for iJet in range(min(len(tree.JetsAK8),5)):
			if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
				AK8_ptParts[-iJet-1] += tree.GenParticles[iPart].Pt()
				AK8_nParts[-iJet-1] += 1
		if isFromHVQuark[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
			for iJet in range(min(len(tree.JetsAK8),5)):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart]:
					AK8jetsWithHVDecendants[-iJet-1] = "1"
					AK8_ptHVParts[-iJet-1] += tree.GenParticles[iPart].Pt()
					AK8_nHVParts[-iJet-1] += 1


	AK8ptFrac = [0,0,0,0,0]
	jetPassesptFracCut = ["0","0","0","0","0"]
	ptFracCut = [0.0,0.0,0.22,0.0,0.0]
	for iJet in range(len(AK8jetsWithHVDecendants)):
		if AK8_ptParts[-iJet-1] != 0:
			AK8ptFrac[-iJet-1] = AK8_ptHVParts[-iJet-1]/AK8_ptParts[-iJet-1]
		if AK8ptFrac[-iJet-1] > ptFracCut[-iJet-1]:
			jetPassesptFracCut[-iJet-1] = "1"
	jetCode = int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2)
	jetCodePassesptFracCut = int(
				str(int(AK8jetsWithHVDecendants[0])*int(jetPassesptFracCut[0])) + 
				str(int(AK8jetsWithHVDecendants[1])*int(jetPassesptFracCut[1])) + 
				str(int(AK8jetsWithHVDecendants[2])*int(jetPassesptFracCut[2])) + 
				str(int(AK8jetsWithHVDecendants[3])*int(jetPassesptFracCut[3])) + 
				str(int(AK8jetsWithHVDecendants[4])*int(jetPassesptFracCut[4])) 
				,2)
	mat_JetsHVParts_nJets[nAK8Jets][int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2)] += 1
	hist_jetCode.Fill(jetCode)
	hist_jetCodePassesptFracCut.Fill(jetCodePassesptFracCut)
	if len(tree.JetsAK8) == 3:
		sdVar = [min(tree.JetsAK8[0].Pt(),tree.JetsAK8[1].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[1].Pt()),min(tree.JetsAK8[0].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[0].Pt()+tree.JetsAK8[2].Pt()),min(tree.JetsAK8[1].Pt(),tree.JetsAK8[2].Pt())/(tree.JetsAK8[1].Pt()+tree.JetsAK8[2].Pt())]
		DeltaRList = [tree.JetsAK8[0].DeltaR(tree.JetsAK8[1]),tree.JetsAK8[0].DeltaR(tree.JetsAK8[2]),tree.JetsAK8[1].DeltaR(tree.JetsAK8[2])]

		#print(isFromHVQuark)
		# figureout which order the jets rae in distance from METPhi
		deltaPhiList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
		IdxList = [0,1,2]
		mindPhiIdx = deltaPhiList.index(min(deltaPhiList))
		IdxList.remove(mindPhiIdx)
		maxdPhiIdx = deltaPhiList.index(max(deltaPhiList))
		IdxList.remove(maxdPhiIdx)
		middPhiIdx = IdxList[0]

		jetFourVectorsWithHVDecendants = []
		jetFourVectorsThatPassPtFracCut = []
		for iJet in range(min(len(tree.JetsAK8),5)):
			if AK8jetsWithHVDecendants[-iJet-1] == "1":
				jetFourVectorsWithHVDecendants.append(tree.JetsAK8[iJet])
				if jetPassesptFracCut[-iJet-1] == "1":
					jetFourVectorsThatPassPtFracCut.append(tree.JetsAK8[iJet])

		hist_MT_12.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		hist_MT_123.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
		hist_MT_ideal.Fill(trans_mass_Njet(jetFourVectorsWithHVDecendants, tree.MET, tree.METPhi))
		hist_MT_ptFracCut.Fill(trans_mass_Njet(jetFourVectorsThatPassPtFracCut, tree.MET, tree.METPhi))
		hist_MT_delta.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)-trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		if jetCode == 3:
			hist_MT_perf.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
			hist_MT_delta_12.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)-trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		else:
			hist_MT_perf.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
			hist_MT_delta_123.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi)-trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		
		for iBin in range(hist_sdVar_13_MTReso.GetXaxis().GetNbins()):
			zCut = hist_sdVar_13_MTReso.GetXaxis().GetBinLowEdge(iBin)
			if (sdVar[1] < zCut):
				hist_MT_sdVar13List[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
			elif (sdVar[1] >= zCut):
				hist_MT_sdVar13List[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))

		for iBin in range(hist_jetPt_maxDphi_MTReso.GetXaxis().GetNbins()):
			zCut = hist_jetPt_maxDphi_MTReso.GetXaxis().GetBinLowEdge(iBin)
			if (tree.JetsAK8[maxdPhiIdx].Pt() > zCut):
				hist_MT_jetPtmaxDphiList[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
			elif (tree.JetsAK8[maxdPhiIdx].Pt() <= zCut):
				hist_MT_jetPtmaxDphiList[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
		
		for iBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetNbins()):
			zCut = hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetBinLowEdge(iBin)
			for jBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetYaxis().GetNbins()):
				ptCut = hist2D_sdVar_jetPtMaxdPhi_MTReso.GetYaxis().GetBinLowEdge(jBin)
				if (tree.JetsAK8[maxdPhiIdx].Pt() > ptCut) or (sdVar[1] < zCut):
					hist2D_MT_List[iBin][jBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
				else:
					hist2D_MT_List[iBin][jBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
		
	elif includeAllEvents == True:
		hist_MT_12.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		hist_MT_123.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		hist_MT_ptFracCut.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
		for iBin in range(hist_sdVar_13_MTReso.GetXaxis().GetNbins()):
			hist_MT_sdVar13List[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))

		for iBin in range(hist_jetPt_maxDphi_MTReso.GetXaxis().GetNbins()):
			hist_MT_jetPtmaxDphiList[iBin].Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))

for histo in histList:
	total = histo.GetEntries()
	if total != 0:
		histo.Scale(1/(total))

iBinMin = 100000
minReso = 100000
for iBin in range(hist_sdVar_13_MTReso.GetXaxis().GetNbins()):
	binReso = hist_MT_sdVar13List[iBin].GetRMS()/hist_MT_sdVar13List[iBin].GetMean()
	hist_sdVar_13_MTReso.SetBinContent(iBin+1, binReso)
	if minReso > binReso:
		iBinMin = iBin
		minReso = binReso
		print(iBin, hist_sdVar_13_MTReso.GetXaxis().GetBinLowEdge(iBin),hist_MT_sdVar13List[iBin].GetMean(), hist_MT_sdVar13List[iBin].GetRMS(),binReso)
print("JETPT CUTS NOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
iBinMinjetPt = 100000
minReso = 100000
for iBin in range(hist_jetPt_maxDphi_MTReso.GetXaxis().GetNbins()):
	binReso = hist_MT_jetPtmaxDphiList[iBin].GetRMS()/hist_MT_jetPtmaxDphiList[iBin].GetMean()
	hist_jetPt_maxDphi_MTReso.SetBinContent(iBin+1, binReso)
	if minReso > binReso:
		iBinMinjetPt = iBin
		minReso = binReso
		print(iBin, hist_jetPt_maxDphi_MTReso.GetXaxis().GetBinLowEdge(iBin),hist_MT_jetPtmaxDphiList[iBin].GetMean(), hist_MT_jetPtmaxDphiList[iBin].GetRMS(),binReso)


print("DOUBLE CUTS NOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
iBinMinsdVar = 100000
jBinMinjetPt = 100000
minReso = 100000


for iBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetNbins()):
	for jBin in range(hist2D_sdVar_jetPtMaxdPhi_MTReso.GetYaxis().GetNbins()):
		binReso = hist2D_MT_List[iBin][jBin].GetRMS()/hist2D_MT_List[iBin][jBin].GetMean()
		hist2D_sdVar_jetPtMaxdPhi_MTReso.SetBinContent(iBin+1, jBin+1, binReso)
		if minReso > binReso:
			iBinMinsdVar = iBin
			jBinMinjetPt = jBin
			minReso = binReso
			print(iBin, jBin, hist2D_sdVar_jetPtMaxdPhi_MTReso.GetXaxis().GetBinLowEdge(iBin), hist2D_sdVar_jetPtMaxdPhi_MTReso.GetYaxis().GetBinLowEdge(jBin),binReso)

binLabels = ["0", "1", "2","12", "3","13","23","123"]

hist_jetCode.SetStats(0)
hist_jetCodePassesptFracCut.SetStats(0)
for x in range(1,9):
	hist_jetCode.GetXaxis().SetBinLabel(x,binLabels[x-1])
	hist_jetCodePassesptFracCut.GetXaxis().SetBinLabel(x,binLabels[x-1])


c1 = rt.TCanvas('c1','c1',1500,1000)
c1.SetRightMargin(0.13)
hist2D_sdVar_jetPtMaxdPhi_MTReso.Draw("colz")
c1.SaveAs("2D_MTResolutionZoom6.png")

drawHistos([hist_sdVar_13_MTReso],"MTReso_sdVar13")
drawHistos([hist_jetPt_maxDphi_MTReso],"MTReso_jetPt_dPhimax")

drawHistos([hist_MT_ptFracCut,hist_MT_sdVar13List[iBinMin],hist_MT_jetPtmaxDphiList[iBinMinjetPt],hist2D_MT_List[iBinMinsdVar][jBinMinjetPt]],"MT_distos_cuts")
drawHistos([hist_MT_12,hist_MT_123,hist_MT_ptFracCut],"MT_distos")
drawHistos([hist_MT_delta,hist_MT_delta_12,hist_MT_delta_123],"MT_delta")
drawHistos([hist_jetCode,hist_jetCodePassesptFracCut],"jetCode")
print("HistName | Mean | StdDev | frac.StdDev")
print(hist_MT_12.GetName() + " | " + str(hist_MT_12.GetMean()) + " | " + str(hist_MT_12.GetRMS()) + " | " + str(hist_MT_12.GetRMS()/hist_MT_12.GetMean()))
print(hist_MT_123.GetName() + " | " + str(hist_MT_123.GetMean()) + " | " + str(hist_MT_123.GetRMS()) + " | " + str(hist_MT_123.GetRMS()/hist_MT_123.GetMean()))
print(hist_MT_perf.GetName() + " | " + str(hist_MT_perf.GetMean()) + " | " + str(hist_MT_perf.GetRMS()) + " | " + str(hist_MT_perf.GetRMS()/hist_MT_perf.GetMean()))
print(hist_MT_ideal.GetName() + " | " + str(hist_MT_ideal.GetMean()) + " | " + str(hist_MT_ideal.GetRMS()) + " | " + str(hist_MT_ideal.GetRMS()/hist_MT_ideal.GetMean()))
print(hist_MT_ptFracCut.GetName() + " | " + str(hist_MT_ptFracCut.GetMean()) + " | " + str(hist_MT_ptFracCut.GetRMS()) + " | " + str(hist_MT_ptFracCut.GetRMS()/hist_MT_ptFracCut.GetMean()))
print("Num Events with 2 or more jets = " + str(nEventsWith2OrMoreJets))
print("Num Events where leading two jet Pt are both above 170 GeV = " + str(nEventsWithLeadingTwoJetsHavingPtGreaterThan170))
print("Num Events with MET to MT ratio above .15 = " + str(nEventsWithMETMTRatioGreaterThanp15))
print("Num Events without leptons = " + str(nEventsPassedPreSelection))

inFile.Close()
