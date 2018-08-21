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
alpha = '0p5' #['0p1', '0p2', '0p5', '1']
f_rvis = 1.0 - 0.3

inFile = rt.TFile.Open("root://cmsxrootd.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV14/"+createInFileName(mZprime, mDark,rinv,alpha))
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
nEventsWith2OrMoreJets = 0
nEventsWithLeadingTwoJetsHavingPtGreaterThan170 = 0
nEventsWithMETMTRatioGreaterThanp15 = 0
nEventsPassedPreSelection = 0


# initilize histograms

hist_AK8DiJetMt = rt.TH1F("AK8DiJetMt",";M;count/a.u.",100,0,5000)
#hist_CA11DiJetMt = rt.TH1F("CA11DiJetMt",";M;count/a.u.",100,0,5000)
hist_HVQuarksMt = rt.TH1F("HVQuarksMt",";M;count/a.u.",100,0,5000)
hist_ZprimeMt = rt.TH1F("ZprimeMt",";M;count/a.u.",100,0,5000)
hist_AK8DiJetMjj = rt.TH1F("AK8DiJetMjj",";M;count/a.u.",100,0,5000)
#hist_CA11DiJetMjj = rt.TH1F("CA11DiJetMjj",";M;count/a.u.",100,0,5000)

hist_AK8AllJetsMt = rt.TH1F("AK8AllJetMt",";M;count/a.u.",100,0,5000)
#hist_CA11AllJetsMt = rt.TH1F("CA11AllJetMt",";M;count/a.u.",100,0,5000)

hist_AK8JetsWithHVParts = rt.TH1F("AK8JetsWithHVParts",";Jets With HV Particles;count/a.u.", 33, 0, 33)
#hist_CA11JetsWithHVParts = rt.TH1F("CA11JetsWithHVParts",";Jets With HV Particles;count/a.u.", 33, 0, 33)

hist_AK8JetsWithHVPartsvsNjets = rt.TH2F("2D_HVJetsvsnAK8Jets",";HVJets;nAK8Jets", 33, 0, 33, 5, 2, 7)
#hist_CA11JetsWithHVPartsvsNjets = rt.TH2F("2D_HVJetsvsnCA11Jets",";HVJets;nCA11Jets", 33, 0, 33, 6, 1, 7)

#hist_nAK8_vs_nCA11 = rt.TH2F("2D_AK8vsCA11",";nAK8;nCA11",5,2,7,6,1,7)
hist_nAK8 = rt.TH1F("nAK8",";nJets;",6,1,7)
#hist_nCA11 = rt.TH1F("nCA11",";nJets;",6,1,7)


nBinsInFracPlots = 20
#hist2d_CA11_nfracHV = rt.TH2F("nCA11JetsvsnFracPartsHV","Fractional Number of Particles per Jet;CA11 Jet Number;% of HV-daughters",5, 0, 5,nBinsInFracPlots,0,1.01)
#hist2d_CA11_nfracHV.SetStats(0)
hist2d_AK8_nfracHV = rt.TH2F("nAK8JetsvsnFracPartsHV","Fractional Number of Particles per Jet;AK8 Jet Number;% of HV-daughters", 3, 0, 3,nBinsInFracPlots,0,1.01)
hist2d_AK8_nfracHV.SetStats(0)

#hist2d_CA11_ptfracHV = rt.TH2F("nCA11JetsvsotFracPartsHV","Fractional pT per Jet;CA11 Jet Number;% of pT from HV-daughters",5, 0, 5,nBinsInFracPlots,0,1.01)
#hist2d_CA11_ptfracHV.SetStats(0)
hist2d_AK8_ptfracHV = rt.TH2F("nAK8JetsvsotFracPartsHV","Fractional pT per Jet (alpha = " + alpha+");AK8 Jet Number;% of pT from HV-daughters", 3, 0, 3,nBinsInFracPlots,0,1.01)
hist2d_AK8_ptfracHV.SetStats(0)

hist2d_AK8_jetPt_ptSMPartsFromHVparts = rt.TH2F("AK8_jetPt_ptSMPartsFromHVparts",";jetPt;SMPartsPt",100,0,2000,100,0,2500)


histList = []

nPlotsMade = 0


mat_JetsHVParts_nJets = [[0 for x in range(32)] for y in range(10)]

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))


	jetCollection = tree.JetsAK8
	nAK8Jets = len(jetCollection)
	#nCA11Jets = len(tree.JetsCA11)
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


	numberOfDaughtersAParticleHas = [0 for x in range(len(tree.GenParticles))]
	for iPart in range(len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if iParent != -1: 
			numberOfDaughtersAParticleHas[iParent] += 1

	AK8jetsWithHVDecendants = ["0","0","0","0","0"]
	#CA11jetsWithHVDecendants = ["0","0","0","0","0"]
	AK8_nHVParts = [0,0,0,0,0]
	#CA11_nHVParts = [0,0,0,0,0]
	AK8_nParts = [0,0,0,0,0]
	#CA11_nParts = [0,0,0,0,0]
	AK8_ptHVParts = [0,0,0,0,0]
	#CA11_ptHVParts = [0,0,0,0,0]
	AK8_ptParts = [0,0,0,0,0]
	#CA11_ptParts = [0,0,0,0,0]
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
		if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or (isFromHVQuark[iParent]):
			isFromHVQuark[iPart] = 1
		#finding what Jet a particle is in, only if it decends from a HVQuark:
		for iJet in range(min(len(tree.JetsAK8),5)):
			if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
				AK8_ptParts[-iJet-1] += tree.GenParticles[iPart].Pt()
				AK8_nParts[-iJet-1] += 1
		#for iJet in range(min(len(tree.JetsCA11),5)):
		#	if tree.JetsCA11[iJet].DeltaR(tree.GenParticles[iPart]) < 1.1 and not numberOfDaughtersAParticleHas[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
		#		CA11_ptParts[-iJet-1] += tree.GenParticles[iPart].Pt()
		#		CA11_nParts[-iJet-1] += 1
		if isFromHVQuark[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
			for iJet in range(min(len(tree.JetsAK8),5)):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart]:
					AK8jetsWithHVDecendants[-iJet-1] = "1"
					AK8_ptHVParts[-iJet-1] += tree.GenParticles[iPart].Pt()
					AK8_nHVParts[-iJet-1] += 1
		#	for iJet in range(min(len(tree.JetsCA11),5)):
		#		if tree.JetsCA11[iJet].DeltaR(tree.GenParticles[iPart]) < 1.1 and not numberOfDaughtersAParticleHas[iPart]:
		#			CA11jetsWithHVDecendants[-iJet-1] = "1"
		#			CA11_ptHVParts[-iJet-1] += tree.GenParticles[iPart].Pt()
		#			CA11_nHVParts[-iJet-1] += 1
	mat_JetsHVParts_nJets[nAK8Jets][int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2)] += 1
	hist_AK8JetsWithHVParts.Fill(int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2))
	hist_AK8JetsWithHVPartsvsNjets.Fill(int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2), nAK8Jets)
	#hist_CA11JetsWithHVParts.Fill(int(CA11jetsWithHVDecendants[0]+CA11jetsWithHVDecendants[1]+CA11jetsWithHVDecendants[2]+CA11jetsWithHVDecendants[3]+CA11jetsWithHVDecendants[4],2))
	#hist_CA11JetsWithHVPartsvsNjets.Fill(int(CA11jetsWithHVDecendants[0]+CA11jetsWithHVDecendants[1]+CA11jetsWithHVDecendants[2]+CA11jetsWithHVDecendants[3]+CA11jetsWithHVDecendants[4],2), nCA11Jets)

	for x in range(5):
		if len(tree.JetsAK8) > x:
			hist2d_AK8_jetPt_ptSMPartsFromHVparts.Fill(tree.JetsAK8[x].Pt(),AK8_ptParts[-x-1])
		#if CA11_nParts[-x-1] != 0:
		#	hist2d_CA11_nfracHV.Fill(x, float(CA11_nHVParts[-x-1])/float(CA11_nParts[-x-1]))
		if AK8_nParts[-x-1] != 0:
			hist2d_AK8_nfracHV.Fill(x, float(AK8_nHVParts[-x-1])/float(AK8_nParts[-x-1]))
		#if CA11_ptParts[-x-1] != 0:
		#	hist2d_CA11_ptfracHV.Fill(x, float(CA11_ptHVParts[-x-1])/float(CA11_ptParts[-x-1]))
		if AK8_ptParts[-x-1] != 0:
			hist2d_AK8_ptfracHV.Fill(x, float(AK8_ptHVParts[-x-1])/float(AK8_ptParts[-x-1]))
	#print(isFromHVQuark)
	# figureout which order the jets rae in distance from METPhi
	deltaPhiList = [tree.DeltaPhi1, tree.DeltaPhi2, tree.DeltaPhi3]
	IdxList = [0,1,2]
	mindPhiIdx = deltaPhiList.index(min(deltaPhiList))
	IdxList.remove(mindPhiIdx)
	maxdPhiIdx = deltaPhiList.index(max(deltaPhiList))
	IdxList.remove(maxdPhiIdx)
	middPhiIdx = IdxList[0]


	if nPlotsMade < 20 and len(tree.JetsAK8) >= 3:
		nPlotsMade += 1
		drawEventEtaPhiPlot(tree.JetsAK8, [], tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, isFromHVQuark, iEvent)
	
	#if len(tree.JetsCA11) == 5:
	#	drawEventEtaPhiPlot(tree.JetsAK8, tree.JetsCA11, tree.GenParticles, tree.GenParticles_PdgId, tree.METPhi, isFromHVQuark, iEvent)
		
	if len(listOfHVQuarks) == 2:
		hist_HVQuarksMt.Fill(trans_mass_Njet(listOfHVQuarks, tree.MET, tree.METPhi))
	hist_AK8DiJetMt.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
	hist_AK8DiJetMjj.Fill((tree.JetsAK8[0]+tree.JetsAK8[1]).M())
	#print((tree.JetsAK8[0]+tree.JetsAK8[1]).M())
	#if len(tree.JetsCA11) >= 2:
	#	hist_CA11DiJetMt.Fill(trans_mass_Njet([tree.JetsCA11[0],tree.JetsCA11[1]], tree.MET, tree.METPhi))
	#	hist_CA11DiJetMjj.Fill((tree.JetsCA11[0]+tree.JetsCA11[1]).M())

	#hist_nAK8_vs_nCA11.Fill(len(tree.JetsAK8), len(tree.JetsCA11))
	hist_nAK8.Fill(len(tree.JetsAK8))
	#hist_nCA11.Fill(len(tree.JetsCA11))
	hist_AK8AllJetsMt.Fill(trans_mass_Njet(tree.JetsAK8, tree.MET, tree.METPhi))
	#hist_CA11AllJetsMt.Fill(trans_mass_Njet(tree.JetsCA11, tree.MET, tree.METPhi))

for histo in histList:
	total = histo.GetEntries()
	if total != 0:
		histo.Scale(1/(total))

#setting bin labels for the 'binary' plot
#binLabels = ["0", "1", "2","12", "3","13","23","123", "4","14","24","124","34","134","234","1234", "5","15","25","125","35","135","235","1235","45","145","245","1245","345","1345","2345","12345"]
binLabels = ["0", "", "","12", "","","","123", "","","","","","","","1234", "","","","","","","","","","","","","","","","12345"]

hist_AK8JetsWithHVParts.SetStats(0)
for x in range(1,32):
	hist_AK8JetsWithHVParts.GetXaxis().SetBinLabel(x,binLabels[x-1])
c1 = rt.TCanvas()

#hist_nAK8_vs_nCA11.Draw("colz")
#c1.SaveAs("nAK8vsnCA11.png")
#hist2d_CA11_nfracHV.Draw("colz")
#c1.SaveAs("CA11_nfracHV.png")
hist2d_AK8_nfracHV.Draw("colz")
c1.SaveAs("AK8_nfracHV.png")
#hist2d_CA11_ptfracHV.Draw("colz")
#c1.SaveAs("CA11_ptfracHV.png")
hist2d_AK8_ptfracHV.Draw("colz")
c1.SaveAs("AK8_ptfracHV.png")
hist2d_AK8_jetPt_ptSMPartsFromHVparts.Draw("colz")
xyLine = rt.TLine(0,0,2000,2000)
xyLine.SetLineColor(rt.kRed)
xyLine.Draw("same")
c1.SaveAs("AK8_jetPt_ptSMPartsFromHVparts.png")
for x in range(32):
	print(str(mat_JetsHVParts_nJets[0][x]) + " " +str(mat_JetsHVParts_nJets[1][x]) + " " +str(mat_JetsHVParts_nJets[2][x]) + " " +str(mat_JetsHVParts_nJets[3][x]) + " " +str(mat_JetsHVParts_nJets[4][x]) + " " +str(mat_JetsHVParts_nJets[5][x]) + " " +str(mat_JetsHVParts_nJets[6][x]) + " " +str(mat_JetsHVParts_nJets[7][x]) + " " +str(mat_JetsHVParts_nJets[8][x]) + " " +str(mat_JetsHVParts_nJets[9][x]))

drawHistos([hist_nAK8],"nJets")
#drawHistos([hist_SDMass_withHV,hist_SDMass_withoutHV],"SoftDropMassLOG", True)
#drawHistos([hist_JetsWithHVParts],"JetsWithHVParts")
#drawHistos([hist_SDMass_withHV,hist_SDMass_withoutHV],"SoftDropMass")
drawHistos([hist_AK8JetsWithHVParts],"JetsWithHVPartsLOG", True)
#drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt, hist_HVQuarksMt], "DiJetMt")
#drawHistos([hist_AK8DiJetMt, hist_AK8DiJetMjj,hist_AK8AllJetsMt],"AK8")
#drawHistos([hist_CA11DiJetMt, hist_CA11DiJetMjj,hist_CA11AllJetsMt],"CA11")
#drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt,hist_AK8AllJetsMt,hist_CA11AllJetsMt], "allJetsLOG", True)
#drawHistos([hist_AK8DiJetMt, hist_CA11DiJetMt,hist_AK8AllJetsMt,hist_CA11AllJetsMt], "allJets")



print("Num Events with 2 or more jets = " + str(nEventsWith2OrMoreJets))
print("Num Events where leading two jet Pt are both above 170 GeV = " + str(nEventsWithLeadingTwoJetsHavingPtGreaterThan170))
print("Num Events with MET to MT ratio above .15 = " + str(nEventsWithMETMTRatioGreaterThanp15))
print("Num Events without leptons = " + str(nEventsPassedPreSelection))

inFile.Close()
