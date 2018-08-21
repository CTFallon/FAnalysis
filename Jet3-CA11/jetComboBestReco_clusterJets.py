import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutations_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+"_3jets.png"

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
nEventsPassedPreSelection = 0
nHVParticles = 0


# initilize histograms

RList = [0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9]

histList_pt_ak8Jets = []
histList_pt_RCJets = []
histList_numHVParticlesInAK8Jets = []
histList_numHVParticlesInRCJets = []
histList_mt12_RCJets = []
histList_jetMergedMT = []

for x in range(7):
	histList_pt_ak8Jets.append(rt.TH1F("ak8Jets_mass_jet"+str(x),"ak8Jets_mass_jet"+str(x)+";mass;count/a.u.",100,0,2500))
	histList_numHVParticlesInAK8Jets.append(rt.TH1I("numInvisInAK8Jet"+str(x),"NumInvisInAK8Jet"+str(x)+";numInv;count/a.u.",100,0,100))

for irVal in range(len(RList)):
	histList_pt_RCJets.append([])
	histList_numHVParticlesInRCJets.append([])
	histList_mt12_RCJets.append(rt.TH1F("RCJets_mt12_jet"+str(RList[irVal]),"RCJets_mt12_jet"+str(RList[irVal])+";mt;count/a.u.",100,0,6000))
	histList_jetMergedMT.append(rt.TH1F("MergedJetMT"+str(RList[irVal]), "mergedJetMT"+str(RList[irVal])+";mt;count",100,0,5000))

	for x in range(7):
		histList_pt_RCJets[irVal].append(rt.TH1F("RCJets_mass_jet"+str(x)+str(RList[irVal]),"RCJets_mass_jet"+str(x)+str(RList[irVal])+";mass;count/a.u.",100,0,2500))
		histList_numHVParticlesInRCJets[irVal].append(rt.TH1I("numInvisInRCJet"+str(x)+str(RList[irVal]),"NumInvisInRCJet"+str(x)+str(RList[irVal])+";numInv;count/a.u.",100,0,100))

hist_mt12_ak8Jets = rt.TH1F("ak8Jets_mt12_jet","ak8Jets_mt12_jet;mt;count/a.u.",100,0,6000)

hist_HVinEvent = rt.TH1F("HVinEvent","HVinEvent;HV;count",100,0,100)

#jetFlowMat = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	jetCollection = tree.JetsAK8
	#print(type(jetCollection[0]))
	nAK8Jets = len(jetCollection)

	#some calculations to determine which histos to fill
	# PreSelection Cuts
	if not (len(jetCollection)>=2): #moreThan2Jets
		continue
	if not ((jetCollection[0].Pt() > 170.0) and (jetCollection[1].Pt() > 170.0)): #leadingJetsPtOver170
		continue
	if not (tree.MET/tree.MT_AK8 > 0.15): #METoverMTgreaterThan0p15
		continue
	if not ((len(tree.Electrons) + len(tree.Muons)) == 0): #leptonNotPresent
		continue
	nEventsPassedPreSelection += 1
	nHVinEvent = 0
	for iPart in range(len(tree.GenParticles)):
		if abs(tree.GenParticles_PdgId[iPart]) >= 4000000:
			nHVParticles += 1
			nHVinEvent += 1
	hist_HVinEvent.Fill(nHVinEvent)
	"""
	for iJet in range(nAK8Jets):
		nHVParts = 0
		for iPart in range(2,len(tree.GenParticles)):
			if tree.GenParticles[iPart].DeltaR(jetCollection[iJet]) < 0.8 and abs(tree.GenParticles_PdgId[iPart]) > 4000000:
				nHVParts += 1
		histList_numHVParticlesInAK8Jets[iJet].Fill(nHVParts)
	"""
	hist_mt12_ak8Jets.Fill(trans_mass_Njet([jetCollection[0],jetCollection[1]],tree.MET,tree.METPhi))


	for x in range(nAK8Jets):
		histList_pt_ak8Jets[x].Fill(jetCollection[x].Mag())
	for irVal in range(len(RList)):
		ak8JetList = []
		ak8JetList.append(jetCollection[0])
		ak8JetList.append(jetCollection[1])
		if nAK8Jets > 2:
			ak8JetList.append(jetCollection[2])
		if nAK8Jets > 3:
			ak8JetList.append(jetCollection[3])
		if nAK8Jets > 4:
			ak8JetList.append(jetCollection[4])
		if nAK8Jets > 5:
			ak8JetList.append(jetCollection[5])
		newJets = []
		jetsMergedTo12 = [0 for x in range(len(ak8JetList))]
		jetsMergedTo12[0] = 1
		jetsMergedTo12[1] = 1
		for i in range(2,len(ak8JetList)):
			dList = []
			for j in range(len(ak8JetList)):
				dTemp = min(1/(ak8JetList[i].Pt())**2,1/(ak8JetList[j].Pt())**2)
				if i != j:
					dTemp *= ((ak8JetList[i].DeltaR(ak8JetList[j]))/RList[irVal])
				dList.append(dTemp)
			if dList[0] == min(dList) or dList[1] == min(dList):
				jetsMergedTo12[i] = 1
		jetMTList = []
		for x in range(len(jetsMergedTo12)):
			if jetsMergedTo12[x] == 1:
				jetMTList.append(ak8JetList[x])
		histList_jetMergedMT[irVal].Fill(trans_mass_Njet(jetMTList,tree.MET,tree.METPhi))
				
		
		"""
		newJets = pseudoAK(ak8JetList,RList[irVal])
		newJets.sort(key=lambda x: x.Pt(), reverse=True)
		nJetsAfterNewCluster = len(newJets)
		for x in range(nJetsAfterNewCluster):
			histList_pt_RCJets[irVal][x].Fill(newJets[x].Mag())
		
		for iJet in range(nJetsAfterNewCluster):
			nHVParts = 0
			for iPart in range(2,len(tree.GenParticles)):
				if tree.GenParticles[iPart].DeltaR(newJets[iJet]) < RList[irVal] and abs(tree.GenParticles_PdgId[iPart]) > 4000000:
					nHVParts += 1
			histList_numHVParticlesInRCJets[irVal][iJet].Fill(nHVParts)
		
		if nJetsAfterNewCluster > 1:
			histList_mt12_RCJets[irVal].Fill(trans_mass_Njet([newJets[0],newJets[1]],tree.MET,tree.METPhi))
		"""	

#for X in jetFlowMat:
#		print(X)
print("Total Num Events Passed Preselection = " + str(nEventsPassedPreSelection))
print("Total Num of Genparticles that are HV = " + str(nHVParticles))
print("Average number of HV Particles per Event = " + str(float(nHVParticles)/float(nEventsPassedPreSelection)))
drawHistos(histList_pt_ak8Jets,"ak8Jets_pt",True)
drawHistos(histList_numHVParticlesInAK8Jets,"numHVinak8Jets",True)

for x in range(len(RList)):
	drawHistos(histList_pt_RCJets[x],"RCJets_mass"+str(RList[x]),True)
	drawHistos(histList_numHVParticlesInRCJets[x],"numHVinRCJets"+str(RList[x]),True)
drawHistos(histList_jetMergedMT,"MergedJetsMT", True)
drawHistos([hist_mt12_ak8Jets]+histList_mt12_RCJets,"MT12_ak8andRC",True)
drawHistos([hist_HVinEvent],"HVinEvent",True)

inFile.Close()
