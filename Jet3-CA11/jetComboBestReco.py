import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutations_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

def mt_calc(jet1, jet2, met, metPhi):
	jetSum = jet1 + jet2
	mjj2 = jetSum.M2()
	term1 = rt.TMath.Sqrt(mjj2 + jetSum.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-jetSum.Phi())*jetSum.Pt()*met
	return rt.TMath.Sqrt(mjj2 + 2*(term1 - term2))

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

def isNear(val1, val2, delta):
	return (abs(val1-val2)<delta)

rt.gStyle.SetOptStat(0)

def absDphi(phi1, phi2):
	dphi = phi1-phi2
	if dphi > rt.TMath.Pi():
		dphi -= 2*rt.TMath.Pi()
	elif dphi < -rt.TMath.Pi():
		dphi += 2*rt.TMath.Pi()
	dphi = rt.TMath.Abs(dphi)
	return dphi

#step1, know what data set we're looking at

mZprime = 3000
mDark = '20'
rinv = '0p3'
alpha = '0p2'
f_rvis = 1.0 - 0.3

inFile = rt.TFile("../"+createInFileName(mZprime, mDark,rinv,alpha),"read")
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
jet1Aligned = 0
jet2Aligned = 0
jet3Aligned = 0
flag=0
nEventsPassedPreSelection = 0
jetTagIsGood = 0
goodTagAndMET1 = 0
goodTagAndMET2 = 0
goodTagAndMET3 = 0
jetAntiTagIsGood = 0
goodAntiTagAndMET1 = 0
goodAntiTagAndMET2 = 0
goodAntiTagAndMET3 = 0
QuarksWithSameJet = 0
badTagEvents = 0
badTagAndMET1 = 0
badTagAndMET2 = 0
badTagAndMET3 = 0

# initilize histograms

hist_MT_HVJets = rt.TH1F("MT_HVQuarkJets","HVQuarkJets;MT;a.u.",100,0,2*mZprime)
hist_MT_pseudoInvisibleHVQuarks = rt.TH1F("pseudoDecayHVQuarks","pseudoHVQuarks;MT;a.u.",100,0,2*mZprime)
hist_MT_bestReco = rt.TH1F("MT_bestReco","bestReco;MT;a.u.",100,0,2*mZprime)

hist_MT12_All = rt.TH1F("MT_12_All", "12_All;MT;a.u.",100,0,2*mZprime)
hist_MT13_All = rt.TH1F("MT_13_All", "13_All;MT;a.u.",100,0,2*mZprime)
hist_MT14_All = rt.TH1F("MT_14_All", "14_All;MT;a.u.",100,0,2*mZprime)
hist_MT23_All = rt.TH1F("MT_23_All", "23_All;MT;a.u.",100,0,2*mZprime)
hist_MT24_All = rt.TH1F("MT_24_All", "24_All;MT;a.u.",100,0,2*mZprime)
hist_MT34_All = rt.TH1F("MT_34_All", "34_All;MT;a.u.",100,0,2*mZprime)

hist_MT123_All = rt.TH1F("MT_123_All", "123_All;MT;a.u.",100,0,2*mZprime)
hist_MT124_All = rt.TH1F("MT_124_All", "124_All;MT;a.u.",100,0,2*mZprime)
hist_MT134_All = rt.TH1F("MT_134_All", "134_All;MT;a.u.",100,0,2*mZprime)
hist_MT234_All = rt.TH1F("MT_234_All", "234_All;MT;a.u.",100,0,2*mZprime)

hist_MT1234_All = rt.TH1F("MT_1234_All", "1234_All;MT;a.u.",100,0,2*mZprime)

histList = []
histList.append(hist_MT12_All)
histList.append(hist_MT13_All)
histList.append(hist_MT14_All)
histList.append(hist_MT23_All)
histList.append(hist_MT24_All)
histList.append(hist_MT34_All)
histList.append(hist_MT123_All)
histList.append(hist_MT124_All)
histList.append(hist_MT134_All)
histList.append(hist_MT234_All)
histList.append(hist_MT1234_All)
histList.append(hist_MT_HVJets)
histList.append(hist_MT_pseudoInvisibleHVQuarks)
histList.append(hist_MT_bestReco)

hist_MT12_All.SetLineColor(2)
hist_MT13_All.SetLineColor(3)
hist_MT14_All.SetLineColor(4)
hist_MT23_All.SetLineColor(5)
hist_MT24_All.SetLineColor(6)
hist_MT34_All.SetLineColor(7)
hist_MT123_All.SetLineColor(8)
hist_MT124_All.SetLineColor(9)
hist_MT134_All.SetLineColor(11)
hist_MT234_All.SetLineColor(12)
hist_MT1234_All.SetLineColor(13)
hist_MT_HVJets.SetLineColor(1)
hist_MT_pseudoInvisibleHVQuarks.SetLineColor(14)
hist_MT_bestReco.SetLineColor(15)
hist_MT_HVJets.SetLineStyle(3)
hist_MT_bestReco.SetLineStyle(2)

jetComboMassDict = {'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboBestMassDict = {'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	#jetCollection = tree.GenJetsAK8
	jetCollection = tree.JetsAK8

	
	#some calculations to determine which histos to fill
	# PreSelection Cuts
	if not (len(tree.JetsAK8)>=2): #moreThan2Jets
		continue
	if not ((tree.JetsAK8[0].Pt() > 170.0) and (tree.JetsAK8[1].Pt() > 170.0)): #leadingJetsPtOver170
		continue
	if not (tree.MET/tree.MT_AK8 > 0.15): #METoverMTgreaterThan0p15
		continue
	if not ((len(tree.Electrons) + len(tree.Muons)) == 0): #leptonNotPresent
		continue
	nEventsPassedPreSelection += 1

	DeltaPhi3 = 99
	DeltaPhi4 = 99
	if len(jetCollection) == 3:
		DeltaPhi3 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
	if len(jetCollection) > 3:
		DeltaPhi4 = rt.TMath.Abs(jetCollection[3].Phi()-tree.METPhi)

	DeltaPhiMin = min(DeltaPhi4, DeltaPhi3, tree.DeltaPhi1_AK8, tree.DeltaPhi2_AK8)

	#if not (DeltaPhiMin == DeltaPhi4):
	#	continue

	#DeltaPhi3_AK8 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
	#DeltaPhiMin_3jet = min(DeltaPhi3_AK8, tree.DeltaPhi1_AK8, tree.DeltaPhi2_AK8)
	#find the HV quarks, find which particles have the HVQuarks as ((N-)grand)parents, and figure out which AK8 jet they belong too
	nHVquark = 0
	nbarHVquark = 0
	jetsWithHVQuarks = [0,0,0,0]
	for iPart in range(len(tree.GenParticles)):
		#print(iPart, tree.GenParticles_PdgId[iPart], tree.GenParticles_ParentIdx[iPart])
		inAJet = 0
		if tree.GenParticles_PdgId[iPart] == 4900101:
			hvQuark = tree.GenParticles[iPart]
			nHVquark += 1
		elif tree.GenParticles_PdgId[iPart] == -4900101:
			barhvQuark = tree.GenParticles[iPart]
			nbarHVquark += 1
		else:
			iPartParentage = []
			parentId = tree.GenParticles_ParentIdx[iPart]
			while not (parentId == -1):
				iPartParentage.append(tree.GenParticles_PdgId[parentId])
				parentId = tree.GenParticles_ParentIdx[parentId]
		for iJet in range(min(len(tree.JetsAK8),4)):
			if (not jetsWithHVQuarks[iJet]) and (iPart > 2) and (tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) > 0.8) and ((4900101 in iPartParentage) or (-4900101 in iPartParentage)):
				jetsWithHVQuarks[iJet] = 1
	if (nHVquark != 1) or (nbarHVquark != 1): # skip events that don't have only 2 unstable HV quarks
		continue


	

	#figure out which particles are decendants of HVquakrs (+/- 4900101)
	#figure out which jet these particles belong too
	listOfJetsWithParticlesFromHVQuarks = []
	for x in range(4):
		if jetsWithHVQuarks[x]:
			listOfJetsWithParticlesFromHVQuarks.append(tree.JetsAK8[x])
	MTofPesudoDecayedHVQuarks = trans_mass_Njet([hvQuark*f_rvis, barhvQuark*f_rvis], tree.MET,tree.METPhi)
	hist_MT_pseudoInvisibleHVQuarks.Fill(MTofPesudoDecayedHVQuarks)
	MTofHVQuarksProducts = trans_mass_Njet(listOfJetsWithParticlesFromHVQuarks,tree.MET,tree.METPhi)
	hist_MT_HVJets.Fill(MTofHVQuarksProducts)
	#create all of the MT masses (should be 11 of them, plus HVJets means 12 total)
	jetComboMassDict['12'] = trans_mass_Njet([jetCollection[0], jetCollection[1]],tree.MET,tree.METPhi)
	hist_MT12_All.Fill(jetComboMassDict['12'])
	if len(jetCollection) >= 3:
		jetComboMassDict['13'] = trans_mass_Njet([jetCollection[0], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['23'] = trans_mass_Njet([jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['123'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		hist_MT13_All.Fill(jetComboMassDict['13'])
		hist_MT23_All.Fill(jetComboMassDict['23'])
		hist_MT123_All.Fill(jetComboMassDict['123'])
	if len(jetCollection) >= 4:
		jetComboMassDict['14'] = trans_mass_Njet([jetCollection[0], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['24'] = trans_mass_Njet([jetCollection[1], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['34'] = trans_mass_Njet([jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['124'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['134'] = trans_mass_Njet([jetCollection[0], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['234'] = trans_mass_Njet([jetCollection[1], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['1234'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		hist_MT14_All.Fill(jetComboMassDict['14']) 
		hist_MT24_All.Fill(jetComboMassDict['24'])
		hist_MT34_All.Fill(jetComboMassDict['34'])
		hist_MT124_All.Fill(jetComboMassDict['124'])
		hist_MT134_All.Fill(jetComboMassDict['134'])
		hist_MT234_All.Fill(jetComboMassDict['234'])
		hist_MT1234_All.Fill(jetComboMassDict['1234'])
	massDiffList = {key : abs(x - MTofHVQuarksProducts) for (key,x) in jetComboMassDict.items()}
	jetComboBestMassDict[min(massDiffList,key=massDiffList.get)] += 1
	hist_MT_bestReco.Fill(jetComboMassDict[min(massDiffList,key=massDiffList.get)])




#step 4, normalize all histograms

for histo in histList:
	area = histo.GetEntries()
	try:
		histo.Scale(1/area)
	except ZeroDivisionError:
		#print("Empty histogram: " + histo.GetName())
		x = 5
	print(histo.GetName() + " " + str(area))
	#print("StdDev of " + histo.GetName() + " is " + str(histo.GetStdDev()))
print("Preselection " +str(nEventsPassedPreSelection))
print("Total " + str(nEvents))
for (key,value) in jetComboBestMassDict.items():
	print((key + " " + str(value)) )
drawHistos([hist_MT12_All,hist_MT_HVJets,hist_MT_pseudoInvisibleHVQuarks],"MT_HVQuarksJetsPseudoDecayed",True)
drawHistos([hist_MT_HVJets,hist_MT12_All,hist_MT13_All,hist_MT14_All,hist_MT23_All,hist_MT24_All,hist_MT34_All],"MT_DiJets",True)
drawHistos([hist_MT_HVJets,hist_MT123_All,hist_MT124_All,hist_MT134_All,hist_MT234_All],"MT_TriJets",True)
drawHistos([hist_MT_HVJets,hist_MT1234_All],"MT_QuadJets",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT13_All, hist_MT14_All, hist_MT123_All, hist_MT124_All, hist_MT134_All, hist_MT1234_All],"MT_Jet1",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT23_All, hist_MT24_All, hist_MT123_All, hist_MT124_All, hist_MT234_All, hist_MT1234_All],"MT_Jet2",True)
drawHistos([hist_MT_HVJets, hist_MT13_All, hist_MT23_All, hist_MT34_All, hist_MT123_All, hist_MT134_All, hist_MT234_All, hist_MT1234_All],"MT_Jet3",True)
drawHistos([hist_MT_HVJets, hist_MT14_All, hist_MT24_All, hist_MT34_All, hist_MT124_All, hist_MT134_All, hist_MT234_All, hist_MT1234_All],"MT_Jet4",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT123_All, hist_MT1234_All,hist_MT_bestReco,hist_MT_pseudoInvisibleHVQuarks],"MT_bestReco",True)

# Draw same jets Construction Here

#draw same METX plots here



inFile.Close()
