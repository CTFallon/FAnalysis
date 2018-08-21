import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

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

def isNear(val1, val2, delta):
	return (abs(val1-val2)<delta)

rt.gStyle.SetOptStat(0)

def absDphi(phi1, phi2):
	dphi = phi1-phi2
	while dphi >= rt.TMath.Pi():
		dphi -= 2*rt.TMath.Pi()
	while dphi < -rt.TMath.Pi():
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
nEventsPassedPreSelection = 0

# initilize histograms
outputFile = rt.TFile("nTuplewithVariousMTCalculations.root","recreate")

jetComboBestMassDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboHVQuarksDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}

jetMETXDict = {'MET1':0,'MET2':0,'MET3':0,'MET4':0}

IndexToMET = {0:'MET1',1:'MET2',2:'MET3',3:'MET4'}

for iEvent in range(nEvents):
	jetComboMassDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
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
	
	#calculat the deltaphi min of the first four leading pT jets, and tag which jet is most cloesly aligned with the MET
	DeltaPhi3 = 99
	DeltaPhi4 = 99
	if len(tree.JetsAK8) >= 3:
		DeltaPhi3 = absDphi(tree.JetsAK8[2].Phi(),tree.METPhi)
	if len(tree.JetsAK8) > 3:
		DeltaPhi4 = absDphi(tree.JetsAK8[3].Phi(),tree.METPhi)
	DeltaPhiList = [tree.DeltaPhi1_AK8, tree.DeltaPhi2_AK8, DeltaPhi3,DeltaPhi4]
	DeltaPhiMin = min(DeltaPhiList)
	#if not (DeltaPhiMin < 0.6): # final selection cut. Maybe the deltaPhiX_METY distributions will differentiate in these events...
		#continue
	METXkey = IndexToMET[DeltaPhiList.index(DeltaPhiMin)]
	#if not (METXkey == 'MET3'): # for tseeing the effects METX has on these plots. Make sure to chance the createImagename function above
	#	continue

	#find the HV quarks, find which particles have the HVQuarks as ((N-)grand)parents, and figure out which AK8 jet they belong too
	nHVquark = 0
	nbarHVquark = 0
	jetsWithHVQuarks = [0,0,0,0]
	for iPart in range(len(tree.GenParticles)):
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
		if iPart > 2  and ((4900101 in iPartParentage) or (-4900101 in iPartParentage)):
			for iJet in range(min(len(tree.JetsAK8),4)):
				if (not jetsWithHVQuarks[iJet]) and (tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8):
					jetsWithHVQuarks[iJet] = 1
	if (nHVquark != 1) or (nbarHVquark != 1): # skip events that don't have only 2 unstable HV quarks
		continue
	jetMETXDict[METXkey] += 1


	

	#figure out which particles are decendants of HVquakrs (+/- 4900101)
	#figure out which jet these particles belong too
	listOfJetsWithParticlesFromHVQuarks = []
	jetsHVQuarksKey = ''
	for x in range(4):
		if jetsWithHVQuarks[x]:
			listOfJetsWithParticlesFromHVQuarks.append(tree.JetsAK8[x])
			jetsHVQuarksKey += str(x+1)
	jetComboHVQuarksDict[jetsHVQuarksKey] += 1
	MTofPesudoDecayedHVQuarks = trans_mass_Njet([hvQuark*f_rvis, barhvQuark*f_rvis], tree.MET,tree.METPhi)
	MTofHVQuarksProducts = trans_mass_Njet(listOfJetsWithParticlesFromHVQuarks,tree.MET,tree.METPhi)
	#create all of the MT masses (should be 16 of them, plus HVJets means 17 total)
	jetComboMassDict['12'] = trans_mass_Njet([jetCollection[0], jetCollection[1]],tree.MET,tree.METPhi)
	jetComboMassDict[''] = trans_mass_Njet([],tree.MET,tree.METPhi)
	jetComboMassDict['1'] = trans_mass_Njet([jetCollection[0]],tree.MET,tree.METPhi)
	jetComboMassDict['2'] = trans_mass_Njet([jetCollection[1]],tree.MET,tree.METPhi)
	if len(jetCollection) >= 3:
		jetComboMassDict['3'] = trans_mass_Njet([jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['13'] = trans_mass_Njet([jetCollection[0], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['23'] = trans_mass_Njet([jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['123'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
	if len(jetCollection) >= 4:
		jetComboMassDict['4'] = trans_mass_Njet([jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['14'] = trans_mass_Njet([jetCollection[0], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['24'] = trans_mass_Njet([jetCollection[1], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['34'] = trans_mass_Njet([jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['124'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['134'] = trans_mass_Njet([jetCollection[0], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['234'] = trans_mass_Njet([jetCollection[1], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
		jetComboMassDict['1234'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2], jetCollection[3]],tree.MET,tree.METPhi)
	massDiffList = {key : abs(x - MTofHVQuarksProducts) for (key,x) in jetComboMassDict.items()}
	BestMassKey = min(massDiffList,key=massDiffList.get)
	jetComboBestMassDict[BestMassKey] += 1
	jetComboBestMassAndMETXDict[BestMassKey][METXkey] += 1
	jetComboBestMassAndNJetsDict[BestMassKey][nJetKey] += 1
	jetComboHVQuarksAndNJetsDict[jetsHVQuarksKey][nJetKey] += 1


outputFile.Write()
outputFile.Close()


inFile.Close()
