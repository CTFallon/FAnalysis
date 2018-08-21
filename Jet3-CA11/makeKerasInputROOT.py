import ROOT as rt
import sys
import numpy as np
from array import array
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
# even tough this is called "makeKerasInputROOT.py" I've modified it to instead create a root file for each type of jet construction we're considering
# specifically, reconstruction using just jet 1, jets 12, jets 13, jets 123, or jets 1234
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
treeOld = inFile.Get("TreeMaker2/PreSelection")
treeOld.SetBranchStatus("*",0)
treeOld.SetBranchStatus("*AK8*",1)
treeOld.SetBranchStatus("MET",1)
treeOld.SetBranchStatus("METPhi",1)
treeOld.SetBranchStatus("Electrons",1)
treeOld.SetBranchStatus("Muons",1)
treeOld.SetBranchStatus("GenParticles",1)
treeOld.SetBranchStatus("GenParticles_PdgId",1)
treeOld.SetBranchStatus("GenParticles_ParentIdx",1)

tree = treeOld.CloneTree(0)

DR12 = array('d',[0.0])
DR13 = array('d',[0.0])
DR14 = array('d',[0.0])
DR23 = array('d',[0.0])
DR24 = array('d',[0.0])
DR34 = array('d',[0.0])

tree.Branch('deltaR12', DR12, 'deltaR12/D')
tree.Branch('deltaR13', DR13, 'deltaR13/D')
tree.Branch('deltaR14', DR14, 'deltaR14/D')
tree.Branch('deltaR23', DR23, 'deltaR23/D')
tree.Branch('deltaR24', DR24, 'deltaR24/D')
tree.Branch('deltaR34', DR34, 'deltaR34/D')

outputFile = rt.TFile("TMVATrainingInput.root","recreate")

tree_2jets_Jet1 = tree.CloneTree(0)
tree_2jets_Jet1.SetName("2jets_jet1")
tree_2jets_Jet12 = tree.CloneTree(0)
tree_2jets_Jet12.SetName("2jets_jet12")
tree_2jets_JetOther = tree.CloneTree(0)
tree_2jets_JetOther.SetName("2jets_jetOther")
treeDict2jets = {'0':tree_2jets_JetOther,
'1':tree_2jets_Jet1,
'2':tree_2jets_JetOther,
'3':tree_2jets_JetOther,
'4':tree_2jets_JetOther,
'12':tree_2jets_Jet12,
'13':tree_2jets_JetOther,
'14':tree_2jets_JetOther,
'23':tree_2jets_JetOther,
'24':tree_2jets_JetOther,
'34':tree_2jets_JetOther,
'123':tree_2jets_JetOther,
'124':tree_2jets_JetOther,
'134':tree_2jets_JetOther,
'234':tree_2jets_JetOther,
'1234':tree_2jets_JetOther}

tree_3jets_Jet12 = tree.CloneTree(0)
tree_3jets_Jet12.SetName("3jets_jet12")
tree_3jets_Jet13 = tree.CloneTree(0)
tree_3jets_Jet13.SetName("3jets_jet13")
tree_3jets_Jet123 = tree.CloneTree(0)
tree_3jets_Jet123.SetName("3jets_jet123")
tree_3jets_JetOther = tree.CloneTree(0)
tree_3jets_JetOther.SetName("3jets_jetOther")
treeDict3jets = {'0':tree_3jets_JetOther,
'1':tree_3jets_JetOther,
'2':tree_3jets_JetOther,
'3':tree_3jets_JetOther,
'4':tree_3jets_JetOther,
'12':tree_3jets_Jet12,
'13':tree_3jets_Jet13,
'14':tree_3jets_JetOther,
'23':tree_3jets_JetOther,
'24':tree_3jets_JetOther,
'34':tree_3jets_JetOther,
'123':tree_3jets_Jet123,
'124':tree_3jets_JetOther,
'134':tree_3jets_JetOther,
'234':tree_3jets_JetOther,
'1234':tree_3jets_JetOther}

tree_4jets_Jet12 = tree.CloneTree(0)
tree_4jets_Jet12.SetName("4jets_jet12")
tree_4jets_Jet13 = tree.CloneTree(0)
tree_4jets_Jet13.SetName("4jets_jet13")
tree_4jets_Jet14 = tree.CloneTree(0)
tree_4jets_Jet14.SetName("4jets_jet14")
tree_4jets_Jet123 = tree.CloneTree(0)
tree_4jets_Jet123.SetName("4jets_jet123")
tree_4jets_Jet124 = tree.CloneTree(0)
tree_4jets_Jet124.SetName("4jets_jet124")
tree_4jets_Jet134 = tree.CloneTree(0)
tree_4jets_Jet134.SetName("4jets_jet134")
tree_4jets_Jet234 = tree.CloneTree(0)
tree_4jets_Jet234.SetName("4jets_jet234")
tree_4jets_Jet1234 = tree.CloneTree(0)
tree_4jets_Jet1234.SetName("4jets_jet1234")
tree_4jets_JetOther = tree.CloneTree(0)
tree_4jets_JetOther.SetName("4jets_jetOther")
treeDict4jets = {'0':tree_4jets_JetOther,
'1':tree_4jets_JetOther,
'2':tree_4jets_JetOther,
'3':tree_4jets_JetOther,
'4':tree_4jets_JetOther,
'12':tree_4jets_Jet12,
'13':tree_4jets_Jet13,
'14':tree_4jets_Jet14,
'23':tree_4jets_JetOther,
'24':tree_4jets_JetOther,
'34':tree_4jets_JetOther,
'123':tree_4jets_Jet123,
'124':tree_4jets_Jet124,
'134':tree_4jets_Jet134,
'234':tree_4jets_Jet234,
'1234':tree_4jets_Jet1234}

nEvents = treeOld.GetEntries()
print("Total Events: " + str(nEvents))
nEventsPassedPreSelection = 0



jetComboBestMassDict = {'0':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboHVQuarksDict = {'0':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetMETXDict = {'MET1':0,'MET2':0,'MET3':0,'MET4':0}
jetComboBestMassAndMETXDict = {
'0':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'1':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'2':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'3':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'4':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'12':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'13':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'14':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'23':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'24':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'34':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'123':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'124':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'134':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'234':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
'1234':{'MET1':0,'MET2':0,'MET3':0,'MET4':0}
}
NJetDict = {'2':0,'3':0,'4+':0}
jetComboBestMassAndNJetsDict = {
'0':{'2':0,'3':0,'4+':0},
'1':{'2':0,'3':0,'4+':0},
'2':{'2':0,'3':0,'4+':0},
'3':{'2':0,'3':0,'4+':0},
'4':{'2':0,'3':0,'4+':0},
'12':{'2':0,'3':0,'4+':0},
'13':{'2':0,'3':0,'4+':0},
'14':{'2':0,'3':0,'4+':0},
'23':{'2':0,'3':0,'4+':0},
'24':{'2':0,'3':0,'4+':0},
'34':{'2':0,'3':0,'4+':0},
'123':{'2':0,'3':0,'4+':0},
'124':{'2':0,'3':0,'4+':0},
'134':{'2':0,'3':0,'4+':0},
'234':{'2':0,'3':0,'4+':0},
'1234':{'2':0,'3':0,'4+':0}
}
jetComboHVQuarksAndNJetsDict = {
'0':{'2':0,'3':0,'4+':0},
'1':{'2':0,'3':0,'4+':0},
'2':{'2':0,'3':0,'4+':0},
'3':{'2':0,'3':0,'4+':0},
'4':{'2':0,'3':0,'4+':0},
'12':{'2':0,'3':0,'4+':0},
'13':{'2':0,'3':0,'4+':0},
'14':{'2':0,'3':0,'4+':0},
'23':{'2':0,'3':0,'4+':0},
'24':{'2':0,'3':0,'4+':0},
'34':{'2':0,'3':0,'4+':0},
'123':{'2':0,'3':0,'4+':0},
'124':{'2':0,'3':0,'4+':0},
'134':{'2':0,'3':0,'4+':0},
'234':{'2':0,'3':0,'4+':0},
'1234':{'2':0,'3':0,'4+':0}
}

IndexToMET = {0:'MET1',1:'MET2',2:'MET3',3:'MET4'}
jetComboToIndex = {'0':0,'1':1,'2':2,'3':3,'4':4,'12':5,'13':6,'14':7,'23':8,'24':9,'34':10,'123':11,'124':12,'134':13,'234':14,'1234':15}

for iEvent in range(nEvents):
	
	jet0 = np.zeros(16, dtype=float)
	#print(jet0)

	jetComboMassDict = {'0':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}

	DR12[0] = 0.
	DR13[0] = 0.
	DR14[0] = 0.
	DR23[0] = 0.
	DR24[0] = 0.
	DR34[0] = 0.

	treeOld.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	#jetCollection = treeOld.GenJetsAK8
	jetCollection = treeOld.JetsAK8

	
	# PreSelection Cuts
	if not (len(treeOld.JetsAK8)>=2): #moreThan2Jets
		continue
	if not ((treeOld.JetsAK8[0].Pt() > 170.0) and (treeOld.JetsAK8[1].Pt() > 170.0)): #leadingJetsPtOver170
		continue
	if not (treeOld.MET/treeOld.MT_AK8 > 0.15): #METoverMTgreaterThan0p15
		continue
	if not ((len(treeOld.Electrons) + len(treeOld.Muons)) == 0): #leptonNotPresent
		continue
	nEventsPassedPreSelection += 1
	# tag how many jets this event has.
	nJetKey = "Hello"
	if len(treeOld.JetsAK8) == 2:
		nJetKey = '2'
	elif len(treeOld.JetsAK8) == 3:
		nJetKey = '3'
	elif len(treeOld.JetsAK8) > 3:
		nJetKey = '4+'
	else:
		print("Number of jets is not 2, 3, 4, or greater!")
	
	#calculat the deltaphi min of the first four leading pT jets, and tag which jet is most cloesly aligned with the MET
	DeltaPhi3 = 99
	DeltaPhi4 = 99
	if len(treeOld.JetsAK8) >= 3:
		DeltaPhi3 = absDphi(treeOld.JetsAK8[2].Phi(),treeOld.METPhi)
	if len(treeOld.JetsAK8) > 3:
		DeltaPhi4 = absDphi(treeOld.JetsAK8[3].Phi(),treeOld.METPhi)
	DeltaPhiList = [treeOld.DeltaPhi1_AK8, treeOld.DeltaPhi2_AK8, DeltaPhi3,DeltaPhi4]
	DeltaPhiMin = min(DeltaPhiList)
	METXkey = IndexToMET[DeltaPhiList.index(DeltaPhiMin)]
	#find the HV quarks, find which particles have the HVQuarks as ((N-)grand)parents, and figure out which AK8 jet they belong too
	nHVquark = 0
	nbarHVquark = 0
	jetsWithHVQuarks = [0,0,0,0]
	for iPart in range(len(treeOld.GenParticles)):
		if treeOld.GenParticles_PdgId[iPart] == 4900101:
			hvQuark = treeOld.GenParticles[iPart]
			nHVquark += 1
		elif treeOld.GenParticles_PdgId[iPart] == -4900101:
			barhvQuark = treeOld.GenParticles[iPart]
			nbarHVquark += 1
		else:
			iPartParentage = []
			parentId = treeOld.GenParticles_ParentIdx[iPart]
			while not (parentId == -1):
				iPartParentage.append(treeOld.GenParticles_PdgId[parentId])
				parentId = treeOld.GenParticles_ParentIdx[parentId]
		if iPart > 2  and ((4900101 in iPartParentage) or (-4900101 in iPartParentage)):
			for iJet in range(min(len(treeOld.JetsAK8),4)):
				if (not jetsWithHVQuarks[iJet]) and (treeOld.JetsAK8[iJet].DeltaR(treeOld.GenParticles[iPart]) < 0.8):
					jetsWithHVQuarks[iJet] = 1
	if (nHVquark != 1) or (nbarHVquark != 1): # skip events that don't have only 2 unstable HV quarks
		continue
	NJetDict[nJetKey] += 1
	jetMETXDict[METXkey] += 1


	

	#figure out which particles are decendants of HVquakrs (+/- 4900101)
	#figure out which jet these particles belong too
	listOfJetsWithParticlesFromHVQuarks = []
	jetsHVQuarksKey = ''
	for x in range(4):
		if jetsWithHVQuarks[x]:
			listOfJetsWithParticlesFromHVQuarks.append(treeOld.JetsAK8[x])
			jetsHVQuarksKey += str(x+1)
	if jetsHVQuarksKey == '':
		print("No jets!")
		jetsHVQuarksKey = '0'
		print(jetsHVQuarksKey)
	jetComboHVQuarksDict[jetsHVQuarksKey] += 1
	MTofPesudoDecayedHVQuarks = trans_mass_Njet([hvQuark*f_rvis, barhvQuark*f_rvis], treeOld.MET,treeOld.METPhi)
	MTofHVQuarksProducts = trans_mass_Njet(listOfJetsWithParticlesFromHVQuarks,treeOld.MET,treeOld.METPhi)
	#create all of the MT masses (should be 11 of them, plus HVJets means 12 total)
	jetComboMassDict['0'] = trans_mass_Njet([],treeOld.MET,treeOld.METPhi)
	jetComboMassDict['1'] = trans_mass_Njet([jetCollection[0]],treeOld.MET,treeOld.METPhi)
	jetComboMassDict['2'] = trans_mass_Njet([jetCollection[1]],treeOld.MET,treeOld.METPhi)
	jetComboMassDict['12'] = trans_mass_Njet([jetCollection[0], jetCollection[1]],treeOld.MET,treeOld.METPhi)
	DR12[0] = jetCollection[0].DeltaR(jetCollection[1])
	if len(jetCollection) >= 3:
		jetComboMassDict['3'] = trans_mass_Njet([jetCollection[2]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['13'] = trans_mass_Njet([jetCollection[0], jetCollection[2]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['23'] = trans_mass_Njet([jetCollection[1], jetCollection[2]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['123'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2]],treeOld.MET,treeOld.METPhi)
		DR13[0] = jetCollection[0].DeltaR(jetCollection[2])
		DR23[0] = jetCollection[1].DeltaR(jetCollection[2])
	if len(jetCollection) >= 4:
		jetComboMassDict['4'] = trans_mass_Njet([jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['14'] = trans_mass_Njet([jetCollection[0], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['24'] = trans_mass_Njet([jetCollection[1], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['34'] = trans_mass_Njet([jetCollection[2], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['124'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['134'] = trans_mass_Njet([jetCollection[0], jetCollection[2], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['234'] = trans_mass_Njet([jetCollection[1], jetCollection[2], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		jetComboMassDict['1234'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2], jetCollection[3]],treeOld.MET,treeOld.METPhi)
		DR14[0] = jetCollection[0].DeltaR(jetCollection[3])
		DR24[0] = jetCollection[1].DeltaR(jetCollection[3])
		DR34[0] = jetCollection[2].DeltaR(jetCollection[3])

	massDiffList = {key : abs(x - MTofHVQuarksProducts) for (key,x) in jetComboMassDict.items()}
	BestMassKey = min(massDiffList,key=massDiffList.get)
	jetComboBestMassDict[BestMassKey] += 1
	jetComboBestMassAndMETXDict[BestMassKey][METXkey] += 1
	jetComboBestMassAndNJetsDict[BestMassKey][nJetKey] += 1
	jetComboHVQuarksAndNJetsDict[jetsHVQuarksKey][nJetKey] += 1
	if len(treeOld.JetsAK8) == 2:
		treeDict2jets[BestMassKey].Fill()
	elif len(treeOld.JetsAK8) == 3:
		treeDict3jets[BestMassKey].Fill()
	elif len(treeOld.JetsAK8) >= 4:
		treeDict4jets[BestMassKey].Fill()
	else:
		print("Num jets doesnt make sense!")


print("Preselection " +str(nEventsPassedPreSelection))
print("Total " + str(nEvents))

print("BestJetCombo nEventsinTree")
for (key,value) in treeDict2jets.items():
	print(key + " " + str(value.GetEntries()))
for (key,value) in treeDict3jets.items():
	print(key + " " + str(value.GetEntries()))
for (key,value) in treeDict4jets.items():
	print(key + " " + str(value.GetEntries()))

outputFile.Write()
outputFile.Close()


inFile.Close()
