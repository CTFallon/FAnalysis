import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutations_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+"_2jets.png"

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
outputFile = rt.TFile("output.root","recreate")
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

hist_MET1_deltaPhi1 = rt.TH1F("MET1_deltaPhi1", "MET1_deltaPhi1;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi1_12 = rt.TH1F("MET1_deltaPhi1_12", "MET1_deltaPhi1_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi1_123 = rt.TH1F("MET1_deltaPhi1_123", "MET1_deltaPhi1_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi2 = rt.TH1F("MET1_deltaPhi2", "MET1_deltaPhi2;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi2_12 = rt.TH1F("MET1_deltaPhi2_12", "MET1_deltaPhi2_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi2_123 = rt.TH1F("MET1_deltaPhi2_123", "MET1_deltaPhi2_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi3 = rt.TH1F("MET1_deltaPhi3", "MET1_deltaPhi3;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi3_12 = rt.TH1F("MET1_deltaPhi3_12", "MET1_deltaPhi3_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi3_123 = rt.TH1F("MET1_deltaPhi3_123", "MET1_deltaPhi3_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi4 = rt.TH1F("MET1_deltaPhi4", "MET1_deltaPhi4;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET2_deltaPhi1 = rt.TH1F("MET2_deltaPhi1", "MET2_deltaPhi1;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi1_12 = rt.TH1F("MET2_deltaPhi1_12", "MET2_deltaPhi1_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi1_123 = rt.TH1F("MET2_deltaPhi1_123", "MET2_deltaPhi1_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi2 = rt.TH1F("MET2_deltaPhi2", "MET2_deltaPhi2;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi2_12 = rt.TH1F("MET2_deltaPhi2_12", "MET2_deltaPhi2_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi2_123 = rt.TH1F("MET2_deltaPhi2_123", "MET2_deltaPhi2_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi3 = rt.TH1F("MET2_deltaPhi3", "MET2_deltaPhi3;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi3_12 = rt.TH1F("MET2_deltaPhi3_12", "MET2_deltaPhi3_12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi3_123 = rt.TH1F("MET2_deltaPhi3_123", "MET2_deltaPhi3_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi4 = rt.TH1F("MET2_deltaPhi4", "MET2_deltaPhi4;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET3_deltaPhi1 = rt.TH1F("MET3_deltaPhi1", "MET3_deltaPhi1;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi1_13 = rt.TH1F("MET3_deltaPhi1_13", "MET3_deltaPhi1_13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi1_123 = rt.TH1F("MET3_deltaPhi1_123", "MET3_deltaPhi1_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi2 = rt.TH1F("MET3_deltaPhi2", "MET3_deltaPhi2;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi2_13 = rt.TH1F("MET3_deltaPhi2_13", "MET3_deltaPhi2_13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi2_123 = rt.TH1F("MET3_deltaPhi2_123", "MET3_deltaPhi2_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi3 = rt.TH1F("MET3_deltaPhi3", "MET3_deltaPhi3;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi3_13 = rt.TH1F("MET3_deltaPhi3_13", "MET3_deltaPhi3_13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi3_123 = rt.TH1F("MET3_deltaPhi3_123", "MET3_deltaPhi3_123;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi4 = rt.TH1F("MET3_deltaPhi4", "MET3_deltaPhi4;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())


hist_MET1_deltaPhi12 = rt.TH1F("MET1_deltaPhi12", "MET1_deltaPhi12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi13 = rt.TH1F("MET1_deltaPhi13", "MET1_deltaPhi13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi14 = rt.TH1F("MET1_deltaPhi14", "MET1_deltaPhi14;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi23 = rt.TH1F("MET1_deltaPhi23", "MET1_deltaPhi23;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi24 = rt.TH1F("MET1_deltaPhi24", "MET1_deltaPhi24;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET1_deltaPhi34 = rt.TH1F("MET1_deltaPhi34", "MET1_deltaPhi34;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET2_deltaPhi12 = rt.TH1F("MET2_deltaPhi12", "MET2_deltaPhi12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi13 = rt.TH1F("MET2_deltaPhi13", "MET2_deltaPhi13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi14 = rt.TH1F("MET2_deltaPhi14", "MET2_deltaPhi14;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi23 = rt.TH1F("MET2_deltaPhi23", "MET2_deltaPhi23;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi24 = rt.TH1F("MET2_deltaPhi24", "MET2_deltaPhi24;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET2_deltaPhi34 = rt.TH1F("MET2_deltaPhi34", "MET2_deltaPhi34;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET3_deltaPhi12 = rt.TH1F("MET3_deltaPhi12", "MET3_deltaPhi12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi13 = rt.TH1F("MET3_deltaPhi13", "MET3_deltaPhi13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi14 = rt.TH1F("MET3_deltaPhi14", "MET3_deltaPhi14;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi23 = rt.TH1F("MET3_deltaPhi23", "MET3_deltaPhi23;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi24 = rt.TH1F("MET3_deltaPhi24", "MET3_deltaPhi24;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET3_deltaPhi34 = rt.TH1F("MET3_deltaPhi34", "MET3_deltaPhi34;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET4_deltaPhi12 = rt.TH1F("MET4_deltaPhi12", "MET4_deltaPhi12;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi13 = rt.TH1F("MET4_deltaPhi13", "MET4_deltaPhi13;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi14 = rt.TH1F("MET4_deltaPhi14", "MET4_deltaPhi14;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi23 = rt.TH1F("MET4_deltaPhi23", "MET4_deltaPhi23;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi24 = rt.TH1F("MET4_deltaPhi24", "MET4_deltaPhi24;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi34 = rt.TH1F("MET4_deltaPhi34", "MET4_deltaPhi34;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())

hist_MET4_deltaPhi1 = rt.TH1F("MET4_deltaPhi1", "MET4_deltaPhi1;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi2 = rt.TH1F("MET4_deltaPhi2", "MET4_deltaPhi2;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi3 = rt.TH1F("MET4_deltaPhi3", "MET4_deltaPhi3;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())
hist_MET4_deltaPhi4 = rt.TH1F("MET4_deltaPhi4", "MET4_deltaPhi4;\Delta \Phi;a.u.",100,0,rt.TMath.Pi())





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
histList.append(hist_MET1_deltaPhi1)
histList.append(hist_MET1_deltaPhi2)
histList.append(hist_MET1_deltaPhi3)
histList.append(hist_MET1_deltaPhi4)
histList.append(hist_MET2_deltaPhi1)
histList.append(hist_MET2_deltaPhi2)
histList.append(hist_MET2_deltaPhi3)
histList.append(hist_MET2_deltaPhi4)
histList.append(hist_MET3_deltaPhi1)
histList.append(hist_MET3_deltaPhi2)
histList.append(hist_MET3_deltaPhi3)
histList.append(hist_MET3_deltaPhi4)
histList.append(hist_MET4_deltaPhi1)
histList.append(hist_MET4_deltaPhi2)
histList.append(hist_MET4_deltaPhi3)
histList.append(hist_MET4_deltaPhi4)
histList.append(hist_MET3_deltaPhi2_13)
histList.append(hist_MET3_deltaPhi2_123)
histList.append(hist_MET3_deltaPhi1_13)
histList.append(hist_MET3_deltaPhi1_123)
histList.append(hist_MET3_deltaPhi3_13)
histList.append(hist_MET3_deltaPhi3_123)
histList.append(hist_MET1_deltaPhi2_12)
histList.append(hist_MET1_deltaPhi2_123)
histList.append(hist_MET1_deltaPhi1_12)
histList.append(hist_MET1_deltaPhi1_123)
histList.append(hist_MET1_deltaPhi3_12)
histList.append(hist_MET1_deltaPhi3_123)
histList.append(hist_MET2_deltaPhi2_12)
histList.append(hist_MET2_deltaPhi2_123)
histList.append(hist_MET2_deltaPhi1_12)
histList.append(hist_MET2_deltaPhi1_123)
histList.append(hist_MET2_deltaPhi3_12)
histList.append(hist_MET2_deltaPhi3_123)
histList.append(hist_MET1_deltaPhi12)
histList.append(hist_MET1_deltaPhi13)
histList.append(hist_MET1_deltaPhi14)
histList.append(hist_MET1_deltaPhi23)
histList.append(hist_MET1_deltaPhi24)
histList.append(hist_MET1_deltaPhi34)
histList.append(hist_MET2_deltaPhi12)
histList.append(hist_MET2_deltaPhi13)
histList.append(hist_MET2_deltaPhi14)
histList.append(hist_MET2_deltaPhi23)
histList.append(hist_MET2_deltaPhi24)
histList.append(hist_MET2_deltaPhi34)
histList.append(hist_MET3_deltaPhi12)
histList.append(hist_MET3_deltaPhi13)
histList.append(hist_MET3_deltaPhi14)
histList.append(hist_MET3_deltaPhi23)
histList.append(hist_MET3_deltaPhi24)
histList.append(hist_MET3_deltaPhi34)
histList.append(hist_MET4_deltaPhi12)
histList.append(hist_MET4_deltaPhi13)
histList.append(hist_MET4_deltaPhi14)
histList.append(hist_MET4_deltaPhi23)
histList.append(hist_MET4_deltaPhi24)
histList.append(hist_MET4_deltaPhi34)

hist_MT_HVJets.SetLineStyle(3)
hist_MT_bestReco.SetLineStyle(2)


jetComboBestMassDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboHVQuarksDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetMETXDict = {'MET1':0,'MET2':0,'MET3':0,'MET4':0}
jetComboBestMassAndMETXDict = {
'':{'MET1':0,'MET2':0,'MET3':0,'MET4':0},
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
'':{'2':0,'3':0,'4+':0},
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
'':{'2':0,'3':0,'4+':0},
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
	# tag how many jets this event has.
	nJetKey = "Hello"
	if len(tree.JetsAK8) == 2:
		nJetKey = '2'
	elif len(tree.JetsAK8) == 3:
		nJetKey = '3'
	elif len(tree.JetsAK8) > 3:
		nJetKey = '4+'
	else:
		print("Number of jets is not 2, 3, 4, or greater!")
	
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
	if not (len(tree.JetsAK8) == 2): # for tseeing the effects  has on these plots. Make sure to chance the createImagename function above
		continue
	if METXkey == 'MET1':
		hist_MET1_deltaPhi1.Fill(DeltaPhiList[0])
		hist_MET1_deltaPhi2.Fill(DeltaPhiList[1])
		hist_MET1_deltaPhi3.Fill(DeltaPhiList[2])
		hist_MET1_deltaPhi4.Fill(DeltaPhiList[3])
	elif METXkey =='MET2':
		hist_MET2_deltaPhi1.Fill(DeltaPhiList[0])
		hist_MET2_deltaPhi2.Fill(DeltaPhiList[1])
		hist_MET2_deltaPhi3.Fill(DeltaPhiList[2])
		hist_MET2_deltaPhi4.Fill(DeltaPhiList[3])
	elif METXkey == 'MET3':
		hist_MET3_deltaPhi1.Fill(DeltaPhiList[0])
		hist_MET3_deltaPhi2.Fill(DeltaPhiList[1])
		hist_MET3_deltaPhi3.Fill(DeltaPhiList[2])
		hist_MET3_deltaPhi4.Fill(DeltaPhiList[3])
	elif METXkey == 'MET4':
		hist_MET4_deltaPhi1.Fill(DeltaPhiList[0])
		hist_MET4_deltaPhi2.Fill(DeltaPhiList[1])
		hist_MET4_deltaPhi3.Fill(DeltaPhiList[2])
		hist_MET4_deltaPhi4.Fill(DeltaPhiList[3])
	else:
		print("METX is not proper!")

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
	NJetDict[nJetKey] += 1
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
	hist_MT_pseudoInvisibleHVQuarks.Fill(MTofPesudoDecayedHVQuarks)
	MTofHVQuarksProducts = trans_mass_Njet(listOfJetsWithParticlesFromHVQuarks,tree.MET,tree.METPhi)
	hist_MT_HVJets.Fill(MTofHVQuarksProducts)
	#create all of the MT masses (should be 16 of them, plus HVJets means 17 total)
	jetComboMassDict['12'] = trans_mass_Njet([jetCollection[0], jetCollection[1]],tree.MET,tree.METPhi)
	jetComboMassDict[''] = trans_mass_Njet([],tree.MET,tree.METPhi)
	jetComboMassDict['1'] = trans_mass_Njet([jetCollection[0]],tree.MET,tree.METPhi)
	jetComboMassDict['2'] = trans_mass_Njet([jetCollection[1]],tree.MET,tree.METPhi)
	hist_MT12_All.Fill(jetComboMassDict['12'])
	if len(jetCollection) >= 3:
		jetComboMassDict['3'] = trans_mass_Njet([jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['13'] = trans_mass_Njet([jetCollection[0], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['23'] = trans_mass_Njet([jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['123'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		hist_MT13_All.Fill(jetComboMassDict['13'])
		hist_MT23_All.Fill(jetComboMassDict['23'])
		hist_MT123_All.Fill(jetComboMassDict['123'])
	if len(jetCollection) >= 4:
		jetComboMassDict['4'] = trans_mass_Njet([jetCollection[3]],tree.MET,tree.METPhi)
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
	BestMassKey = min(massDiffList,key=massDiffList.get)
	jetComboBestMassDict[BestMassKey] += 1
	jetComboBestMassAndMETXDict[BestMassKey][METXkey] += 1
	jetComboBestMassAndNJetsDict[BestMassKey][nJetKey] += 1
	jetComboHVQuarksAndNJetsDict[jetsHVQuarksKey][nJetKey] += 1
	hist_MT_bestReco.Fill(jetComboMassDict[min(massDiffList,key=massDiffList.get)])
	hist_massDiff_best12.Fill(jetComboMassDict[min(massDiffList,key=massDiffList.get)] - jetComboMassDict['12'])

	if METXkey == 'MET1':
		hist_MET1_deltaPhi12.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()))
		if len(tree.JetsAK8) > 2:
			hist_MET1_deltaPhi13.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[2].Phi()))
			hist_MET1_deltaPhi23.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[2].Phi()))
		if len(tree.JetsAK8) > 3:
			hist_MET1_deltaPhi14.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET1_deltaPhi24.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET1_deltaPhi34.Fill(absDphi(tree.JetsAK8[2].Phi(),tree.JetsAK8[3].Phi()))
	if METXkey == 'MET2':
		hist_MET2_deltaPhi12.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()))
		if len(tree.JetsAK8) > 2:
			hist_MET2_deltaPhi13.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[2].Phi()))
			hist_MET2_deltaPhi23.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[2].Phi()))
		if len(tree.JetsAK8) > 3:
			hist_MET2_deltaPhi14.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET2_deltaPhi24.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET2_deltaPhi34.Fill(absDphi(tree.JetsAK8[2].Phi(),tree.JetsAK8[3].Phi()))
	if METXkey == 'MET3':
		hist_MET3_deltaPhi12.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()))
		if len(tree.JetsAK8) > 2:
			hist_MET3_deltaPhi13.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[2].Phi()))
			hist_MET3_deltaPhi23.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[2].Phi()))
		if len(tree.JetsAK8) > 3:
			hist_MET3_deltaPhi14.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET3_deltaPhi24.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET3_deltaPhi34.Fill(absDphi(tree.JetsAK8[2].Phi(),tree.JetsAK8[3].Phi()))
	if METXkey == 'MET4':
		hist_MET4_deltaPhi12.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[1].Phi()))
		if len(tree.JetsAK8) > 2:
			hist_MET4_deltaPhi13.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[2].Phi()))
			hist_MET4_deltaPhi23.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[2].Phi()))
		if len(tree.JetsAK8) > 3:
			hist_MET4_deltaPhi14.Fill(absDphi(tree.JetsAK8[0].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET4_deltaPhi24.Fill(absDphi(tree.JetsAK8[1].Phi(),tree.JetsAK8[3].Phi()))
			hist_MET4_deltaPhi34.Fill(absDphi(tree.JetsAK8[2].Phi(),tree.JetsAK8[3].Phi()))
	if BestMassKey == '13' and METXkey == 'MET3':
		hist_MET3_deltaPhi1_13.Fill(DeltaPhiList[0])
		hist_MET3_deltaPhi2_13.Fill(DeltaPhiList[1])
		hist_MET3_deltaPhi3_13.Fill(DeltaPhiList[2])
	if BestMassKey == '123' and METXkey == 'MET3':
		hist_MET3_deltaPhi1_123.Fill(DeltaPhiList[0])
		hist_MET3_deltaPhi2_123.Fill(DeltaPhiList[1])
		hist_MET3_deltaPhi3_123.Fill(DeltaPhiList[2])
	if BestMassKey == '12' and METXkey == 'MET1':
		hist_MET1_deltaPhi1_12.Fill(DeltaPhiList[0])
		hist_MET1_deltaPhi2_12.Fill(DeltaPhiList[1])
		hist_MET1_deltaPhi3_12.Fill(DeltaPhiList[2])
	if BestMassKey == '123' and METXkey == 'MET1':
		hist_MET1_deltaPhi1_123.Fill(DeltaPhiList[0])
		hist_MET1_deltaPhi2_123.Fill(DeltaPhiList[1])
		hist_MET1_deltaPhi3_123.Fill(DeltaPhiList[2])
	if BestMassKey == '12' and METXkey == 'MET2':
		hist_MET2_deltaPhi1_12.Fill(DeltaPhiList[0])
		hist_MET2_deltaPhi2_12.Fill(DeltaPhiList[1])
		hist_MET2_deltaPhi3_12.Fill(DeltaPhiList[2])
	if BestMassKey == '123' and METXkey == 'MET2':
		hist_MET2_deltaPhi1_123.Fill(DeltaPhiList[0])
		hist_MET2_deltaPhi2_123.Fill(DeltaPhiList[1])
		hist_MET2_deltaPhi3_123.Fill(DeltaPhiList[2])


#step 4, normalize all histograms

#for histo in histList:
	#area = histo.GetEntries()
	#try:
	#	histo.Scale(1/area)
	#except ZeroDivisionError:
		#print("Empty histogram: " + histo.GetName())
	#	x = 5
	#print(histo.GetName() + " " + str(area))
	#print("StdDev of " + histo.GetName() + " is " + str(histo.GetStdDev()))
print("Preselection " +str(nEventsPassedPreSelection))
print("Total " + str(nEvents))
"""
print("JetComboHVQuarks nEvents")
for (key,value) in jetComboHVQuarksDict.items():
	print(key + " " + str(value))

print("JetCombo nEvents")
for (key,value) in jetComboBestMassDict.items():
	print(key + " " + str(value))

print("METX nEvents")
for (key, value) in jetMETXDict.items():
	print(key + " " + str(value))

print("nJets nEvents")
for (key, value) in NJetDict.items():
	print(key + " " + str(value))

print("JetCombo MET1 MET2 MET3 MET4")
for (key1, value1) in jetComboBestMassAndMETXDict.items():
	print(key1 + " " + str(value1["MET1"]) + " " + str(value1["MET2"]) + " " + str(value1["MET3"]) + " " + str(value1["MET4"]))

print("JetCombo 2jets 3jets 4+jets")
for (key1, value1) in jetComboBestMassAndNJetsDict.items():
	print(key1 + " " + str(value1["2"]) + " " + str(value1["3"]) + " " + str(value1["4+"]))

print("JetComboHVQuarks 2jets 3jets 4+jets")
for (key1, value1) in jetComboHVQuarksAndNJetsDict.items():
	print(key1 + " " + str(value1["2"]) + " " + str(value1["3"]) + " " + str(value1["4+"]))
"""
drawHistos([hist_MT12_All,hist_MT_HVJets,hist_MT_pseudoInvisibleHVQuarks],"MT_HVQuarksJetsPseudoDecayed",True)
drawHistos([hist_MT_HVJets,hist_MT12_All,hist_MT13_All,hist_MT14_All,hist_MT23_All,hist_MT24_All,hist_MT34_All],"MT_DiJets",True)
drawHistos([hist_MT_HVJets,hist_MT123_All,hist_MT124_All,hist_MT134_All,hist_MT234_All],"MT_TriJets",True)
drawHistos([hist_MT_HVJets,hist_MT1234_All],"MT_QuadJets",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT13_All, hist_MT14_All, hist_MT123_All, hist_MT124_All, hist_MT134_All, hist_MT1234_All],"MT_Jet1",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT23_All, hist_MT24_All, hist_MT123_All, hist_MT124_All, hist_MT234_All, hist_MT1234_All],"MT_Jet2",True)
drawHistos([hist_MT_HVJets, hist_MT13_All, hist_MT23_All, hist_MT34_All, hist_MT123_All, hist_MT134_All, hist_MT234_All, hist_MT1234_All],"MT_Jet3",True)
drawHistos([hist_MT_HVJets, hist_MT14_All, hist_MT24_All, hist_MT34_All, hist_MT124_All, hist_MT134_All, hist_MT234_All, hist_MT1234_All],"MT_Jet4",True)
drawHistos([hist_MT_HVJets, hist_MT12_All, hist_MT123_All, hist_MT1234_All,hist_MT_bestReco,hist_MT_pseudoInvisibleHVQuarks],"MT_bestReco",True)
drawHistos([hist_MET1_deltaPhi1,hist_MET1_deltaPhi2,hist_MET1_deltaPhi3,hist_MET1_deltaPhi4],"deltaPhi_MET1",True)
drawHistos([hist_MET2_deltaPhi1,hist_MET2_deltaPhi2,hist_MET2_deltaPhi3,hist_MET2_deltaPhi4],"deltaPhi_MET2",True)
drawHistos([hist_MET3_deltaPhi1,hist_MET3_deltaPhi2,hist_MET3_deltaPhi3,hist_MET3_deltaPhi4],"deltaPhi_MET3",True)
drawHistos([hist_MET3_deltaPhi2,hist_MET3_deltaPhi2_13,hist_MET3_deltaPhi2_123],"deltaPhi2_MET3",True)
drawHistos([hist_MET3_deltaPhi1,hist_MET3_deltaPhi1_13,hist_MET3_deltaPhi1_123],"deltaPhi1_MET3",True)
drawHistos([hist_MET3_deltaPhi3,hist_MET3_deltaPhi3_13,hist_MET3_deltaPhi3_123],"deltaPhi3_MET3",True)

drawHistos([hist_MET1_deltaPhi2,hist_MET1_deltaPhi2_12,hist_MET1_deltaPhi2_123],"deltaPhi2_MET1",True)
drawHistos([hist_MET1_deltaPhi1,hist_MET1_deltaPhi1_12,hist_MET1_deltaPhi1_123],"deltaPhi1_MET1",True)
drawHistos([hist_MET1_deltaPhi3,hist_MET1_deltaPhi3_12,hist_MET1_deltaPhi3_123],"deltaPhi3_MET1",True)

drawHistos([hist_MET2_deltaPhi2,hist_MET2_deltaPhi2_12,hist_MET2_deltaPhi2_123],"deltaPhi2_MET2",True)
drawHistos([hist_MET2_deltaPhi1,hist_MET2_deltaPhi1_12,hist_MET2_deltaPhi1_123],"deltaPhi1_MET2",True)
drawHistos([hist_MET2_deltaPhi3,hist_MET2_deltaPhi3_12,hist_MET2_deltaPhi3_123],"deltaPhi3_MET2",True)

drawHistos([hist_MET4_deltaPhi1,hist_MET4_deltaPhi2,hist_MET4_deltaPhi3,hist_MET4_deltaPhi4],"deltaPhi_MET4",True)

drawHistos([hist_MET1_deltaPhi12,hist_MET1_deltaPhi13,hist_MET1_deltaPhi14,hist_MET1_deltaPhi23,hist_MET1_deltaPhi24,hist_MET1_deltaPhi34],"DeltaPhiJets_MET1",True)
drawHistos([hist_MET2_deltaPhi12,hist_MET2_deltaPhi13,hist_MET2_deltaPhi14,hist_MET2_deltaPhi23,hist_MET2_deltaPhi24,hist_MET2_deltaPhi34],"DeltaPhiJets_MET2",True)
drawHistos([hist_MET3_deltaPhi12,hist_MET3_deltaPhi13,hist_MET3_deltaPhi14,hist_MET3_deltaPhi23,hist_MET3_deltaPhi24,hist_MET3_deltaPhi34],"DeltaPhiJets_MET3",True)
drawHistos([hist_MET4_deltaPhi12,hist_MET4_deltaPhi13,hist_MET4_deltaPhi14,hist_MET4_deltaPhi23,hist_MET4_deltaPhi24,hist_MET4_deltaPhi34],"DeltaPhiJets_MET4",True)

outputFile.Write()
outputFile.Close()


inFile.Close()
