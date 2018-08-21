import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(0)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetsAK8_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

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

def drawHistos(histo, plotName, logY = False):
	canvas = rt.TCanvas()
	canvas.SetCanvasSize(900,600)
	histStack = rt.THStack()
	for hist in histo:
		histStack.Add(hist)
	if logY == True:
		canvas.SetLogy()
	histStack.Draw("nostack")
	canvas.BuildLegend(0.75,0.75,0.9,0.9,"")
	canvas.SetTitle(plotName)
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
alpha = '0p5'

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
hist_goodTag_All = rt.TH1F("gt_All",";MT;a.u.",100,0,1.3*mZprime)
hist_goodTag_MET1 = rt.TH1F("gt_MET1",";MT;a.u.",100,0,1.3*mZprime)
hist_goodTag_MET2 = rt.TH1F("gt_MET2",";MT;a.u.",100,0,1.3*mZprime)
hist_goodTag_MET3 = rt.TH1F("gt_MET3",";MT;a.u.",100,0,1.3*mZprime)
hist_goodAntiTag_All = rt.TH1F("gat_All",";MT;a.u.",100,0,1.3*mZprime)
hist_goodAntiTag_MET1 = rt.TH1F("gat_MET1",";MT;a.u.",100,0,1.3*mZprime)
hist_goodAntiTag_MET2 = rt.TH1F("gat_MET2",";MT;a.u.",100,0,1.3*mZprime)
hist_goodAntiTag_MET3 = rt.TH1F("gat_MET3",";MT;a.u.",100,0,1.3*mZprime)
hist_badTag_All = rt.TH1F("bt_All",";MT;a.u.",100,0,1.3*mZprime)
hist_badTag_MET1 = rt.TH1F("bt_MET1",";MT;a.u.",100,0,1.3*mZprime)
hist_badTag_MET2 = rt.TH1F("bt_MET2",";MT;a.u.",100,0,1.3*mZprime)
hist_badTag_MET3 = rt.TH1F("bt_MET3",";MT;a.u.",100,0,1.3*mZprime)
hist_ref_MT = rt.TH1F("ref_All", "Reference;MT;a.u.",100,0,1.3*mZprime)

hist_triJetMT_All = rt.TH1F("triJetMT_All","All;;a.u.",100,0,1.5*mZprime)
hist_triJetMT_MET1 = rt.TH1F("triJetMT_MET1","MET1;MT;a.u.",100,0,1.5*mZprime)
hist_triJetMT_MET2 = rt.TH1F("triJetMT_MET2","MET2;MT;a.u.",100,0,1.5*mZprime)
hist_triJetMT_MET3 = rt.TH1F("triJetMT_MET3","MET3;MT;a.u.",100,0,1.5*mZprime)

histList = []
histList.append(hist_goodTag_All)
histList.append(hist_goodTag_MET1)
histList.append(hist_goodTag_MET2)
histList.append(hist_goodTag_MET3)
histList.append(hist_goodAntiTag_All)
histList.append(hist_goodAntiTag_MET1)
histList.append(hist_goodAntiTag_MET2)
histList.append(hist_goodAntiTag_MET3)
histList.append(hist_badTag_All)
histList.append(hist_badTag_MET1)
histList.append(hist_badTag_MET2)
histList.append(hist_badTag_MET3)
histList.append(hist_ref_MT)
histList.append(hist_triJetMT_All)
histList.append(hist_triJetMT_MET1)
histList.append(hist_triJetMT_MET2)
histList.append(hist_triJetMT_MET3)

for histo in histList:
	if "All" in histo.GetName():
		histo.SetLineColor(rt.kBlack)
	elif "MET1" in histo.GetName():
		histo.SetLineColor(rt.kBlue)
	elif "MET2" in histo.GetName():
		histo.SetLineColor(rt.kGreen)
	elif "MET3" in histo.GetName():
		histo.SetLineColor(rt.kRed)
	if "gt" in histo.GetName():
		histo.SetLineStyle(1)
	elif "gat" in histo.GetName():
		histo.SetLineStyle(2)
	elif "bt" in histo.GetName():
		histo.SetLineStyle(3)
	elif "tri" in histo.GetName():
		histo.SetLineStyle(4)

hist_ref_MT.SetLineColor(rt.kMagenta)


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
	#additional selection, we need 3+ jets
	if not (len(tree.JetsAK8)>=3):
		continue
	DeltaPhi3_AK8 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
	DeltaPhiMin_3jet = min(DeltaPhi3_AK8, tree.DeltaPhi1_AK8, tree.DeltaPhi2_AK8)
	#figure out which jet alignes with which HVquark
	nHVquark = 0
	nbarHVquark = 0
	for iPart in range(len(tree.GenParticles)):
		if tree.GenParticles_PdgId[iPart] == 4900101:
			hvQuark = tree.GenParticles[iPart]
			nHVquark += 1
		elif tree.GenParticles_PdgId[iPart] == -4900101:
			barhvQuark = tree.GenParticles[iPart]
			nbarHVquark += 1
		elif tree.GenParticles_PdgId[iPart] == 4900023:
			ZpIdx = iPart

	if (nHVquark != 1) or (nbarHVquark != 1):
	#	print("quark problem, EventNumber: "+str(iEvent))
	#	print("Zprime daughters' PDGID's: ")
	#	for iPart in range(len(tree.GenParticles)):
	#		if tree.GenParticles_ParentIdx[iPart] == ZpIdx:
	#			print(tree.GenParticles_PdgId[iPart])
		continue

	

	listOfZpDaughters = []
	for iPart in range(len(tree.GenParticles)):
		if tree.GenParticles_ParentIdx[iPart] == ZpIdx:
			listOfZpDaughters.append(tree.GenParticles_PdgId[iPart])
	if len(listOfZpDaughters) != 2:
		print("Zprime has " + str(len(listOfZpDaughters)) + " daughters. Their PDG IDs are:")
		for x in listOfZpDaughters:
			print(str(x))

	
	deltaRjetHVquark = [hvQuark.DeltaR(jetCollection[0]),hvQuark.DeltaR(jetCollection[1]),hvQuark.DeltaR(jetCollection[2])]
	deltaRjetbarHVquark = [barhvQuark.DeltaR(jetCollection[0]),barhvQuark.DeltaR(jetCollection[1]),barhvQuark.DeltaR(jetCollection[2])]
	
	deltaRjethvQuarkMinIndex = deltaRjetHVquark.index(min(deltaRjetHVquark))
	deltaRjetbarhvQuarkMinIndex = deltaRjetbarHVquark.index(min(deltaRjetbarHVquark))

	if deltaRjethvQuarkMinIndex == deltaRjetbarhvQuarkMinIndex:
		QuarksWithSameJet += 1
		continue


	if not isNear(trans_mass_Njet(jetCollection[0:2], tree.MET,tree.METPhi),mt_calc(jetCollection[0],jetCollection[1], tree.MET,tree.METPhi),0.001):
		print("new MT calc doesn't match with old one!")
	hist_triJetMT_All.Fill(trans_mass_Njet(jetCollection[0:3], tree.MET,tree.METPhi))
	if isNear(DeltaPhiMin_3jet, tree.DeltaPhi1_AK8, 0.01):
		jet1Aligned += 1
		alignedJetIndex = 0
		hist_triJetMT_MET1.Fill(trans_mass_Njet(jetCollection[0:3], tree.MET,tree.METPhi))
	elif isNear(DeltaPhiMin_3jet, tree.DeltaPhi2_AK8, 0.01):
		jet2Aligned += 1
		alignedJetIndex = 1
		hist_triJetMT_MET2.Fill(trans_mass_Njet(jetCollection[0:3], tree.MET,tree.METPhi))
	elif isNear(DeltaPhiMin_3jet, DeltaPhi3_AK8, 0.01):
		jet3Aligned += 1
		alignedJetIndex = 2
		hist_triJetMT_MET3.Fill(trans_mass_Njet(jetCollection[0:3], tree.MET,tree.METPhi))
	else:
		print("DeltaPhiMin doesn't match with any of the jets!")
		alignedJetIndex = False
	# Jet tagging algorithm
	# first pass: take MET jet as primary jet, and highest pt of remaining jets as secondary
	
	if (type(alignedJetIndex) == int):
		iPrimJet = alignedJetIndex
		if iPrimJet == 0:
			iSecoJet = 1
			iExtrJet = 2
		elif iPrimJet == 1:
			iSecoJet = 0
			iExtrJet = 2
		elif iPrimJet == 2:
			iSecoJet = 0
			iExtrJet = 1

		#define tag types
		# Good Tag is when Prim and Seco are aligned with HVquarks
		# Good Anti-Tag is when Prim and Extr are aligned with HVquarks
		# bad tag is when Seco and Extr are aligned with HVquarks
		
		goodTag = ((iPrimJet == deltaRjethvQuarkMinIndex) or (iPrimJet == deltaRjetbarhvQuarkMinIndex)) and ((iSecoJet == deltaRjethvQuarkMinIndex) or (iSecoJet == deltaRjetbarhvQuarkMinIndex))
		goodAntiTag = ((iPrimJet == deltaRjethvQuarkMinIndex) or (iPrimJet == deltaRjetbarhvQuarkMinIndex)) and ((iExtrJet == deltaRjethvQuarkMinIndex) or (iExtrJet == deltaRjetbarhvQuarkMinIndex))
		badTag = ((iSecoJet == deltaRjethvQuarkMinIndex) or (iSecoJet == deltaRjetbarhvQuarkMinIndex)) and ((iExtrJet == deltaRjethvQuarkMinIndex) or (iExtrJet == deltaRjetbarhvQuarkMinIndex))
		hist_ref_MT.Fill((hvQuark+barhvQuark).M())
		#hist_ref_MT.Fill(tree.GenParticles[ZpIdx].Mt())
		if (goodTag):
			jetTagIsGood += 1
			hist_goodTag_All.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			if alignedJetIndex == 0:
				goodTagAndMET1 += 1
				hist_goodTag_MET1.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 1:
				goodTagAndMET2 += 1
				hist_goodTag_MET2.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 2:
				goodTagAndMET3 += 1
				hist_goodTag_MET3.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
		if (goodAntiTag): #anti-tag
			jetAntiTagIsGood += 1
			hist_goodAntiTag_All.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iExtrJet],tree.MET,tree.METPhi))
			if alignedJetIndex == 0:
				goodAntiTagAndMET1 += 1
				hist_goodAntiTag_MET1.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iExtrJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 1:
				goodAntiTagAndMET2 += 1
				hist_goodAntiTag_MET2.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iExtrJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 2:
				goodAntiTagAndMET3 += 1
				hist_goodAntiTag_MET3.Fill(mt_calc(jetCollection[iPrimJet],jetCollection[iExtrJet],tree.MET,tree.METPhi))
		if (badTag):
			badTagEvents+=1
			hist_badTag_All.Fill(mt_calc(jetCollection[iExtrJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			if alignedJetIndex == 0:
				badTagAndMET1 += 1
				hist_badTag_MET1.Fill(mt_calc(jetCollection[iExtrJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 1:
				badTagAndMET2 += 1
				hist_badTag_MET2.Fill(mt_calc(jetCollection[iExtrJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))
			elif alignedJetIndex == 2:
				badTagAndMET3 += 1
				hist_badTag_MET3.Fill(mt_calc(jetCollection[iExtrJet],jetCollection[iSecoJet],tree.MET,tree.METPhi))



#step 4, normalize all histograms

for histo in histList:
	area = histo.GetEntries()
	try:
		histo.Scale(1/area)
	except ZeroDivisionError:
		print("Empty histogram: " + histo.GetName())

#4a: print events in each region
print("Total Events: " + str(nEvents))
print("Events that passed preSeelction: "+ str(nEventsPassedPreSelection))
print("Total Events with atleast 3jets: " + str(jet1Aligned+jet2Aligned+jet3Aligned))
print("MET aligned with jet 1: " + str(jet1Aligned))
print("MET aligned with jet 2: " + str(jet2Aligned))
print("MET aligned with jet 3: " + str(jet3Aligned))
print("GoodTag Events: " + str(jetTagIsGood))
print("GoodTag with MET1: "+ str(goodTagAndMET1))
print("GoodTag with MET2: "+ str(goodTagAndMET2))
print("GoodTag with MET3: "+ str(goodTagAndMET3))
print("GoodAntiTag Events: " + str(jetAntiTagIsGood))
print("GoodAntiTag with MET1: "+ str(goodAntiTagAndMET1))
print("GoodAntiTag with MET2: "+ str(goodAntiTagAndMET2))
print("GoodAntiTag with MET3: "+ str(goodAntiTagAndMET3))
print("badTag Events: " + str(badTagEvents))
print("badTag with MET1: "+ str(badTagAndMET1))
print("badTag with MET2: "+ str(badTagAndMET2))
print("badTag with MET3: "+ str(badTagAndMET3))
print("Quarks Aligned with same Jet: " + str(QuarksWithSameJet))

for hist in histList:
	print("StdDev of " + hist.GetName() + " is " + str(hist.GetStdDev()))

drawHistos([hist_ref_MT,hist_goodTag_All,hist_goodTag_MET1,hist_goodTag_MET2,hist_goodTag_MET3],"goodTags_MT.png",True)
drawHistos([hist_ref_MT,hist_goodAntiTag_All,hist_goodAntiTag_MET1,hist_goodAntiTag_MET2,hist_goodAntiTag_MET3],"goodAntiTags_MT.png",True)
drawHistos([hist_ref_MT,hist_badTag_All,hist_badTag_MET1,hist_badTag_MET2,hist_badTag_MET3],"badTags_MT.png",True)

drawHistos([hist_ref_MT,hist_goodTag_All,hist_goodAntiTag_All,hist_badTag_All],"All_Events_MT.png",True)
drawHistos([hist_ref_MT,hist_goodTag_MET1,hist_goodAntiTag_MET1,hist_badTag_MET1,hist_triJetMT_MET1],"MET1_Events_MT.png",True)
drawHistos([hist_ref_MT,hist_goodTag_MET2,hist_goodAntiTag_MET2,hist_badTag_MET2,hist_triJetMT_MET2],"MET2_Events_MT.png",True)
drawHistos([hist_ref_MT,hist_goodTag_MET3,hist_goodAntiTag_MET3,hist_badTag_MET3,hist_triJetMT_MET3],"MET3_Events_MT.png",True)
drawHistos([hist_triJetMT_All,hist_triJetMT_MET1,hist_triJetMT_MET2,hist_triJetMT_MET3,hist_ref_MT],"TriJetMT",True)

inFile.Close()
