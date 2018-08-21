import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutationsCA11_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

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
	canvas = rt.TCanvas()
	canvas.SetCanvasSize(900,600)
	histStack = rt.THStack()
	histStack.SetTitle(histo[0].GetTitle())
	for hist in histo:
		histStack.Add(hist)
	if logY == True:
		canvas.SetLogy()
	histStack.SetTitle(plotName)
	if not stack:
		histStack.Draw("nostack")
	else:
		histStack.Draw()
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
hist_ref_MT = rt.TH1F("MT_quarks", "HVQuarks;MT;a.u.",100,0,1.5*mZprime)
hist_ref2_MT = rt.TH1F("MT_event", "Event;MT;a.u.",100,0,1.5*mZprime)

hist_MT12_All = rt.TH1F("MT_12_All", "12_All;MT;a.u.",100,0,1.5*mZprime)
hist_MT13_All = rt.TH1F("MT_13_All", "13_All;MT;a.u.",100,0,1.5*mZprime)
hist_MT23_All = rt.TH1F("MT_23_All", "23_All;MT;a.u.",100,0,1.5*mZprime)
hist_MTtri_All = rt.TH1F("MT_tri_All","tri_All;;a.u.",100,0,1.5*mZprime)
hist_MTmax_All = rt.TH1F("MT_max_All","max_All;;a.u.",100,0,1.5*mZprime)
hist_deltaPhi1_All = rt.TH1F("dPhi1_All", "dPhi1_All;dPhi;a.u.",100,0,3.2)
hist_deltaPhi2_All = rt.TH1F("dPhi2_All", "dPhi2_All;dPhi;a.u.",100,0,3.2)
hist_deltaPhi3_All = rt.TH1F("dPhi3_All", "dPhi3_All;dPhi;a.u.",100,0,3.2)

hist_MT12_MET1 = rt.TH1F("MT_12_MET1", "12_MET1;MT;a.u.",100,0,1.5*mZprime)
hist_MT13_MET1 = rt.TH1F("MT_13_MET1", "13_MET1;MT;a.u.",100,0,1.5*mZprime)
hist_MT23_MET1 = rt.TH1F("MT_23_MET1", "23_MET1;MT;a.u.",100,0,1.5*mZprime)
hist_MTtri_MET1 = rt.TH1F("MT_tri_MET1","tri_MET1;MT;a.u.",100,0,1.5*mZprime)
hist_MTmax_MET1 = rt.TH1F("MT_max_MET1","max_MET1;;a.u.",100,0,1.5*mZprime)
hist_deltaPhi1_MET1 = rt.TH1F("dPhi1_MET1", "dPhi1_MET1;dPhi;a.u.",100,0,3.2)
hist_deltaPhi2_MET1 = rt.TH1F("dPhi2_MET1", "dPhi2_MET1;dPhi;a.u.",100,0,3.2)
hist_deltaPhi3_MET1 = rt.TH1F("dPhi3_MET1", "dPhi3_MET1;dPhi;a.u.",100,0,3.2)

hist_MT12_MET2 = rt.TH1F("MT_12_MET2", "12_MET2;MT;a.u.",100,0,1.5*mZprime)
hist_MT13_MET2 = rt.TH1F("MT_13_MET2", "13_MET2;MT;a.u.",100,0,1.5*mZprime)
hist_MT23_MET2 = rt.TH1F("MT_23_MET2", "23_MET2;MT;a.u.",100,0,1.5*mZprime)
hist_MTtri_MET2 = rt.TH1F("MT_tri_MET2","tri_MET2;MT;a.u.",100,0,1.5*mZprime)
hist_MTmax_MET2 = rt.TH1F("MT_max_MET2","max_MET2;;a.u.",100,0,1.5*mZprime)
hist_deltaPhi1_MET2 = rt.TH1F("dPhi1_MET2", "dPhi1_MET2;dPhi;a.u.",100,0,3.2)
hist_deltaPhi2_MET2 = rt.TH1F("dPhi2_MET2", "dPhi2_MET2;dPhi;a.u.",100,0,3.2)
hist_deltaPhi3_MET2 = rt.TH1F("dPhi3_MET2", "dPhi3_MET2;dPhi;a.u.",100,0,3.2)

hist_MT12_MET3 = rt.TH1F("MT_12_MET3", "12_MET3;MT;a.u.",100,0,1.5*mZprime)
hist_MT13_MET3 = rt.TH1F("MT_13_MET3", "13_MET3;MT;a.u.",100,0,1.5*mZprime)
hist_MT23_MET3 = rt.TH1F("MT_23_MET3", "23_MET3;MT;a.u.",100,0,1.5*mZprime)
hist_MTtri_MET3 = rt.TH1F("MT_tri_MET3","tri_MET3;MT;a.u.",100,0,1.5*mZprime)
hist_MTmax_MET3 = rt.TH1F("MT_max_MET3","max_MET3;;a.u.",100,0,1.5*mZprime)
hist_deltaPhi1_MET3 = rt.TH1F("dPhi1_MET3", "dPhi1_MET3;dPhi;a.u.",100,0,3.2)
hist_deltaPhi2_MET3 = rt.TH1F("dPhi2_MET3", "dPhi2_MET3;dPhi;a.u.",100,0,3.2)
hist_deltaPhi3_MET3 = rt.TH1F("dPhi3_MET3", "dPhi3_MET3;dPhi;a.u.",100,0,3.2)

hist_deltaR_MET1 = rt.TH1F("deltaR_MET1","deltaR_MET1;deltaR;a.u.",100,0,6)
hist_deltaR_MET2 = rt.TH1F("deltaR_MET2","deltaR_MET2;deltaR;a.u.",100,0,6)
hist_deltaR_MET3 = rt.TH1F("deltaR_MET3","deltaR_MET3;deltaR;a.u.",100,0,6)


hist_MT_1X = rt.TH1F("MT_1X","1X;MT;a.u.",100,0,1.5*mZprime)
hist_MT_HVquarks = rt.TH1F("MT_HVQuarksJets","HVQuarksJets;MT;a.u.",100,0,1.5*mZprime)

hist_diff_plotMax = rt.TH1F("diff_plotMax","Max-Mid;\Delta MT;a.u.", 100,0, mZprime)
hist_diff_plotMin = rt.TH1F("diff_plotMin","Mid-Min;\Delta MT;a.u.", 100,0, mZprime)
hist_diff_plotMid = rt.TH1F("diff_plotMid","abs(Max-Min)/2;\Delta MT;a.u.", 100,0, mZprime)



histList = []
histList.append(hist_ref_MT)
histList.append(hist_ref2_MT)
histList.append(hist_MT12_All)
histList.append(hist_MT13_All)
histList.append(hist_MT23_All)
histList.append(hist_MTtri_All)
histList.append(hist_MT12_MET1)
histList.append(hist_MT13_MET1)
histList.append(hist_MT23_MET1)
histList.append(hist_MTtri_MET1)
histList.append(hist_MT12_MET2)
histList.append(hist_MT13_MET2)
histList.append(hist_MT23_MET2)
histList.append(hist_MTtri_MET2)
histList.append(hist_MT12_MET3)
histList.append(hist_MT13_MET3)
histList.append(hist_MT23_MET3)
histList.append(hist_MTtri_MET3)
histList.append(hist_MT_1X)
histList.append(hist_MT_HVquarks)
histList.append(hist_MTmax_All)
histList.append(hist_MTmax_MET1)
histList.append(hist_MTmax_MET2)
histList.append(hist_MTmax_MET3)
histList.append(hist_diff_plotMax)
histList.append(hist_diff_plotMin)
histList.append(hist_diff_plotMid)
histList.append(hist_deltaPhi1_All)
histList.append(hist_deltaPhi2_All)
histList.append(hist_deltaPhi3_All)
histList.append(hist_deltaPhi1_MET1)
histList.append(hist_deltaPhi2_MET1)
histList.append(hist_deltaPhi3_MET1)
histList.append(hist_deltaPhi1_MET2)
histList.append(hist_deltaPhi2_MET2)
histList.append(hist_deltaPhi3_MET2)
histList.append(hist_deltaPhi1_MET3)
histList.append(hist_deltaPhi2_MET3)
histList.append(hist_deltaPhi3_MET3)
histList.append(hist_deltaR_MET1)
histList.append(hist_deltaR_MET2)
histList.append(hist_deltaR_MET3)

hist_ref_MT.SetLineColor(rt.kMagenta)
hist_ref_MT.SetLineStyle(2)
hist_ref2_MT.SetLineColor(rt.kBlack)
hist_ref2_MT.SetLineStyle(2)
hist_MT_1X.SetLineColor(rt.kRed)
hist_MT_HVquarks.SetLineColor(rt.kBlue)


for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	#jetCollection = tree.GenJetsAK8
	jetCollection = tree.JetsCA11

	
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
	#additional selection, we need 3+ jets in our jet collection
	if not (len(tree.JetsAK8)>=3):
		continue
	DeltaPhi3_CA11 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
	DeltaPhiMin_3jet = min(DeltaPhi3_CA11, tree.DeltaPhi1_CA11, tree.DeltaPhi2_CA11)
	#find the HV quarks
	nHVquark = 0
	nbarHVquark = 0
	for iPart in range(len(tree.GenParticles)):
		if tree.GenParticles_PdgId[iPart] == 4900101:
			hvQuark = tree.GenParticles[iPart]
			nHVquark += 1
		elif tree.GenParticles_PdgId[iPart] == -4900101:
			barhvQuark = tree.GenParticles[iPart]
			nbarHVquark += 1
	if (nHVquark != 1) or (nbarHVquark != 1): # skip events that don't have only 2 unstable HV quarks
		continue
	


	deltaRjetHVquark = [hvQuark.DeltaR(jetCollection[0]),hvQuark.DeltaR(jetCollection[1]),hvQuark.DeltaR(jetCollection[2])]
	deltaRjetbarHVquark = [barhvQuark.DeltaR(jetCollection[0]),barhvQuark.DeltaR(jetCollection[1]),barhvQuark.DeltaR(jetCollection[2])]
	
	deltaRjethvQuarkMinIndex = deltaRjetHVquark.index(min(deltaRjetHVquark))
	deltaRjetbarhvQuarkMinIndex = deltaRjetbarHVquark.index(min(deltaRjetbarHVquark))

	hist_ref_MT.Fill((hvQuark+barhvQuark).M())
	hist_ref2_MT.Fill(tree.MT_CA11)
	if not (deltaRjethvQuarkMinIndex == deltaRjetbarhvQuarkMinIndex):
		MT_HVQuarks_Jets = trans_mass_Njet([jetCollection[deltaRjethvQuarkMinIndex], jetCollection[deltaRjetbarhvQuarkMinIndex]],tree.MET,tree.METPhi)
		hist_MT_HVquarks.Fill(MT_HVQuarks_Jets)

	jet12Mt = trans_mass_Njet([jetCollection[0], jetCollection[1]],tree.MET,tree.METPhi)
	jet13Mt = trans_mass_Njet([jetCollection[0], jetCollection[2]],tree.MET,tree.METPhi)
	jet23Mt = trans_mass_Njet([jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
	trijetMt = trans_mass_Njet([jetCollection[0], jetCollection[1],jetCollection[2]],tree.MET,tree.METPhi)
	mtList = [jet12Mt, jet13Mt, jet23Mt]
	jetMtmax = max(mtList)
	jetMtmin = min(mtList)
	mtList.remove(jetMtmax)
	mtList.remove(jetMtmin)
	jetMtmid = mtList[0]

	hist_MT12_All.Fill(jet12Mt)
	hist_MT13_All.Fill(jet13Mt)
	hist_MT23_All.Fill(jet23Mt)
	hist_MTmax_All.Fill(jetMtmax)
	hist_MTtri_All.Fill(trijetMt)
	hist_deltaPhi1_All.Fill(tree.DeltaPhi1_CA11)
	hist_deltaPhi2_All.Fill(tree.DeltaPhi2_CA11)
	hist_deltaPhi3_All.Fill(DeltaPhi3_CA11)
	if (tree.DeltaPhi1_CA11 == DeltaPhiMin_3jet):
		hist_MT12_MET1.Fill(jet12Mt)
		hist_MT13_MET1.Fill(jet13Mt)
		hist_MT23_MET1.Fill(jet23Mt)
		hist_MTmax_MET1.Fill(jetMtmax)
		hist_MTtri_MET1.Fill(trijetMt)
		hist_deltaPhi1_MET1.Fill(tree.DeltaPhi1_CA11)
		hist_deltaPhi2_MET1.Fill(tree.DeltaPhi2_CA11)
		hist_deltaPhi3_MET1.Fill(DeltaPhi3_CA11)
		hist_deltaR_MET1.Fill(jetCollection[1].DeltaR(jetCollection[2]))
	elif (tree.DeltaPhi2_CA11 == DeltaPhiMin_3jet):
		hist_MT12_MET2.Fill(jet12Mt)
		hist_MT13_MET2.Fill(jet13Mt)
		hist_MT23_MET2.Fill(jet23Mt)
		hist_MTmax_MET2.Fill(jetMtmax)
		hist_MTtri_MET2.Fill(trijetMt)
		hist_deltaPhi1_MET2.Fill(tree.DeltaPhi1_CA11)
		hist_deltaPhi2_MET2.Fill(tree.DeltaPhi2_CA11)
		hist_deltaPhi3_MET2.Fill(DeltaPhi3_CA11)
		hist_MT_1X.Fill(jet12Mt)
		hist_deltaR_MET2.Fill(jetCollection[0].DeltaR(jetCollection[2]))
	elif (DeltaPhi3_CA11 == DeltaPhiMin_3jet):
		hist_diff_plotMax.Fill(jetMtmax-jetMtmid)
		hist_diff_plotMin.Fill(jetMtmid-jetMtmin)
		hist_diff_plotMid.Fill((jetMtmax-jetMtmin)/2)
		hist_MT12_MET3.Fill(jet12Mt)
		hist_MT13_MET3.Fill(jet13Mt)
		hist_MT23_MET3.Fill(jet23Mt)
		hist_MTmax_MET3.Fill(jetMtmax)
		hist_MTtri_MET3.Fill(trijetMt)
		hist_deltaPhi1_MET3.Fill(tree.DeltaPhi1_CA11)
		hist_deltaPhi2_MET3.Fill(tree.DeltaPhi2_CA11)
		hist_deltaPhi3_MET3.Fill(DeltaPhi3_CA11)
		hist_MT_1X.Fill(jet13Mt)
		hist_deltaR_MET3.Fill(jetCollection[0].DeltaR(jetCollection[1]))


for histo in histList:
	if "All" in histo.GetName():
		histo.SetLineColor(rt.kBlack)
	elif "MET1" in histo.GetName():
		histo.SetLineColor(rt.kBlue)
	elif "MET2" in histo.GetName():
		histo.SetLineColor(rt.kGreen)
	elif "MET3" in histo.GetName():
		histo.SetLineColor(rt.kRed)
drawHistos([hist_deltaPhi1_MET1,hist_deltaPhi1_MET2,hist_deltaPhi1_MET3],"dPhi1Stack",False,True)
drawHistos([hist_deltaPhi2_MET1,hist_deltaPhi2_MET2,hist_deltaPhi2_MET3],"dPhi2Stack",False,True)
drawHistos([hist_deltaPhi3_MET1,hist_deltaPhi3_MET2,hist_deltaPhi3_MET3],"dPhi3Stack",False,True)
drawHistos([hist_deltaPhi1_MET1,hist_deltaPhi1_MET2,hist_deltaPhi1_MET3],"dPhi1StackLogy",True,True)
drawHistos([hist_deltaPhi2_MET1,hist_deltaPhi2_MET2,hist_deltaPhi2_MET3],"dPhi2StackLogy",True,True)
drawHistos([hist_deltaPhi3_MET1,hist_deltaPhi3_MET2,hist_deltaPhi3_MET3],"dPhi3StackLogy",True,True)


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
print("Events that have 3+ jets and preSeelction: " + str(hist_ref_MT.GetEntries()))

for hist in histList:
	print("StdDev of " + hist.GetName() + " is " + str(hist.GetStdDev()))



# Draw same jets Construction Here
drawHistos([hist_ref_MT, hist_MT_1X, hist_MT_HVquarks,hist_ref2_MT], "MT_1X_HVQuarks", True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT12_All,hist_MT12_MET1,hist_MT12_MET2,hist_MT12_MET3],"MT_12",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT13_All,hist_MT13_MET1,hist_MT13_MET2,hist_MT13_MET3],"MT_13",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT23_All,hist_MT23_MET1,hist_MT23_MET2,hist_MT23_MET3],"MT_23",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MTtri_All,hist_MTtri_MET1,hist_MTtri_MET2,hist_MTtri_MET3],"MT_tri",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MTmax_All, hist_MTmax_MET1, hist_MTmax_MET2, hist_MTmax_MET3], "MT_max",True)
drawHistos([hist_deltaPhi1_All,hist_deltaPhi1_MET1,hist_deltaPhi1_MET2,hist_deltaPhi1_MET3],"dPhi1",True)
drawHistos([hist_deltaPhi2_All,hist_deltaPhi2_MET1,hist_deltaPhi2_MET2,hist_deltaPhi2_MET3],"dPhi2",True)
drawHistos([hist_deltaPhi3_All,hist_deltaPhi3_MET1,hist_deltaPhi3_MET2,hist_deltaPhi3_MET3],"dPhi3",True)

hist_diff_plotMax.SetLineColor(rt.kRed)
hist_diff_plotMid.SetLineColor(rt.kMagenta)
hist_diff_plotMin.SetLineColor(rt.kBlue)
drawHistos([hist_diff_plotMax,hist_diff_plotMin,hist_diff_plotMid],"difference_plot",True)
drawHistos([hist_deltaR_MET1,hist_deltaR_MET2,hist_deltaR_MET3],"deltaR")
for histo in histList:
	name = histo.GetName()
	if "12" in name or "Phi1" in name:
		histo.SetLineColor(rt.kBlue)
	elif "13" in name or "Phi2" in name:
		histo.SetLineColor(rt.kGreen)
	elif "23" in name or "Phi3" in name:
		histo.SetLineColor(rt.kRed)
	elif "tri" in name:
		histo.SetLineColor(rt.kBlack)
	elif 'max' in name:
		histo.SetLineColor(rt.kCyan)

#draw same METX plots here
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT12_All,hist_MT13_All,hist_MT23_All,hist_MTtri_All],"MT_All",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT12_MET1,hist_MT13_MET1,hist_MT23_MET1,hist_MTtri_MET1],"MT_MET1",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT12_MET2,hist_MT13_MET2,hist_MT23_MET2,hist_MTtri_MET2],"MT_MET2",True)
drawHistos([hist_ref_MT,hist_ref2_MT, hist_MT12_MET3,hist_MT13_MET3,hist_MT23_MET3,hist_MTtri_MET3],"MT_MET3",True)
drawHistos([hist_deltaPhi1_All,hist_deltaPhi2_All,hist_deltaPhi3_All],"deltaPhi_All",True)
drawHistos([hist_deltaPhi1_MET1,hist_deltaPhi2_MET1,hist_deltaPhi3_MET1],"deltaPhi_MET1",True)
drawHistos([hist_deltaPhi1_MET2,hist_deltaPhi2_MET2,hist_deltaPhi3_MET2],"deltaPhi_MET2",True)
drawHistos([hist_deltaPhi1_MET3,hist_deltaPhi2_MET3,hist_deltaPhi3_MET3],"deltaPhi_MET3",True)



inFile.Close()
