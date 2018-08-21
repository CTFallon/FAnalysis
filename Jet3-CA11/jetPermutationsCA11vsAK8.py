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
nEventsPassedPreSelection = 0

# initilize histograms
hist_ref_MT = rt.TH1F("MT_quarks", "HVQuarks;MT;a.u.",100,0,1.5*mZprime)
#hist_ref2_MT = rt.TH1F("MT_event", "Event;MT;a.u.",100,0,1.5*mZprime)

hist_MT_AK8_12    = rt.TH1F("MT_AK8_12"   , "AK8_12;MT;a.u."   ,100,0,2*mZprime)
hist_MT_AK8_123   = rt.TH1F("MT_AK8_123"  , "AK8_123;MT;a.u."  ,100,0,2*mZprime)
hist_MT_AK8_1234  = rt.TH1F("MT_AK8_1234" , "AK8_1234;MT;a.u." ,100,0,2*mZprime)
hist_MT_CA11_12   = rt.TH1F("MT_CA11_12"  , "CA11_12;MT;a.u."  ,100,0,2*mZprime)
hist_MT_CA11_123  = rt.TH1F("MT_CA11_123" , "CA11_123;MT;a.u." ,100,0,2*mZprime)
hist_MT_CA11_1234 = rt.TH1F("MT_CA11_1234", "CA11_1234;MT;a.u.",100,0,2*mZprime)



histList = []
histList.append(hist_ref_MT)
#histList.append(hist_ref2_MT)
histList.append(hist_MT_AK8_12)
histList.append(hist_MT_AK8_123)
histList.append(hist_MT_AK8_1234)
histList.append(hist_MT_CA11_12)
histList.append(hist_MT_CA11_123)
histList.append(hist_MT_CA11_1234)

hist_ref_MT.SetLineColor(rt.kMagenta)
hist_ref_MT.SetLineStyle(2)
#hist_ref2_MT.SetLineColor(rt.kBlack)
#hist_ref2_MT.SetLineStyle(2)

hist_MT_AK8_12.SetLineColor(2)
hist_MT_AK8_123.SetLineColor(3)
hist_MT_AK8_1234.SetLineColor(4)
hist_MT_CA11_12.SetLineColor(5)
hist_MT_CA11_123.SetLineColor(6)
hist_MT_CA11_1234.SetLineColor(7)


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
	#additional selection, we need 3+ jets in our jet collection
	#if not (len(tree.JetsAK8)>=3):
	#	continue
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
	
	hist_ref_MT.Fill((hvQuark+barhvQuark).M())
	#hist_ref2_MT.Fill(tree.MT_AK8)
	
	hist_MT_AK8_12.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1]], tree.MET, tree.METPhi))
	if len(tree.JetsAK8) >= 3:
		hist_MT_AK8_123.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]], tree.MET, tree.METPhi))
	if len(tree.JetsAK8) >= 4:
		hist_MT_AK8_1234.Fill(trans_mass_Njet([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2],tree.JetsAK8[3]], tree.MET, tree.METPhi))
	if len(tree.JetsCA11) >= 2:
		hist_MT_CA11_12.Fill(trans_mass_Njet([tree.JetsCA11[0],tree.JetsCA11[1]], tree.MET, tree.METPhi))
	if len(tree.JetsCA11) >= 3:
		hist_MT_CA11_123.Fill(trans_mass_Njet([tree.JetsCA11[0],tree.JetsCA11[1],tree.JetsCA11[2]], tree.MET, tree.METPhi))
	if len(tree.JetsCA11) >= 4:
		hist_MT_CA11_1234.Fill(trans_mass_Njet([tree.JetsCA11[0],tree.JetsCA11[1],tree.JetsCA11[2],tree.JetsCA11[3]], tree.MET, tree.METPhi))

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

for hist in histList:
	print("StdDev of " + hist.GetName() + " is " + str(hist.GetStdDev()))

drawHistos([hist_ref_MT,hist_MT_AK8_12,hist_MT_AK8_123,hist_MT_AK8_1234,hist_MT_CA11_12,hist_MT_CA11_123,hist_MT_CA11_1234],"AllJets",True)
drawHistos([hist_ref_MT,hist_MT_AK8_12,hist_MT_AK8_123,hist_MT_AK8_1234],"AK8Jets",True)
drawHistos([hist_ref_MT,hist_MT_CA11_12,hist_MT_CA11_123,hist_MT_CA11_1234],"CA11Jets",True)
drawHistos([hist_ref_MT,hist_MT_AK8_12,hist_MT_AK8_123,hist_MT_CA11_12],"triJetAK8vsDiJetCA11",True)



inFile.Close()
