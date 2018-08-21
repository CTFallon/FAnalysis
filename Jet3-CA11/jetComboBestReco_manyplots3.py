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

def getSigBkgEffiency(histSig, histListBkg,textName):
	print("Saveing effieciecies for " + textName + " to file")
	outText = open(textName+".txt","w")
	nBins = histSig.GetNbinsX()
	#outText.write("Number of Bins is " + str(nBins) + "\n")
	outText.write("SigEff BkgEff"+ "\n")
	sigIntegral = histSig.GetEntries()
	bkgIntegral = 0
	for bkgHist in histListBkg:
		bkgIntegral += bkgHist.GetEntries()
	for startBin in range(nBins):
		sigVal = histSig.Integral(startBin,nBins)
		bkgVal = 0
		for bkgHist in histListBkg:
			bkgVal += bkgHist.Integral(startBin,nBins)
		#outText.write(str(float(sigVal)/float(sigIntegral)) + " " + str(float(bkgVal)/float(bkgIntegral))+ "\n")
		outText.write(str(sigVal) + " " + str(bkgVal)+ "\n")
	outText.close()
	print("Done")


def pseudoAK(jetList, R):
	endJetList = []
	while len(jetList) > 1:
		nJets = len(jetList)
		dMat = [[0 for i in range(nJets)] for j in range(nJets)]
		dMatMin, iMin, jMin = 100, 10, 10
		for i in range(nJets):
			for j in range(nJets):
				if i == j:
					dMat[i][j] = 1/jetList[i].Pt()
				else:
					dMat[i][j] = min(1/jetList[i],1/jetList[j])*jetList[i].DeltaR(jetList[j])/R
				if dMatMin > dMat[i][j]:
					dMatMin, iMin, jMin = dMat[i][j], i, j
		if iMin != jMin:
			newJet = jetList[i] + jetList[j]
			jetList.pop(iMin)
			jetList.pop(jMin)
			jetList.append(newJet)
		if iMin == jMin:
			endJetList.append(jetList.pop(iMin))
	else:
		endJetList.append(jetList.pop())
	return endJetList
		
			

def FOM(nSig, nBkg):
	try:
		return float(nSig)/rt.TMath.Sqrt(float(nBkg))
	except ZeroDivisionError:
		return 0.0

def makeFOMPlot(histSig, histListBkg, Name):
	histFOM = rt.TH1F(Name+"sig", histSig.GetTitle(), histSig.GetNbinsX(), histSig.GetXaxis().GetBinLowEdge(1), histSig.GetXaxis().GetBinUpEdge(histSig.GetNbinsX()))
	histBkg = rt.TH1F(Name+"bkg", histSig.GetTitle(), histSig.GetNbinsX(), histSig.GetXaxis().GetBinLowEdge(1), histSig.GetXaxis().GetBinUpEdge(histSig.GetNbinsX()))
	for hist in histListBkg:
		histBkg.Add(hist,1)
	for iBin in range(histFOM.GetNbinsX()):
		histFOM.SetBinContent(iBin, FOM(histSig.Integral(iBin, histSig.GetNbinsX()), histBkg.Integral(iBin, histBkg.GetNbinsX())))
	canvas = rt.TCanvas("canvas","canvas",900,600)
	histFOM.Draw()
	canvas.SaveAs(Name+".png")
	print(Name + " " + str(histFOM.GetMaximum()) + " " + str(histFOM.GetXaxis().GetBinLowEdge(histFOM.GetMaximumBin())))
	histFOM.Delete()
	histBkg.Delete()



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

inFile = rt.TFile("../../"+createInFileName(mZprime, mDark,rinv,alpha),"read")
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
nEventsPassedPreSelection = 0
nEventsInGroup = 0

# initilize histograms
outputFile = rt.TFile("output.root","recreate")
hist_MT_bestReco = rt.TH1F("MT_bestReco","bestReco;MT;a.u.",100,0,2*mZprime)


hist_DR_12 = rt.TH1F("DR_12","DR_12;\Delta R;count/a.u.",100,0,6)
hist_DR_12_bm12 = rt.TH1F("DR_12_12","DR_12_12;\Delta R;count/a.u.",100,0,6)
hist_DR_12_bm123 = rt.TH1F("DR_12_123","DR_12_123;\Delta R;count/a.u.",100,0,6)

hist_DR_13 = rt.TH1F("DR_13","DR_13;\Delta R;count/a.u.",100,0,6)
hist_DR_13_bm12 = rt.TH1F("DR_13_12","DR_13_12;\Delta R;count/a.u.",100,0,6)
hist_DR_13_bm123 = rt.TH1F("DR_13_123","DR_13_123;\Delta R;count/a.u.",100,0,6)

hist_DR_23 = rt.TH1F("DR_23","DR_23;\Delta R;count/a.u.",100,0,6)
hist_DR_23_bm12 = rt.TH1F("DR_23_12","DR_23_12;\Delta R;count/a.u.",100,0,6)
hist_DR_23_bm123 = rt.TH1F("DR_23_123","DR_23_123;\Delta R;count/a.u.",100,0,6)

hist_MT12 = rt.TH1F("MT12","MT12;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT12_bm12 = rt.TH1F("MT12_bm12","MT12_bm12;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT12_bm123 = rt.TH1F("MT12_bm123","MT12_bm123;MT;count/a.u.",100,0,1.6*mZprime)

hist_MT13 = rt.TH1F("MT13","MT13;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT13_bm12 = rt.TH1F("MT13_bm12","MT13_bm12;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT13_bm123 = rt.TH1F("MT13_bm123","MT13_bm123;MT;count/a.u.",100,0,1.6*mZprime)

hist_MT123 = rt.TH1F("MT123","MT123;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT123_bm12 = rt.TH1F("MT123_bm12","MT123_bm12;MT;count/a.u.",100,0,1.6*mZprime)
hist_MT123_bm123 = rt.TH1F("MT123_bm123","MT123_bm123;MT;count/a.u.",100,0,1.6*mZprime)

hist_MET = rt.TH1F("MET","MET;MET;count/a.u.",100,0,0.6*mZprime)
hist_MET_bm12 = rt.TH1F("MET_bm12","MET_bm12;MET;count/a.u.",100,0,0.6*mZprime)
hist_MET_bm123 = rt.TH1F("MET_bm123","MET_bm123;MET;count/a.u.",100,0,0.6*mZprime)

hist_dPhi1 = rt.TH1F("\Delta\Phi1","\Delta\Phi1;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi1_bm12 = rt.TH1F("\Delta\Phi1_bm12","\Delta\Phi1_bm12;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi1_bm123 = rt.TH1F("\Delta\Phi1_bm123","\Delta\Phi1_bm123;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())

hist_dPhi2 = rt.TH1F("\Delta\Phi2","\Delta\Phi2;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi2_bm12 = rt.TH1F("\Delta\Phi2_bm12","\Delta\Phi2_bm12;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi2_bm123 = rt.TH1F("\Delta\Phi2_bm123","\Delta\Phi2_bm123;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())

hist_dPhi3 = rt.TH1F("\Delta\Phi3","\Delta\Phi3;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi3_bm12 = rt.TH1F("\Delta\Phi3_bm12","\Delta\Phi3_bm12;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhi3_bm123 = rt.TH1F("\Delta\Phi3_bm123","\Delta\Phi3_bm123;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())

hist_dPhiMin = rt.TH1F("\Delta\PhiMin","\Delta\PhiMin;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhiMin_bm12 = rt.TH1F("\Delta\PhiMin_bm12","\Delta\PhiMin_bm12;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())
hist_dPhiMin_bm123 = rt.TH1F("\Delta\PhiMin_bm123","\Delta\PhiMin_bm123;\Delta\Phi;count/a.u.",100,0,rt.TMath.Pi())

hist_METSignificance = rt.TH1F("METSignificance","METSignificance;METSignificance;count/a.u.",100,0,4000)
hist_METSignificance_bm12 = rt.TH1F("METSignificance_bm12","METSignificance_bm12;METSignificance;count/a.u.",100,0,4000)
hist_METSignificance_bm123 = rt.TH1F("METSignificance_bm123","METSignificance_bm123;METSignificance;count/a.u.",100,0,4000)

hist_jet3BoostedMomentum = rt.TH1F("jet3BoostedP","jet3BoostedP;Momentum;count/a.u.", 100,0,4000)
hist_jet3BoostedMomentum_bm12 = rt.TH1F("jet3BoostedP_bm12","jet3BoostedP_bm12;Momentum;count/a.u.", 100,0,4000)
hist_jet3BoostedMomentum_bm123 = rt.TH1F("jet3BoostedP_bm123","jet3BoostedP_bm123;Momentum;count/a.u.", 100,0,4000)

hist_jet3BoostedPt = rt.TH1F("jet3BoostedPt","jet3BoostedPt;Pt;count/a.u.", 100,0,4000)
hist_jet3BoostedPt_bm12 = rt.TH1F("jet3BoostedPt_bm12","jet3BoostedPt_bm12;Pt;count/a.u.", 100,0,4000)
hist_jet3BoostedPt_bm123 = rt.TH1F("jet3BoostedPt_bm123","jet3BoostedPt_bm123;Pt;count/a.u.", 100,0,4000)

hist_jet3BoostedEta = rt.TH1F("jet3BoostedEta","jet3BoostedEta;Eta;count/a.u.", 100, -5, 5)
hist_jet3BoostedEta_bm12 = rt.TH1F("jet3BoostedEta_bm12","jet3BoostedEta_bm12;Eta;count/a.u.", 100, -5, 5)
hist_jet3BoostedEta_bm123 = rt.TH1F("jet3BoostedEta_bm123","jet3BoostedEta_bm123;Eta;count/a.u.", 100, -5, 5)

hist_jet3BoostedPhi = rt.TH1F("jet3BoostedPhi","jet3BoostedPhi;Phi;count/a.u.", 100,-rt.TMath.Pi(),rt.TMath.Pi())
hist_jet3BoostedPhi_bm12 = rt.TH1F("jet3BoostedPhi_bm12","jet3BoostedPhi_bm12;Phi;count/a.u.", 100,-rt.TMath.Pi(),rt.TMath.Pi())
hist_jet3BoostedPhi_bm123 = rt.TH1F("jet3BoostedPhi_bm123","jet3BoostedPhi_bm123;Phi;count/a.u.", 100,-rt.TMath.Pi(),rt.TMath.Pi())

hist_jet3BoosteddEta = rt.TH1F("jet3BoosteddEta","jet3BoosteddEta;dEta;count/a.u.", 100, 0, 6)
hist_jet3BoosteddEta_bm12 = rt.TH1F("jet3BoosteddEta_bm12","jet3BoosteddEta_bm12;dEta;count/a.u.", 100, 0, 6)
hist_jet3BoosteddEta_bm123 = rt.TH1F("jet3BoosteddEta_bm123","jet3BoosteddEta_bm123;dEta;count/a.u.", 100, 0, 6)

hist_jet3BoosteddPhi = rt.TH1F("jet3BoosteddPhi","jet3BoosteddPhi;dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_jet3BoosteddPhi_bm12 = rt.TH1F("jet3BoosteddPhi_bm12","jet3BoosteddPhi_bm12;dPhi;count/a.u.", 100,0,rt.TMath.Pi())
hist_jet3BoosteddPhi_bm123 = rt.TH1F("jet3BoosteddPhi_bm123","jet3BoosteddPhi_bm123;dPhi;count/a.u.", 100,0,rt.TMath.Pi())

hist_jet3BoosteddR = rt.TH1F("jet3BoosteddR","jet3BoosteddR;dR;count/a.u.", 100, 0,7)
hist_jet3BoosteddR_bm12 = rt.TH1F("jet3BoosteddR_bm12","jet3BoosteddR_bm12;dER;count/a.u.", 100, 0,7)
hist_jet3BoosteddR_bm123 = rt.TH1F("jet3BoosteddR_bm123","jet3BoosteddR_bm123;dR;count/a.u.", 100, 0,7)

hist_BoostVecMag = rt.TH1F("BoostVecMag","BoostVecMag;Mag;count/a.u.",100,0,1)
hist_BoostVecMag_bm12 = rt.TH1F("BoostVecMag_bm12","BoostVecMag_bm12;Mag;count/a.u.",100,0,1)
hist_BoostVecMag_bm123 = rt.TH1F("BoostVecMag_bm123","BoostVecMag_bm123;Mag;count/a.u.",100,0,1)

hist_jet12boosted = rt.TH1F("testBoosted", "pt;pt;count/a.u.", 1000, -50, 4050)


#2d plots. varios DR's for each best mass case, and for all...

hist_2DR_12v13_All = rt.TH2F("2DR_12v13_All","2DR_12v13_All;\Delta R 12;\Delta R 13", 100,0,6,100,0,6)
hist_2DR_12v13_bm12 = rt.TH2F("2DR_12v13_bm12","2DR_12v13_bm12;\Delta R 12;\Delta R 13", 100,0,6,100,0,6)
hist_2DR_12v13_bm123 = rt.TH2F("2DR_12v13_bm123","2DR_12v13_bm123;\Delta R 12;\Delta R 13", 100,0,6,100,0,6)

hist_2DR_12v23_All = rt.TH2F("2DR_12v23_All","2DR_12v23_All;\Delta R 12;\Delta R 23", 100,0,6,100,0,6)
hist_2DR_12v23_bm12 = rt.TH2F("2DR_12v23_bm12","2DR_12v23_bm12;\Delta R 12;\Delta R 23", 100,0,6,100,0,6)
hist_2DR_12v23_bm123 = rt.TH2F("2DR_12v23_bm123","2DR_12v23_bm123;\Delta R 12;\Delta R 23", 100,0,6,100,0,6)

hist_2DR_13v23_All = rt.TH2F("2DR_13v23_All","2DR_13v23_All;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)
hist_2DR_13v23_bm12 = rt.TH2F("2DR_13v23_bm12","2DR_13v23_bm12;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)
hist_2DR_13v23_bm123 = rt.TH2F("2DR_13v23_bm123","2DR_13v23_bm123;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)

hist_2DRboost_12v13 = rt.TH2F("2DRboost_12v13","2DRboost_12v13;\Delta R 12;\Delta R 13", 100,3,6,100,0,6)
hist_2DRboost_12v13_bm12 = rt.TH2F("2DRboost_12v13_bm12","2DRboost_12v13_bm12;\Delta R 12;\Delta R 13", 100,3,6,100,0,6)
hist_2DRboost_12v13_bm123 = rt.TH2F("2DRboost_12v13_bm123","2DRboost_12v13_bm123;\Delta R 12;\Delta R 13", 100,3,6,100,0,6)

hist_2DRboost_12v23 = rt.TH2F("2DRboost_12v23","2DRboost_12v23;\Delta R 12;\Delta R 23", 100,3,6,100,0,6)
hist_2DRboost_12v23_bm12 = rt.TH2F("2DRboost_12v23_bm12","2DRboost_12v23_bm12;\Delta R 12;\Delta R 23", 100,3,6,100,0,6)
hist_2DRboost_12v23_bm123 = rt.TH2F("2DRboost_12v23_bm123","2DRboost_12v23_bm123;\Delta R 12;\Delta R 23", 100,3,6,100,0,6)


hist_2DRboost_13v23 = rt.TH2F("2DRboost_13v23","2DRboost_13v23;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)
hist_2DRboost_13v23_bm12 = rt.TH2F("2DRboost_13v23_bm12","2DRboost_13v23_bm12;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)
hist_2DRboost_13v23_bm123 = rt.TH2F("2DRboost_13v23_bm123","2DRboost_13v23_bm123;\Delta R 13;\Delta R 23", 100,0,6,100,0,6)

#2d plots fo deltaphi
hist_2dPhi_dPhi1vdPhi2_All = rt.TH2F("2dPhi_dPhi1vdPhi2_All","2dPhi_dPhi1vdPhi2_All;\Delta\Phi1;\Delta\Phi2",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi1vdPhi2_bm12 = rt.TH2F("2dPhi_dPhi1vdPhi2_bm12","2dPhi_dPhi1vdPhi2_bm12;\Delta\Phi1;\Delta\Phi2",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi1vdPhi2_bm123 = rt.TH2F("2dPhi_dPhi1vdPhi2_bm123","2dPhi_dPhi1vdPhi2_bm123;\Delta\Phi1;\Delta\Phi2",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())

hist_2dPhi_dPhi1vdPhi3_All = rt.TH2F("2dPhi_dPhi1vdPhi3_All","2dPhi_dPhi1vdPhi3_All;\Delta\Phi1;\Delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi1vdPhi3_bm12 = rt.TH2F("2dPhi_dPhi1vdPhi3_bm12","2dPhi_dPhi1vdPhi3_bm12;\Delta\Phi1;\Delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi1vdPhi3_bm123 = rt.TH2F("2dPhi_dPhi1vdPhi3_bm123","2dPhi_dPhi1vdPhi3_bm123;\Delta\Phi1;\Delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())

hist_2dPhi_dPhi2vdPhi3_All = rt.TH2F("2dPhi_dPhi2vdPhi3_All","2dPhi_dPhi2vdPhi3_All;\Delta\Phi2;\delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi2vdPhi3_bm12 = rt.TH2F("2dPhi_dPhi2vdPhi3_bm12","2dPhi_dPhi2vdPhi3_bm12;\Delta\Phi2;\delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())
hist_2dPhi_dPhi2vdPhi3_bm123 = rt.TH2F("2dPhi_dPhi2vdPhi3_bm123","2dPhi_dPhi2vdPhi3_bm123;\Delta\Phi2;\delta\Phi3",100,0,rt.TMath.Pi(),100,0,rt.TMath.Pi())


histList = []
histList.append(hist_MT_bestReco)
histList.append(hist_DR_12)

histList2d = []
histList2d.append(hist_2DR_12v13_All)
histList2d.append(hist_2DR_12v13_bm12)
histList2d.append(hist_2DR_12v13_bm123)
histList2d.append(hist_2DR_12v23_All)
histList2d.append(hist_2DR_12v23_bm12)
histList2d.append(hist_2DR_12v23_bm123)
histList2d.append(hist_2DR_13v23_All)
histList2d.append(hist_2DR_13v23_bm12)
histList2d.append(hist_2DR_13v23_bm123)
histList2d.append(hist_2dPhi_dPhi1vdPhi2_All)
histList2d.append(hist_2dPhi_dPhi1vdPhi2_bm12)
histList2d.append(hist_2dPhi_dPhi1vdPhi2_bm123)
histList2d.append(hist_2dPhi_dPhi1vdPhi3_All)
histList2d.append(hist_2dPhi_dPhi1vdPhi3_bm12)
histList2d.append(hist_2dPhi_dPhi1vdPhi3_bm123)
histList2d.append(hist_2dPhi_dPhi2vdPhi3_All)
histList2d.append(hist_2dPhi_dPhi2vdPhi3_bm12)
histList2d.append(hist_2dPhi_dPhi2vdPhi3_bm123)
histList2d.append(hist_2DRboost_12v13)
histList2d.append(hist_2DRboost_12v13_bm12)
histList2d.append(hist_2DRboost_12v13_bm123)
histList2d.append(hist_2DRboost_12v23)
histList2d.append(hist_2DRboost_12v23_bm12)
histList2d.append(hist_2DRboost_12v23_bm123)
histList2d.append(hist_2DRboost_13v23)
histList2d.append(hist_2DRboost_13v23_bm12)
histList2d.append(hist_2DRboost_13v23_bm123)





hist_MT_bestReco.SetLineStyle(2)


jetComboBestMassDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboDeltaRTag = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetComboHVQuarksDict = {'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}
jetMETXDict = {'MET1':0,'MET2':0,'MET3':0,'MET4':0}
IndexToMET = {0:'MET1',1:'MET2',2:'MET3',3:'MET4'}


jetComboDeltaRCutBMTruth = {
'111':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'110':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'101':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'011':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'100':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'010':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'001':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0},
'000':{'':0,'1':0,'2':0,'3':0,'4':0,'12':0,'13':0,'14':0,'23':0,'24':0,'34':0,'123':0,'124':0,'134':0,'234':0,'1234':0}}

nEventsAK83CA112 = 0
nEventsCA11MTisAK83MT = 0

for iEvent in range(nEvents):
	jetComboMassDict = {'':0,'1':0,'2':0,'12':0}
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	jetCollection = tree.JetsAK8

	
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
	
	if not (len(jetCollection) == 3): # for tseeing the effects  has on these plots. Make sure to chance the createImagename function above
		continue
	
	#calculat the deltaphi min of the first four leading pT jets, and tag which jet is most cloesly aligned with the MET
	DeltaPhi1 = tree.DeltaPhi1_AK8
	DeltaPhi2 = tree.DeltaPhi2_AK8
	DeltaPhi3 = 99
	DeltaPhi4 = 99
	if len(jetCollection) >= 3:
		DeltaPhi3 = absDphi(jetCollection[2].Phi(),tree.METPhi)
	if len(jetCollection) > 3:
		DeltaPhi4 = absDphi(jetCollection[3].Phi(),tree.METPhi)
	DeltaPhiList = [DeltaPhi1, DeltaPhi2, DeltaPhi3,DeltaPhi4]
	DeltaPhiMin = min(DeltaPhiList)
	METXkey = IndexToMET[DeltaPhiList.index(DeltaPhiMin)]

	#if not (METXkey == 'MET3'):
	#	continue
	nEventsInGroup += 1
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
			for iJet in range(min(len(jetCollection),4)):
				if (not jetsWithHVQuarks[iJet]) and (jetCollection[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8):
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
			listOfJetsWithParticlesFromHVQuarks.append(jetCollection[x])
			jetsHVQuarksKey += str(x+1)
	jetComboHVQuarksDict[jetsHVQuarksKey] += 1
	MTofHVQuarksProducts = trans_mass_Njet(listOfJetsWithParticlesFromHVQuarks,tree.MET,tree.METPhi)
	#create all of the MT masses
	jetComboMassDict[''] = trans_mass_Njet([],tree.MET,tree.METPhi)
	jetComboMassDict['1'] = trans_mass_Njet([jetCollection[0]],tree.MET,tree.METPhi)
	jetComboMassDict['2'] = trans_mass_Njet([jetCollection[1]],tree.MET,tree.METPhi)
	jetComboMassDict['12'] = trans_mass_Njet([jetCollection[0], jetCollection[1]],tree.MET,tree.METPhi)
	if len(jetCollection) >= 3:
		jetComboMassDict['3'] = trans_mass_Njet([jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['13'] = trans_mass_Njet([jetCollection[0], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['23'] = trans_mass_Njet([jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
		jetComboMassDict['123'] = trans_mass_Njet([jetCollection[0], jetCollection[1], jetCollection[2]],tree.MET,tree.METPhi)
	massDiffList = {key : abs(x - MTofHVQuarksProducts) for (key,x) in jetComboMassDict.items()}
	BestMassKey = min(massDiffList,key=massDiffList.get)
	jetComboBestMassDict[BestMassKey] += 1
	hist_MT_bestReco.Fill(jetComboMassDict[min(massDiffList,key=massDiffList.get)])
	if not (BestMassKey == jetsHVQuarksKey):
		print("Best Mass is not HV Quarks! " + BestMassKey + " " + jetsHVQuarksKey)

	dR12 = jetCollection[0].DeltaR(jetCollection[1])
	if len(jetCollection) >= 3:
		dR13 = jetCollection[0].DeltaR(jetCollection[2])
		dR23 = jetCollection[1].DeltaR(jetCollection[2])

	#boost the 3rd jet into the CoM frame of the leading two jets
	if len(jetCollection) >= 3:
		boostedJet3 = rt.TLorentzVector(jetCollection[2].X(),jetCollection[2].Y(),jetCollection[2].Z(),jetCollection[2].T())
		boostedJet2 = rt.TLorentzVector(jetCollection[1].X(),jetCollection[1].Y(),jetCollection[1].Z(),jetCollection[1].T())
		boostedJet1 = rt.TLorentzVector(jetCollection[0].X(),jetCollection[0].Y(),jetCollection[0].Z(),jetCollection[0].T())
		jet12vector = jetCollection[0]+jetCollection[1]
		boost = -(jet12vector).BoostVector()
		boostedJet3.Boost(boost)
		boostedJet2.Boost(boost)
		boostedJet1.Boost(boost)
		jet12vector.Boost(boost)
		
	hist_DR_12.Fill(dR12)
	hist_DR_13.Fill(dR13)
	hist_DR_23.Fill(jetCollection[1].DeltaR(jetCollection[2]))
	hist_2DR_12v13_All.Fill(dR12,dR13)
	hist_2DR_12v23_All.Fill(dR12,dR23)
	hist_2DR_13v23_All.Fill(dR13,dR23)
	hist_2DRboost_12v13.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet1.DeltaR(boostedJet3))
	hist_2DRboost_12v23.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet2.DeltaR(boostedJet3))
	hist_2DRboost_13v23.Fill(boostedJet1.DeltaR(boostedJet3),boostedJet2.DeltaR(boostedJet3))
	hist_MT12.Fill(jetComboMassDict['12'])
	hist_MT13.Fill(jetComboMassDict['13'])
	hist_MT123.Fill(jetComboMassDict['123'])
	hist_MET.Fill(tree.MET)
	hist_dPhi1.Fill(DeltaPhi1)
	hist_dPhi2.Fill(DeltaPhi2)
	hist_dPhi3.Fill(DeltaPhi3)
	hist_dPhiMin.Fill(DeltaPhiMin)
	hist_METSignificance.Fill(tree.METSignificance)
	hist_2dPhi_dPhi1vdPhi2_All.Fill(DeltaPhi1,DeltaPhi2)
	hist_2dPhi_dPhi1vdPhi3_All.Fill(DeltaPhi1,DeltaPhi3)
	hist_2dPhi_dPhi2vdPhi3_All.Fill(DeltaPhi2,DeltaPhi3)
	hist_jet3BoostedMomentum.Fill(boostedJet3.P())
	hist_jet3BoostedPt.Fill(boostedJet3.Pt())
	hist_jet3BoostedPhi.Fill(boostedJet3.Phi())
	hist_jet3BoostedEta.Fill(boostedJet3.Eta())
	hist_jet3BoosteddPhi.Fill(absDphi(boostedJet3.Phi(),boostedJet1.Phi()))
	hist_jet3BoosteddEta.Fill(rt.TMath.Abs(boostedJet3.Eta()-boostedJet1.Eta()))
	hist_jet3BoosteddR.Fill(boostedJet3.DeltaR(boostedJet1))
	hist_jet12boosted.Fill(jet12vector.P())
	hist_BoostVecMag.Fill(boost.Mag())
	if not (BestMassKey == '123'):
		hist_DR_12_bm12.Fill(dR12)
		hist_DR_13_bm12.Fill(dR13)
		hist_DR_23_bm12.Fill(dR23)
		hist_2DR_12v13_bm12.Fill(dR12,dR13)
		hist_2DR_12v23_bm12.Fill(dR12,dR23)
		hist_2DR_13v23_bm12.Fill(dR13,dR23)
		hist_2DRboost_12v13_bm12.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet1.DeltaR(boostedJet3))
		hist_2DRboost_12v23_bm12.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet2.DeltaR(boostedJet3))
		hist_2DRboost_13v23_bm12.Fill(boostedJet1.DeltaR(boostedJet3),boostedJet2.DeltaR(boostedJet3))
		hist_MT12_bm12.Fill(jetComboMassDict['12'])
		hist_MT123_bm12.Fill(jetComboMassDict['123'])
		hist_MET_bm12.Fill(tree.MET)
		hist_dPhi1_bm12.Fill(DeltaPhi1)
		hist_dPhi2_bm12.Fill(DeltaPhi2)
		hist_dPhi3_bm12.Fill(DeltaPhi3)
		hist_dPhiMin_bm12.Fill(DeltaPhiMin)
		hist_METSignificance_bm12.Fill(tree.METSignificance)
		hist_2dPhi_dPhi1vdPhi2_bm12.Fill(DeltaPhi1,DeltaPhi2)
		hist_2dPhi_dPhi1vdPhi3_bm12.Fill(DeltaPhi1,DeltaPhi3)
		hist_2dPhi_dPhi2vdPhi3_bm12.Fill(DeltaPhi2,DeltaPhi3)
		hist_jet3BoostedMomentum_bm12.Fill(boostedJet3.P())
		hist_jet3BoostedPt_bm12.Fill(boostedJet3.Pt())
		hist_jet3BoostedPhi_bm12.Fill(boostedJet3.Phi())
		hist_jet3BoostedEta_bm12.Fill(boostedJet3.Eta())
		hist_jet3BoosteddPhi_bm12.Fill(absDphi(boostedJet3.Phi(),boostedJet1.Phi()))
		hist_jet3BoosteddEta_bm12.Fill(rt.TMath.Abs(boostedJet3.Eta()-boostedJet1.Eta()))
		hist_jet3BoosteddR_bm12.Fill(boostedJet3.DeltaR(boostedJet1))
		hist_BoostVecMag_bm12.Fill(boost.Mag())
	elif BestMassKey == '123':
		hist_DR_12_bm123.Fill(dR12)
		hist_DR_13_bm123.Fill(dR13)
		hist_DR_23_bm123.Fill(dR23)
		hist_2DR_12v13_bm123.Fill(dR12,dR13)
		hist_2DR_12v23_bm123.Fill(dR12,dR23)
		hist_2DR_13v23_bm123.Fill(dR13,dR23)
		hist_2DRboost_12v13_bm123.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet1.DeltaR(boostedJet3))
		hist_2DRboost_12v23_bm123.Fill(boostedJet1.DeltaR(boostedJet2),boostedJet2.DeltaR(boostedJet3))
		hist_2DRboost_13v23_bm123.Fill(boostedJet1.DeltaR(boostedJet3),boostedJet2.DeltaR(boostedJet3))
		hist_MT12_bm123.Fill(jetComboMassDict['12'])
		hist_MT123_bm123.Fill(jetComboMassDict['123'])
		hist_MET_bm123.Fill(tree.MET)
		hist_dPhi1_bm123.Fill(DeltaPhi1)
		hist_dPhi2_bm123.Fill(DeltaPhi2)
		hist_dPhi3_bm123.Fill(DeltaPhi3)
		hist_dPhiMin_bm123.Fill(DeltaPhiMin)
		hist_METSignificance_bm123.Fill(tree.METSignificance)
		hist_2dPhi_dPhi1vdPhi2_bm123.Fill(DeltaPhi1,DeltaPhi2)
		hist_2dPhi_dPhi1vdPhi3_bm123.Fill(DeltaPhi1,DeltaPhi3)
		hist_2dPhi_dPhi2vdPhi3_bm123.Fill(DeltaPhi2,DeltaPhi3)
		hist_jet3BoostedMomentum_bm123.Fill(boostedJet3.P())
		hist_jet3BoostedPt_bm123.Fill(boostedJet3.Pt())
		hist_jet3BoostedPhi_bm123.Fill(boostedJet3.Phi())
		hist_jet3BoostedEta_bm123.Fill(boostedJet3.Eta())
		hist_jet3BoosteddPhi_bm123.Fill(absDphi(boostedJet3.Phi(),boostedJet1.Phi()))
		hist_jet3BoosteddEta_bm123.Fill(rt.TMath.Abs(boostedJet3.Eta()-boostedJet1.Eta()))
		hist_jet3BoosteddR_bm123.Fill(boostedJet3.DeltaR(boostedJet1))
		hist_BoostVecMag_bm123.Fill(boost.Mag())

	if BestMassKey == '123':
		nEventsAK83CA112 += 1
		if isNear(jetComboMassDict['123'], tree.MT_CA11, 5.0): 
			nEventsCA11MTisAK83MT += 1




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
print("nEvents with 3 AK8 and 2 CA11 = " + str(nEventsAK83CA112))
print("nEvents where MT_CA11[0,1] is near MT_AK8[0,1,2] = "+ str(nEventsCA11MTisAK83MT))
print("Preselection " +str(nEventsPassedPreSelection))
print("Total " + str(nEvents))
nEventsInGroup = hist_DR_12.GetEntries()
print("nEvents in Group = " + str(nEventsInGroup))
print("nEvents with 12 " + str(hist_DR_12_bm12.GetEntries()) + " (" + str(float(hist_DR_12_bm12.GetEntries())/float(nEventsInGroup)) + ")")
print("nEvents with 123 " + str(hist_DR_12_bm123.GetEntries()) + " (" + str(float(hist_DR_12_bm123.GetEntries())/float(nEventsInGroup)) + ")")
print("dRCut 0 1 2 3 4 12 13 14 23 24 34 123 124 234 1234")
for (key, value1) in jetComboDeltaRCutBMTruth.items():
	print(key + " " +str(jetComboDeltaRCutBMTruth[key]['']) + " " +str(jetComboDeltaRCutBMTruth[key]['1']) + " " +str(jetComboDeltaRCutBMTruth[key]['2']) + " " +str(jetComboDeltaRCutBMTruth[key]['3']) + " " +str(jetComboDeltaRCutBMTruth[key]['4']) + " " +str(jetComboDeltaRCutBMTruth[key]['12']) + " " +str(jetComboDeltaRCutBMTruth[key]['13']) + " " +str(jetComboDeltaRCutBMTruth[key]['14']) + " " +str(jetComboDeltaRCutBMTruth[key]['23']) + " " +str(jetComboDeltaRCutBMTruth[key]['24']) + " " +str(jetComboDeltaRCutBMTruth[key]['34']) + " " +str(jetComboDeltaRCutBMTruth[key]['123']) + " " +str(jetComboDeltaRCutBMTruth[key]['124']) + " " +str(jetComboDeltaRCutBMTruth[key]['234']) + " " +str(jetComboDeltaRCutBMTruth[key]['1234']))

drawHistos([hist_DR_12,hist_DR_13,hist_DR_23],"DeltaR_All",True)
drawHistos([hist_DR_12_bm12,hist_DR_13_bm12,hist_DR_23_bm12],"DeltaR_bm12",True)
drawHistos([hist_DR_12_bm123,hist_DR_13_bm123,hist_DR_23_bm123],"DeltaR_bm123",True)

drawHistos([hist_DR_12,hist_DR_12_bm12,hist_DR_12_bm123],"DeltaR_between12",True)
drawHistos([hist_DR_13,hist_DR_13_bm12,hist_DR_13_bm123],"DeltaR_between13",True)
drawHistos([hist_DR_23,hist_DR_23_bm12,hist_DR_23_bm123],"DeltaR_between23",True)

drawHistos([hist_MT12,hist_MT13,hist_MT123],"MT_All",True)
drawHistos([hist_MT12_bm12,hist_MT13_bm12,hist_MT123_bm12],"MT_bm12",True)
drawHistos([hist_MT12_bm123,hist_MT13_bm123,hist_MT123_bm123],"MT_bm123",True)

drawHistos([hist_MT12,hist_MT12_bm12,hist_MT12_bm123],"MT12",True)
drawHistos([hist_MT13,hist_MT13_bm12,hist_MT13_bm123],"MT13",True)
drawHistos([hist_MT123,hist_MT123_bm12,hist_MT123_bm123],"MT123",True)
drawHistos([hist_MET,hist_MET_bm12,hist_MET_bm123],"MET",True)
drawHistos([hist_jet12boosted],'testPlot',True)
drawHistos([hist_dPhi1,hist_dPhi1_bm12,hist_dPhi1_bm123],"DeltaPhi1",True)
drawHistos([hist_dPhi2,hist_dPhi2_bm12,hist_dPhi2_bm123],"DeltaPhi2",True)
drawHistos([hist_dPhi3,hist_dPhi3_bm12,hist_dPhi3_bm123],"DeltaPhi3",True)
drawHistos([hist_dPhiMin,hist_dPhiMin_bm12,hist_dPhiMin_bm123],"DeltaPhiMin",True)

drawHistos([hist_jet3BoostedMomentum,hist_jet3BoostedMomentum_bm12,hist_jet3BoostedMomentum_bm123],"jet3BoostedMomentum",True)
drawHistos([hist_jet3BoostedPt,hist_jet3BoostedPt_bm12,hist_jet3BoostedPt_bm123],"jet3BoostedPt",True)
drawHistos([hist_jet3BoostedPhi,hist_jet3BoostedPhi_bm12,hist_jet3BoostedPhi_bm123],"jet3BoostedPhi",True)
drawHistos([hist_jet3BoostedEta,hist_jet3BoostedEta_bm12,hist_jet3BoostedEta_bm123],"jet3BoostedEta",True)
drawHistos([hist_jet3BoosteddPhi,hist_jet3BoosteddPhi_bm12,hist_jet3BoosteddPhi_bm123],"jet3BoosteddPhi",True)
drawHistos([hist_jet3BoosteddEta,hist_jet3BoosteddEta_bm12,hist_jet3BoosteddEta_bm123],"jet3BoosteddEta",True)
drawHistos([hist_jet3BoosteddR,hist_jet3BoosteddR_bm12,hist_jet3BoosteddR_bm123],"jet3BoosteddR",True)
drawHistos([hist_BoostVecMag,hist_BoostVecMag_bm12,hist_BoostVecMag_bm123],"BoostVecMag",True)

drawHistos([hist_METSignificance,hist_METSignificance_bm123,hist_METSignificance_bm12],"METSignificance",True)

canv1 = rt.TCanvas()
for hist in histList2d:
	hist.Draw('colz')
	canv1.SaveAs(createImageName(mZprime, mDark, rinv, alpha, hist.GetName()))

makeFOMPlot(hist_DR_12_bm123, [hist_DR_12_bm12],"DR12_bm123")
makeFOMPlot(hist_DR_13_bm123, [hist_DR_13_bm12],"DR13_bm123")
makeFOMPlot(hist_DR_23_bm123, [hist_DR_23_bm12],"DR23_bm123")

makeFOMPlot(hist_MT12_bm123, [hist_MT12_bm12],"MT12_bm123")
makeFOMPlot(hist_MT13_bm123, [hist_MT13_bm12],"MT13_bm123")
makeFOMPlot(hist_MT123_bm123, [hist_MT123_bm12],"MT123_bm123")

makeFOMPlot(hist_MET_bm123, [hist_MET_bm12],"MET_bm123")

makeFOMPlot(hist_dPhi1_bm123, [hist_dPhi1_bm12],"dPhi1_bm123")
makeFOMPlot(hist_dPhi2_bm123, [hist_dPhi2_bm12],"dPhi2_bm123")
makeFOMPlot(hist_dPhi3_bm123, [hist_dPhi3_bm12],"dPhi3_bm123")

makeFOMPlot(hist_DR_12_bm12, [hist_DR_12_bm123],"DR12_bm12")
makeFOMPlot(hist_DR_13_bm12, [hist_DR_13_bm123],"DR13_bm12")
makeFOMPlot(hist_DR_23_bm12, [hist_DR_23_bm123],"DR23_bm12")

makeFOMPlot(hist_MT12_bm12, [hist_MT12_bm123],"MT12_bm12")
makeFOMPlot(hist_MT13_bm12, [hist_MT13_bm123],"MT13_bm12")
makeFOMPlot(hist_MT123_bm12, [hist_MT123_bm123],"MT123_bm12")

makeFOMPlot(hist_MET_bm12, [hist_MET_bm123],"MET_bm12")

makeFOMPlot(hist_dPhi1_bm12, [hist_dPhi1_bm123],"dPhi1_bm12")
makeFOMPlot(hist_dPhi2_bm12, [hist_dPhi2_bm123],"dPhi2_bm12")
makeFOMPlot(hist_dPhi3_bm12, [hist_dPhi3_bm123],"dPhi3_bm12")



outputFile.Write()
outputFile.Close()


inFile.Close()
