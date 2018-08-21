import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetMTTests_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

def plotTwo2D(hist1, hist2,outputName):
	canv =rt.TCanvas()
	hist1.SetMarkerStyle(2)
	hist2.SetMarkerStyle(3)
	hist1.SetMarkerColorAlpha(rt.kBlue,0.5)
	hist2.SetMarkerColorAlpha(rt.kGreen,0.5)
	hist1.Draw()
	hist2.Draw("same")
	canv.SaveAs(outputName)

def makeCornerPlot3D(inputHisto, outputName):
	#make canvas, and split into 9 sections:
	canvas = rt.TCanvas("canvas","canvas", 900,900)
	#  x <>  <>   1 2 3
	# yx  y  <>   4 5 6
	# zx zy   z   7 8 9
	canvas.Divide(3,3,0,0)
	T = inputHisto.GetTitle()
	X = inputHisto.GetXaxis().GetTitle()
	Y = inputHisto.GetYaxis().GetTitle()
	Z = inputHisto.GetZaxis().GetTitle()

	hist_X = inputHisto.Project3D("x")
	hist_X.GetXaxis().SetLabelSize(0)
	hist_X.GetYaxis().CenterTitle()
	hist_X.SetTitle(";;"+X)

	#hist_XY = inputHisto.Project3D("xy")
	#hist_XY.GetYaxis().SetLabelSize(0)
	#hist_XY.GetXaxis().SetLabelSize(0)
	#hist_XY.SetTitle(";;")
	
	#hist_XZ = inputHisto.Project3D("xz")
	#hist_XZ.GetYaxis().SetLabelSize(0)
	#hist_XZ.GetXaxis().SetLabelSize(0)
	#hist_XZ.SetTitle(";;")
	
	hist_YX = inputHisto.Project3D("yx")
	hist_YX.GetXaxis().SetLabelSize(0)
	hist_YX.GetYaxis().CenterTitle()
	hist_YX.SetTitle(";"+X+";"+Y)
	
	hist_Y = inputHisto.Project3D("y")
	hist_Y.GetYaxis().SetLabelSize(0)
	hist_Y.GetXaxis().SetLabelSize(0)
	hist_Y.SetTitle(";;")
	
	#hist_YZ = inputHisto.Project3D("yz")
	#hist_YZ.GetXaxis().SetLabelSize(0)
	#hist_YZ.GetYaxis().SetLabelSize(0)
	#hist_YZ.SetTitle(";;")
	
	hist_ZX = inputHisto.Project3D("zx")
	hist_ZX.GetXaxis().CenterTitle()
	hist_ZX.GetYaxis().CenterTitle()
	hist_ZX.SetTitle(";"+X+";"+Z)
	
	hist_ZY = inputHisto.Project3D("zy")
	hist_ZY.GetYaxis().SetLabelSize(0)
	hist_ZY.GetXaxis().CenterTitle()
	hist_ZY.SetTitle(";"+Y+";")
	
	hist_Z = inputHisto.Project3D("z")
	hist_Z.GetYaxis().SetLabelSize(0)
	hist_Z.GetXaxis().CenterTitle()
	hist_Z.SetTitle(";"+Z+";")
	
	canvas.cd(1)
	hist_X.DrawCopy()
	#canvas.cd(2)
	#hist_XY.DrawCopy("cont")
	#canvas.cd(3)
	#hist_XZ.DrawCopy("cont")
	canvas.cd(4)
	hist_YX.DrawCopy("cont")
	canvas.cd(5)
	hist_Y.DrawCopy()
	#canvas.cd(6)
	#hist_YZ.DrawCopy("cont")
	canvas.cd(7)
	hist_ZX.DrawCopy("cont")
	canvas.cd(8)
	hist_ZY.DrawCopy("cont")
	canvas.cd(9)
	hist_Z.DrawCopy()
	canvas.SaveAs(outputName)

def makeCornerPlot2D(inputHisto, outputName):
	#make canvas, and split into 4 sections:
	canvas = rt.TCanvas("canvas","canvas", 900,900)
	#  x <>   1 2
	# yx  y   3 4

	canvas.Divide(2,2,0,0)

	T = inputHisto.GetTitle()
	X = inputHisto.GetXaxis().GetTitle()
	Y = inputHisto.GetYaxis().GetTitle()

	hist_YX = inputHisto
	hist_YX.SetTitle("")
	hist_YX.GetXaxis().CenterTitle()
	hist_YX.GetYaxis().CenterTitle()
	hist_X = inputHisto.ProjectionX()
	hist_X.GetXaxis().SetLabelSize(0)
	hist_X.GetYaxis().CenterTitle()
	hist_Y = inputHisto.ProjectionY()
	hist_Y.GetYaxis().SetLabelSize(0)
	hist_Y.GetXaxis().CenterTitle()

	canvas.cd(3)
	hist_YX.DrawCopy("cont")
	canvas.cd(1)
	hist_X.DrawCopy()
	
	canvas.cd(4)
	hist_Y.DrawCopy("hbar2")
	
	canvas.SaveAs(outputName)

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

# initilize histograms

hist_MT_All = rt.TH1F("All MT", "All;MT;a.u.",100,0,1.3*mZprime)
hist_MT_2jets = rt.TH1F("2jets MT", "nJets==2;MT;a.u.",100,0,1.3*mZprime)
hist_MT_3jets = rt.TH1F("3jets MT", "nJets==3;MT;a.u.",100,0,1.3*mZprime)
hist_MT_4jets = rt.TH1F("4jets MT", "nJets>=4;MT;a.u.",100,0,1.3*mZprime)

hist_MT3jet_3jets = rt.TH1F("3jets MT3jet", "3jet MT;MT",100,0,1.8*mZprime)
hist_MTvsMT3jet_3jets = rt.TH2F("3jets MTvsMT3jet", ";MT;3Jet MT",100,0,1.3*mZprime,100,0,1.8*mZprime)
prof_MTvsMT3jet_3jets = rt.TProfile("3jets MTvsMT3jet profile","MT;3jetMT;a.u.",100,0,1.3*mZprime,0,1.8*mZprime)

hist_MVsEtaVsPhi = rt.TH3F("","T;Mass;Eta;Phi",100,0,30*int(mDark),100,-3,3,100,-rt.TMath.Pi(),rt.TMath.Pi())
hist_MtMt3Jet_dPhiMin2Jet = rt.TH3F("","T;MT;dPhiMin2jet;MT3Jet",100,0,1.3*mZprime,100,0,rt.TMath.Pi(),100,0,1.8*mZprime)
hist_MtMt3Jet_dPhiMin3Jet = rt.TH3F("","T;MT;dPhiMin3jet;MT3Jet",100,0,1.3*mZprime,100,0,rt.TMath.Pi(),100,0,1.8*mZprime)
hist_MtMt3Jet_MET = rt.TH3F("","T;MT;MET;MT3Jet",100,0,1.3*mZprime,100,0,0.7*mZprime,100,0,1.8*mZprime)


histList = []
histList.append(hist_MT_All)
histList.append(hist_MT_2jets)
histList.append(hist_MT_3jets)
histList.append(hist_MT_4jets)

histList.append(hist_MT3jet_3jets)

nEventsPassedPreSelection = 0
nEventsWithOnly2Jets = 0
nEventsWithOnly3Jets = 0
nEventsWith4OrMoreJets = 0

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	#jetCollection = tree.GenJetsAK8
	jetCollection = tree.JetsAK8
	nJets = len(jetCollection)
	
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
	
	#do stuff for all events
	hist_MT_All.Fill(tree.MT_AK8)
	for jet in jetCollection:
		hist_MVsEtaVsPhi.Fill(jet.M(),jet.Eta(),jet.Phi())
	if nJets == 2: #do stuff here for events with 2 jets
		nEventsWithOnly2Jets += 1
		hist_MT_2jets.Fill(tree.MT_AK8)
	elif nJets == 3: #do stuff here for events wth 3 jets
		nEventsWithOnly3Jets += 1
		DeltaPhi3_AK8 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
		DeltaPhiMin_3jet = min(DeltaPhi3_AK8, tree.DeltaPhi1_AK8, tree.DeltaPhi2_AK8)
		hist_MT_3jets.Fill(tree.MT_AK8)
		MT3jet = trans_mass_Njet(jetCollection[0:3],tree.MET,tree.METPhi)
		hist_MT3jet_3jets.Fill(MT3jet)
		hist_MTvsMT3jet_3jets.Fill(tree.MT_AK8,MT3jet)
		hist_MtMt3Jet_dPhiMin2Jet.Fill(tree.MT_AK8,tree.DeltaPhiMin_AK8,MT3jet)
		hist_MtMt3Jet_dPhiMin3Jet.Fill(tree.MT_AK8,DeltaPhiMin_3jet,MT3jet)
		hist_MtMt3Jet_MET.Fill(tree.MT_AK8,tree.MET,MT3jet)
		prof_MTvsMT3jet_3jets.Fill(tree.MT_AK8,MT3jet,1)
	elif nJets >= 4: # do stuff here for events with 4+ jets
		nEventsWith4OrMoreJets += 1
		DeltaPhi3_AK8 = rt.TMath.Abs(jetCollection[2].Phi()-tree.METPhi)
		DeltaPhiMin_3jet = min(DeltaPhi3_AK8, tree.DeltaPhi2_AK8, tree.DeltaPhi1_AK8)
		DeltaPhi4_AK8 = rt.TMath.Abs(jetCollection[3].Phi()-tree.METPhi)
		DeltaPhiMin_4jet = min(DeltaPhi4_AK8, DeltaPhi3_AK8, tree.DeltaPhi2_AK8, tree.DeltaPhi1_AK8)
"""

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

	if not (deltaRjethvQuarkMinIndex == deltaRjetbarhvQuarkMinIndex):
		MT_HVQuarks_Jets = trans_mass_Njet([jetCollection[deltaRjethvQuarkMinIndex], jetCollection[deltaRjetbarhvQuarkMinIndex]],tree.MET,tree.METPhi)

	# fill 'All' histograms here
	if (tree.DeltaPhi1_AK8 == DeltaPhiMin_3jet):
	elif (tree.DeltaPhi2_AK8 == DeltaPhiMin_3jet):
	elif (DeltaPhi3_AK8 == DeltaPhiMin_3jet):
"""

print("Events with 2 jets: " + str(nEventsWithOnly2Jets+nEventsWithOnly3Jets+nEventsWith4OrMoreJets))
print("Events with only 2 jets: " + str(nEventsWithOnly2Jets))
print("Events with 3 jets: " + str(nEventsWithOnly3Jets+nEventsWith4OrMoreJets))
print("Events with only 3 jets: " + str(nEventsWithOnly3Jets))
print("Events with 4 jets: " + str(nEventsWith4OrMoreJets))


#step 4, normalize & color all histograms

for histo in histList:
	name = histo.GetName()
	if "All" in name:
		histo.SetLineColor(rt.kBlack)
	elif "2jets" in name:
		histo.SetLineColor(rt.kBlue)
	elif "3jets" in name:
		histo.SetLineColor(rt.kGreen)
	elif "4jets" in name:
		histo.SetLineColor(rt.kRed)
	else:
		print(name + " doesn't contain 'All', '2jets', '3jets', or '4jets'")
	area = histo.GetEntries()
	#area = 1
	try:
		histo.Scale(1/area)
	except ZeroDivisionError:
		print("Empty histogram: " + histo.GetName())

#5: save histogram plots
hist_MT3jet_3jets.SetLineColor(rt.kBlue)
#hist_MT4jet_4jets.SetLineColor(rt.kBlue)
drawHistos([hist_MT_All,hist_MT_2jets,hist_MT_3jets],"MT",True)
drawHistos([hist_MT_All,hist_MT3jet_3jets,hist_MT_3jets],"MT_3jets",True)

c1 = rt.TCanvas()
hist_MTvsMT3jet_3jets.Draw("colz")
c1.SaveAs(createImageName(mZprime, mDark, rinv, alpha, "Mtvs3jetMT_3jets"))
prof_MTvsMT3jet_3jets.Draw()
c1.SaveAs(createImageName(mZprime, mDark, rinv, alpha, "Mtvs3jetMT_3jets_profile"))

makeCornerPlot3D(hist_MVsEtaVsPhi,createImageName(mZprime, mDark, rinv, alpha, "3dcorner_example"))
makeCornerPlot2D(hist_MTvsMT3jet_3jets,createImageName(mZprime, mDark, rinv, alpha, "2dcorner_3jetsMTvsMT3Jet"))
makeCornerPlot3D(hist_MtMt3Jet_dPhiMin2Jet,createImageName(mZprime, mDark, rinv, alpha, "3dcorner_dPhiMin2jet"))
makeCornerPlot3D(hist_MtMt3Jet_dPhiMin3Jet,createImageName(mZprime, mDark, rinv, alpha, "3dcorner_dPhiMin3jet"))
makeCornerPlot3D(hist_MtMt3Jet_MET,createImageName(mZprime, mDark, rinv, alpha, "3dcorner_MET"))

plotTwo2D(hist_MtMt3Jet_MET.Project3D("xy"),hist_MtMt3Jet_MET.Project3D("zy"), createImageName(mZprime, mDark, rinv, alpha, "double2Dplot_MET"))
plotTwo2D(hist_MtMt3Jet_dPhiMin2Jet.Project3D("xy"),hist_MtMt3Jet_dPhiMin2Jet.Project3D("zy"), createImageName(mZprime, mDark, rinv, alpha, "double2Dplot_dPhiMin2Jet"))
plotTwo2D(hist_MtMt3Jet_dPhiMin3Jet.Project3D("xy"),hist_MtMt3Jet_dPhiMin3Jet.Project3D("zy"), createImageName(mZprime, mDark, rinv, alpha, "double2Dplot_dPhiMin3jet"))

inFile.Close()
