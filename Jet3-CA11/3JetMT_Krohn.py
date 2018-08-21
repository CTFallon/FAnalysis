import ROOT as rt
import sys
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
#Script to compute the MT resolution at several cuts along a variable
# should also include the MT distribution, along with (mean or peak, stdDev) for the all 12, all 123, 'perfect', and cuts.

def createInFileName(mZ, mD, rI, aD):
	return "PrivateSamples.SVJ_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_n-1000_0_RA2AnalysisTree.root"

def createImageName(mZ, mD, rI, aD, plotName):
	return "PS_SVJ_JetPermutations_mZprime-"+str(mZ)+"_mDark-"+mD+"_rinv-"+rI+"_alpha-"+aD+"_"+plotName+".png"

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
	if len(histo) == 1:
		histStack.SetTitle(histo[0].GetTitle())
	else:
		histStack.SetTitle(plotName)
	for hist in histo:
		hist.SetLineColor(histo.index(hist)+1)
		histStack.Add(hist)
	if logY == True:
		canvas.SetLogy()
	if not stack:
		histStack.Draw("NOSTACK")
	else:
		histStack.Draw("")
	histStack.GetXaxis().SetTitle(histo[0].GetXaxis().GetTitle())
	if len(histo) > 1:
		canvas.BuildLegend(0.75,0.75,0.9,0.9,"")
	canvas.SaveAs(createImageName(mZprime, mDark, rinv, alpha, plotName))



#step1, know what data set we're looking at

mZprime = 3000
mDark = '20'
rinv = '0p3'
alpha = '0p2' #['0p1', '0p2', '0p5', '1']
f_rvis = 1.0 - 0.3

includeAllEvents = True

inFile = rt.TFile.Open("root://cmsxrootd.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV14/"+createInFileName(mZprime, mDark,rinv,alpha))
tree = inFile.Get("TreeMaker2/PreSelection")

nEvents = tree.GetEntries()
print("Total Events: " + str(nEvents))
nEventsWith2OrMoreJets = 0
nEventsWithLeadingTwoJetsHavingPtGreaterThan170 = 0
nEventsWithMETMTRatioGreaterThanp15 = 0
nEventsPassedPreSelection = 0


# initilize histograms
histList = []

y = rt.TH1F("rapidity","Absolute Jet Rapidity;y;count//a.u.",100,0,4)
rat_mpt = rt.TH1F("deltaVar","m/p_T;m/p_T;count//a.u.",100,0,1)
rat_pt = rt.TH1F("ptRatio","pT Ratio;max(pt)/min(pt);count//a.u.",100, 0.9, 9)
abs_y = rt.TH1F("deltaRapidity","deltaRapidity;|yi-yj|;count//a.u.",100,-0.05,3.5)
rat_mptmpt = rt.TH1F("deltaVarRatio","\Delta_i//\Delta_j;\Delta_i//\Delta_j;count//a.u.",100,0.95,6)

y_Passed = rt.TH1F("rapidity_passed","Passed;y;count//a.u.",100,0,4)
rat_mpt_Passed = rt.TH1F("deltaVar_passed","Passed;m/p_T;count//a.u.",100,0,1)
rat_pt_Passed = rt.TH1F("ptRatio_passed","Passed;max(pt)/min(pt);count//a.u.",100, 0.9, 9)
abs_y_Passed = rt.TH1F("deltaRapidity_passed","Passed;|yi-yj|;count//a.u.",100,-0.05,3.5)
rat_mptmpt_Passed = rt.TH1F("deltaVarRatio_passed","Passed;\Delta_i//\Delta_j;count//a.u.",100,0.95,6)

y_Failed = rt.TH1F("rapidity_Failed","Failed;y;count//a.u.",100,0,4)
rat_mpt_Failed = rt.TH1F("deltaVar_Failed","Failed;m/p_T;count//a.u.",100,0,1)
rat_pt_Failed = rt.TH1F("ptRatio_Failed","Failed;max(pt)/min(pt);count//a.u.",100, 0.9, 9)
abs_y_Failed = rt.TH1F("deltaRapidity_Failed","Failed;|yi-yj|;count//a.u.",100,-0.05,3.5)
rat_mptmpt_Failed = rt.TH1F("deltaVarRatio_Failed","Failed;\Delta_i//\Delta_j;count//a.u.",100,0.95,6)

num_dJets = rt.TH1F("num_dJets","Number of Distinguished Jets;N(dJets);count//a.u.",8,0,8)
dJetsIdx = rt.TH1F("dJetsIdx","Index of Distinguished Jets;Index;count//a.u.",8,0,8)

nISRTaggedISR = 0
nISRTaggedFSR = 0
nFSRTaggedISR = 0
nFSRTaggedFSR = 0
nTotalJets = 0

for iEvent in range(nEvents):
	tree.GetEvent(iEvent)
	if iEvent%1000==0:
		print("EventNumber:"+str(iEvent))

	jetCollection = tree.JetsAK8
	nJets = len(jetCollection)
	#some calculations to determine which histos to fill
	# PreSelection Cuts
	if not (nJets==3): #moreThan2Jets, temp set to ==3 for 3JetMt classification
		continue
	nEventsWith2OrMoreJets += 1
	if not ((jetCollection[0].Pt() > 170.0) and (jetCollection[1].Pt() > 170.0)): #leadingJetsPtOver170
		continue
	nEventsWithLeadingTwoJetsHavingPtGreaterThan170 += 1
	if not (tree.MET/tree.MT_AK8 > 0.15): #METoverMTgreaterThan0p15
		continue
	nEventsWithMETMTRatioGreaterThanp15 += 1
	if not ((len(tree.Electrons) + len(tree.Muons)) == 0): #leptonNotPresent
		continue
	nEventsPassedPreSelection += 1
	nTotalJets += min(nJets,5)

	# 'truth' level ISR jets

	numberOfDaughtersAParticleHas = [0 for x in range(len(tree.GenParticles))]
	for iPart in range(len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if iParent != -1: 
			numberOfDaughtersAParticleHas[iParent] += 1
	AK8jetsWithHVDecendants = ["0","0","0","0","0"]
	AK8_nHVParts = [0,0,0,0,0]
	AK8_nParts = [0,0,0,0,0]
	AK8_ptHVParts = [0,0,0,0,0]
	AK8_ptParts = [0,0,0,0,0]
	#make vector of length nGenParts that is 1 if the particle came from a HV quark
	isFromHVQuark = [0 for x in range(len(tree.GenParticles))]
	listOfHVQuarks = []
	for iPart in range(len(tree.GenParticles)):
		iParent = tree.GenParticles_ParentIdx[iPart]
		if abs(tree.GenParticles_PdgId[iPart]) == 4900101:
			listOfHVQuarks.append(tree.GenParticles[iPart])
		if iParent >= iPart:
			print("Ut-oh, the parent has a higher index than the child...")
		if (abs(tree.GenParticles_PdgId[iParent]) == 4900101) or (isFromHVQuark[iParent]):
			isFromHVQuark[iPart] = 1
		#finding what Jet a particle is in, only if it decends from a HVQuark:
		for iJet in range(min(len(tree.JetsAK8),5)):
			if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
				AK8_ptParts[-iJet-1] += tree.GenParticles[iPart].Pt()
				AK8_nParts[-iJet-1] += 1
		if isFromHVQuark[iPart] and abs(tree.GenParticles_PdgId[iPart]) < 4900000:
			for iJet in range(min(len(tree.JetsAK8),5)):
				if tree.JetsAK8[iJet].DeltaR(tree.GenParticles[iPart]) < 0.8 and not numberOfDaughtersAParticleHas[iPart]:
					AK8jetsWithHVDecendants[-iJet-1] = "1"
					AK8_ptHVParts[-iJet-1] += tree.GenParticles[iPart].Pt()
					AK8_nHVParts[-iJet-1] += 1


	AK8ptFrac = [0,0,0,0,0]
	jetPassesptFracCut = ["0","0","0","0","0"]
	ptFracCut = [0.0,0.0,0.22,0.0,0.0]
	for iJet in range(len(AK8jetsWithHVDecendants)):
		if AK8_ptParts[-iJet-1] != 0:
			AK8ptFrac[-iJet-1] = AK8_ptHVParts[-iJet-1]/AK8_ptParts[-iJet-1]
		if AK8ptFrac[-iJet-1] > ptFracCut[-iJet-1]:
			jetPassesptFracCut[-iJet-1] = "1"
	jetCode = int(AK8jetsWithHVDecendants[0]+AK8jetsWithHVDecendants[1]+AK8jetsWithHVDecendants[2]+AK8jetsWithHVDecendants[3]+AK8jetsWithHVDecendants[4],2)
	jetCodePassesptFracCut = int(
				str(int(AK8jetsWithHVDecendants[0])*int(jetPassesptFracCut[0])) + 
				str(int(AK8jetsWithHVDecendants[1])*int(jetPassesptFracCut[1])) + 
				str(int(AK8jetsWithHVDecendants[2])*int(jetPassesptFracCut[2])) + 
				str(int(AK8jetsWithHVDecendants[3])*int(jetPassesptFracCut[3])) + 
				str(int(AK8jetsWithHVDecendants[4])*int(jetPassesptFracCut[4])) 
				,2)


	# create and fill jet varaibles
	y_i        = [abs(jetCollection[iJet].Rapidity()) for iJet in range(nJets)]
	rat_mpt_i  = [jetCollection[iJet].M()/jetCollection[iJet].Pt() for iJet in range(nJets)]
	rat_pt_ij  = [[] for iJet in range(nJets)]
	abs_y_ij   = [[] for iJet in range(nJets)]
	rat_mpt_ij = [[] for iJet in range(nJets)]
	for iJet in range(nJets):
		for jJet in range(nJets):
			if iJet != jJet:
				rat_pt_ij[iJet].append(max(jetCollection[iJet].Pt(),jetCollection[jJet].Pt())/min(jetCollection[iJet].Pt(),jetCollection[jJet].Pt()))
				abs_y_ij[iJet].append(abs(y_i[iJet]-y_i[jJet]))
				rat_mpt_ij[iJet].append(max(rat_mpt_i[iJet],rat_mpt_i[jJet])/min(rat_mpt_i[iJet],rat_mpt_i[jJet]))
	# fill jet variable histograms
	for iJet in range(nJets):
		y.Fill(y_i[iJet])
		rat_mpt.Fill(rat_mpt_i[iJet])
		rat_pt.Fill(min(rat_pt_ij[iJet]))
		abs_y.Fill(min(abs_y_ij[iJet]))
		rat_mptmpt.Fill(min(rat_mpt_ij[iJet]))
	#impliment Krohn algorithm

	jetDistinguished = [[0,0,0] for iJet in range(nJets)]
	distinguishedJetsIdx = []
	for iJet in range(nJets): #NOT passing means the jet is presumed to be FSR
		# if a jet Passes at least one test, we consider it "distinguished" (in the terms of Krohn)
		jetPassedptTest = 0
		jetPassedyTest = 0
		jetPasseddeltaTest = 0
		if min(rat_pt_ij[iJet]) > 2:
			jetPassedptTest = 1
		if min(abs_y_ij[iJet]) > 1.5:
			jetPassedyTest = 1
		if min(rat_mpt_ij[iJet]) > 1.5:
			jetPasseddeltaTest = 1
		jetDistinguished[iJet] = [jetPassedptTest, jetPassedyTest, jetPasseddeltaTest]
		if jetPassedptTest + jetPassedyTest + jetPasseddeltaTest >= 1:
			distinguishedJetsIdx.append(iJet)

	for iJet in range(min(nJets,5)):
		if not int(AK8jetsWithHVDecendants[-iJet-1]):
			y_Passed.Fill(y_i[iJet])
			rat_mpt_Passed.Fill(rat_mpt_i[iJet])
			rat_pt_Passed.Fill(min(rat_pt_ij[iJet]))
			abs_y_Passed.Fill(min(abs_y_ij[iJet]))
			rat_mptmpt_Passed.Fill(min(rat_mpt_ij[iJet]))
		else:
			y_Failed.Fill(y_i[iJet])
			rat_mpt_Failed.Fill(rat_mpt_i[iJet])
			rat_pt_Failed.Fill(min(rat_pt_ij[iJet]))
			abs_y_Failed.Fill(min(abs_y_ij[iJet]))
			rat_mptmpt_Failed.Fill(min(rat_mpt_ij[iJet]))

	num_dJets.Fill(len(distinguishedJetsIdx))
	for iJet in distinguishedJetsIdx:
		dJetsIdx.Fill(iJet)

	for iJet in range(min(nJets,5)):
		if int(AK8jetsWithHVDecendants[-iJet-1])*int(jetPassesptFracCut[-iJet-1]) == 0:
			trueISR = 1
		else:
			trueISR = 0
		if iJet in distinguishedJetsIdx:
			taggedISR = 1
		else:
			taggedISR = 0
		if trueISR and taggedISR:
			nISRTaggedISR += 1
		elif trueISR and not taggedISR:
			nISRTaggedFSR += 1
		elif not trueISR and taggedISR:
			nFSRTaggedISR += 1
		elif not trueISR and not taggedISR:
			nFSRTaggedFSR += 1
		else:
			print("Jet is funky...")

print("ISR Tagged ISR : " + str(nISRTaggedISR))
print("ISR Tagged FSR : " + str(nISRTaggedFSR))
print("FSR Tagged FSR : " + str(nFSRTaggedFSR))
print("FSR Tagged ISR : " + str(nFSRTaggedISR))
print("Total Jets     : " + str(nTotalJets))

drawHistos([y,rat_mpt,rat_pt,abs_y,rat_mptmpt],"test")
drawHistos([y,y_Passed,y_Failed],"Rapidity")
drawHistos([rat_mpt,rat_mpt_Passed,rat_mpt_Failed],"Ratio_MiPti")
drawHistos([rat_pt,rat_pt_Passed,rat_pt_Failed],"Ratio_PtiPtj")
drawHistos([abs_y,abs_y_Passed,abs_y_Failed],"DeltaRapidity")
drawHistos([rat_mptmpt,rat_mptmpt_Passed,rat_mptmpt_Failed],"Ratio_DeltaiDeltaj")
drawHistos([num_dJets],"numDJets",True)
drawHistos([dJetsIdx],"dJetsIdx",True)
		
		

	
inFile.Close()
