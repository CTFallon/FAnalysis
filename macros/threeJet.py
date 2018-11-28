from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array
from random import randint

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	# adding friend tree
	ff = rt.TFile.Open(self.extraDir+"outFriend.root")
	friendTree = ff.Get("friend")
	tree.AddFriend(friendTree)
	# added friend tree
	nEvents = tree.GetEntries()
	print("n events = " + str(nEvents))
	
	# Turn off all branches, selective turn on branches
	tree.SetBranchStatus("*", 0)
	tree.SetBranchStatus("RunNum",1)
	tree.SetBranchStatus("EvtNum",1)
	tree.SetBranchStatus("*AK8*", 1)
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)

	# branches from friend
	tree.SetBranchStatus("passedPreSelection",1)
	#tree.SetBranchStatus("numGenParts",1)
	#tree.SetBranchStatus("genParticleInAK8Jet",1)
	tree.SetBranchStatus("genParticleIsFromHVQuark",1)
	#tree.SetBranchStatus("numberOfDaughtersAParticleHas",1)
	#tree.SetBranchStatus("fracPtFromHVQuarks",1)
	tree.SetBranchStatus("numHVPartsInJet",1)
	tree.SetBranchStatus("numSMPartsInJet",1)
	tree.SetBranchStatus("iJetMaxDeltaPhi",1)
	tree.SetBranchStatus("pTMaxDeltaPhi",1)
	tree.SetBranchStatus("dPhiMaxDeltaPhi",1)
	#tree.SetBranchStatus("zPrimept",1)
	#tree.SetBranchStatus("zPrimephi",1)
	tree.SetBranchStatus("pGJ_visible",1)
	tree.SetBranchStatus("pGJ_invis",1)
	tree.SetBranchStatus("pGJ_every",1)
	tree.SetBranchStatus("fracVisPTfromVisHVQ",1)
	tree.SetBranchStatus("fracInvPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromAllHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVisHVQ",1)
	tree.SetBranchStatus("fracTotPTfromInvHVQ",1)
	tree.SetBranchStatus("fracTotPTfromVis",1)
	tree.SetBranchStatus("fracVisHVQtoInvHVQ",1)

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	
	hist_dPhi = self.makeTH1F("hist_dPhi",
		"DeltaPhi Between all Jets and MET",
		100,0,rt.TMath.Pi())
	hist_dPhi_FSR = self.makeTH1F("hist_dPhi_FSR",
		"DeltaPhi Between FSR Jets and MET",
		100,0,rt.TMath.Pi())
	hist_dPhi_ISR = self.makeTH1F("hist_dPhi_ISR",
		"DeltaPhi Between ISR Jets and MET",
		100,0,rt.TMath.Pi())
	
	hist_dR = self.makeTH1F("hist_dR",
		"Delta R between Jet(MaxDPhi) and other Jets;deltaR;",
		100,0,6)
	hist_dR_FSR = self.makeTH1F("hist_dR_FSR",
		"Delta R between Jet(MaxDPhi) and FSR Jets;deltaR;",
		100,0,6)
	hist_dR_ISR = self.makeTH1F("hist_dR_ISR",
		"Delta R between Jet(MaxDPhi) and ISR Jets;deltaR;",
		100,0,6)


	hist_avgDR = self.makeTH1F("hist_avgDR",
		"Average DeltaR between JetMaxDPhi and other jets",
		100,0,6)
	hist_avgDR_2FSR = self.makeTH1F("hist_avgDR_2FSR",
		"Average DeltaR between JetMaxDPhi and other jets (2 FSR in event)",
		100,0,6)
	hist_avgDR_3FSR = self.makeTH1F("hist_avgDR_3FSR",
		"Average DeltaR between JetMaxDPhi and other jets (3 FSR in event)",
		100,0,6)

	hist_tThrust = self.makeTH1F("hist_tThrust",
		"Transverse Thrust",100,0.5,1.1)
	hist_tThrust_2FSR = self.makeTH1F("hist_tThrust_2FSR",
		"Transverse Thrust for evnets with 2FSR",100,0.5,1.1)
	hist_tThrust_3FSR = self.makeTH1F("hist_tThrust_3FSR",
		"Transverse Thrust for events with 3 FSR",100,0.5,1.1)

	hist_tThrustTheta = self.makeTH1F("hist_tThrustTheta",
		"Trancverse Thrust Angle", 100, -rt.TMath.Pi(),3*rt.TMath.Pi())
	hist_tThrustTheta_2FSR = self.makeTH1F("hist_tThrustTheta_2FSR",
		"Trancverse Thrust Angle for evnets with 2FSR", 100, -rt.TMath.Pi(),3*rt.TMath.Pi())
	hist_tThrustTheta_3FSR = self.makeTH1F("hist_tThrustTheta_3FSR",
		"Trancverse Thrust Anglefor events with 3 FSR", 100, -rt.TMath.Pi(),3*rt.TMath.Pi())
	ctr=0
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if tree.passedPreSelection != 1: # skip events that failed preSelection
			continue
		# the above means that each event has at least 2 AK8 jets above 170 GeV
		# and has zero leptons, with MET/MT above 0.15
		nJets = len(tree.JetsAK8)
		if nJets == 2: 
			continue
		nFSRJets = 0
		for i in tree.JetsAK8_isISR:
			if not i:
				nFSRJets += 1
		#if nFSRJets != 3:
		#	continue
		# Find the jet that is furthest the MET
		ctr+=1
		dPhiMetList = []
		ISRJetsIndex = []
		FSRJetsIndex = []
		for i in range(nJets):
			dPhi_temp = deltaPhi(tree.JetsAK8[i].Phi(),tree.METPhi)
			dPhiMetList.append(dPhi_temp)
			hist_dPhi.Fill(dPhi_temp)
			if tree.JetsAK8_isISR[i]:
				ISRJetsIndex.append(i)
				hist_dPhi_ISR.Fill(dPhi_temp)
			else:
				FSRJetsIndex.append(i)
				hist_dPhi_FSR.Fill(dPhi_temp)
		maxDPhi = min(dPhiMetList)
		maxDPhiIndex = dPhiMetList.index(maxDPhi)
		# Now that we know the jet furthest from the met, we need to know
		# if the jet closest to it is FSR or ISR
		drTotal = []
		nOtherJets = 0.
		for i in range(nJets):
			if i == maxDPhiIndex:
				continue
			nOtherJets += 1
			dR_temp = tree.JetsAK8[maxDPhiIndex].DeltaR(tree.JetsAK8[i])
			drTotal.append(dR_temp)
			hist_dR.Fill(dR_temp)
			if i in FSRJetsIndex:
				hist_dR_FSR.Fill(dR_temp)
			else:
				hist_dR_ISR.Fill(dR_temp)

		drTotal = min(drTotal)
		hist_avgDR.Fill(drTotal)
		tt, ttt = transverseThrust([tree.JetsAK8[0],tree.JetsAK8[1],tree.JetsAK8[2]])
		hist_tThrust.Fill(tt)
		hist_tThrustTheta.Fill(ttt)
		if nFSRJets == 2:
			hist_avgDR_2FSR.Fill(drTotal)
			hist_tThrust_2FSR.Fill(tt)
			hist_tThrustTheta_2FSR.Fill(ttt)
		if nFSRJets == 3:
			hist_avgDR_3FSR.Fill(drTotal)
			hist_tThrust_3FSR.Fill(tt)
			hist_tThrustTheta_3FSR.Fill(ttt)
		if ctr%1000 == 0:
			sphericity([tree.JetsAK8[0],tree.JetsAK8[1]])
		#	drawEventEtaPhiPlot(tree.JetsAK8,
		#			tree.JetsAK8_isISR,
		#			tree.GenParticles,
		#			tree.GenParticles_PdgId,
		#			tree.METPhi,
		#			tree.genParticleIsFromHVQuark,
		#			iEvent,
		#			self.extraDir+"/etaPhi/")


	makePlots([hist_dR,hist_dR_FSR,hist_dR_ISR],
		self.extraDir,"DeltaR")
	makePlots([hist_dPhi,hist_dPhi_FSR,hist_dPhi_ISR],
		self.extraDir,"DeltaPhiMet")
	makePlots([hist_avgDR,hist_avgDR_2FSR,hist_avgDR_3FSR],
		self.extraDir,"AverageDeltaR",log=True)
	makePlots([hist_tThrust,hist_tThrust_2FSR,hist_tThrust_3FSR],
		self.extraDir,"TransverseThrust")
	makePlots([hist_tThrustTheta,hist_tThrustTheta_2FSR,hist_tThrustTheta_3FSR],
		self.extraDir,"TransverseThrustTheta")

	

def addLoop():
	baseClass.loop = loop

def deltaPhi(phi1, phi2):
	x = phi1 - phi2
	while x >= rt.TMath.Pi():
		x = x - 2*rt.TMath.Pi()
	while x < -rt.TMath.Pi():
		x = x + 2*rt.TMath.Pi()
	return abs(x)

def makePlots(hList,direct, name, log = False, lC = [0.7,0.7,1.0,1.0]):
	c1 = rt.TCanvas("c1","c1",900,600)
	if log:
		c1.SetLogy()
	for iH in range(len(hList)):
		if iH == 0:
			hList[iH].SetLineColor(1)
			hList[iH].Draw()
		else:
			hList[iH].SetLineColor(iH+1)
			hList[iH].Draw('same')
	c1.BuildLegend(lC[0],lC[1],lC[2],lC[3])
	c1.SaveAs(direct+name+".png")

def drawEventEtaPhiPlot(jetCollectionAK8, jetISR, partCol, particlePDGID,METPhi, isFromHVQuark, plotNumber, edir):
	canv = rt.TCanvas("canv","canv",1600,800)
	histAxis = rt.TH2F("axisHsito", ";\eta;\phi",100,-6,6,100,-rt.TMath.Pi(),rt.TMath.Pi())
	histAxis.SetStats(False)
	histAxis.Draw()
	objectList = []
	for iJet in range(len(jetCollectionAK8)):
		objectList.append(rt.TEllipse(jetCollectionAK8[iJet].Eta(), jetCollectionAK8[iJet].Phi(), 0.8, 0.8))
		objectList[-1].SetLineColor(iJet+1)
		objectList[-1].SetLineWidth(2)
		if jetISR[iJet]:
			objectList[-1].SetLineStyle(2)
	for iPart in range(len(partCol)):
		objectList.append(rt.TMarker(partCol[iPart].Eta(),partCol[iPart].Phi(),2))
		objectList[-1].SetMarkerSize(2)
		if abs(particlePDGID[iPart]) == 4900101: # HV Quarks are
			objectList[-1].SetMarkerColor(3) # Green
			objectList[-1].SetMarkerStyle(20) # Solid Triangles
			objectList[-1].SetMarkerSize(3) # Large
		elif abs(particlePDGID[iPart] == 4900023): # Z' is
			objectList[-1].SetMarkerStyle(34) #  thick cross
			objectList[-1].SetMarkerSize(3) # large
		elif abs(particlePDGID[iPart]) > 4900000: # HV Particles are
			objectList[-1].SetMarkerColor(2) # Red
		if isFromHVQuark[iPart]: # decendants of HV quarks are 
			objectList[-1].SetMarkerStyle(5) # X's
	objectList.append(rt.TLine(-6,METPhi,6,METPhi))
	objectList[-1].SetLineStyle(2)
	objectList[-1].SetLineColor(4)
	objectList[-1].SetLineWidth(2)
	for thing in objectList:
		thing.Draw()
	canv.SaveAs(edir+"etaPhi_"+str(plotNumber)+".png")



def transverseThrust(jets):
	tThrust = rt.TF1("tThrust",
	"(abs(cos(x)*[0]+sin(x)*[1])+abs(cos(x)*[2]+sin(x)*[3])+abs(cos(x)*[4]+sin(x)*[5]))/(sqrt(pow([0],2)+pow([1],2))+sqrt(pow([2],2)+pow([3],2))+sqrt(pow([4],2)+pow([5],2)))",0,2*rt.TMath.Pi())
	tThrust.SetParameters(jets[0].Px(),jets[0].Py(),
							jets[1].Px(),jets[1].Py(),
							jets[2].Px(),jets[2].Py(),)
	return tThrust.GetMaximum(), tThrust.GetMaximumX()

def sphericity(jets): # jets is a list of n TLorentzVectors
	#S_T == 2(l_2)/(l_2+l_1)
	# where l_i is the ith eigenvalue of the following matrix:
	# 1/(sum[pT_i]) * sum[(1/pT_i)[[MATRIX]]]
	# where MATRIX is the 2x2 matrix of 
	# [   px_i^2  px_i*py_i ]
	# [ py_i*px_i    py_i^2 ]  
	# tends towards 0 for pencil-like events
	# tends towards 1 for isotropic events
	norminv = 0.
	mat1 = 0.
	mat2 = 0.
	mat3 = 0.
	for jet in jets:
		norminv += jet.Pt()
		mat1 += jet.Px()**2/jet.Pt()
		mat2 += jet.Px()*jet.Py()/jet.Pt()
		mat3 += jet.Py()**2/jet.Pt()
	sMat = rt.TMatrixD(2,2)
	sMat[0][0] = mat1/norminv
	sMat[1][0] = mat2/norminv
	sMat[0][1] = mat2/norminv
	sMat[1][1] = mat3/norminv
	eigSys = rt.TMatrixDEigen(sMat)
	eigVal = eigSys.GetEigenValuesRe()
	print(eigVal[0],eigVal[1],2*eigVal[1]/(eigVal[0]+eigVal[1]))
	return 2*eigVal[1]/(eigVal[0]+eigVal[1])

	




























