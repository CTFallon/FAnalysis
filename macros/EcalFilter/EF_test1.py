from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array
from EventListFilter import *

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

def loop(self):
	# set up trees or chains
	#f = rt.TFile.Open(self.inputFileList[0])
	#tree = f.Get(self.treeNameList[0])
	tree = self.getChain(self.treeNameList[0])
	# added friend tree
	nEvents = tree.GetEntries()
	

	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	print("n events = " + str(nEvents))

	hist_LeadPhi_allEvents = self.makeTH1F("hist_LeadPhi_allEvents","Lead AK8 Jet Phi of all Events;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_SubPhi_allEvents = self.makeTH1F("hist_SubPhi_allEvents","Subead AK8 Jet Phi of all Events;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_LeadPhi_passBoundry = self.makeTH1F("hist_LeadPhi_passBoundry","Lead AK8 Jet Phi of Events that pass Boundry;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_SubPhi_passBoundry = self.makeTH1F("hist_SubPhi_passBoundry","Subead AK8 Jet Phi of Events that pass Boundry;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_LeadPhi_passTrigPrim = self.makeTH1F("hist_LeadPhi_passTrigPrim","Lead AK8 Jet Phi of Events that pass TrigPrim;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_SubPhi_passTrigPrim = self.makeTH1F("hist_SubPhi_passTrigPrim","Subead AK8 Jet Phi of Events that pass TrigPrim;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_LeadPhi_passBoth = self.makeTH1F("hist_LeadPhi_passBoth","Lead AK8 Jet Phi of Events that pass Both;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())
	hist_SubPhi_passBoth = self.makeTH1F("hist_SubPhi_passBoth","Subead AK8 Jet Phi of Events that pass Both;Phi;Events",50,-rt.TMath.Pi(),rt.TMath.Pi())

	hist_Ak4_subDistToMask = self.makeTH1F("hist_Ak4_subDistToMask","Minimum Distance of SubAK4 jet to Masked Cells",100,0,2)
	listOfMaskedPoints = [[-1.2,-1.19],[-0.912,2.03],[-0.912,3.01],[-0.816,-1.75],[-0.72,-2.17],[-0.72,-0.77],[-0.528,2.73],[-0.432,2.73],[-0.336,0.21],[-0.24,0.07],[-0.24,0.21],[-0.144,-2.59],[-0.144,0.77],[-0.048,0.91],[0.144,1.75],[0.912,1.75],[0.912,2.87],[1.008,0.63],[1.296,-0.49],[-1.584,0.63],[-0.816,1.47],[-0.72,-2.31],[-0.144,0.07],[-0.048,-2.59],[-0.048,0.77],[0.048,0.91],[1.104,-3.15],[1.488,2.73]]

	#create eventfilter list from Kevin's code:
	# need to know if Data or MC
	if "Data" in self.fileID:
		elf = EventListFilter(("EcalFilterLists/Run2016B_ver2-Nano1June2019_ver2-v2_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016C-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016D-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016E-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016F-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016G-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Run2016H-Nano1June2019-v1_JetHT_EcalDeadCellBoundaryEnergyFilterList.txt")) 
	elif "QCD" in self.fileID:
		elf = EventListFilterMC(("EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8_ext1_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8_EcalDeadCellBoundaryEnergyFilterList.txt","EcalFilterLists/Summer16v5-Nano1June2019_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8_ext2_EcalDeadCellBoundaryEnergyFilterList.txt",))
	else:
		print("Dun goofed!")
	print()

	
	# Turn off all branches, selective turn on branches
	tree.SetBranchStatus("*", 0)
	tree.SetBranchStatus("RunNum",1)
	tree.SetBranchStatus("LumiBlockNum",1)
	tree.SetBranchStatus("EvtNum",1)
	tree.SetBranchStatus("*AK8*", 1)
	tree.SetBranchStatus("Jets", 1)
	tree.SetBranchStatus("DeltaPhi*", 1)
	tree.SetBranchStatus("MET",1)
	tree.SetBranchStatus("METPhi",1)
	tree.SetBranchStatus("EcalDeadCellTriggerPrimitiveFilter",1)
	tree.SetBranchStatus("GenJets",1)
	#tree.SetBranchStatus("GenParticles*",1)
	#tree.SetBranchStatus("Electrons",1)
	#tree.SetBranchStatus("Muons",1)
	if (("Jets" in self.fileID) or ("QCD" in self.fileID)): # only need to do this for MC bkg
		if "16" in self.fileID:
			lumi = 35921.036
			print("2016 Lumi")
		elif "17" in self.fileID:
			lumi = 41521.331
			print("2017 Lumi")
		elif "18PRE" in self.fileID:
			lumi = 21071.460
			print("2018 pre Lumi")
		elif "18POST" in self.fileID:
			lumi = 38621.232
			print("2018 post Lumi")
		elif "18" in self.fileID:
			lumi = 59692.692
			print("2018 full lumi")
		else:
			print("Dont know what total lumi to use. default to 40 fb-1")
			lumi = 40000.
	else:
		lumi = 1

	nPassBoundry = 0
	nPassTrigPrim = 0
	nPassBoth = 0
	nTotal = 0.0
	#boundryStr = ""
	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()
	#run:lumi:event
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET, self.fileID):
				continue
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)):
			weight = tree.Weight*tree.puWeightNew
			if weight == 0.0:
				print("Weights: {} {}".format(tree.Weight, tree.puWeightNew))
		else:
			weight = 1.
		pass_Boundry = False
		pass_TrigPrim = False
		if "Data" in self.fileID:
			pass_Boundry = elf.CheckEvent(tree.RunNum, tree.LumiBlockNum, tree.EvtNum)
		elif "QCD" in self.fileID:
			pass_Boundry = elf.CheckEvent(tree.RunNum, tree.LumiBlockNum, tree.EvtNum, int(tree.GenJets[0].Pt()))
		#boundryStr += str(int(pass_Boundry))
	
		pass_TrigPrim = bool(tree.EcalDeadCellTriggerPrimitiveFilter)

		nPassBoundry += int(pass_Boundry)
		nPassTrigPrim += int(pass_TrigPrim)
		nPassBoth += int(pass_Boundry and pass_TrigPrim)
		nTotal += 1

	
		hist_LeadPhi_allEvents.Fill(tree.JetsAK8[0].Phi(),lumi*weight)
		hist_SubPhi_allEvents.Fill(tree.JetsAK8[1].Phi(),lumi*weight)
		if pass_Boundry:
			hist_LeadPhi_passBoundry.Fill(tree.JetsAK8[0].Phi(),lumi*weight)
			hist_SubPhi_passBoundry.Fill(tree.JetsAK8[1].Phi(),lumi*weight)
		if pass_TrigPrim:
			hist_LeadPhi_passTrigPrim.Fill(tree.JetsAK8[0].Phi(),lumi*weight)
			hist_SubPhi_passTrigPrim.Fill(tree.JetsAK8[1].Phi(),lumi*weight)
		if pass_Boundry and pass_TrigPrim:
			hist_LeadPhi_passBoth.Fill(tree.JetsAK8[0].Phi(),lumi*weight)
			hist_SubPhi_passBoth.Fill(tree.JetsAK8[1].Phi(),lumi*weight)

		minDist = 100
		for pair in listOfMaskedPoints:
			eta, phi = pair[0], pair[1]
			distFromJetToPair = rt.TMath.Sqrt((tree.Jets[1].Eta() - eta)**2 + (tree.Jets[1].Phi() - phi)**2)
			if distFromJetToPair < minDist:
				minDist = distFromJetToPair
		hist_Ak4_subDistToMask.Fill(minDist,lumi*weight)


	print("Total Events.......... {} {}".format(int(nTotal),nTotal/nTotal))
	print("Events Pass TrigPrim.. {} {}".format(nPassTrigPrim,nPassTrigPrim/nTotal))
	print("Events Pass Boundry... {} {}".format(nPassBoundry,nPassBoundry/nTotal))
	print("Events Pass Both...... {} {}".format(nPassBoth,nPassBoth/nTotal))
	#print(boundryStr)
	

def addLoop():
	baseClass.loop = loop


