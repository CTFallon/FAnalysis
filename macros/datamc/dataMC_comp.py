from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

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
	print("n events = " + str(nEvents))


	# initalize histograms to be made, or create Friend tree to be filled
	self.outRootFile.cd()

	# list of branches to plot
	plotDict = {#key = var name, value = [nBins, binLow, binHigh, title]
	'JetsAK8[0].Pt()':[100,0,6000,"Leading Jet Pt;Jet Pt;Events"],
	'MET':[100,0,6000,"MET; MET; Events"],
	'JetsAK8_girth[0]':[100,0,1,"Leading Jet AK8 Girth; Girth; Events"],
	}
	branchList = tree.GetListOfBranches()
	branchListNames = []
	for thing in branchList:
		branchListNames.append(thing.GetName())
	tree.SetBranchStatus("*",0)
	tree.SetBranchStatus("RunNum" ,1)
	tree.SetBranchStatus("HEMOptVetoFilter" ,1)
	tree.SetBranchStatus("madHT",1)
	tree.SetBranchStatus("GenElectrons",1)
	tree.SetBranchStatus("GenMuons",1)
	tree.SetBranchStatus("GenTaus",1)
	tree.SetBranchStatus("GenMET",1)
	for plotVar in plotDict.keys():
		if plotVar.split("[")[0] in branchListNames:
			tree.SetBranchStatus(plotVar.split("[")[0],1)
		else:
			print(plotVar.split("[")[0] +" not in tree")

	if (("Jets" in self.fileID) or ("QCD" in self.fileID)): # only need to do this for MC bkg
		tree.SetBranchStatus("Weight",1)
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
	histDict = {}

	for plotVar, histSpecs in plotDict.items():
		histDict[plotVar] = self.makeTH1F(plotVar+"_"+self.fileID,histSpecs[3],histSpecs[0],histSpecs[1],histSpecs[2]) 
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET):
				continue		
		
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)): # Bkg MC get tree weight, data and signal MC get weight == 1
			weight = tree.Weight
		else: 
			weight = 1.
		
		if "Data18" in self.fileID:
			print(tree.RunNum, tree.HEMOptVetoFilter)
			if tree.RunNum >= 319077:
				if tree.HEMOptVetoFilter == 0:
					continue

		for plotVar in plotDict.keys():
			plotVarList = plotVar.replace("]","[").split("[")
			if '' in plotVarList: plotVarList.remove('')
			if len(plotVarList) == 1:
				histDict[plotVar].Fill(getattr(tree,plotVar),weight*lumi)
			elif len(plotVarList) == 2:
				histDict[plotVar].Fill(getattr(tree,plotVarList[0])[int(plotVarList[1])],weight*lumi)
			elif len(plotVarList) == 3:	
				histDict[plotVar].Fill(getattr(getattr(tree,plotVarList[0])[int(plotVarList[1])],plotVarList[2][1:-2])(),weight*lumi)
				# getattr is funky for methods. for a public variable of a class, getattr(obj, attr) works.
				# for a public function, needs getattr(obj,func)(args of func)
				# since JetsAK8[0].Pt() is a public function, the proper way to get that value using getattr is
				# getattr(JetsAK8[0],Pt)(), which we do here, but JetsAK8[0] is replaced by getattr(tree,JetsAK8)[0]

def addLoop():
	baseClass.loop = loop


