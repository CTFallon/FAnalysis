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

	nBins = 1000
	# list of branches to plot
	plotDict = {#key = var name, value = [varType, nBins, binLow, binHigh, title]
				# varType can be "s" - single value (ie 'MET')
				#				 "sF" - single value, but a function of (not sure if this exists, but just in case)
				#				 "vA" - vector, all values (ie 'JetsAK8_girth')
				#				 "vI" - vector, only index value (ie 'JetsAK8_girth[0]')
				#				 "vAF" - vector, all values but function (ie 'JetsAK8.Pt()' - Pt of all ak8 Jets)
				#				 "vIF" - vector, indexed function (ie 'JetsAK8[0].Pt()', only Pt of leading AK8 Jet)
				#				 "vR", "vRF
	'MET':["s",nBins,200,2000,self.fileID+";MET; Events"],
	'METPhi':["s",nBins,-rt.TMath.Pi(),rt.TMath.Pi(),self.fileID+";MET#Phi; Events"],
	'MHT':["s",nBins,200,2000,self.fileID+";MHT; Events"],
	'MT_AK8':["s",nBins,1500,6000,self.fileID+";MT;Events"],
	'MJJ_AK8':["s",nBins,0,5000,self.fileID+";Mjj;Events"],
	'JetsAK8[0].Pt()':["vIF",nBins,0,3000,self.fileID+";Jet Pt;Events"],
	'JetsAK8[1].Pt()':["vIF",nBins,0,1800,self.fileID+";Jet Pt;Events"],
	'JetsAK8[0].Eta()':["vIF",nBins,-2.4,2.4,self.fileID+";Jet Eta;Events"],
	'JetsAK8[1].Eta()':["vIF",nBins,-2.4,2.4,self.fileID+";Jet Eta;Events"],
	'JetsAK8[0].Phi()':["vIF",nBins,-rt.TMath.Pi(),rt.TMath.Pi(),self.fileID+";Jet Phi;Events"],
	'JetsAK8[1].Phi()':["vIF",nBins,-rt.TMath.Pi(),rt.TMath.Pi(),self.fileID+";Jet Phi;Events"],
	'JetsAK8_girth[0]':["vI",nBins,0,0.5,self.fileID+"; Girth; Events"],
	'JetsAK8_girth[1]':["vI",nBins,0,0.5,self.fileID+"; Girth; Events"],
	'JetsAK8_softDropMass[0]':["vI",nBins,0,600,self.fileID+"; SoftDrop Mass; Events"],
	'JetsAK8_softDropMass[1]':["vI",nBins,0,450,self.fileID+"; SoftDrop Mass; Events"],
	'JetsAK8_axismajor[0]':["vI",nBins,0,0.5,self.fileID+"; Major Axis; Events"],
	'JetsAK8_axismajor[1]':["vI",nBins,0,0.5,self.fileID+"; Major Axis; Events"],
	'JetsAK8_axisminor[0]':["vI",nBins,0,0.3,self.fileID+"; Minor Axis; Events"],
	'JetsAK8_axisminor[1]':["vI",nBins,0,0.3,self.fileID+"; Minor Axis; Events"],
	'JetsAK8_ptdrlog[0]':["vI",nBins,0,400,self.fileID+"; ptdrlog; Events"],
	'JetsAK8_ptdrlog[1]':["vI",nBins,0,400,self.fileID+"; ptdrlog; Events"],
	'JetsAK8_ptD[0]':["vI",nBins,0.1,1,self.fileID+"; ptD; Events"],
	'JetsAK8_ptD[1]':["vI",nBins,0.1,1,self.fileID+"; ptD; Events"],
	'JetsAK8_maxBvsAll[0]':["vI",nBins,0,1,self.fileID+"; maxBvsAll; Events"],
	'JetsAK8_maxBvsAll[1]':["vI",nBins,0,1,self.fileID+"; maxBvsAll; Events"],
	'JetsAK8_ecfN2b1[0]':["vI",nBins,0,0.5,self.fileID+"; ecfN2b1; Events"],
	'JetsAK8_ecfN2b1[1]':["vI",nBins,0,0.5,self.fileID+"; ecfN2b1; Events"],
	'JetsAK8_ecfN3b1[0]':["vI",nBins,0,4,self.fileID+"; ecfN3b1; Events"],
	'JetsAK8_ecfN3b1[1]':["vI",nBins,0,4,self.fileID+"; ecfN3b1; Events"],
	'JetsAK8_chargedHadronEnergyFraction[0]':["vI",nBins,0,1,self.fileID+"; fChgHad; Events"],
	'JetsAK8_chargedHadronEnergyFraction[1]':["vI",nBins,0,1,self.fileID+"; fChgHad; Events"],
	'JetsAK8_neutralHadronEnergyFraction[0]':["vI",nBins,0,1,self.fileID+"; fNeuHad; Events"],
	'JetsAK8_neutralHadronEnergyFraction[1]':["vI",nBins,0,1,self.fileID+"; fNeuHad; Events"],
	'JetsAK8_electronEnergyFraction[0]':["vI",nBins,0,1,self.fileID+"; fEle; Events"],
	'JetsAK8_electronEnergyFraction[1]':["vI",nBins,0,1,self.fileID+"; fEle; Events"],
	'JetsAK8_muonEnergyFraction[0]':["vI",nBins,0,1,self.fileID+"; fMu; Events"],
	'JetsAK8_muonEnergyFraction[1]':["vI",nBins,0,1,self.fileID+"; fMu; Events"],
	'JetsAK8_photonEnergyFraction[0]':["vI",nBins,0,1,self.fileID+"; f#gamma; Events"],
	'JetsAK8_photonEnergyFraction[1]':["vI",nBins,0,1,self.fileID+"; f#gamma; Events"],
	'deltaR12':["spec",nBins,0.8,3.5,self.fileID+";#Delta R(j_{1},j_{2});Events"],
	'metR':["spec",nBins,0.15,0.7,self.fileID+";MET/m_{T};Events"],
	'nJetsAK8':["spec",7,1,8,self.fileID+";nAK8 Jets;Events"],
	'nJetsAK4':["spec",40,0,40,self.fileID+";nAK4 Jets;Events"],
	'tau32_lead':["spec",nBins,0,1,self.fileID+";Leading Jet #tau_{32};Events"],
	'tau21_lead':["spec",nBins,0,1,self.fileID+";Subleading Jet #tau_{21};Events"],
	'tau32_sub':["spec",nBins,0,1,self.fileID+";Leading Jet #tau_{32};Events"],
	'tau21_sub':["spec",nBins,0,1,self.fileID+";Subeading Jet #tau_{21};Events"],
	'JetsAK8_bdtSVJtag[0]':["vI",nBins,0,1,self.fileID+"; SVJ BDT Output; Events"],
	'JetsAK8_bdtSVJtag[1]':["vI",nBins,0,1,self.fileID+"; SVJ BDT Output; Events"],
	'DeltaPhi1':["s",nBins,0,rt.TMath.Pi(),self.fileID+"; #Delta#phi(j_{1}, MET); Events"],
	'DeltaPhi2':["s",nBins,0,rt.TMath.Pi(),self.fileID+"; #Delta#phi(j_{2}, MET); Events"]
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
	tree.SetBranchStatus("JetsAK8_NsubjettinessTau1",1)
	tree.SetBranchStatus("JetsAK8_NsubjettinessTau2",1)
	tree.SetBranchStatus("JetsAK8_NsubjettinessTau3",1)
	tree.SetBranchStatus("MT_AK8",1)
	tree.SetBranchStatus("Jets",1)
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
	histDictFilter = {}
#	histDictFilter2 = {}

	for plotVar, histSpecs in plotDict.items():
		histDict[plotVar] = self.makeTH1F(plotVar+"_"+self.fileID,histSpecs[4],histSpecs[1],histSpecs[2],histSpecs[3]) 
#		histDictFilter[plotVar] = self.makeTH1F(plotVar+"_"+self.fileID+"_filter",histSpecs[4],histSpecs[1],histSpecs[2],histSpecs[3])
#		histDictFilter2[plotVar] = self.makeTH1F(plotVar+"_"+self.fileID+"_filter2",histSpecs[4],histSpecs[1],histSpecs[2],histSpecs[3])
	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET, self.fileID):
				continue		
		
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)): #or ("ST1" in self.fileID)): # Bkg MC get tree weight, data, SingleTop,  and signal MC get weight == 1
			weight = tree.Weight
		else: 
			weight = 1.
		
		if "18" in self.fileID:
			#print("2018 data")
			if (("PRE" in self.fileID) and ("Data" in self.fileID) and (tree.RunNum >= 319077)):
				continue
			elif "POST" in self.fileID:
				if ((("Data" in self.fileID) and (tree.RunNum < 319077)) or (tree.HEMOptVetoFilter == 0)):
					continue

		for plotVar in plotDict.keys():
			if plotDict[plotVar][0] == "s": # branch
				#print(plotVar)
				histDict[plotVar].Fill(getattr(tree,plotVar),weight*lumi)
#				histDictFilter[plotVar].Fill(getattr(tree,plotVar),weight*lumi*tree.ecalBadCalibReducedFilter)
#				histDictFilter2[plotVar].Fill(getattr(tree,plotVar),weight*lumi*tree.ecalBadCalibReducedExtraFilter)
			elif plotDict[plotVar][0] == "sF": # branch.func()
				varStrHelper = plotVar.split(".")
				bName, bFunc = varStrHelper[0],varStrHelper[1][:-2]
				#print(bname, bFunc)
				histDict[plotVar].Fill(getattr(getattr(tree, bName), bFunc)(), weight*lumi)
#				histDictFilter[plotVar].Fill(getattr(getattr(tree, bName), bFunc)(), weight*lumi*tree.ecalBadCalibReducedFilter)
#				histDictFilter2[plotVar].Fill(getattr(getattr(tree, bName), bFunc)(), weight*lumi*tree.ecalBadCalibReducedExtraFilter)
			elif plotDict[plotVar][0] == "vA": # branch
				#print(plotVar)
				for value in getattr(tree,plotVar):
					histDict[plotVar].Fill(value,weight*lumi)
#					histDictFilter[plotVar].Fill(value,weight*lumi*tree.ecalBadCalibReducedFilter)
#					histDictFilter2[plotVar].Fill(value,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
			elif plotDict[plotVar][0] == "vI": # branch[index]
				varStrHelper = plotVar.replace("]","[").split("[")
				bName, bIndex = varStrHelper[0],varStrHelper[1]
				#print(bName, bIndex)
				histDict[plotVar].Fill(getattr(tree,bName)[int(bIndex)],weight*lumi)
#				histDictFilter[plotVar].Fill(getattr(tree,bName)[int(bIndex)],weight*lumi*tree.ecalBadCalibReducedFilter)
#				histDictFilter2[plotVar].Fill(getattr(tree,bName)[int(bIndex)],weight*lumi*tree.ecalBadCalibReducedExtraFilter)
			elif plotDict[plotVar][0] == "vAF":# branch.func()
				bName, bFunc = plotVar.split(".")[0],plotVar.split(".")[1][:-2]
				#print(bName, bFunc)
				for value in getattr(tree,bName):
					histDict[plotVar].Fill(getattr(value, bFunc)(), weight*lumi)
#					histDictFilter[plotVar].Fill(getattr(value, bFunc)(), weight*lumi*tree.ecalBadCalibReducedFilter)
#					histDictFilter2[plotVar].Fill(getattr(value, bFunc)(), weight*lumi*tree.ecalBadCalibReducedExtraFilter)
			elif plotDict[plotVar][0] == "vIF":# branch[index].func()
				varStrHelper = plotVar.replace("]","[").replace("[",".").split(".")
				bName, bIndex, bFunc = varStrHelper[0],varStrHelper[1],varStrHelper[3][0:-2]
				histDict[plotVar].Fill(getattr(getattr(tree,bName)[int(bIndex)],bFunc)(),weight*lumi)
#				histDictFilter[plotVar].Fill(getattr(getattr(tree,bName)[int(bIndex)],bFunc)(),weight*lumi*tree.ecalBadCalibReducedFilter)
#				histDictFilter2[plotVar].Fill(getattr(getattr(tree,bName)[int(bIndex)],bFunc)(),weight*lumi*tree.ecalBadCalibReducedExtraFilter)
		#special ones that don't fit the normal stuff
		detlaR12 = tree.JetsAK8[0].DeltaR(tree.JetsAK8[1])
		metR = tree.MET/tree.MT_AK8
		nJetsAK8 = len(tree.JetsAK8)
		nJetsAK4 = len(tree.Jets)
		try:
			tau23_lead = tree.JetsAK8_NsubjettinessTau3[0]/tree.JetsAK8_NsubjettinessTau2[0]
		except ZeroDivisionError:
			tau23_lead = 0
		try:
			tau12_lead = tree.JetsAK8_NsubjettinessTau2[0]/tree.JetsAK8_NsubjettinessTau1[0]
		except ZeroDivisionError:
			tau12_lead = 0
		try:
			tau23_sub = tree.JetsAK8_NsubjettinessTau3[1]/tree.JetsAK8_NsubjettinessTau2[1]
		except ZeroDivisionError:
			tau23_sub = 0
		try:
			tau12_sub = tree.JetsAK8_NsubjettinessTau2[1]/tree.JetsAK8_NsubjettinessTau1[1]
		except ZeroDivisionError:
			tau12_sub = 0
		histDict['deltaR12'].Fill(detlaR12,weight*lumi)
		histDict['metR'].Fill(metR,weight*lumi)
		histDict['nJetsAK8'].Fill(nJetsAK8,weight*lumi)
		histDict['nJetsAK4'].Fill(nJetsAK4,weight*lumi)
		histDict['tau32_lead'].Fill(tau23_lead,weight*lumi)
		histDict['tau21_lead'].Fill(tau12_lead,weight*lumi)
		histDict['tau32_sub'].Fill(tau23_sub,weight*lumi)
		histDict['tau21_sub'].Fill(tau12_sub,weight*lumi)

#		histDictFilter['deltaR12'].Fill(detlaR12,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['metR'].Fill(metR,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['nJetsAK8'].Fill(nJetsAK8,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['nJetsAK4'].Fill(nJetsAK4,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['tau32_lead'].Fill(tau23_lead,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['tau21_lead'].Fill(tau12_lead,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['tau32_sub'].Fill(tau23_sub,weight*lumi*tree.ecalBadCalibReducedFilter)
#		histDictFilter['tau21_sub'].Fill(tau12_sub,weight*lumi*tree.ecalBadCalibReducedFilter)

#		histDictFilter2['deltaR12'].Fill(detlaR12,weight*lumi*tree.ecalBadCalibReducedExtraFilter*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['metR'].Fill(metR,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['nJetsAK8'].Fill(nJetsAK8,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['nJetsAK4'].Fill(nJetsAK4,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['tau32_lead'].Fill(tau23_lead,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['tau21_lead'].Fill(tau12_lead,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['tau32_sub'].Fill(tau23_sub,weight*lumi*tree.ecalBadCalibReducedExtraFilter)
#		histDictFilter2['tau21_sub'].Fill(tau12_sub,weight*lumi*tree.ecalBadCalibReducedExtraFilter)

				
			
				# getattr is funky for methods. for a public variable of a class, getattr(obj, attr) works.
				# for a public function, needs getattr(obj,func)(args of func)
				# since JetsAK8[0].Pt() is a public function, the proper way to get that value using getattr is
				# getattr(JetsAK8[0],Pt)(), which we do here, but JetsAK8[0] is replaced by getattr(tree,JetsAK8)[0]

def addLoop():
	baseClass.loop = loop


