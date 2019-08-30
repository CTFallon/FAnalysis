from analysisBase import baseClass
import ROOT as rt
import tdrstyle
from array import array

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(1)
tdrstyle.setTDRStyle()

# want two histograms for each set (leadjet, subjet, bothjets)
# both jetPt
# one of only jets that pass BDT WP (numer)
# one of all jets in each jetset (denom)

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

	nBins = 2800
	nBinsVarSize = 19
	#binEdges = array('f',[200.,570.,620.,650.,680.,700.,730.,750.,770.,810.,830.,870.,900.,930.,970.,1020.,1110.,1270.,3000.])
	binEdges = array('f',[200.,283.,326.,359.,388.,415.,444.,475.,511.,554.,602.,647.,685.,720.,752.,786.,822.,870.,941.,1079])

	# list of branches to plot
	plotDict = {#key = var name, value = [varType, nBins, binLow, binHigh, title]
				# varType can be "s" - single value (ie 'MET')
				#				 "sF" - single value, but a function of (not sure if this exists, but just in case)
				#				 "vA" - vector, all values (ie 'JetsAK8_girth')
				#				 "vI" - vector, only index value (ie 'JetsAK8_girth[0]')
				#				 "vAF" - vector, all values but function (ie 'JetsAK8.Pt()' - Pt of all ak8 Jets)
				#				 "vIF" - vector, indexed function (ie 'JetsAK8[0].Pt()', only Pt of leading AK8 Jet)
				#				 "vR", "vRF
	'JetsAK8[0].Pt()':["vIF",nBins,200,3000,self.fileID+"; Jet Pt, lead jet; Events"],
	'JetsAK8[1].Pt()':["vIF",nBins,200,3000,self.fileID+"; Jet Pt, sub jet; Events"],
	'JetsAK8.Pt()':["vAF",nBins,200,3000,self.fileID+"; Jet Pt, both jets; Events"]
	}

	histDict_numer = {}
	histDict_denom = {}
	
	# define the WP
	WP = 0.6
	for plotVar, histSpecs in plotDict.items():
		#histDict_numer[plotVar] = self.makeTH1F(plotVar+"_numer_"+self.fileID,histSpecs[4],histSpecs[1],histSpecs[2],histSpecs[3]) 
		#histDict_denom[plotVar] = self.makeTH1F(plotVar+"_denom_"+self.fileID,histSpecs[4],histSpecs[1],histSpecs[2],histSpecs[3])
		histDict_numer[plotVar] = self.makeTH1Fvarbins(plotVar+"_numer_"+self.fileID,histSpecs[4],nBinsVarSize,binEdges) 
		histDict_denom[plotVar] = self.makeTH1Fvarbins(plotVar+"_denom_"+self.fileID,histSpecs[4],nBinsVarSize,binEdges)
	
	listOfBranchesOneMightNeed = ["RunNum","HEMOptVetoFilter" ,"madHT","GenElectrons","GenMuons","GenTaus","GenMET","MT_AK8","JetsAK8_bdtSVJtag","Weight","puWeight"]
	listOfBranches = tree.GetListOfBranches()
	for branch in listOfBranches:
		if branch.GetName() in listOfBranchesOneMightNeed:
			branch.SetStatus(1)
		else:
			branch.SetStatus(0)

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

	for iEvent in range(nEvents):
		if iEvent%1000 == 0:
			print("Event: " + str(iEvent) + "/" + str(nEvents))
		tree.GetEvent(iEvent)
		if "TT" in self.fileID:
			if self.stitchTT(tree.GetFile().GetName().split("/")[-1], tree.madHT, len(tree.GenElectrons),len(tree.GenMuons), len(tree.GenTaus), tree.GenMET, self.fileID):
				continue
		if (("Jets" in self.fileID) or ("QCD" in self.fileID)):
			#print("Weights: {} {}".format(tree.Weight, tree.puWeight))
			weight = tree.Weight*tree.puWeight
			if weight == 0.0:
				print("Weights: {} {}".format(tree.Weight, tree.puWeight))
		else:
			weight = 1.
		if ("18" in self.fileID):
			# four cases:
			# Data18PRE only use events from runs < 319077
			# mc18PRE use eveything
			# Data18POST only use envets from runs >= 319077 and pass HEMO
			# mc18POST only use evnets that pass HEMO
			if "PRE" in self.fileID:
				if (("Data" in self.fileID) and (tree.RunNum >= 319077)):
					continue
				else:
					pass
			elif ("POST" in self.fileID):
				if tree.HEMOptVetoFilter == 0:
					continue
				else:
					if (("Data" in self.fileID) and (tree.RunNum < 319077)):
						continue
				
				# getattr is funky for methods. for a public variable of a class, getattr(obj, attr) works.
				# for a public function, needs getattr(obj,func)(args of func)
				# since JetsAK8[0].Pt() is a public function, the proper way to get that value using getattr is
				# getattr(JetsAK8[0],Pt)(), which we do here, but JetsAK8[0] is replaced by getattr(tree,JetsAK8)[0]
		
		for plotVar in plotDict.keys():
			if plotDict[plotVar][0] == "vAF":# branch.func()
				bName, bFunc = plotVar.split(".")[0],plotVar.split(".")[1][:-2]
				#print(bName, bFunc)
				for index, value in enumerate(getattr(tree,bName)):
					if index < 2:
						histDict_denom[plotVar].Fill(getattr(value, bFunc)(), weight*lumi)
						if tree.JetsAK8_bdtSVJtag[index] > WP:
							histDict_numer[plotVar].Fill(getattr(value, bFunc)(), weight*lumi)
			elif plotDict[plotVar][0] == "vIF":# branch[index].func()
				varStrHelper = plotVar.replace("]","[").replace("[",".").split(".")
				bName, bIndex, bFunc = varStrHelper[0],varStrHelper[1],varStrHelper[3][0:-2]
				histDict_denom[plotVar].Fill(getattr(getattr(tree,bName)[int(bIndex)],bFunc)(),weight*lumi)
				if tree.JetsAK8_bdtSVJtag[int(bIndex)] > WP:
					histDict_numer[plotVar].Fill(getattr(getattr(tree,bName)[int(bIndex)],bFunc)(),weight*lumi)



def addLoop():
	baseClass.loop = loop


