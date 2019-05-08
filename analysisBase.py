# base class for FAnalysis, copied from Lucien Lo's PFG HcalTupleAnalyzer, adapted to python
import ROOT as rt
import datetime
class baseClass:
	def __init__(self, fileID, treeList, rootDir, outFileName):
		#self.fileList = fileList # .py file of dictoary
		self.fileID = fileID
		self.treeList = treeList # list of Tree names in each ROOT file
		self.outFileName = outFileName
		self.objects = []
		self.loadFileList()
		self.loadTreeList()
		self.extraDir = rootDir+'/'
		self.loadOutFile()

	def loadFileList(self):
		from input_conf.inputRoot import fileDict
		from input_conf.inputRoot import plotDict
		try:
			self.inputFileList = fileDict[self.fileID]
		except KeyError:
			print("Key does not exist in fileDict. Switching to Plotting Mode.")
			try:
				self.inputFileList = plotDict[self.fileID]
			except KeyError:
				exit("Key does not exist in plotDict. Exiting.")
			
	def loadTreeList(self):
		treeFile = open("input_conf/"+self.treeList)
		self.treeNameList = []
		for line in treeFile:
			self.treeNameList.append(line[:-1])#ignore new line character
		treeFile.close()

	def loadOutFile(self):
		self.outRootFile = rt.TFile.Open(
			self.extraDir+self.outFileName,
			'RECREATE'
			)
		
	def getChain(self, tree_name):
		chain = rt.TChain(tree_name)
		for fileName in self.inputFileList:
			chain.Add(fileName)
		return chain

	def stitchTT(self, eventFileName, madHT, GenMET, nLeptons):
		# return True means skip event
		#print(eventFileName)
		name = eventFileName.split("/")[-1].split("_")
		if "201" in name[2]:
			if not ((madHT < 600) and (nLeptons != 0)):
				return True
		elif "HT" in name[2]:
			if not (madHT >= 600):
				return True
		elif "201" in name[3]:
			if not (madHT < 600 and GenMET < 150):
				return True
		elif "201" in name[4]:
			if not (madHT < 600 and GenMET >= 150):
				return True
		else:
			print("TTBar stitiching error!")
			print(name)
			return True
		return False
	
	def makeTH1F(self, name, nbinsx, xmin, xmax):
		hist = rt.TH1F(name, name, nbinsx, xmin, xmax)
		self.objects.append(hist)
		return hist

	def makeTH1F(self, name, title, nbinsx, xmin, xmax):
		hist = rt.TH1F(name, title, nbinsx, xmin, xmax)
		self.objects.append(hist)
		return hist

	def makeTH2F(self, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax):
		hist = rt.TH2F(name, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
		self.objects.append(hist)
		return hist

	def makeTH2F(self, name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax):
		hist = rt.TH2F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
		self.objects.append(hist)
		return hist

	def makeTH3F(self, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax):
		hist = rt.TH3F(name, name, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax)
		self.objects.append(hist)
		return hist

	def makeTH3F(self, name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax):
		hist = rt.TH3F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax)
		self.objects.append(hist)
		return hist
	
	def makeTGraph(self, n, x, y):
		graph = rt.TGraph(n,x,y)
		self.objects.append(graph)
		return graph
	
	def makeTGraph(self):
		graph = rt.TGraph()
		self.objects.append(graph)
		return graph

	def makePng(self, LoHi, name, doLeg = True, log = False, doCum = False):
		c1 = rt.TCanvas("c1","c1",1200,900)
		if log:
			c1.SetLogy()
		stack = rt.THStack()
		LoH = []
		for h in LoHi:
			LoH.append(h.Clone())
		for i in range(len(LoH)):
			LoH[i].SetLineColor(i+1)
			stack.Add(LoH[i])
		stack.Draw("nostack hist")
		stack.GetXaxis().SetTitle(LoH[0].GetXaxis().GetTitle())
		stack.GetYaxis().SetTitle("Count")
		stack.SetTitle(name)
		c1.Modified()
		if doLeg:
			c1.BuildLegend(0.7,0.7,0.9,0.9)
		c1.SaveAs(self.extraDir+name+".png")
		if doCum:
			sCum = rt.THStack()
			for i in range(len(LoH)):
				sCum.Add(LoH[i].GetCumulative())
			sCum.Draw("nostack")
			sCum.GetYaxis().SetTitle("Cumulative Fraction")
			sCum.Modified()
			if doLeg:
				c1.BuildLegend(0.7,0.7,0.9,0.9)
			c1.SaveAs(self.extraDir+name+"_cum.png")
		for i in range(len(LoH)):
			if LoH[i].Integral() != 0:
				LoH[i].Scale(1/LoH[i].Integral())
		stack.Draw("nostack hist")
		stack.GetYaxis().SetTitle("a.u.")
		c1.Modified()
		if doLeg:
			c1.BuildLegend(0.7,0.7,0.9,0.9)
		c1.SaveAs(self.extraDir+name+"_norm.png")
		if doCum:
			sCum = rt.THStack()
			for i in range(len(LoH)):
				sCum.Add(LoH[i].GetCumulative())
			sCum.Draw("nostack")
			sCum.GetYaxis().SetTitle("Normalized Cumulative Fraction")
			sCum.Modified()
			if doLeg:
				c1.BuildLegend(0.7,0.7,0.9,0.9)
			c1.SaveAs(self.extraDir+name+"_normcum.png")

	def makeRatio(self, h2, h1 ,name, doLeg = True, log = False):
		# C++ Author: Olivier Couet, adapted to python by Colin Fallon
		# Define the Canvas
		c = rt.TCanvas("c", "canvas", 800, 800)

		# Upper plot will be in pad1
		pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
		pad1.SetBottomMargin(0) # Upper and lower plot are joined
		pad1.SetGridx()         # Vertical grid
		pad1.Draw()             # Draw the upper pad: pad1
		pad1.cd()               # pad1 becomes the current pad
		stack = rt.THStack()	# use THStack so that top of histograms dont get cut off
		h1.SetStats(0)       # No statistics on upper plot
		stack.Add(h1)			
		stack.Add(h2)
		stack.Draw("nostack hist")         # Draw h2 on top of h1
		stack.GetYaxis().SetTitle(h1.GetTitle().split(" ")[0])
		if doLeg:
			pad1.BuildLegend(0.75,0.75,0.95,0.95)

		# Do not draw the Y axis label on the upper plot and redraw a small
		# axis instead, in order to avoid the first label (0) to be clipped.
		#h1.GetYaxis().SetLabelSize(0.)
		#axis = rt.TGaxis( -5, 20, -5, 220, 20,220,510,"")
		#axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
		#axis.SetLabelSize(15)
		#axis.Draw()

		# lower plot will be in pad
		c.cd()          # Go back to the main canvas before defining pad2
		pad2 = rt.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
		pad2.SetTopMargin(0)
		pad2.SetBottomMargin(0.2)
		pad2.SetGridx() # vertical grid
		pad2.Draw()
		pad2.cd()       # pad2 becomes the current pad

		# Define the ratio plot
		h3 = h1.Clone("h3")
		h3.SetLineColor(rt.kBlack)
		h3.SetMinimum(0.00)  # Define Y ..
		h3.SetMaximum(2.00) # .. range
		h3.Sumw2()
		h3.SetStats(0)      # No statistics on lower plot
		h3.Divide(h2)
		h3.SetMarkerStyle(21)
		h3.Draw("ep")       # Draw the ratio plot

		# h1 settings
		h1.SetLineColor(rt.kBlue+1)
		h1.SetLineWidth(2)

		# Y axis h1 plot settings
		h1.GetYaxis().SetTitleSize(20)
		h1.GetYaxis().SetTitleFont(43)
		h1.GetYaxis().SetTitleOffset(1.55)

		# h2 settings
		h2.SetLineColor(rt.kRed)
		h2.SetLineWidth(2)

		# Ratio plot (h3) settings
		h3.SetTitle("") # Remove the ratio title

		# Y axis ratio plot settings
		h3.GetYaxis().SetTitle(h1.GetTitle().split(" ")[0]+"/"+h2.GetTitle().split(" ")[0])
		h3.GetYaxis().SetNdivisions(505)
		h3.GetYaxis().SetTitleSize(20)
		h3.GetYaxis().SetTitleFont(43)
		h3.GetYaxis().SetTitleOffset(1.55)
		h3.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
		h3.GetYaxis().SetLabelSize(15)

		# X axis ratio plot settings
		h3.GetXaxis().SetTitleSize(20)
		h3.GetXaxis().SetTitleFont(43)
		h3.GetXaxis().SetTitleOffset(4.)
		h3.GetXaxis().SetLabelFont(43); # Absolute font size in pixel (precision 3)
		h3.GetXaxis().SetLabelSize(15)

		#save as .png
		c.SaveAs(self.extraDir+name+"_ratio.png")

	def makeRatioStack(self, stack, data , name, doLeg = True, log = False):
		# C++ Author: Olivier Couet, adapted to python by Colin Fallon
		# Define the Canvas
		c = rt.TCanvas("c", "canvas", 800, 800)

		# Upper plot will be in pad1
		pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
		pad1.SetBottomMargin(0) # Upper and lower plot are joined
		pad1.SetGridx()         # Vertical grid
		pad1.Draw()             # Draw the upper pad: pad1
		pad1.cd()               # pad1 becomes the current pad
		data.SetStats(0)       # No statistics on upper plot
		stack.Draw("hist")
		data.Draw("same E1")
		if doLeg:
			pad1.BuildLegend(0.75,0.75,0.95,0.95)

		# Do not draw the Y axis label on the upper plot and redraw a small
		# axis instead, in order to avoid the first label (0) to be clipped.
		#h1.GetYaxis().SetLabelSize(0.)
		#axis = rt.TGaxis( -5, 20, -5, 220, 20,220,510,"")
		#axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
		#axis.SetLabelSize(15)
		#axis.Draw()

		# lower plot will be in pad
		c.cd()          # Go back to the main canvas before defining pad2
		pad2 = rt.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
		pad2.SetTopMargin(0)
		pad2.SetBottomMargin(0.2)
		pad2.SetGridx() # vertical grid
		pad2.Draw()
		pad2.cd()       # pad2 becomes the current pad

		# Define the ratio plot
		h3 = data.Clone("h3")
		h3.SetLineColor(rt.kBlack)
		h3.SetMinimum(0.00)  # Define Y ..
		h3.SetMaximum(2.00) # .. range
		h3.Sumw2()
		h3.SetStats(0)      # No statistics on lower plot
		h3.Divide(stack.GetStack().Last())
		h3.SetMarkerStyle(21)
		h3.Draw("ep")       # Draw the ratio plot

		# h1 settings
		data.SetLineColor(rt.kBlue+1)
		data.SetLineWidth(2)

		# Y axis h1 plot settings
		data.GetYaxis().SetTitleSize(20)
		data.GetYaxis().SetTitleFont(43)
		data.GetYaxis().SetTitleOffset(1.55)

		# Ratio plot (h3) settings
		h3.SetTitle("") # Remove the ratio title

		# Y axis ratio plot settings
		h3.GetYaxis().SetTitle("Data/Background")
		h3.GetYaxis().SetNdivisions(505)
		h3.GetYaxis().SetTitleSize(20)
		h3.GetYaxis().SetTitleFont(43)
		h3.GetYaxis().SetTitleOffset(1.55)
		h3.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
		h3.GetYaxis().SetLabelSize(15)

		# X axis ratio plot settings
		h3.GetXaxis().SetTitleSize(20)
		h3.GetXaxis().SetTitleFont(43)
		h3.GetXaxis().SetTitleOffset(4.)
		h3.GetXaxis().SetLabelFont(43); # Absolute font size in pixel (precision 3)
		h3.GetXaxis().SetLabelSize(15)

		#save as .png
		c.SaveAs(self.extraDir+name+"_ratio.png")

	def passedPreselection(self, eventJetsAK8, eventMET, eventMT_AK8, eventNElectrons, eventNMuons):
		nJets = len(eventJetsAK8)
		if nJets < 2:
			return False
		jet1Pt = eventJetsAK8[0].Pt()
		jet2Pt = eventJetsAK8[1].Pt()
		jet1Eta = eventJetsAK8[0].Eta()
		jet2Eta = eventJetsAK8[1].Eta()
		deltaEta = abs(jet1Eta-jet2Eta)
		METoverMT = eventMET/eventMT_AK8
		nLep = eventNElectrons+eventNMuons
		#missing MET filters?
		
		test1 = (nJets >= 2)
		test2 = (jet1Pt > 200) and (jet2Pt > 200)
		test3 = ((abs(jet1Eta) < 2.4) and (abs(jet2Eta) < 2.4))
		test4 = (deltaEta < 1.5)

		test5 = (METoverMT > 0.15)
		test6 = (eventMT_AK8 > 1500)

		test7 = (nLep == 0)

		testJets = test1 and test2 and test3 and test4
		testMT = test5 and test6
		testLeps = test7
		if testJets and testLeps and testMT:
			return True
		else:
			return False
			
		

	def make2dPng(self, hist, name):
		c1 = rt.TCanvas("c1","c1",1200,900)
		hist.Draw("colz")
		c1.SaveAs(self.extraDir+name+".png")

	def write(self):
		self.outRootFile.cd()
		for thing in self.objects:
			try:
				thing.Write()
				print("Object {} has been written".format(thing.GetName()))
			except AttributeError:
				print("Object {} cannot be written!".format(thing))

	def selfprint(self):
		print(str(datetime.datetime.now()))
		print("------------------------")
		print("Analysing these files:")
		for line in self.inputFileList:
			print(line)
		print("Using these trees:")
		for line in self.treeNameList:
			print(line)
		print("And saving output here:")
		print(self.extraDir+self.outFileName)
		print("------------------------")
	
	def run(self):
		self.selfprint()
		self.loop()
		self.write()













































