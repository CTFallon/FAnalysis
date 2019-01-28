# base class for FAnalysis, copied from Lucien Lo's PFG HcalTupleAnalyzer, adapted to python
import ROOT as rt
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
		self.inputFileList = []
		try:
			self.inputFileList.append(fileDict[self.fileID])#ignore new line character
		except KeyError:
			exit("Key does not exist in dictonary")
			
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
	
	def makeTGraph(self, n, x, y):
		graph = rt.TGraph(n,x,y)
		self.objects.append(graph)
		return graph
	
	def makeTGraph(self):
		graph = rt.TGraph()
		self.objects.append(graph)
		return graph

	def makePng(self, LoH, name, doLeg = True):
		c = rt.TCanvas("c1","c1",1200,900)
		stack = rt.THStack()
		for i in range(len(LoH)):
			LoH[i].SetLineColor(i+1)
			stack.Add(LoH[i])
		stack.Draw("nostack")
		if doLeg:
			c.BuildLegend(0.6,0.2,0.8,0.4)
		c.SaveAs(self.extraDir+name+".png")
		for i in range(len(LoH)):
			if LoH[i].Integral() != 0:
				LoH[i].Scale(1/LoH[i].Integral())
		stack.Draw("nostack")
		c.SaveAs(self.extraDir+name+"_norm.png")

	def write(self):
		self.outRootFile.cd()
		for thing in self.objects:
			try:
				thing.Write()
				print("Object {} has been written".format(thing.GetName()))
			except AttributeError:
				print("Object {} cannot be written!".format(thing))

	def selfprint(self):
		print("------------------------")
		print("Analysing these files:")
		for line in self.inputFileList:
			print(line)
		print("Using these trees:")
		for line in self.treeNameList:
			print(line)
		print("And saving output here:")
		print(self.outFileName)
		print("------------------------")
	
	def run(self):
		self.selfprint()
		self.loop()
		self.write()













































