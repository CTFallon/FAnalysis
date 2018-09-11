# base class for FAnalysis, copied from Lucien Lo's PFG HcalTupleAnalyzer, adapted to python
import ROOT as rt
class baseClass:
	def __init__(self, fileList, treeList, rootDirct, outFileName):
		self.fileList = fileList # list of ROOT files to read
		self.treeList = treeList # list of Tree names in each ROOT file
		self.outFileName = outFileName
		self.objects = []
		self.loadFileList()
		self.loadTreeList()
		self.extraDir = rootDirt+'/'
		self.loadOutFile()

	def loadFileList(self):
		inFile = open(self.fileList)
		self.inputFileList = []
		for line in inFile:
			self.inputFileList.append(line[:-1])#ignore new line character
		inFile.close()

	def loadTreeList(self):
		treeFile = open(self.treeList)
		self.treeNameList = []
		for line in treeFile:
			self.treeNameList.append(line[:-1])#ignore new line character
		treeFile.close()

	def loadOutFile(self):
		self.outRootFile = rt.TFile.Open(self.extraDir+self.outFileName, 'RECREATE')
		
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

	def write(self):
		self.outRootFile.cd()
		for thing in self.objects:
			thing.Write()

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













































