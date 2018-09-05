from analysisBase import baseClass

def loop(self):
	# set up trees
	f = rt.TFile.Open(self.inputFileList[0])
	tree = f.Get(self.treeNameList[0])
	
	nEvents = tree.GetEntries()
	
	for iEvent in range(nEvents):
		tree.GetEvent(iEvent)
		print(len(tree.JetsAK8))
	

def addLoop():
	baseClass.loop = loop


