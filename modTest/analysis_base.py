import ROOT as rt



class analysisBase:
	def __init__(self):
		self.objList = []

	def createStack(self,stackName):
		stack = rt.THStack("h_"+stackName+"_stack",stackName+"_stack")
		self.objList.append(stack)
		return stack

	def createTH1F(self,histName,nBins,xLo,xHi):
		hist = rt.TH1F("h_"+histName, histName, nBins, xLo, xHi)
		self.objList.append(hist)
		return hist

	def writeObjects(self):
		for obj in self.objList:
			obj.Write()

	def run(self,config):
		self.loop(config)
		self.writeObjects()
