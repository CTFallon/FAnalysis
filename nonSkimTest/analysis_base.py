import ROOT as rt



class analysisBase:
	def __init__(self):
		self.stackList = []
		self.histList = []

	def createStack(self,stackName):
		stack = rt.THStack("h_"+stackName+"_stack",stackName+"_stack")
		self.stackList.append(stack)
		return stack

	def createTH1F(self,histName,nBins,xLo,xHi):
		hist = rt.TH1F("h_"+histName, histName, nBins, xLo, xHi)
		self.histList.append(hist)
		return hist

	def writeHists(self):
		for hist in self.histList:
			hist.Write()

	def writeStacks(self):
		for stack in self.stackList:
			stack.Write()
