import ConfigParser as cp


class Fparser:
	def __init__(self, inputFile):
		self.config = cp.SafeConfigParser()
		self.config.read(inputFile)

	def getbkgList(self):
		return self.config.sections()

	def getNumberOfFiles(self,bkgType):
		return self.config.getint(bkgType,'nfile')

	def getFileName(self,bkgType,n):
		return self.config.get(bkgType,'file'+str(n)+".Name")

	def getColor(self,bkgType):
		return self.config.getint(bkgType,'color')

