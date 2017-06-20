from analysis_base import *

class analysisClass(analysisBase):
	def __init__(self):
		analysisBase.__init__(self)

	def loop(self, config): # Loop runs once
		h_GenPart_Mag_stack = self.createStack("GenPart_Mag")
		h_GenPart_3P_stack = self.createStack("GenPart_3P")
		h_GenPart_Pt_stack = self.createStack("GenPart_Pt")

		#chain all the relavant background files together
		for bkgType in config.getbkgList(): # from here to the end runs for each section in the .ini file
			tempChain = rt.TChain(bkgType)
			for iFile in range(config.getNumberOfFiles(bkgType)):
				fileName = config.getFileName(bkgType,iFile)
				tempChain.Add(fileName)

			#creation of histograms
			h_Mag_temp = self.createTH1F("GenPart_Mag_"+bkgType,1000,-1,45)
			h_3P_temp = self.createTH1F("GenPart_3P_"+bkgType,1000,0,1000)
			h_Pt_temp = self.createTH1F("GenPart_Pt_"+bkgType,1000,0,1000)
			h_Mag_temp.SetFillColor(config.getColor(bkgType))
			h_3P_temp.SetFillColor(config.getColor(bkgType))
			h_Pt_temp.SetFillColor(config.getColor(bkgType))


			nEvents = tempChain.GetEntries()
			for iEvent in range(nEvents): # from here to the .Add lines runs once for every event in the TChain created above
				tempChain.GetEvent(iEvent)
				if iEvent == 0:
					print "Number of Events is " + str(nEvents) + " and weight is " +str(tempChain.Weight)+"."
				if iEvent%1000 == 0:
					print "Processing Event Number:", iEvent, "in background", bkgType
				if tempChain.BTags == 3:
					for particle in tempChain.GenParticles:
						h_Mag_temp.Fill(particle.Mag(),tempChain.Weight)
						h_3P_temp.Fill(particle.P(),tempChain.Weight)
						h_Pt_temp.Fill(particle.Pt(),tempChain.Weight)	

			h_GenPart_Mag_stack.Add(h_Mag_temp,"hist")
			h_GenPart_3P_stack.Add(h_3P_temp,"hist")
			h_GenPart_Pt_stack.Add(h_Pt_temp,"hist")

			del h_Mag_temp
			del h_3P_temp
			del h_Pt_temp
			del tempChain
