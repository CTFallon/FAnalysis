from analysis_base import * # this clas inherits from the analysisBase class

# first, it creates any and all root objects that are needed outside the scope of a specific background

# then, it starts looping over each background
#			for each background it:
#					creates a TChain object of all the background files
#					creates all the local root objects needed for the background
#					loops over all events in the TChain object and fills the local (and if needed, global) histograms
#					deletes the background root objects from memory 
#								(as long as you use the analysis_base 'create____' functions, they will be saved later)
class analysisClass(analysisBase): 
	def __init__(self): 
		analysisBase.__init__(self)

	def loop(self, config): # only change things within this function
		#define 'global' histograms
		JetsAK8_Mag_stack = self.createStack("JetsAK8_Mag")
		Electrons_Mag_stack = self.createStack("Electrons_Mag")
		#chain all the relavant background files together
		# this shouldn't need to be modified if there is only one section in the input, thats fine
		for bkgType in config.getbkgList():
			tempChain = rt.TChain(bkgType)
			for iFile in range(config.getNumberOfFiles(bkgType)):
				fileName = config.getFileName(bkgType,iFile)
				tempChain.Add(fileName)

			#creation of local histograms
			# this is where you create the background-specific histograms
			# dont forget the 'create____' functions from analysis_base
			h_JetsAK8_Mag_temp = self.createTH1F("JetsAK8_Mag_"+bkgType,1000,-10,400)
			h_JetsAK8_Mag_temp.SetFillColor(config.getColor(bkgType))
			h_Electrons_Mag_temp = self.createTH1F("Electrons_Mag_"+bkgType,1000,-.5,.5)
			h_Electrons_Mag_temp.SetFillColor(config.getColor(bkgType))


			nEvents = tempChain.GetEntries()
			print "There are " +str(nEvents) + " events in "+ bkgType+"."
			for iEvent in range(nEvents): #loop over all events....
				tempChain.GetEvent(iEvent)
				#This is where your cuts should go. 
				if iEvent % 10000 == 0:
					print "Procesing event " + str(iEvent)
				for particle in tempChain.JetsAK8:
					h_JetsAK8_Mag_temp.Fill(particle.Mag(),1)
				for particle in tempChain.Electrons:
					h_Electrons_Mag_temp.Fill(particle.Mag(),1)
				

			
			JetsAK8_Mag_stack.Add(h_JetsAK8_Mag_temp)
			Electrons_Mag_stack.Add(h_Electrons_Mag_temp)
			del h_JetsAK8_Mag_temp
			del h_Electrons_Mag_temp # delete the temp objects so that they're freed up for next loop
			del tempChain # as long as you used the create functions that store the objects in 
			# the analysis object's lists, they will be written later on.
