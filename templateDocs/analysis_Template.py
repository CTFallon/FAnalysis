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
		# generally, histograms that are not background-specific
		# don't forget to use the 'create___' functions defined (or need to be defined) in analysis_base
		# example: 
		h_GenPart_Mag_stack = self.createStack("GenPart_Mag")

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
			h_Mag_temp = self.createTH1F("GenPart_Mag_"+bkgType,1000,-1,45)
			h_3P_temp.SetFillColor(config.getColor(bkgType))


			nEvents = tempChain.GetEntries()
			for iEvent in range(nEvents): #loop over all events....
				tempChain.GetEvent(iEvent)
				if iEvent%1000 == 0: #not needed, but nice to know the program is working...
					print "Processing Event Number:", iEvent, "in background", bkgType
				#This is where your cuts should go. In this example, we're cutting out all
				# events that have 0, 1, 2, or more than 3 b quark tags
				if tempChain.BTags == 3:
					for particle in tempChain.GenParticles: # this is because GenParticles is a vector, some histos won't need this
						h_Mag_temp.Fill(particle.Mag(),tempChain.Weight)
			h_GenPart_Mag_stack.Add(h_Mag_temp,"hist")

			del h_Mag_temp # delete the temp objects so that they're freed up for next loop
			del tempChain # as long as you used the create functions that store the objects in 
			# the analysis object's lists, they will be written later on.
