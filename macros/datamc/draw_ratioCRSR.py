import ROOT as rt
rt.gROOT.SetBatch(True)

#run from FAnalysis folder with python macros/datamc/draw_ratioCRSR.py


list16 = ["QCD16","TTJets16","WJets16","ZJets16","Data16"]
list17 = ["QCD17","TTJets17","WJets17","ZJets17","Data17"]
list18PRE = ["QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE","Data18PRE"]
list18POST = ["QCD18POST","TTJets18POST","WJets18POST","ZJets18POST","Data18POST"]

bkgList = [list16, list17, list18PRE, list18POST]

colorDict = {"QC":602, "TT":798,"WJ":801,"ZJ":881,"ba":2}

def makeRatio(histCR, histSR , name, direc, doLeg = True, log = False, doNorm = False):
	# C++ Author: Olivier Couet, adapted to python by Colin Fallon
	# Define the Canvas
	c = rt.TCanvas("c", "canvas", 800, 800)

	# Upper plot will be in pad1
	pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0) # Upper and lower plot are joined
	pad1.SetGridx()         # Vertical grid
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	if log:
		pad1.SetLogy()
	histCR.SetStats(0)       # No statistics on upper plot
	if histSR.GetNbinsX() == 1000:
		histSR.Rebin(20)
		histCR.Rebin(20)
	histSR.SetTitle("Signal Region")
	histCR.SetTitle("Control Region")
	histSR.SetLineColor(1)
	histCR.SetLineColor(2)
	if doNorm:
		histCR.Scale(1/histCR.Integral())
		histSR.Scale(1/histSR.Integral())
	stack = rt.THStack()
	stack.Add(histCR)
	stack.Add(histSR)
	stack.SetTitle(name)
	stack.Draw("nostack hist")
	#histCR.Draw("hist")
	#histSR.Draw("same hist")
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
	h3 = histSR.Clone("h3")
	h3.SetLineColor(rt.kBlack)
	h3.SetMinimum(0.00)  # Define Y ..
	if doNorm:
		h3.SetMaximum(2.00)
	else:
		h3.SetMaximum(1.50) # .. range
	h3.Sumw2()
	h3.SetStats(0)      # No statistics on lower plot
	h3.Divide(histCR)
	h3.SetMarkerStyle(21)
	h3.Draw("ep")       # Draw the ratio plot

	# Y axis h1 plot settings
	histSR.GetYaxis().SetTitleSize(20)
	histSR.GetYaxis().SetTitleFont(43)
	histSR.GetYaxis().SetTitleOffset(1.55)

	# Ratio plot (h3) settings
	h3.SetTitle("") # Remove the ratio title

	one = rt.TF1("one", '1', h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax())
	one.SetLineStyle(2)
	one.SetLineColor(rt.kBlack)
	one.Draw("same")

	# Y axis ratio plot settings
	h3.GetYaxis().SetTitle("SR/CR")
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
	
	c.Update()
	#save as .png
	c.SaveAs(direc+name)

directory = "outputs/dmc/plotsSRCR_CRde_2p2/"
for listYear in bkgList:
	leadSRstack = rt.THStack()
	leadCRstack = rt.THStack()
	subSRstack = rt.THStack()
	subCRstack = rt.THStack()
	allSRstack = rt.THStack()
	allCRstack = rt.THStack()
	for dataSet in listYear:
		print(dataSet)
		listOfSR = []
		listOfCR = []
		_file1 = rt.TFile.Open("outputs/dmc/"+dataSet+"/dataMC_comp_nomf.root","read") # signal region
		#print("JetsAK8_bdtSVJtag[0]_"+dataSet,"JetsAK8_bdtSVJtag[1]_"+dataSet,"JetsAK8_bdtSVJtag_"+dataSet)
		listOfSR.append(_file1.Get("JetsAK8_bdtSVJtag[0]_"+dataSet).Clone(dataSet+"_SRleadBDT"))
		listOfSR[-1].SetDirectory(0)
		leadSRstack.Add(listOfSR[-1])
		listOfSR.append(_file1.Get("JetsAK8_bdtSVJtag[1]_"+dataSet).Clone(dataSet+"_SRsubBDT"))
		listOfSR[-1].SetDirectory(0)
		subSRstack.Add(listOfSR[-1])
		listOfSR.append(_file1.Get("JetsAK8_bdtSVJtag_"+dataSet).Clone(dataSet+"_SRallBDT"))
		listOfSR[-1].SetDirectory(0)
		allSRstack.Add(listOfSR[-1])
		_file1.Close()
		_file2 = rt.TFile.Open("outputs/dmc/"+dataSet+"/dataMC_comp_CRde.root","read") # control region
		listOfCR.append(_file2.Get("JetsAK8_bdtSVJtag[0]_"+dataSet).Clone(dataSet+"_CRleadBDT"))
		listOfCR[-1].SetDirectory(0)
		leadCRstack.Add(listOfCR[-1])
		listOfCR.append(_file2.Get("JetsAK8_bdtSVJtag[1]_"+dataSet).Clone(dataSet+"_CRsubBDT"))
		listOfCR[-1].SetDirectory(0)
		subCRstack.Add(listOfCR[-1])
		listOfCR.append(_file2.Get("JetsAK8_bdtSVJtag_"+dataSet).Clone(dataSet+"_CRallBDT"))
		listOfCR[-1].SetDirectory(0)
		allCRstack.Add(listOfCR[-1])
		_file2.Close()
		#makeRatio(histCR, histSR , name, direc, doLeg = True, log = False)
		makeRatio(listOfCR[0], listOfSR[0] , dataSet+"_leadBDT_SRCRde.png", directory, doLeg = True, log = True)
		makeRatio(listOfCR[1], listOfSR[1] , dataSet+"_subBDT_SRCRde.png", directory, doLeg = True, log = True)
		makeRatio(listOfCR[2], listOfSR[2] , dataSet+"_allBDT_SRCRde.png", directory, doLeg = True, log = True)
		makeRatio(listOfCR[0], listOfSR[0] , dataSet+"_leadBDTnorm_SRCRde.png", directory, doLeg = True, log = True, doNorm = True)
		makeRatio(listOfCR[1], listOfSR[1] , dataSet+"_subBDTnorm_SRCRde.png", directory, doLeg = True, log = True, doNorm = True)
		makeRatio(listOfCR[2], listOfSR[2] , dataSet+"_allBDTnorm_SRCRde.png", directory, doLeg = True, log = True, doNorm = True)
	makeRatio(leadCRstack.GetStack()[3], leadSRstack.GetStack()[3] , listYear[0][3:]+"_leadBDT_SRCRdesum.png", directory, doLeg = True, log = True)
	makeRatio(subCRstack.GetStack()[3], subSRstack.GetStack()[3] , listYear[0][3:]+"_subBDT_SRCRdesum.png", directory, doLeg = True, log = True)
	makeRatio(allCRstack.GetStack()[3], allSRstack.GetStack()[3] , listYear[0][3:]+"_allBDT_SRCRdesum.png", directory, doLeg = True, log = True)
	makeRatio(leadCRstack.GetStack()[3], leadSRstack.GetStack()[3] , listYear[0][3:]+"_leadBDTnorm_SRCRdesum.png", directory, doLeg = True, log = True, doNorm = True)
	makeRatio(subCRstack.GetStack()[3], subSRstack.GetStack()[3] , listYear[0][3:]+"_subBDTnorm_SRCRdesum.png", directory, doLeg = True, log = True, doNorm = True)
	makeRatio(allCRstack.GetStack()[3], allSRstack.GetStack()[3] , listYear[0][3:]+"_allBDTnorm_SRCRdesum.png", directory, doLeg = True, log = True, doNorm = True)

