import ROOT as rt
rt.gROOT.SetBatch(True)

#run from FAnalysis folder with python macros/bdtOpt/draw_bdtPlots_threeWP.py


list16 = ["QCD16","TTJets16","WJets16","ZJets16", "base"]
list17 = ["QCD17","TTJets17","WJets17","ZJets17", "base"]
list18PRE = ["QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE", "base"]
list18POST = ["QCD18POST","TTJets18POST","WJets18POST","ZJets18POST", "base"]

listZ = ["z15", "z20", "z25", "base", "z35", "z40","z45"]
listD = ["d10","base","d40","d60","d80","d00"]
listR = ["r02","base","r05","r07","r08"]
listA = ["adl", "base", "adh"]

colorDict = {"QC":602, "TT":798,"WJ":801,"ZJ":881,"ba":2}

def makePng(LoHi, name, direc, doLeg = True, log = False, nBinsX = 50, noStack = True):
	c1 = rt.TCanvas("c1","c1",1200,900)
	if log:
		c1.SetLogy()
	stack = rt.THStack()
	LoH = []
	for h in LoHi[:-1]:
		LoH.append(h.Clone())
	LoH.sort(key=lambda x: x.Integral())
	if noStack:
		for i in range(len(LoH)):
			LoH[i].SetLineColor(colorDict[LoH[i].GetName()[:2]])
			if LoH[i].GetNbinsX() != nBinsX:
				LoH[i].Rebin(LoH[i].GetNbinsX()/nBinsX)
			stack.Add(LoH[i])
		base = LoHi[-1].Clone()
		if base.GetNbinsX() != nBinsX:
			base.Rebin(base.GetNbinsX()/nBinsX)
		base.SetLineColor(rt.kRed)
		stack.Draw("nostack hist")
		base.Draw("same")
	else:
		for i in range(len(LoH)):
			LoH[i].SetLineColor(colorDict[LoH[i].GetName()[:2]])
			LoH[i].SetFillColor(colorDict[LoH[i].GetName()[:2]])
			LoH[i].SetFillStyle(1001)
			if LoH[i].GetNbinsX() != nBinsX:
				LoH[i].Rebin(LoH[i].GetNbinsX()/nBinsX)
			stack.Add(LoH[i])
		base = LoHi[-1].Clone()
		if base.GetNbinsX() != nBinsX:
			base.Rebin(base.GetNbinsX()/nBinsX)
		base.SetLineColor(rt.kRed)
		stack.Draw("hist")
		base.Draw("same")
	stack.GetXaxis().SetTitle(LoH[0].GetXaxis().GetTitle())
	stack.GetYaxis().SetTitle("Count")
	stack.SetTitle(name)
	stack.SetMaximum(max(stack.GetStack().Last().GetMaximum(),base.GetMaximum())*1.1)
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.6,0.6,0.9)
	c1.SaveAs(direc+name+".png")
	nEventsStack = stack.GetStack().Last().Integral()
	nEventsBase = base.Integral()
	print("{}, {}, {}".format(name, nEventsStack,nEventsBase))
	for i in range(len(LoH)):
		if LoH[i].Integral() != 0:
			LoH[i].Scale(1/LoH[i].Integral())
	stack.Draw("nostack hist")
	base.Draw("same")
	stack.GetYaxis().SetTitle("a.u.")
	stack.SetMaximum(max(stack.GetStack().Last().GetMaximum(),base.GetMaximum())*1.1)
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.6,0.6,0.9)
	c1.SaveAs(direc+name+"_norm.png")

def makeRatioStack(stack, data , name, direc, doLeg = True, log = False):
	# C++ Author: Olivier Couet, adapted to python by Colin Fallon
	# Define the Canvas
	c = rt.TCanvas("c", "canvas", 800, 800)

	# Upper plot will be in pad1
	pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0) # Upper and lower plot are joined
	pad1.SetGridx()         # Vertical grid
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	data.SetStats(0)       # No statistics on upper plot
	stack.Draw("hist")
	data.Draw("same E1")
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
	h3 = data.Clone("h3")
	h3.SetLineColor(rt.kBlack)
	h3.SetMinimum(0.00)  # Define Y ..
	h3.SetMaximum(2.00) # .. range
	h3.Sumw2()
	h3.SetStats(0)      # No statistics on lower plot
	h3.Divide(stack.GetStack().Last())
	h3.SetMarkerStyle(21)
	h3.Draw("ep")       # Draw the ratio plot

	# h1 settings
	data.SetLineColor(rt.kBlue+1)
	data.SetLineWidth(2)

	# Y axis h1 plot settings
	data.GetYaxis().SetTitleSize(20)
	data.GetYaxis().SetTitleFont(43)
	data.GetYaxis().SetTitleOffset(1.55)

	# Ratio plot (h3) settings
	h3.SetTitle("") # Remove the ratio title

	# Y axis ratio plot settings
	h3.GetYaxis().SetTitle("Data/Background")
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

	#save as .png
	c.SaveAs(direc+name+"_ratio.png")

#bkgList = [list16, list17, list18PRE, list18POST]
bkgList = [list17]
sigList = [listZ, listD, listR, listA]

#varaibles = ['MET','METPhi','metR','JetsAK8_bdtSVJtag[0]','JetsAK8_bdtSVJtag[1]','JetsAK8_bdtSVJtag']
varaibles = ['JetsAK8_bdtSVJtag[0]']
#wpLabel = "tight"
#wpLabel = "medium"
#wpLabel = "loose"
for wpLabel in ["loose","medium","tight"]:
	for listYear in bkgList:
		for var in varaibles:
			histList_svj0 = []
			histList_svj1 = []
			histList_svj2 = []
			for dataSet in listYear:
				#print(dataSet)
				_file1 = rt.TFile.Open("outputs/bdtOpt/"+dataSet+"/bdt_threeWP_"+wpLabel+".root","read")
				#_file1.ls()
				histList_svj0.append(_file1.Get(var+"_svj0_"+dataSet).Clone(dataSet+"_svj0"))
				histList_svj1.append(_file1.Get(var+"_svj1_"+dataSet).Clone(dataSet+"_svj1"))
				histList_svj2.append(_file1.Get(var+"_svj2_"+dataSet).Clone(dataSet+"_svj2"))
				histList_svj0[-1].SetDirectory(0)
				histList_svj1[-1].SetDirectory(0)
				histList_svj2[-1].SetDirectory(0)

			directory = "outputs/bdtOpt/Plots_threeWP/"
			makePng(histList_svj0, "svj0_stack_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=False)
			makePng(histList_svj1, "svj1_stack_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=False)
			makePng(histList_svj2, "svj2_stack_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=False)
			makePng(histList_svj0, "svj0_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=True)
			makePng(histList_svj1, "svj1_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=True)
			makePng(histList_svj2, "svj2_"+listYear[0][3:]+"_"+wpLabel+"_"+var, directory, doLeg = True, log = True, noStack=True)




















