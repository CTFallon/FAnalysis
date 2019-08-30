import ROOT as rt
rt.gROOT.SetBatch(True)

#run from FAnalysis folder with python macros/bdtOpt/draw_bdtPlots.py


list16 = ["QCD16","TTJets16","WJets16","ZJets16", "base"]
list17 = ["QCD17","TTJets17","WJets17","ZJets17"]#, "base"]
list18PRE = ["QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE"]#, "base"]
list18POST = ["QCD18POST","TTJets18POST","WJets18POST","ZJets18POST"]#, "base"]

listZ = ["z15", "z20", "z25", "base", "z35", "z40","z45"]
listD = ["d10","base","d40","d60","d80","d00"]
listR = ["r02","base","r05","r07","r08"]
listA = ["adl", "base", "adh"]

colorDict = {"QC":602, "TT":798,"WJ":801,"ZJ":881,"ba":2}

def makePng(LoHi, name, direc, doLeg = True, log = False, nBinsX = 50, stackLines = True):
	c1 = rt.TCanvas("c1","c1",1200,900)
	if log:
		c1.SetLogy()
	stack = rt.THStack()
	LoH = []
	for h in LoHi[:-1]:
		LoH.append(h.Clone())
	LoH.sort(key=lambda x: x.Integral())
	if stackLines:
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
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.7,0.6,0.9)
	c1.SaveAs(direc+name+".png")
	for i in range(len(LoH)):
		if LoH[i].Integral() != 0:
			LoH[i].Scale(1/LoH[i].Integral())
	stack.Draw("nostack hist")
	stack.GetYaxis().SetTitle("a.u.")
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.7,0.6,0.9)
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

bkgList = [list16, list17, list18PRE, list18POST]
sigList = [listZ, listD, listR, listA]

#cuts = [  0, 50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
#cuts = [500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000]
#cuts = [900,905,910,915,920,925,930,935,940,945,950,955,960,965,970,975,980,985,990,995,1000]
#cuts = [980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000]
cuts = [x for x in range(0,1000,1)]

nLeadPassDict = {}
nLeadTotDict = {}
nSubPassDict = {}
nSubTotDict = {}
nAllPassDict = {}
nAllTotDict = {}

for key in ["QCD16","TTJets16","WJets16","ZJets16","QCD17","TTJets17","WJets17","ZJets17","QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE","QCD18POST","TTJets18POST","WJets18POST","ZJets18POST","base"]:
	nLeadPassDict[key] = {}
	nLeadTotDict[key] = {}
	nSubPassDict[key] = {}
	nSubTotDict[key] = {}
	nAllPassDict[key] = {}
	nAllTotDict[key] = {}
	for cutVal in cuts:
		cutKey = str(cutVal)
		nLeadPassDict[key][cutKey] = 0
		nLeadTotDict[key][cutKey] = 0
		nSubPassDict[key][cutKey] = 0
		nSubTotDict[key][cutKey] = 0
		nAllPassDict[key][cutKey] = 0
		nAllTotDict[key][cutKey] = 0

for listYear in bkgList:
	histList_leadBDT = []
	histList_subBDT = []
	histList_allBDT = []
	for dataSet in listYear:
		print(dataSet)
		_file1 = rt.TFile.Open("outputs/bdtOpt/"+dataSet+"/bdt_opti.root","read")
		histList_leadBDT.append(_file1.Get("hist_leadBDT").Clone(dataSet+"_leadBDT"))
		histList_subBDT.append(_file1.Get("hist_subBDT").Clone(dataSet+"_subBDT"))
		histList_allBDT.append(_file1.Get("hist_allBDT").Clone(dataSet+"_allBDT"))
		histList_leadBDT[-1].SetDirectory(0)
		histList_subBDT[-1].SetDirectory(0)
		histList_allBDT[-1].SetDirectory(0)
		for i in cuts:
			nLeadPassDict[dataSet][str(i)] = histList_leadBDT[-1].Integral(i, 1000)
			nLeadTotDict[dataSet][str(i)] = histList_leadBDT[-1].Integral(0, 1000)
			nSubPassDict[dataSet][str(i)] = histList_subBDT[-1].Integral(i, 1000)
			nSubTotDict[dataSet][str(i)] = histList_subBDT[-1].Integral(0, 1000)
			nAllPassDict[dataSet][str(i)] = histList_allBDT[-1].Integral(i, 1000)
			nAllTotDict[dataSet][str(i)] = histList_allBDT[-1].Integral(0, 1000)
	
#	stack_leadBDT = rt.THStack()
#	stack_leadBDT.Add(histList_leadBDT[0])
#	stack_leadBDT.Add(histList_leadBDT[1])
#	stack_leadBDT.Add(histList_leadBDT[2])
#	stack_leadBDT.Add(histList_leadBDT[3])
#	stack_subBDT = rt.THStack()
#	stack_subBDT.Add(histList_subBDT[0])
#	stack_subBDT.Add(histList_subBDT[1])
#	stack_subBDT.Add(histList_subBDT[2])
#	stack_subBDT.Add(histList_subBDT[3])
#	stack_allBDT = rt.THStack()
#	stack_allBDT.Add(histList_allBDT[0])
#	stack_allBDT.Add(histList_allBDT[1])
#	stack_allBDT.Add(histList_allBDT[2])
#	stack_allBDT.Add(histList_allBDT[3])
	directory = "outputs/bdtOpt/Plots2/"
	makePng(histList_leadBDT, "leadBDT_stack_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=False)
	makePng(histList_subBDT, "subBDT_stack_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=False)
	makePng(histList_allBDT, "allBDT_stack_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=False)
	makePng(histList_leadBDT, "leadBDT_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=True)
	makePng(histList_subBDT, "subBDT_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=True)
	makePng(histList_allBDT, "allBDT_"+listYear[0][3:], directory, doLeg = True, log = True, stackLines=True)
#	makeRatioStack(stack_leadBDT, histList_leadBDT[4] , "leadBDT_"+listYear[0][3:], directory, doLeg = True, log = True)
#	makeRatioStack(stack_subBDT, histList_subBDT[4] , "subBDT_"+listYear[0][3:], directory, doLeg = True, log = True)
#	makeRatioStack(stack_allBDT, histList_allBDT[4] , "allBDT_"+listYear[0][3:], directory, doLeg = True, log = True)



print("AllJets 16Pass 16tot 17Pass 17Tot 18PREpass 18PREtot 18POSTpass 18POSTtot basePass baseTotal")
for wp in cuts:
	bkgPass16 = nAllPassDict["QCD16"][str(wp)]+nAllPassDict["TTJets16"][str(wp)]+nAllPassDict["WJets16"][str(wp)]+nAllPassDict["ZJets16"][str(wp)]
	bkgTot16 = nAllTotDict["QCD16"][str(wp)]+nAllTotDict["TTJets16"][str(wp)]+nAllTotDict["WJets16"][str(wp)]+nAllTotDict["ZJets16"][str(wp)]
	bkgPass17 = nAllPassDict["QCD17"][str(wp)]+nAllPassDict["TTJets17"][str(wp)]+nAllPassDict["WJets17"][str(wp)]+nAllPassDict["ZJets17"][str(wp)]
	bkgTot17 = nAllTotDict["QCD17"][str(wp)]+nAllTotDict["TTJets17"][str(wp)]+nAllTotDict["WJets17"][str(wp)]+nAllTotDict["ZJets17"][str(wp)]
	bkgPass18PRE = nAllPassDict["QCD18PRE"][str(wp)]+nAllPassDict["TTJets18PRE"][str(wp)]+nAllPassDict["WJets18PRE"][str(wp)]+nAllPassDict["ZJets18PRE"][str(wp)]
	bkgTot18PRE = nAllTotDict["QCD18PRE"][str(wp)]+nAllTotDict["TTJets18PRE"][str(wp)]+nAllTotDict["WJets18PRE"][str(wp)]+nAllTotDict["ZJets18PRE"][str(wp)]
	bkgPass18POST = nAllPassDict["QCD18POST"][str(wp)]+nAllPassDict["TTJets18POST"][str(wp)]+nAllPassDict["WJets18POST"][str(wp)]+nAllPassDict["ZJets18POST"][str(wp)]
	bkgTot18POST = nAllTotDict["QCD18POST"][str(wp)]+nAllTotDict["TTJets18POST"][str(wp)]+nAllTotDict["WJets18POST"][str(wp)]+nAllTotDict["ZJets18POST"][str(wp)]
	print("{} {} {} {} {} {} {} {} {} {} {}".format(wp*0.001,
		bkgPass16, bkgTot16,
		bkgPass17, bkgTot17,
		bkgPass18PRE, bkgTot18PRE,
		bkgPass18POST, bkgTot18POST,
		nAllPassDict["base"][str(wp)], nAllTotDict["base"][str(wp)]))

print("--")
print("LeadJet 16Pass 16tot 17Pass 17Tot 18PREpass 18PREtot 18POSTpass 18POSTtot basePass baseTotal")
for wp in cuts:
	bkgPass16 = nLeadPassDict["QCD16"][str(wp)]+nLeadPassDict["TTJets16"][str(wp)]+nLeadPassDict["WJets16"][str(wp)]+nLeadPassDict["ZJets16"][str(wp)]
	bkgTot16 = nLeadTotDict["QCD16"][str(wp)]+nLeadTotDict["TTJets16"][str(wp)]+nLeadTotDict["WJets16"][str(wp)]+nLeadTotDict["ZJets16"][str(wp)]
	bkgPass17 = nLeadPassDict["QCD17"][str(wp)]+nLeadPassDict["TTJets17"][str(wp)]+nLeadPassDict["WJets17"][str(wp)]+nLeadPassDict["ZJets17"][str(wp)]
	bkgTot17 = nLeadTotDict["QCD17"][str(wp)]+nLeadTotDict["TTJets17"][str(wp)]+nLeadTotDict["WJets17"][str(wp)]+nLeadTotDict["ZJets17"][str(wp)]
	bkgPass18PRE = nLeadPassDict["QCD18PRE"][str(wp)]+nLeadPassDict["TTJets18PRE"][str(wp)]+nLeadPassDict["WJets18PRE"][str(wp)]+nLeadPassDict["ZJets18PRE"][str(wp)]
	bkgTot18PRE = nLeadTotDict["QCD18PRE"][str(wp)]+nLeadTotDict["TTJets18PRE"][str(wp)]+nLeadTotDict["WJets18PRE"][str(wp)]+nLeadTotDict["ZJets18PRE"][str(wp)]
	bkgPass18POST = nLeadPassDict["QCD18POST"][str(wp)]+nLeadPassDict["TTJets18POST"][str(wp)]+nLeadPassDict["WJets18POST"][str(wp)]+nLeadPassDict["ZJets18POST"][str(wp)]
	bkgTot18POST = nLeadTotDict["QCD18POST"][str(wp)]+nLeadTotDict["TTJets18POST"][str(wp)]+nLeadTotDict["WJets18POST"][str(wp)]+nLeadTotDict["ZJets18POST"][str(wp)]
	print("{} {} {} {} {} {} {} {} {} {} {}".format(wp*0.001,
		bkgPass16, bkgTot16,
		bkgPass17, bkgTot17,
		bkgPass18PRE, bkgTot18PRE,
		bkgPass18POST, bkgTot18POST,
		nLeadPassDict["base"][str(wp)], nLeadTotDict["base"][str(wp)]))

print("--")
print("SubJets 16Pass 16tot 17Pass 17Tot 18PREpass 18PREtot 18POSTpass 18POSTtot basePass baseTotal")
for wp in cuts:
	bkgPass16 = nSubPassDict["QCD16"][str(wp)]+nSubPassDict["TTJets16"][str(wp)]+nSubPassDict["WJets16"][str(wp)]+nSubPassDict["ZJets16"][str(wp)]
	bkgTot16 = nSubTotDict["QCD16"][str(wp)]+nSubTotDict["TTJets16"][str(wp)]+nSubTotDict["WJets16"][str(wp)]+nSubTotDict["ZJets16"][str(wp)]
	bkgPass17 = nSubPassDict["QCD17"][str(wp)]+nSubPassDict["TTJets17"][str(wp)]+nSubPassDict["WJets17"][str(wp)]+nSubPassDict["ZJets17"][str(wp)]
	bkgTot17 = nSubTotDict["QCD17"][str(wp)]+nSubTotDict["TTJets17"][str(wp)]+nSubTotDict["WJets17"][str(wp)]+nSubTotDict["ZJets17"][str(wp)]
	bkgPass18PRE = nSubPassDict["QCD18PRE"][str(wp)]+nSubPassDict["TTJets18PRE"][str(wp)]+nSubPassDict["WJets18PRE"][str(wp)]+nSubPassDict["ZJets18PRE"][str(wp)]
	bkgTot18PRE = nSubTotDict["QCD18PRE"][str(wp)]+nSubTotDict["TTJets18PRE"][str(wp)]+nSubTotDict["WJets18PRE"][str(wp)]+nSubTotDict["ZJets18PRE"][str(wp)]
	bkgPass18POST = nSubPassDict["QCD18POST"][str(wp)]+nSubPassDict["TTJets18POST"][str(wp)]+nSubPassDict["WJets18POST"][str(wp)]+nSubPassDict["ZJets18POST"][str(wp)]
	bkgTot18POST = nSubTotDict["QCD18POST"][str(wp)]+nSubTotDict["TTJets18POST"][str(wp)]+nSubTotDict["WJets18POST"][str(wp)]+nSubTotDict["ZJets18POST"][str(wp)]
	print("{} {} {} {} {} {} {} {} {} {} {}".format(wp*0.001,
		bkgPass16, bkgTot16,
		bkgPass17, bkgTot17,
		bkgPass18PRE, bkgTot18PRE,
		bkgPass18POST, bkgTot18POST,
		nSubPassDict["base"][str(wp)], nSubTotDict["base"][str(wp)]))


#for listYear in sigList:
#	histList_leadBDT = []
#	histList_subBDT = []
#	histList_allBDT = []
#	for dataSet in listYear:
#		_file1 = rt.TFile.Open("outputs/bdtOpt/"+dataSet+"/bdt_opti.root","read")
#		histList_leadBDT.append(_file1.Get("hist_leadBDT").Clone("leadBDT_"+dataSet))
#		histList_subBDT.append(_file1.Get("hist_subBDT").Clone("subBDT_"+dataSet))
#		histList_allBDT.append(_file1.Get("hist_allBDT").Clone("allBDT_"+dataSet))
#		histList_leadBDT[-1].SetDirectory(0)
#		histList_subBDT[-1].SetDirectory(0)
#		histList_allBDT[-1].SetDirectory(0)
##		for i in cuts:
##			nPassLead = histList_leadBDT[-1].Integral(i, 1000)
##			nPassSub = histList_subBDT[-1].Integral(i, 1000)
##			nPassAll = histList_allBDT[-1].Integral(i, 1000)
##			nTotLead = histList_leadBDT[-1].Integral(0, 1000)
##			nTotSub = histList_subBDT[-1].Integral(0, 1000)
##			nTotAll = histList_allBDT[-1].Integral(0, 1000)
##			print("{} {} {} {} {} {} {} {}".format(dataSet, float(i)/1000., nPassLead, nTotLead, nPassSub, nTotSub, nPassAll, nTotAll))
#	directory = "outputs/bdtOpt/Plots/"
#	makePng(histList_leadBDT, "leadBDT_"+listYear[0][0], directory, doLeg = True, log = True)
#	makePng(histList_subBDT, "subBDT_"+listYear[0][0], directory, doLeg = True, log = True)
#	makePng(histList_allBDT, "allBDT_"+listYear[0][0], directory, doLeg = True, log = True)




















