import ROOT as rt
from array import array
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(0)
#run from FAnalysis folder with python macros/scalefact/draw_sfPlots.py


list16 = ["QCD16","TTJets16","WJets16","ZJets16", "Data16"]
list17 = ["QCD17","TTJets17","WJets17","ZJets17", "Data17"]
list18PRE = ["QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE", "Data18PRE"]
list18POST = ["QCD18POST","TTJets18POST","WJets18POST","ZJets18POST", "Data18POST"]

listZ = ["z15", "z20", "z25", "base", "z35", "z40","z45"]
listD = ["d10","base","d40","d60","d80","d00"]
listR = ["r02","base","r05","r07","r08"]
listA = ["adl", "base", "adh"]

colorDict = {"QC":602, "TT":798,"WJ":801,"ZJ":881,"Bk":2,"Da":1}

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

def makeRatio(denomR, numerR, name, direc, doLeg = True, log = False):
	# C++ Author: Olivier Couet, adapted to python by Colin Fallon
	# Define the Canvas
	c = rt.TCanvas("c", "canvas", 800, 800)

	# Upper plot will be in pad1
	pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0) # Upper and lower plot are joined
	pad1.SetGridx()         # Vertical grid
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	#pad1.SetLogy()

	#denom.Rebin(denom.GetNbinsX()/50)
	#numer.Rebin(numer.GetNbinsX()/50)
	if denomR.GetNbinsX() > 50:
		# C++ from kPedro, adapted to python by cFallon
		# https://github.com/kpedro88/Analysis/blob/ea74a1bc3a410ae0ef8e05bc831ab4eb32817b91/scripts/plotTrigEff.C#L146-L158
		minEvents = denom.Integral()/20.
		newbins = array("d",[])
		#print(newbins, type(newbins))
		smallSum = 0
		newbins.append(denomR.GetBinLowEdge(1))
		for iBin in range(denomR.GetNbinsX()):
			smallSum += denomR.GetBinContent(iBin)
			if(smallSum>minEvents or iBin==denomR.GetNbinsX()):
				newbins.append(denomR.GetBinLowEdge(iBin+1))
				smallSum = 0
		#for binedge in newbins:
		#	print(binedge)
		#return
		denom = denomR.Rebin(len(newbins)-1,denomR.GetName()+"_rebin",newbins)
		numer = numerR.Rebin(len(newbins)-1,numerR.GetName()+"_rebin",newbins)
		print(name, newbins)
	else:
		denom = denomR.Clone()
		numer = numerR.Clone()
	denom.SetTitle("All Jets (denom)")
	numer.SetTitle("Passing Jets (numer)")
	denom.SetStats(0)       # No statistics on upper plot
	denom.SetMinimum(numer.GetMinimum()*0.8)
	denom.GetXaxis().SetRangeUser(200,3000)
	numer.GetXaxis().SetRangeUser(200,3000)
	denom.Draw("E1 X+")
	numer.Draw("same E1 X+")
	if doLeg:
		pad1.BuildLegend(0.75,0.75,0.95,0.95)
	denom.SetTitle("")
	pad1.Update()

	# Do not draw the Y axis label on the upper plot and redraw a small
	# axis instead, in order to avoid the first label (0) to be clipped.
	#numer.GetYaxis().SetLabelSize(0.)
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
	#h3 = numer.Clone("h3")
	h3 = rt.TGraphAsymmErrors(numer, denom)
	h3.SetLineColor(rt.kBlack)
	h3.SetMinimum(0.00)  # Define Y ..
	h3.SetMaximum(0.30) # .. range
	#h3.Sumw2()
	#h3.SetStats(0)      # No statistics on lower plot
	#h3.Divide(denom)
	h3.SetMarkerStyle(21)
	h3.Draw("")       # Draw the ratio plot

	# h1 settings
	numer.SetLineColor(rt.kBlue+1)
	numer.SetLineWidth(2)

	# Y axis h1 plot settings
	numer.GetYaxis().SetTitleSize(20)
	numer.GetYaxis().SetTitleFont(43)
	numer.GetYaxis().SetTitleOffset(1.55)

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
	h3.GetXaxis().SetTitle("Jet Pt")
	h3.GetXaxis().SetTitleSize(20)
	h3.GetXaxis().SetTitleFont(43)
	h3.GetXaxis().SetTitleOffset(4.)
	h3.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
	h3.GetXaxis().SetLabelSize(15)
	#h3.GetXaxis().SetLimits(200,3000)
	x = rt.Double()
	y = rt.Double()
#	print(name)
#	print(denomR.GetNbinsX(), "denomR bins")
#	print(numerR.GetNbinsX(), "numerR bins")
#	print("iBin x y HistoLow")
#	for iBin in range(h3.GetN()):
#		h3.GetPoint(iBin, x, y)
#		print("{} {} {} {}".format(iBin, x,y, int(denomR.GetXaxis().GetBinLowEdge(iBin))))
	pad2.Update()

	#save as .png
	c.SaveAs(direc+name+"_ratio.png")
	
bkgList = [list16, list17, list18PRE, list18POST]
sigList = [listZ, listD, listR, listA]

histDict_numer_pt0 = {}
histDict_numer_pt1 = {}
histDict_numer_ptb = {}
histDict_denom_pt0 = {}
histDict_denom_pt1 = {}
histDict_denom_ptb = {}


for yearList in bkgList:
	for dataSet in yearList:
		#print(dataSet)
		_file1 = rt.TFile.Open("outputs/scalefact/"+dataSet+"/scalefact_jetpt_CRDE.root","read")
		histDict_numer_pt0[dataSet] = _file1.Get("JetsAK8[0].Pt()_numer_"+dataSet)
		histDict_numer_pt0[dataSet].SetDirectory(0)
		histDict_numer_pt1[dataSet] = _file1.Get("JetsAK8[1].Pt()_numer_"+dataSet)
		histDict_numer_pt1[dataSet].SetDirectory(0)
		histDict_numer_ptb[dataSet] = _file1.Get("JetsAK8.Pt()_numer_"+dataSet)
		histDict_numer_ptb[dataSet].SetDirectory(0)
		
		histDict_denom_pt0[dataSet] = _file1.Get("JetsAK8[0].Pt()_denom_"+dataSet)
		histDict_denom_pt0[dataSet].SetDirectory(0)
		histDict_denom_pt1[dataSet] = _file1.Get("JetsAK8[1].Pt()_denom_"+dataSet)
		histDict_denom_pt1[dataSet].SetDirectory(0)
		histDict_denom_ptb[dataSet] = _file1.Get("JetsAK8.Pt()_denom_"+dataSet)
		histDict_denom_ptb[dataSet].SetDirectory(0)
	
#		direct = "outputs/scalefact/plots_test/"
#		makePng([listOfNumer[0],listOfDenom[0]], "leadJetPt", direct, doLeg = True, log = False, nBinsX = 50, stackLines = True)
#		makePng([listOfNumer[1],listOfDenom[1]], "subJetPt", direct, doLeg = True, log = False, nBinsX = 50, stackLines = True)
#		makePng([listOfNumer[2],listOfDenom[2]], "bothJetPt", direct, doLeg = True, log = False, nBinsX = 50, stackLines = True)
#		makeRatio(listOfDenom[0],listOfNumer[0], dataSet+"_leadJetPt", direct, doLeg = True, log = False)
#		makeRatio(listOfDenom[1],listOfNumer[1], dataSet+"_subJetPt", direct, doLeg = True, log = False)
#		makeRatio(listOfDenom[2],listOfNumer[2], dataSet+"_bothJetPt", direct, doLeg = True, log = False)

histDict_numer_sumBkgs_pt0 = {"16":histDict_numer_pt0["QCD16"].Clone("Bkg16_numer_pt0"),"17":histDict_numer_pt0["QCD17"].Clone("Bkg17_numer_pt0"),"18PRE":histDict_numer_pt0["QCD18PRE"].Clone("Bkg18PRE_numer_pt0"),"18POST":histDict_numer_pt0["QCD18POST"].Clone("Bkg18POST_numer_pt0")}
histDict_numer_sumBkgs_pt1 = {"16":histDict_numer_pt1["QCD16"].Clone("Bkg16_numer_pt1"),"17":histDict_numer_pt1["QCD17"].Clone("Bkg17_numer_pt1"),"18PRE":histDict_numer_pt1["QCD18PRE"].Clone("Bkg18PRE_numer_pt1"),"18POST":histDict_numer_pt1["QCD18POST"].Clone("Bkg18POST_numer_pt1")}
histDict_numer_sumBkgs_ptb = {"16":histDict_numer_ptb["QCD16"].Clone("Bkg16_numer_ptb"),"17":histDict_numer_ptb["QCD17"].Clone("Bkg17_numer_ptb"),"18PRE":histDict_numer_ptb["QCD18PRE"].Clone("Bkg18PRE_numer_ptb"),"18POST":histDict_numer_ptb["QCD18POST"].Clone("Bkg18POST_numer_ptb")}
histDict_denom_sumBkgs_pt0 = {"16":histDict_denom_pt0["QCD16"].Clone("Bkg16_denom_pt0"),"17":histDict_denom_pt0["QCD17"].Clone("Bkg17_denom_pt0"),"18PRE":histDict_denom_pt0["QCD18PRE"].Clone("Bkg18PRE_denom_pt0"),"18POST":histDict_denom_pt0["QCD18POST"].Clone("Bkg18POST_denom_pt0")}
histDict_denom_sumBkgs_pt1 = {"16":histDict_denom_pt1["QCD16"].Clone("Bkg16_denom_pt1"),"17":histDict_denom_pt1["QCD17"].Clone("Bkg17_denom_pt1"),"18PRE":histDict_denom_pt1["QCD18PRE"].Clone("Bkg18PRE_denom_pt1"),"18POST":histDict_denom_pt1["QCD18POST"].Clone("Bkg18POST_denom_pt1")}
histDict_denom_sumBkgs_ptb = {"16":histDict_denom_ptb["QCD16"].Clone("Bkg16_denom_ptb"),"17":histDict_denom_ptb["QCD17"].Clone("Bkg17_denom_ptb"),"18PRE":histDict_denom_ptb["QCD18PRE"].Clone("Bkg18PRE_denom_ptb"),"18POST":histDict_denom_ptb["QCD18POST"].Clone("Bkg18POST_denom_ptb")}

dictOfRatios_pt0 = {}
dictOfRatios_pt1 = {}
dictOfRatios_ptb = {}

for yearKey in ["16","17","18PRE","18POST"]:
	histDict_numer_sumBkgs_pt0[yearKey].SetTitle("Bkg"+yearKey)
	histDict_numer_sumBkgs_pt1[yearKey].SetTitle("Bkg"+yearKey)
	histDict_numer_sumBkgs_ptb[yearKey].SetTitle("Bkg"+yearKey)
	histDict_denom_sumBkgs_pt0[yearKey].SetTitle("Bkg"+yearKey)
	histDict_denom_sumBkgs_pt1[yearKey].SetTitle("Bkg"+yearKey)
	histDict_denom_sumBkgs_ptb[yearKey].SetTitle("Bkg"+yearKey)
	for bkgKey in ["TTJets","WJets","ZJets"]:
		histDict_numer_sumBkgs_pt0[yearKey].Add(histDict_numer_pt0[bkgKey+yearKey])
		histDict_numer_sumBkgs_pt1[yearKey].Add(histDict_numer_pt1[bkgKey+yearKey])
		histDict_numer_sumBkgs_ptb[yearKey].Add(histDict_numer_ptb[bkgKey+yearKey])
		histDict_denom_sumBkgs_pt0[yearKey].Add(histDict_denom_pt0[bkgKey+yearKey])
		histDict_denom_sumBkgs_pt1[yearKey].Add(histDict_denom_pt1[bkgKey+yearKey])
		histDict_denom_sumBkgs_ptb[yearKey].Add(histDict_denom_ptb[bkgKey+yearKey])
	direct = "outputs/scalefact/plots_test2/"
	#makeRatio(histDict_denom_sumBkgs_pt0[yearKey],histDict_numer_sumBkgs_pt0[yearKey], yearKey+"_leadJetPt", direct, doLeg = True, log = False)
	#makeRatio(histDict_denom_sumBkgs_pt1[yearKey],histDict_numer_sumBkgs_pt1[yearKey], yearKey+"_subJetPt", direct, doLeg = True, log = False)
	#makeRatio(histDict_denom_sumBkgs_ptb[yearKey],histDict_numer_sumBkgs_ptb[yearKey], yearKey+"_bothJetPt", direct, doLeg = True, log = False)
	#print("{} {} {} {} {}".format(histDict_numer_pt0["QCD"+yearKey].Integral(),histDict_numer_pt0["TTJets"+yearKey].Integral(),histDict_numer_pt0["WJets"+yearKey].Integral(),histDict_numer_pt0["ZJets"+yearKey].Integral(),histDict_numer_sumBkgs_pt0[yearKey].Integral()))
	dictOfRatios_pt0["Data"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_pt0["Data"+yearKey], histDict_denom_pt0["Data"+yearKey])
	dictOfRatios_pt1["Data"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_pt1["Data"+yearKey], histDict_denom_pt1["Data"+yearKey])
	dictOfRatios_ptb["Data"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_ptb["Data"+yearKey], histDict_denom_ptb["Data"+yearKey])
	dictOfRatios_pt0["Bkg"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_sumBkgs_pt0[yearKey], histDict_denom_sumBkgs_pt0[yearKey])
	dictOfRatios_pt1["Bkg"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_sumBkgs_pt1[yearKey], histDict_denom_sumBkgs_pt1[yearKey])
	dictOfRatios_ptb["Bkg"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_sumBkgs_ptb[yearKey], histDict_denom_sumBkgs_ptb[yearKey])
	for bkgKey in ["QCD","TTJets","WJets","ZJets"]:
		dictOfRatios_pt0[bkgKey+yearKey] = rt.TGraphAsymmErrors(histDict_numer_pt0[bkgKey+yearKey], histDict_denom_pt0[bkgKey+yearKey])
		dictOfRatios_pt1[bkgKey+yearKey] = rt.TGraphAsymmErrors(histDict_numer_pt1[bkgKey+yearKey], histDict_denom_pt1[bkgKey+yearKey])
		dictOfRatios_ptb[bkgKey+yearKey] = rt.TGraphAsymmErrors(histDict_numer_ptb[bkgKey+yearKey], histDict_denom_ptb[bkgKey+yearKey])


dictOfSF_pt0 = {}
dictOfSF_pt1 = {}
dictOfSF_ptb = {}

for yearKey in ["16", "17","18PRE","18POST"]:
	dictOfSF_pt0[yearKey] = dictOfRatios_pt0["Data"+yearKey].Clone("ScaleFactors_pt0_"+yearKey)
	dictOfSF_pt1[yearKey] = dictOfRatios_pt1["Data"+yearKey].Clone("ScaleFactors_pt1_"+yearKey)
	dictOfSF_ptb[yearKey] = dictOfRatios_ptb["Data"+yearKey].Clone("ScaleFactors_ptb_"+yearKey)
	for iPoint in range(dictOfSF_pt0[yearKey].GetN()):
		xT_pt0 = rt.Double()
		yT_pt0 = rt.Double()
		xB_pt0 = rt.Double()
		yB_pt0 = rt.Double()
		xT_pt1 = rt.Double()
		yT_pt1 = rt.Double()
		xB_pt1 = rt.Double()
		yB_pt1 = rt.Double()
		xT_ptb = rt.Double()
		yT_ptb = rt.Double()
		xB_ptb = rt.Double()
		yB_ptb = rt.Double()
		dictOfRatios_pt0["Data"+yearKey].GetPoint(iPoint, xT_pt0, yT_pt0)
		dictOfRatios_pt0["Bkg"+yearKey].GetPoint(iPoint, xB_pt0, yB_pt0)
		dictOfRatios_pt1["Data"+yearKey].GetPoint(iPoint, xT_pt1, yT_pt1)
		dictOfRatios_pt1["Bkg"+yearKey].GetPoint(iPoint, xB_pt1, yB_pt1)
		dictOfRatios_ptb["Data"+yearKey].GetPoint(iPoint, xT_ptb, yT_ptb)
		dictOfRatios_ptb["Bkg"+yearKey].GetPoint(iPoint, xB_ptb, yB_ptb)
		try:
			dictOfSF_pt0[yearKey].SetPoint(iPoint,xT_pt0, yT_pt0/yB_pt0)
		except ZeroDivisionError:
			dictOfSF_pt0[yearKey].SetPoint(iPoint,xT_pt0, 1.)
		try:
			dictOfSF_pt1[yearKey].SetPoint(iPoint,xT_pt1, yT_pt1/yB_pt1)
		except ZeroDivisionError:
			dictOfSF_pt1[yearKey].SetPoint(iPoint,xT_pt1, 1.)
		try:
			dictOfSF_ptb[yearKey].SetPoint(iPoint,xT_ptb, yT_ptb/yB_ptb)
		except ZeroDivisionError:
			dictOfSF_ptb[yearKey].SetPoint(iPoint,xT_ptb, 1.)
		print(xT_pt0, xB_pt0, xT_pt1, xB_pt1, xT_ptb, xB_ptb)

for key in dictOfRatios_pt0:
	print(key)
	dictOfRatios_pt0[key].SetLineColor(colorDict[key[:2]])
	dictOfRatios_pt1[key].SetLineColor(colorDict[key[:2]])
	dictOfRatios_ptb[key].SetLineColor(colorDict[key[:2]])
c = rt.TCanvas()
for yearKey in ["16","17","18PRE","18POST"]:
	dictOfRatios_pt0["QCD"+yearKey].SetMinimum(0)
	dictOfRatios_pt0["QCD"+yearKey].SetMaximum(1)
	dictOfRatios_pt0["QCD"+yearKey].Draw("")
	for setKey in ["TTJets","WJets","ZJets","Bkg","Data"]:
		dictOfRatios_pt0[setKey+yearKey].Draw("same")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" Leading Jet Pt")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_pt0_ratio.png")
	c.Clear()

for yearKey in ["16","17","18PRE","18POST"]:
	dictOfRatios_pt1["QCD"+yearKey].SetMinimum(0)
	dictOfRatios_pt1["QCD"+yearKey].SetMaximum(1)
	dictOfRatios_pt1["QCD"+yearKey].Draw("")
	for setKey in ["TTJets","WJets","ZJets","Bkg","Data"]:
		dictOfRatios_pt1[setKey+yearKey].Draw("same")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" Subleading Jet Pt")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_pt1_ratio.png")
	c.Clear()

for yearKey in ["16","17","18PRE","18POST"]:
	dictOfRatios_ptb["QCD"+yearKey].SetMinimum(0)
	dictOfRatios_ptb["QCD"+yearKey].SetMaximum(1)
	dictOfRatios_ptb["QCD"+yearKey].Draw("")
	for setKey in ["TTJets","WJets","ZJets","Bkg","Data"]:
		dictOfRatios_ptb[setKey+yearKey].Draw("same")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" Both Leading Jets' Pt")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_ptb_ratio.png")
	c.Clear()

for yearKey in ["16","17","18PRE","18POST"]:
	dictOfSF_pt0[yearKey].SetMinimum(-0.1)
	dictOfSF_pt0[yearKey].SetMaximum(2)
	dictOfSF_pt0[yearKey].SetLineColor(1)
	dictOfSF_pt0[yearKey].SetTitle("Leading Jet")
	dictOfSF_pt0[yearKey].Draw("")
	dictOfSF_pt1[yearKey].SetLineColor(2)
	dictOfSF_pt1[yearKey].SetTitle("Subleading Jet")
	dictOfSF_pt1[yearKey].Draw("same")
	dictOfSF_ptb[yearKey].SetLineColor(3)
	dictOfSF_ptb[yearKey].SetTitle("Both Jets")
	dictOfSF_ptb[yearKey].Draw("same")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" ScaleFactors")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_SF_ratio.png")
	c.Clear()

