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

histDict_numer = {}
histDict_denom = {}


for yearList in bkgList:
	for dataSet in yearList:
		#print(dataSet)
		_file1 = rt.TFile.Open("outputs/scalefact/"+dataSet+"/scalefact_deltaeta_CRDE.root","read")
		histDict_numer[dataSet] = _file1.Get("deltaEta_numer_"+dataSet)
		histDict_numer[dataSet].SetDirectory(0)
		
		histDict_denom[dataSet] = _file1.Get("deltaEta_denom_"+dataSet)
		histDict_denom[dataSet].SetDirectory(0)

histDict_numer_sumBkgs = {"16":histDict_numer["QCD16"].Clone("Bkg16_numer"),"17":histDict_numer["QCD17"].Clone("Bkg17_numer"),"18PRE":histDict_numer["QCD18PRE"].Clone("Bkg18PRE_numer"),"18POST":histDict_numer["QCD18POST"].Clone("Bkg18POST_numer")}
histDict_denom_sumBkgs = {"16":histDict_denom["QCD16"].Clone("Bkg16_denom"),"17":histDict_denom["QCD17"].Clone("Bkg17_denom"),"18PRE":histDict_denom["QCD18PRE"].Clone("Bkg18PRE_denom"),"18POST":histDict_denom["QCD18POST"].Clone("Bkg18POST_denom")}

dictOfRatios = {}

for yearKey in ["16","17","18PRE","18POST"]:
	histDict_numer_sumBkgs[yearKey].SetTitle("Bkg"+yearKey)
	histDict_denom_sumBkgs[yearKey].SetTitle("Bkg"+yearKey)
	for bkgKey in ["TTJets","WJets","ZJets"]:
		histDict_numer_sumBkgs[yearKey].Add(histDict_numer[bkgKey+yearKey])
		histDict_denom_sumBkgs[yearKey].Add(histDict_denom[bkgKey+yearKey])
	direct = "outputs/scalefact/plots_deltaEta/"
	dictOfRatios["Data"+yearKey] = rt.TGraphAsymmErrors(histDict_numer["Data"+yearKey], histDict_denom["Data"+yearKey])
	dictOfRatios["Bkg"+yearKey] = rt.TGraphAsymmErrors(histDict_numer_sumBkgs[yearKey], histDict_denom_sumBkgs[yearKey])
	for bkgKey in ["QCD","TTJets","WJets","ZJets"]:
		dictOfRatios[bkgKey+yearKey] = rt.TGraphAsymmErrors(histDict_numer[bkgKey+yearKey], histDict_denom[bkgKey+yearKey])


dictOfSF = {}

for yearKey in ["16", "17","18PRE","18POST"]:
	dictOfSF[yearKey] = dictOfRatios["Data"+yearKey].Clone("ScaleFactors_"+yearKey)
	for iPoint in range(dictOfSF[yearKey].GetN()):
		xT = rt.Double()
		yT = rt.Double()
		xB = rt.Double()
		yB = rt.Double()
		dictOfRatios["Data"+yearKey].GetPoint(iPoint, xT, yT)
		dictOfRatios["Bkg"+yearKey].GetPoint(iPoint, xB, yB)
		try:
			dictOfSF[yearKey].SetPoint(iPoint,xT, yT/yB)
		except ZeroDivisionError:
			dictOfSF[yearKey].SetPoint(iPoint,xT, 1.)
		print(xT, xB)

for key in dictOfRatios:
	print(key)
	dictOfRatios[key].SetLineColor(colorDict[key[:2]])
c = rt.TCanvas()
for yearKey in ["16","17","18PRE","18POST"]:
	dictOfRatios["QCD"+yearKey].SetMinimum(0)
	dictOfRatios["QCD"+yearKey].SetMaximum(0.2)
	dictOfRatios["QCD"+yearKey].Draw("")
	for setKey in ["TTJets","WJets","ZJets","Bkg","Data"]:
		dictOfRatios[setKey+yearKey].Draw("same")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" #Delta#eta(j_{1}, j_{2})")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_ratio.png")
	c.Clear()

for yearKey in ["16","17","18PRE","18POST"]:
	dictOfSF[yearKey].SetMinimum(0.8)
	dictOfSF[yearKey].SetMaximum(2.0)
	dictOfSF[yearKey].SetLineColor(1)
	dictOfSF[yearKey].SetTitle("#Delta#eta({j_1, j_2}")
	dictOfSF[yearKey].Draw("")
	title = rt.TPaveText(0.3,0.91,0.7,0.99,"brNDC")
	title.AddText(yearKey+" ScaleFactors")
	title.Draw()
	c.BuildLegend()
	c.SaveAs(direct+yearKey+"_SF_ratio.png")
	c.Clear()

