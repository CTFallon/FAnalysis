import ROOT as rt
from array import array
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptTitle(0)
#run from FAnalysis folder with python macros/sfABCD/draw_sfABCD.py


list16 = ["QCD16","TTJets16","WJets16","ZJets16"]
list17 = ["QCD17","TTJets17","WJets17","ZJets17"]
list18PRE = ["QCD18PRE","TTJets18PRE","WJets18PRE","ZJets18PRE"]
list18POST = ["QCD18POST","TTJets18POST","WJets18POST","ZJets18POST"]

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
	for h in LoHi:
		print(h.GetName(),h.Integral())
		LoH.append(h.Clone())
	LoH.sort(key=lambda x: x.Integral())
	if stackLines:
		for i in range(len(LoH)):
			LoH[i].SetLineColor(colorDict[LoH[i].GetName()[:2]])
			#if LoH[i].GetNbinsX() != nBinsX:
			#	LoH[i].Rebin(LoH[i].GetNbinsX()/nBinsX)
			stack.Add(LoH[i])
		#base = LoHi[-1].Clone()
		#if base.GetNbinsX() != nBinsX:
		#	base.Rebin(base.GetNbinsX()/nBinsX)
		#base.SetLineColor(rt.kRed)
		stack.Draw("nostack hist")
		#base.Draw("same")
	else:
		for i in range(len(LoH)):
			LoH[i].SetLineColor(colorDict[LoH[i].GetName()[:2]])
			LoH[i].SetFillColor(colorDict[LoH[i].GetName()[:2]])
			LoH[i].SetFillStyle(1001)
			#if LoH[i].GetNbinsX() != nBinsX:
			#	LoH[i].Rebin(LoH[i].GetNbinsX()/nBinsX)
			stack.Add(LoH[i])
		#base = LoHi[-1].Clone()
		#if base.GetNbinsX() != nBinsX:
		#	base.Rebin(base.GetNbinsX()/nBinsX)
		#base.SetLineColor(rt.kRed)
		stack.Draw("hist")
		#base.Draw("same")
	stack.GetXaxis().SetTitle(LoH[0].GetXaxis().GetTitle())
	stack.GetYaxis().SetTitle("Count")
	stack.SetTitle(name)
	stack.SetMaximum(10**5)
	stack.SetMinimum(10**-2)
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.7,0.6,0.9)
	c1.SaveAs(direc+name+".png")
	normStack = rt.THStack()
	for i in range(len(LoH)):
		if LoH[i].Integral() != 0:
			LoH[i].Scale(1/LoH[i].Integral())
		normStack.Add(LoH[i])
	normStack.Draw("hist")
	normStack.GetXaxis().SetTitle(LoH[0].GetXaxis().GetTitle())
	normStack.GetYaxis().SetTitle("a.u.")
	c1.Modified()
	if doLeg:
		c1.BuildLegend(0.4,0.7,0.6,0.9)
	c1.SaveAs(direc+name+"_norm.png")

def makeRatio(histCR, histSR , name, direc, doLeg = True, log = False, doNorm = True):
	# C++ Author: Olivier Couet, adapted to python by Colin Fallon
	# Define the Canvas
	c = rt.TCanvas("c", "canvas", 800, 800)
	print(histCR.GetName(), histCR.Integral())
	print(histSR.GetName(), histSR.Integral())
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
	histSR.SetTitle("Low #Delta#eta")
	histCR.SetTitle("High #Delta#eta")
	histSR.SetLineColor(1)
	histCR.SetLineColor(2)
	if doNorm:
		try:
			histCR.Scale(1/histCR.Integral())
			histSR.Scale(1/histSR.Integral())
		except ZeroDivisionError:
			pass
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
	c.SaveAs(direc+name+"_ratio.png")
	
bkgList = [list16, list17, list18PRE, list18POST]
sigList = [listZ, listD, listR, listA]


# Getting a 1D histgoram with the correct events:
# The 3D histos are X:Y:Z -> nSVJtags:RT:variable
#		nSVJtags has 3 bins with edges [0,1,2,3]
#		RT has 100 same-width bins with range [0,1] #this might change to 0.15 to 0.8
#		ingeneral, varaible has 1k bins with various ranges
# so, to plot the variable for Region A1:
# 	take the histgram (hA1) from _dEtaLow for the varaible you want
#	and call hA1.ProjectionZ("_A1",2,2,25,100)
# 		where 	"_A1" is appended to the name
#				2,2 is the limits of the X axis (nSVJ) (in bin index)
#				25,100 is the limits of the Y axis (RT) (in bin index)

# Regions unprime and prime are _dEtaLow and _dEtaHigh, respectivly.
# regions A and B are RT > .25 (iBin 26) and 0.15 < RT < 0.25 (iBins 16 and 25)
# regions 0, 1, and 2 are numSVJtags == 0, == 1, == 2, (iBins 1, 2, 3)

#Lets start with just getting nEvents from various places...
direc = "outputs/sfABCD/plots_sameMTcut_correct/"

dictOfBackgroundTotals = {}

#	'LeadJetPt'
#	'SubLeadJetPt'
#	'BothJetPt'
#	'MT'
#	'MET'
#	'METPhi'
#	'DeltaEta'
#	'nSVJTag'
#	'RT'

#for var1 in ["LeadJetPt", "SubLeadJetPt","BothJetPt","MT","MET","METPhi","DeltaEta","nSVJTag","RT"]:
for var1 in ["METPhi"]: # use code JFURND to find lines for tesing METPhi problem
#	var1 = "MT"
#	stackDict_nSVJ = {}
#	stackDict_RT = {}
#	stackDict_dEta = {}
	stackDict_var1 = {}
	listOfRegions = ["A0","A1","A2","B0","B1","B2"]
	for region in listOfRegions:
#		stackDict_nSVJ[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_nSVJ",region+"_nSVJ;nSVJ;events")
#		stackDict_RT[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_RT",region+"_RT;RT;events")
#		stackDict_dEta[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_dEta",region+"_dEta;dEta;events")
		stackDict_var1[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_"+var1,region+"_"+var1+";"+var1+";events")
		dictOfBackgroundTotals["high_"+var1+"_"+region] = rt.TH1F("high_METPhi_"+region,"high_METPhi_"+region,20, -rt.TMath.Pi(),rt.TMath.Pi())

	for bkg in ["ZJets16","WJets16","TTJets16","QCD16"]:
		for region in listOfRegions:
			_file1 = rt.TFile.Open("outputs/sfABCD/"+bkg+"/sfABCD_test2_dEtaHigh.root","read")

#			temp_nSVJ =  _file1.Get(bkg+"_nSVJTag_"+region)
#			temp_RT = _file1.Get(bkg+"_RT_"+region)
#			temp_dEta = _file1.Get(bkg+"_DeltaEta_"+region)
			temp_var1 = _file1.Get(bkg+"_"+var1+"_"+region)

#			temp_nSVJ.SetDirectory(0)
#			temp_RT.SetDirectory(0)
#			temp_dEta.SetDirectory(0)
			temp_var1.SetDirectory(0)
			print(var1, region, bkg, temp_var1.GetName(), temp_var1.Integral()) #JFURND
		
			_file1.Close()
			if temp_var1.GetNbinsX() == 1000:
				temp_var1.Rebin(50)
	

#			stackDict_nSVJ[region][bkg] = temp_nSVJ.Clone()
#			stackDict_RT[region][bkg] = temp_RT.Clone()
#			stackDict_dEta[region][bkg] = temp_dEta.Clone()
			stackDict_var1[region][bkg] = temp_var1.Clone()
			dictOfBackgroundTotals["high_"+var1+"_"+region].Add(stackDict_var1[region][bkg])

#			stackDict_nSVJ[region][bkg].SetTitle(region[0]+"'"+region[1]+stackDict_nSVJ[region][bkg].GetTitle()[2:])
#			stackDict_RT[region][bkg].SetTitle(region[0]+"'"+region[1]+stackDict_RT[region][bkg].GetTitle()[2:])
#			stackDict_dEta[region][bkg].SetTitle(region[0]+"'"+region[1]+stackDict_dEta[region][bkg].GetTitle()[2:])
			stackDict_var1[region][bkg].SetTitle(region[0]+"'"+region[1]+stackDict_var1[region][bkg].GetTitle()[2:])

	for region in listOfRegions:
#		makePng([stackDict_nSVJ[region]["QCD16"],stackDict_nSVJ[region]["TTJets16"],stackDict_nSVJ[region]["WJets16"],stackDict_nSVJ[region]["ZJets16"]], region+"_dEtaHigh_nSVJ", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
#		makePng([stackDict_RT[region]["QCD16"],stackDict_RT[region]["TTJets16"],stackDict_RT[region]["WJets16"],stackDict_RT[region]["ZJets16"]], region+"_dEtaHigh_RT", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
#		makePng([stackDict_dEta[region]["QCD16"],stackDict_dEta[region]["TTJets16"],stackDict_dEta[region]["WJets16"],stackDict_dEta[region]["ZJets16"]], region+"_dEtaHigh_dEta", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
		makePng([stackDict_var1[region]["QCD16"],stackDict_var1[region]["TTJets16"],stackDict_var1[region]["WJets16"],stackDict_var1[region]["ZJets16"]], region+"_dEtaHigh_"+var1, direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)	

#	stackDict_nSVJ = {}
#	stackDict_RT = {}
#	stackDict_dEta = {}
	stackDict_var1l = {}
	for region in listOfRegions:
#		stackDict_nSVJ[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_nSVJ",region+"_nSVJ;nSVJ;events")
#		stackDict_RT[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_RT",region+"_RT;RT;events")
#		stackDict_dEta[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_dEta",region+"_dEta;dEta;events")
		stackDict_var1l[region] = {"ZJets16":rt.TH1D(), "WJets16":rt.TH1D(), "TTJets16":rt.TH1D(), "QCD16":rt.TH1D()}#rt.THStack(region+"_"+var1,region+"_"+var1+";"+var1+";events")
		dictOfBackgroundTotals["low_"+var1+"_"+region] = rt.TH1F("low_METPhi_"+region,"low_METPhi_"+region,20, -rt.TMath.Pi(),rt.TMath.Pi())

	for bkg in ["ZJets16","WJets16","TTJets16","QCD16"]:
		for region in listOfRegions:
			_file1 = rt.TFile.Open("outputs/sfABCD/"+bkg+"/sfABCD_test2_lowdEta_mt1850.root","read")

#			temp_nSVJ =  _file1.Get(bkg+"_nSVJTag_"+region)
#			temp_RT = _file1.Get(bkg+"_RT_"+region)
#			temp_dEta = _file1.Get(bkg+"_DeltaEta_"+region)
			temp_var1 = _file1.Get(bkg+"_"+var1+"_"+region)

#			temp_nSVJ.SetDirectory(0)
#			temp_RT.SetDirectory(0)
#			temp_dEta.SetDirectory(0)
			temp_var1.SetDirectory(0)
			print(var1, region, bkg, temp_var1.GetName(), temp_var1.Integral()) #JFURND
		
			_file1.Close()
			if temp_var1.GetNbinsX() == 1000:
				temp_var1.Rebin(50)
		

#			stackDict_nSVJ[region][bkg] = temp_nSVJ.Clone()
#			stackDict_RT[region][bkg] = temp_RT.Clone()
#			stackDict_dEta[region][bkg] = temp_dEta.Clone()
			stackDict_var1l[region][bkg] = temp_var1.Clone()
			dictOfBackgroundTotals["low_"+var1+"_"+region].Add(stackDict_var1l[region][bkg])

	for region in listOfRegions:
#		makePng([stackDict_nSVJ[region]["QCD16"],stackDict_nSVJ[region]["TTJets16"],stackDict_nSVJ[region]["WJets16"],stackDict_nSVJ[region]["ZJets16"]], region+"_dEtaLow_nSVJ", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
#		makePng([stackDict_RT[region]["QCD16"],stackDict_RT[region]["TTJets16"],stackDict_RT[region]["WJets16"],stackDict_RT[region]["ZJets16"]], region+"_dEtaLow_RT", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
#		makePng([stackDict_dEta[region]["QCD16"],stackDict_dEta[region]["TTJets16"],stackDict_dEta[region]["WJets16"],stackDict_dEta[region]["ZJets16"]], region+"_dEtaLow_dEta", direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)
		makePng([stackDict_var1l[region]["QCD16"],stackDict_var1l[region]["TTJets16"],stackDict_var1l[region]["WJets16"],stackDict_var1l[region]["ZJets16"]], region+"_dEtaLow_"+var1, direc, doLeg = True, log = True, nBinsX = 50, stackLines = False)	
		makeRatio(dictOfBackgroundTotals["high_"+var1+"_"+region], dictOfBackgroundTotals["low_"+var1+"_"+region], region+"_"+var1, direc, doLeg = True, log = True)

