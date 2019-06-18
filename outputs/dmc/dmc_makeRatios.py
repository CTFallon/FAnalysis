import ROOT as rt
rt.gROOT.SetBatch(True)

def makeRatioBkgStack(bkgList, data , title, xlabel, ylabel, name, doLeg = True, log = False):
	# bkgList is list of bkgHistos
	# data should be a TH1 of data
	# title is the histogram title
	# xlabel is histogram's x axis title
	# ylabel is histogram's y axis title
	# name is what the .png file will be named.

	# C++ Author: Olivier Couet, adapted to python by Colin Fallon
	# Define the Canvas
	c = rt.TCanvas("c", "canvas", 800, 800)
	# define stack of bkgHistos
	stack = rt.THStack()
	for histo in bkgList:
		stack.Add(histo)
	# Upper plot will be in pad1
	pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0) # Upper and lower plot are joined
	pad1.SetGridx()         # Vertical grid
	pad1.SetGridy()         # Horizontal grid
	if log:
		pad1.SetLogy()
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	#stack.SetStats(0)       # No statistics on upper plot
	#data.SetStats(0)       # No statistics on upper plot
	ks = stack.GetStack().Last().KolmogorovTest(data)
	if data.GetNbinsX() == 1000:
		for thing in stack.GetStack():
			thing.Rebin(10)
		data.Rebin(10)
	chi2 = stack.GetStack().Last().Chi2Test(data,"CHI2/NDF")
	if data.GetNbinsX() == 100:
		for thing in stack.GetStack():
			thing.Rebin(2)
		data.Rebin(2)

	stack.SetMinimum(1)
	stack.SetMaximum(max(stack.GetStack().Last().GetMaximum(),data.GetMaximum())*1.5)
	stack.Draw("hist")        
	data.Draw("E1 same")
	stack.GetYaxis().SetTitle(ylabel)
	if doLeg:
		if (("Phi" in xlabel) or ("(MET)" in xlabel) or ("tau" in xlabel) or ("DeltaR" in xlabel) or ("Output" in xlabel) or ("Num." in xlabel) or ("#eta" in xlabel) or ("ptD" in xlabel) or ("BvsAll" in xlabel) or ("ecf" in xlabel) or ("chgH" in xlabel)):
			leg = rt.TLegend(0.3,0.0,0.7,0.3,title,"brNDC") # pos bot mid
		elif (("#Delta#phi" in xlabel)):
			leg = rt.TLegend(0.3,0.6,0.7,0.9,title,"brNDC") # pos top mid
		elif (("Filter" in xlabel)):
			leg = rt.TLegend(0.5,0.0,0.9,0.3,title,"brNDC") # pos bot right
		else:		
			leg = rt.TLegend(0.5,0.6,0.9,0.9,title,"brNDC") # pos top right
		leg.SetNColumns(2)
		leg.AddEntry(data,"Data","EP")
		leg.AddEntry(0,"","")
		for hist in bkgList[::-1]:
			histID = hist.GetName().split("_")[2]
			leg.AddEntry(hist, histID, "F" )
		leg.AddEntry(0,"KS = {:.2f}, #chi^2/ndf = {:.2f}".format(ks, chi2),"")
		leg.Draw()
	pad1.Update()

	# Do not draw the Y axis label on the upper plot and redraw a small
	# axis instead, in order to avoid the first label (0) to be clipped.
	#h1.GetYaxis().SetLabelSize(0.)
	#axis = rt.TGaxis( -5, 20, -5, 220, 20,220,510,"")
	#axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
	#axis.SetLabelSize(15)
	#axis.Draw()

	# lower plot will be in pad2
	c.cd()          # Go back to the main canvas before defining pad2
	pad2 = rt.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.2)
	pad2.SetGridx() # vertical grid
	pad2.SetGridy()         # Horizontal grid
	pad2.Draw()
	pad2.cd()       # pad2 becomes the current pad

	# Define the ratio plot
	totBkg = stack.GetStack().Last()
	h3 = data.Clone("h3")
	h3.SetLineColor(rt.kBlack)
	h3.SetMinimum(0.50)  # Define Y ..
	h3.SetMaximum(1.50) # .. range
	h3.Sumw2()
	h3.SetStats(0)      # No statistics on lower plot
	h3.Divide(totBkg)
	h3.SetMarkerStyle(21)
	h3.Draw("ep")       # Draw the ratio plot

	# data settings
	data.SetMarkerColor(rt.kBlack)

	# Y axis h1 plot settings
	data.GetYaxis().SetTitleSize(20)
	data.GetYaxis().SetTitleFont(43)
	data.GetYaxis().SetTitleOffset(1.55)

	# Ratio plot (h3) settings
	h3.SetTitle("") # Remove the ratio title

	# Y axis ratio plot settings
	h3.GetYaxis().SetTitle("Data/Bkg")
	h3.GetYaxis().SetNdivisions(505)
	h3.GetYaxis().SetTitleSize(20)
	h3.GetYaxis().SetTitleFont(43)
	h3.GetYaxis().SetTitleOffset(1.55)
	h3.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
	h3.GetYaxis().SetLabelSize(15)

	# X axis ratio plot settings
	h3.GetXaxis().SetTitle(xlabel)
	h3.GetXaxis().SetTitleSize(20)
	h3.GetXaxis().SetTitleFont(43)
	h3.GetXaxis().SetTitleOffset(4.)
	h3.GetXaxis().SetLabelFont(43); # Absolute font size in pixel (precision 3)
	h3.GetXaxis().SetLabelSize(15)

	pad1.Update()
	pad2.Update()
	c.Update()
	const = rt.TF1("const", '[0]', h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax())
	line = rt.TF1("line", '[0]+[1]*x', h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax())
	one = rt.TF1("one", '1', h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax())
	one.SetLineStyle(2)
	#h3.Fit("const","Q")
	one.SetLineColor(rt.kBlack)
	one.Draw("same")
	#h3.Fit("line","Q+")
	#save as .png
	c.SaveAs("plots5/"+name)
	print(name.split("_")[0] + ' ' + name.split("_")[1] + ' ' + name.split("_")[2]+ " {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} ".format(
bkgList[0].Integral(),bkgList[1].Integral(),bkgList[2].Integral(),bkgList[3].Integral(),totBkg.Integral(),data.Integral(),
bkgList[0].GetEntries(),bkgList[1].GetEntries(),bkgList[2].GetEntries(),bkgList[3].GetEntries(),totBkg.GetEntries(),data.GetEntries(), const.GetParameter(0), line.GetParameter(0), line.GetParameter(1)))
	c.Delete()

bkgList = ['TTJets',  'WJets', 'ZJets','QCD']#, "ST"]

# make plots comparing each data year to bkg year 
varList = {
	'MET':['MET [GeV]'],
	'METPhi':['#Phi(MET)'],
	'MHT':['MHT [GeV]'],
	'JetsAK8[0].Pt()':["Leading Jet Pt [GeV]"],
	'JetsAK8[1].Pt()':["Subleading Jet Pt [GeV]"],
	'JetsAK8[0].Eta()':["Leading Jet #eta"],
	'JetsAK8[1].Eta()':["Subleading Jet #eta"],
	'JetsAK8[0].Phi()':["Leading Jet #phi"],
	'JetsAK8[1].Phi()':["Subleading Jet #phi"],
	'JetsAK8_girth[0]':["Leading Jet Girth"],
	'JetsAK8_girth[1]':["Subleading Jet Girth"],
	'JetsAK8_softDropMass[0]':["Leading Jet m_{SD} [GeV]"],
	'JetsAK8_softDropMass[1]':["Subleading Jet m_{SD} [GeV]"],
	'JetsAK8_axismajor[0]':["Leading Jet Major Axis"],
	'JetsAK8_axismajor[1]':["Subleading Jet Major Axis"],
	'JetsAK8_axisminor[0]':["Leading Jet Minor Axis"],
	'JetsAK8_axisminor[1]':["Subleading Jet Minor Axis"],
	'JetsAK8_ptdrlog[0]':["Leading Jet ptdrlog"],
	'JetsAK8_ptdrlog[1]':["Subleading Jet ptdrlog"],
	'JetsAK8_ptD[0]':["Leading Jet ptD"],
	'JetsAK8_ptD[1]':["Subleading Jet ptD"],
	'JetsAK8_maxBvsAll[0]':["Leading Jet maxBvsAll"],
	'JetsAK8_maxBvsAll[1]':["Subleading Jet maxBvsAll"],
	'JetsAK8_ecfN2b1[0]':["Leading Jet ecfN2b1"],
	'JetsAK8_ecfN2b1[1]':["Subleading Jet ecfN2b1"],
	'JetsAK8_ecfN3b1[0]':["Leading Jet ecfN3b1"],
	'JetsAK8_ecfN3b1[1]':["Subleading Jet ecfN3b1"],
	'JetsAK8_chargedHadronEnergyFraction[0]':["Leading Jet chgHadEnrgyFrctn"],
	'JetsAK8_chargedHadronEnergyFraction[1]':["Subleading Jet chgHadEnrgyFrctn"],
	'JetsAK8_neutralHadronEnergyFraction[0]':["Leading Jet neuHadEnrgyFrctn"],
	'JetsAK8_neutralHadronEnergyFraction[1]':["Subleading Jet neuHadEnrgyFrctn"],
	'JetsAK8_electronEnergyFraction[0]':["Leading Jet eleEnrgyFrctn"],
	'JetsAK8_electronEnergyFraction[1]':["Subleading Jet eleEnrgyFrctn"],
	'JetsAK8_muonEnergyFraction[0]':["Leading Jet muEnrgyFrctn"],
	'JetsAK8_muonEnergyFraction[1]':["Subleading Jet muEnrgyFrctn"],
	'JetsAK8_photonEnergyFraction[0]':["Leading Jet #gammaEnrgyFrctn"],
	'JetsAK8_photonEnergyFraction[1]':["Subleading Jet #gammaEnrgyFrctn"],
	'deltaR12':["#DeltaR(j_{1,2})"],
	'metR':["MET/m_{T}"],
	'nJetsAK8':['Num. AK8 Jets'],
	'nJetsAK4':['Num. AK4 Jets'],
	'tau32_lead':['Leading Jet #tau_{32}'],
	'tau21_lead':['Leading Jet #tau_{21}'],
	'tau32_sub':['Subleading Jet #tau_{32}'],
	'tau21_sub':['Subleading Jet #tau_{21}'],
	'JetsAK8_bdtSVJtag[0]':["Leading Jet SVJ BDT Output"],
	'JetsAK8_bdtSVJtag[1]':["Subleading Jet SVJ BDT Output"],
	'DeltaPhi1':["#Delta#phi(j_{1}, MET)"],
	'DeltaPhi2':["#Delta#phi(j_{2}, MET)"],
	'ecalBadCalibReducedFilter':['ECAL Bad Calibration Reduced Filter'],
	'ecalBadCalibReducedExtraFilter':['ECAL Bad Calibration Reduced Extra Filter']
}

# <var>_<sample>

bkgColorDict = {"QCD":41,"TTJets":42,"WJets":43,"ZJets":44}#, "ST":45}

#for each year
	# for each filter
		# for each bkgground
		# for data
	# draw

print("var yr Filter yTT yW yZ yQCD yBkg yData nTT nW nZ nQCD nBkg nData avgRatio")
for year in ['16','17','18PRE', '18POST']:
	for filt in [""]:#,"_filter","_filter2"]:
		for var in varList.keys():
			showLog = True
#			tD = rt.TH1F()
#			tQ = rt.TH1F()
#			tT = rt.TH1F()
#			tZ = rt.TH1F()
#			tW = rt.TH1F()
#			tS = rt.TH1F()
#			tempHist = {"Data":tD,"QCD":tQ,"TTJets":tT,"WJets":tW,"ZJets":tZ,"ST":tS}

#			for bkgGroup in bkgList:
#				bkgName = bkgGroup+year
#				_file = rt.TFile.Open(bkgName+"/dataMC_comp.root","READ")
#				_file.GetObject(var+"_"+bkgName+filt,tempHist[bkgGroup])
#				tempHist[bkgGroup] = tempHist[bkgGroup].Clone(tempHist[bkgGroup].GetName()+"_")
#				tempHist[bkgGroup].SetDirectory(0)
#				tempHist[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
#				tempHist[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
#				tempHist[bkgGroup].SetFillStyle(1001)
#				_file.Close()
#			_file2 = rt.TFile.Open("Data"+year+"/dataMC_comp.root","READ")
#			_file2.GetObject(var+"_Data"+year+filt,tempHist["Data"])
#			tempHist["Data"] = tempHist["Data"].Clone(tempHist["Data"].GetName()+"_")
#			tempHist["Data"].SetDirectory(0)
#			_file2.Close()

#			if "Phi" in var:
#				if not ("Delta" in var):
#					showLog = False
#			makeRatioBkgStack(
#				[tempHist["TTJets"],tempHist["WJets"],tempHist["ZJets"],tempHist["QCD"], tempHist["ST"]], #bkg list
#				tempHist["Data"], # data
#				varList[var][0]+" "+year+filt, # title
#				varList[var][0], # xlabel
#				"Events", # ylabel
#				var.translate(None,"_[]().")+"_"+year+filt+"_ratio_ST.png", #name
#				doLeg = True, # doLegned
#				log = showLog) # doLog
			tD = rt.TH1F()
			tQ = rt.TH1F()
			tT = rt.TH1F()
			tZ = rt.TH1F()
			tW = rt.TH1F()
			tempHist = {"Data":tD,"QCD":tQ,"TTJets":tT,"WJets":tW,"ZJets":tZ}

			for bkgGroup in bkgList:
				bkgName = bkgGroup+year
				_file = rt.TFile.Open(bkgName+"/dataMC_comp.root","READ")
				_file.GetObject(var+"_"+bkgName,tempHist[bkgGroup])
				tempHist[bkgGroup] = tempHist[bkgGroup].Clone(tempHist[bkgGroup].GetName()+"_")
				tempHist[bkgGroup].SetDirectory(0)
				tempHist[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
				tempHist[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
				tempHist[bkgGroup].SetFillStyle(1001)
				nBins = tempHist[bkgGroup].GetNbinsX()
				tempHist[bkgGroup].SetBinContent(1,tempHist[bkgGroup].GetBinContent(0)+tempHist[bkgGroup].GetBinContent(1))
				tempHist[bkgGroup].SetBinContent(nBins,tempHist[bkgGroup].GetBinContent(nBins)+tempHist[bkgGroup].GetBinContent(nBins+1))
				_file.Close()
			_file2 = rt.TFile.Open("Data"+year+"/dataMC_comp.root","READ")
			_file2.GetObject(var+"_Data"+year,tempHist["Data"])
			tempHist["Data"] = tempHist["Data"].Clone(tempHist["Data"].GetName()+"_")
			tempHist["Data"].SetDirectory(0)
			nBins = tempHist["Data"].GetNbinsX()
			tempHist["Data"].SetBinContent(1,tempHist["Data"].GetBinContent(0)+tempHist["Data"].GetBinContent(1))
			tempHist["Data"].SetBinContent(nBins,tempHist["Data"].GetBinContent(nBins)+tempHist["Data"].GetBinContent(nBins+1))
			_file2.Close()
			makeRatioBkgStack(
				[tempHist["TTJets"],tempHist["WJets"],tempHist["ZJets"],tempHist["QCD"]], #bkg list
				tempHist["Data"], # data
				varList[var][0]+" "+year+filt, # title
				varList[var][0], # xlabel
				"Events", # ylabel
				var.translate(None,"_[]().")+"_"+year+filt+"_ratio.png", #name
				doLeg = True, # doLegned
				log = showLog) # doLog
















