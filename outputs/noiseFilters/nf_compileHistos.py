import ROOT as rt
rt.gROOT.SetBatch(True)

toCopyToExcel = open("copyToExcel_2.txt","w")


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
	if log:
		pad1.SetLogy()
	pad1.Draw()             # Draw the upper pad: pad1
	pad1.cd()               # pad1 becomes the current pad
	#stack.SetStats(0)       # No statistics on upper plot
	#data.SetStats(0)       # No statistics on upper plot
	stack.SetMinimum(1)
	stack.Draw("hist")        
	data.Draw("E1 same")
	stack.GetYaxis().SetTitle(ylabel)
	if doLeg:
		leg = rt.TLegend(0.5,0.6,0.9,0.9,"20"+name.split("_")[1]+" "+name.split("_")[2],"brNDC")
		leg.AddEntry(data,"Data","P")
		for hist in bkgList[::-1]:
			histID = hist.GetName().split("_")[2][0:-2]
			leg.AddEntry(hist, histID, "F" )
		leg.Draw()

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
	pad2.Draw()
	pad2.cd()       # pad2 becomes the current pad

	# Define the ratio plot
	totBkg = stack.GetStack().Last()
	h3 = data.Clone("h3")
	h3.SetLineColor(rt.kBlack)
	h3.SetMinimum(0.00)  # Define Y ..
	h3.SetMaximum(2.00) # .. range
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
	h3.Fit("const","Q")
	h3.Fit("line","Q+")
	#save as .png
	#c.SaveAs(name)
	toCopyToExcel.write(name.split("_")[0] + ' ' + name.split("_")[1] + ' ' + name.split("_")[2]+ " {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
bkgList[0].Integral(),bkgList[1].Integral(),bkgList[2].Integral(),bkgList[3].Integral(),totBkg.Integral(),data.Integral(),
bkgList[0].GetEntries(),bkgList[1].GetEntries(),bkgList[2].GetEntries(),bkgList[3].GetEntries(),totBkg.GetEntries(),data.GetEntries(), const.GetParameter(0), line.GetParameter(0), line.GetParameter(1)))
	c.Delete()

def makeSigBkg(bkgList, sigList , title, xlabel, ylabel, name, doLeg = True, log = False):
	c = rt.TCanvas("c", "canvas", 800, 800)
	stack = rt.THStack()
	stackSig = rt.THStack()
	for histo in bkgList:
		stack.Add(histo)
	for i, histo in enumerate(sigList):
		histo.SetLineColor(i+1)
		stackSig.Add(histo)
	if log:
		c.SetLogy()
	stack.SetMinimum(1)
	stack.Draw("hist")
	stackSig.Draw("nostack same hist")
	stack.GetYaxis().SetTitle(ylabel)
	if doLeg:
		leg = rt.TLegend(0.5,0.6,0.9,0.9,"20"+name.split("_")[1]+" "+name.split("_")[2],"brNDC")
		for hist in bkgList[::-1]:
			histID = hist.GetName().split("_")[2][0:-2]
			leg.AddEntry(hist, histID, "F" )
		for hist in sigList:
			histID = hist.GetName().split("_")[2]
			leg.AddEntry(hist, histID, "L" )
		leg.Draw()

	c.Update()
	#save as .png
	#c.SaveAs(name)
	if xlabel == "MT":
		toCopyToExcel.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
name.split("_")[2], name.split("_")[0],name.split("_")[3],name.split("_")[1],
bkgList[0].Integral(),bkgList[0].GetEntries(),
bkgList[1].Integral(),bkgList[1].GetEntries(),
bkgList[2].Integral(),bkgList[2].GetEntries(),
bkgList[3].Integral(),bkgList[3].GetEntries(),
sigList[0].Integral(),
sigList[1].Integral(),
sigList[2].Integral(),
sigList[3].Integral(),
sigList[4].Integral()
))
	c.Delete()

bkgList = ['TTJets',  'WJets', 'ZJets','QCD']

# make plots comparing each filter to data filters 
filterList = [
	"All",
	"None",
	"BadChargedCandidateFilter",
	"BadPFMuonFilter",
	"CSCTightHaloFilter",
	"ecalBadCalibFilter",
	"EcalDeadCellTriggerPrimitiveFilter",
	"eeBadScFilter",
	"globalSuperTightHalo2016Filter",
	"globalTightHalo2016Filter",
	"HBHEIsoNoiseFilter",
	"HBHENoiseFilter",
	"PrimaryVertexFilter",
	"METRatioFilter",
	"MuonJetFilter",
	"EcalNoiseJetFilter",
	"HTRatioFilter",
	"HTRatioTightFilter",
	"HTRatioDPhiFilter",
	"HTRatioDPhiTightFilter",
	"LowNeutralJetFilter",
	"LowNeutralJetTightFilter"]

# hist_<variable>_<fileID>_<filter>

bkgColorDict = {"QCD":41,"TTJets":42,"WJets":43,"ZJets":44}


doData = False
doSignal = True

if doData:
	toCopyToExcel.write("var yr Filter yTT yW yZ yQCD yBkg yData nTT nW nZ nQCD nBkg nData avgRatio\n")
	for year in ['16','17']:
		print(year)
		for filt in filterList:
			print(filt)
			tDMT = rt.TH1F()
			tQMT = rt.TH1F()
			tTMT = rt.TH1F()
			tZMT = rt.TH1F()
			tWMT = rt.TH1F()
			tempHistMT = {"Data":tDMT,"QCD":tQMT,"TTJets":tTMT,"WJets":tWMT,"ZJets":tZMT}
			tDpt0 = rt.TH1F()
			tQpt0 = rt.TH1F()
			tTpt0 = rt.TH1F()
			tZpt0 = rt.TH1F()
			tWpt0 = rt.TH1F()
			tempHistpt0 = {"Data":tDpt0,"QCD":tQpt0,"TTJets":tTpt0,"WJets":tWpt0,"ZJets":tZpt0}
			tDpt1 = rt.TH1F()
			tQpt1 = rt.TH1F()
			tTpt1 = rt.TH1F()
			tZpt1 = rt.TH1F()
			tWpt1 = rt.TH1F()
			tempHistpt1 = {"Data":tDpt1,"QCD":tQpt1,"TTJets":tTpt1,"WJets":tWpt1,"ZJets":tZpt1}
			for bkgGroup in bkgList:
				bkgName = bkgGroup+year
				_file = rt.TFile.Open(bkgName+"/nF_test.root","READ")
				_file.GetObject("hist_MT_"+bkgName+"_"+filt,tempHistMT[bkgGroup])
				tempHistMT[bkgGroup] = tempHistMT[bkgGroup].Clone(tempHistMT[bkgGroup].GetName()+"_")
				tempHistMT[bkgGroup].SetDirectory(0)
				tempHistMT[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
				tempHistMT[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
				tempHistMT[bkgGroup].SetFillStyle(1001)
				_file.GetObject("hist_Jet0Pt_"+bkgName+"_"+filt,tempHistpt0[bkgGroup])
				tempHistpt0[bkgGroup] = tempHistpt0[bkgGroup].Clone(tempHistpt0[bkgGroup].GetName()+"_")
				tempHistpt0[bkgGroup].SetDirectory(0)
				tempHistpt0[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
				tempHistpt0[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
				tempHistpt0[bkgGroup].SetFillStyle(1001)
				_file.GetObject("hist_Jet1Pt_"+bkgName+"_"+filt,tempHistpt1[bkgGroup])
				tempHistpt1[bkgGroup] = tempHistpt1[bkgGroup].Clone(tempHistpt1[bkgGroup].GetName()+"_")
				tempHistpt1[bkgGroup].SetDirectory(0)
				tempHistpt1[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
				tempHistpt1[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
				tempHistpt1[bkgGroup].SetFillStyle(1001)
				_file.Close()
			_file2 = rt.TFile.Open("Data"+year+"/nF_test.root","READ")
			_file2.GetObject("hist_MT_Data"+year+"_"+filt,tempHistMT["Data"])
			tempHistMT["Data"] = tempHistMT["Data"].Clone(tempHistMT["Data"].GetName()+"_")
			tempHistMT["Data"].SetDirectory(0)
			_file2.GetObject("hist_Jet0Pt_Data"+year+"_"+filt,tempHistpt0["Data"])
			tempHistpt0["Data"] = tempHistpt0["Data"].Clone(tempHistpt0["Data"].GetName()+"_")
			tempHistpt0["Data"].SetDirectory(0)
			_file2.GetObject("hist_Jet1Pt_Data"+year+"_"+filt,tempHistpt1["Data"])
			tempHistpt1["Data"] = tempHistpt1["Data"].Clone(tempHistpt1["Data"].GetName()+"_")
			tempHistpt1["Data"].SetDirectory(0)
			_file2.Close()

			makeRatioBkgStack([tempHistMT["TTJets"],tempHistMT["WJets"],tempHistMT["ZJets"],tempHistMT["QCD"]],tempHistMT["Data"],filt,"MT","Events","MT_"+year+"_"+filt+"_ratio.png", doLeg = True, log = True)
			makeRatioBkgStack([tempHistpt0["TTJets"],tempHistpt0["WJets"],tempHistpt0["ZJets"],tempHistpt0["QCD"]],tempHistpt0["Data"],filt,"Jet0Pt","Events","Jet0Pt_"+year+"_"+filt+"_ratio.png", doLeg = True, log = True)
			makeRatioBkgStack([tempHistpt1["TTJets"],tempHistpt1["WJets"],tempHistpt1["ZJets"],tempHistpt1["QCD"]],tempHistpt1["Data"],filt,"Jet1Pt","Events","Jet1Pt_"+year+"_"+filt+"_ratio.png", doLeg = True, log = True)

mZ10 = ['z05','z10','z20','base-17','z40']
mZLow = ['z10','z15','z20','z25','base-17']
mZHigh = ['base-17','z32','z35','z40','z45']

mDLow = ['d1','d5','d10','base-17','d30']
mDHigh = ['d60','d70','d80','d90','d100']
mDEven = ['base-17','d40','d60','d80','d100']

rILow = ['r01','r02','base-17','r04','r05',]
rIHigh = ['r06','r07','r08','r09','r10']

aD = ['al','base-17','ah', 'base-16','a2']

coreList = ['base-16', 'base-17','a2','d80','z40']

mZScan1 = ['z05', 'z10', 'z15', 'z20', 'z25']
mZScan2 = ['z25', 'base-17', 'z35', 'z40', 'z45']

filterList = [
	"All",
	"None",
	"BadChargedCandidateFilter",
	"BadPFMuonFilter",
	#"CSCTightHaloFilter",
	#"ecalBadCalibFilter",
	"EcalDeadCellTriggerPrimitiveFilter",
	"eeBadScFilter",
	"globalSuperTightHalo2016Filter",
	#"globalTightHalo2016Filter",
	"HBHEIsoNoiseFilter",
	"HBHENoiseFilter",
	#"PrimaryVertexFilter",
	"METRatioFilter",
	"MuonJetFilter",
	"EcalNoiseJetFilter",
	#"HTRatioFilter",
	#"HTRatioTightFilter",
	#"HTRatioDPhiFilter",
	"HTRatioDPhiTightFilter",
	"LowNeutralJetFilter"]#,
	#"LowNeutralJetTightFilter"]

if doSignal:
	for year in ["16",'17']:
		for j, signalList in enumerate([mZ10, mZLow, mZHigh, mDLow, mDHigh,mDEven, rILow, rIHigh, aD, coreList, mZScan1, mZScan2]):
			toCopyToExcel.write("year sigListIdx Filter var nTT yTT nW yW nZ yZ nQCD yQCD n{0} n{1} n{2} n{3} n{4}\n".format(signalList[0],signalList[1],signalList[2],signalList[3],signalList[4],))
			for filt in filterList:
				tS0MT = rt.TH1F()
				tS1MT = rt.TH1F()
				tS2MT = rt.TH1F()
				tS3MT = rt.TH1F()
				tS4MT = rt.TH1F()
				tQMT = rt.TH1F()
				tTMT = rt.TH1F()
				tZMT = rt.TH1F()
				tWMT = rt.TH1F()
				tempHistMT = {signalList[0]:tS0MT,signalList[1]:tS1MT,signalList[2]:tS2MT,signalList[3]:tS3MT,signalList[4]:tS4MT,"QCD":tQMT,"TTJets":tTMT,"WJets":tWMT,"ZJets":tZMT}
				tS0pt0 = rt.TH1F()
				tS1pt0 = rt.TH1F()
				tS2pt0 = rt.TH1F()
				tS3pt0 = rt.TH1F()
				tS4pt0 = rt.TH1F()
				tQpt0 = rt.TH1F()
				tTpt0 = rt.TH1F()
				tZpt0 = rt.TH1F()
				tWpt0 = rt.TH1F()
				tempHistpt0 = {signalList[0]:tS0pt0,signalList[1]:tS1pt0,signalList[2]:tS2pt0,signalList[3]:tS3pt0,signalList[4]:tS4pt0,"QCD":tQpt0,"TTJets":tTpt0,"WJets":tWpt0,"ZJets":tZpt0}
				tS0pt1 = rt.TH1F()
				tS1pt1 = rt.TH1F()
				tS2pt1 = rt.TH1F()
				tS3pt1 = rt.TH1F()
				tS4pt1 = rt.TH1F()
				tQpt1 = rt.TH1F()
				tTpt1 = rt.TH1F()
				tZpt1 = rt.TH1F()
				tWpt1 = rt.TH1F()
				tempHistpt1 = {signalList[0]:tS0pt1,signalList[1]:tS1pt1,signalList[2]:tS2pt1,signalList[3]:tS3pt1,signalList[4]:tS4pt1,"QCD":tQpt1,"TTJets":tTpt1,"WJets":tWpt1,"ZJets":tZpt1}
				for bkgGroup in bkgList:
					bkgName = bkgGroup+year
					_file = rt.TFile.Open(bkgName+"/nF_test.root","READ")
					_file.GetObject("hist_MT_"+bkgName+"_"+filt,tempHistMT[bkgGroup])
					tempHistMT[bkgGroup] = tempHistMT[bkgGroup].Clone(tempHistMT[bkgGroup].GetName()+"_")
					tempHistMT[bkgGroup].SetDirectory(0)
					tempHistMT[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
					tempHistMT[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
					tempHistMT[bkgGroup].SetFillStyle(1001)
					_file.GetObject("hist_Jet0Pt_"+bkgName+"_"+filt,tempHistpt0[bkgGroup])
					tempHistpt0[bkgGroup] = tempHistpt0[bkgGroup].Clone(tempHistpt0[bkgGroup].GetName()+"_")
					tempHistpt0[bkgGroup].SetDirectory(0)
					tempHistpt0[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
					tempHistpt0[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
					tempHistpt0[bkgGroup].SetFillStyle(1001)
					_file.GetObject("hist_Jet1Pt_"+bkgName+"_"+filt,tempHistpt1[bkgGroup])
					tempHistpt1[bkgGroup] = tempHistpt1[bkgGroup].Clone(tempHistpt1[bkgGroup].GetName()+"_")
					tempHistpt1[bkgGroup].SetDirectory(0)
					tempHistpt1[bkgGroup].SetLineColor(bkgColorDict[bkgGroup])
					tempHistpt1[bkgGroup].SetFillColor(bkgColorDict[bkgGroup])
					tempHistpt1[bkgGroup].SetFillStyle(1001)
					_file.Close()
				for sigGroup in signalList:
					sigName = sigGroup
					_file = rt.TFile.Open("signal/"+sigName+"/nF_test.root","READ")
					_file.GetObject("hist_MT_"+sigName+"_"+filt,tempHistMT[sigGroup])
					tempHistMT[sigGroup] = tempHistMT[sigGroup].Clone(tempHistMT[sigGroup].GetName()+"_")
					tempHistMT[sigGroup].SetDirectory(0)
					_file.GetObject("hist_Jet0Pt_"+sigName+"_"+filt,tempHistpt0[sigGroup])
					tempHistpt0[sigGroup] = tempHistpt0[sigGroup].Clone(tempHistpt0[sigGroup].GetName()+"_")
					tempHistpt0[sigGroup].SetDirectory(0)
					_file.GetObject("hist_Jet1Pt_"+sigName+"_"+filt,tempHistpt1[sigGroup])
					tempHistpt1[sigGroup] = tempHistpt1[sigGroup].Clone(tempHistpt1[sigGroup].GetName()+"_")
					tempHistpt1[sigGroup].SetDirectory(0)
					_file.Close()

				makeSigBkg([tempHistMT["TTJets"],tempHistMT["WJets"],tempHistMT["ZJets"],tempHistMT["QCD"]],[tempHistMT[signalList[0]],tempHistMT[signalList[1]],tempHistMT[signalList[2]],tempHistMT[signalList[3]],tempHistMT[signalList[4]]],filt,"MT","Events",str(j)+"_MT_"+year+"_"+filt+"_sig.png", doLeg = True, log = True)
				makeSigBkg([tempHistpt0["TTJets"],tempHistpt0["WJets"],tempHistpt0["ZJets"],tempHistpt0["QCD"]],[tempHistpt0[signalList[0]],tempHistpt0[signalList[1]],tempHistpt0[signalList[2]],tempHistpt0[signalList[3]],tempHistpt0[signalList[4]]],filt,"Jet0Pt","Events",str(j)+"_AK4_Jet0Pt_"+year+"_"+filt+"_sig.png", doLeg = True, log = True)
				makeSigBkg([tempHistpt1["TTJets"],tempHistpt1["WJets"],tempHistpt1["ZJets"],tempHistpt1["QCD"]],[tempHistpt1[signalList[0]],tempHistpt1[signalList[1]],tempHistpt1[signalList[2]],tempHistpt1[signalList[3]],tempHistpt1[signalList[4]]],filt,"Jet1Pt","Events",str(j)+"_AK4_Jet1Pt_"+year+"_"+filt+"_sig.png", doLeg = True, log = True)

toCopyToExcel.close()

