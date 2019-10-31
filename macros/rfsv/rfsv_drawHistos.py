import ROOT as rt
import sys

#rt.gStyle.SetPalette(54,0)
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

listOfDudes = [
["3000","20","0.3","peak", .184, .141],
["2000","20","0.3","peak",.168,.163],
["4000","20","0.3","peak",.147,.163],
["3000","50","0.3","peak",.135,.105],
["3000","100","0.3","peak",.141,.141],
["3000","20","0.5","peak",.189,.180],
["3000","20","0.7","peak",.201,.174],
["3000","20","0.3","low",.163,.138],
["3000","20","0.3","high",.160,.126]]

histLimsDict = {
"eventMT":[24, 1500, 3800],
"jet0Pt":[30, 0, 3000],
"jet1Pt":[30, 0, 3000],
"jet0Eta":[11, -2.4, 2.4],
"jet1Eta":[11, -2.4, 2.4],
'eventRT':[25, 0.0, 1.0],
'eventMET':[32, 100, 1600]}
#print("dSet mZ mD rI a nB nD fD rD iH fH rH")
print("iBin MT Base Double Half")
saveDir = "plots_ANbins_3/"
for mZprime, mDark, rInv, Alpha, errorR, errorF in listOfDudes:
	idString = "mZprime-{}_mDark-{}_rinv-{}_alpha-{}".format(mZprime, mDark, rInv, Alpha)
	
	#take rf, xH and xD
	# plot the suckers ontop of one another.
	# color xH red for half
	# color xD green for Double
	# color rf black for base
	for dSet in ["r","f"]:
		if dSet == "r":
			DSET = "Renormlization"
			error = errorR
		elif dSet == "f":
			DSET = "Factorization"
			error = errorF
		_fileBase = rt.TFile.Open("Skims/rf_"+idString+".root","read")
		_fileHalf = rt.TFile.Open("Skims/"+dSet+"H_"+idString+".root","read")
		_fileDoub = rt.TFile.Open("Skims/"+dSet+"D_"+idString+".root","read")

		treeBase = _fileBase.Get("tree_rf")
		treeHalf = _fileHalf.Get("tree_"+dSet+"H")
		treeDoub = _fileDoub.Get("tree_"+dSet+"D")

		selcJetPt = "(jet0Pt>200)&&(jet1Pt>200)"
		selcJetEta = "(abs(jet0Eta)<2.4)&&(abs(jet1Eta)<2.4)"
		selcJetAbsEta = "(abs(jet0Eta-jet1Eta)<1.5)"
		selcRTMT = "(eventRT>0.15)&&(eventMT>1500)"
		preSelc = selcJetPt+"&&"+selcJetEta+"&&"+selcJetAbsEta+"&&"+selcRTMT
		#preSelc = ""
		for var in ["eventMT"]:#, "jet0Pt", "jet1Pt", "jet0Eta","jet1Eta",'eventRT','eventMET']:
			c1 = rt.TCanvas("c1","c1",900,600)
			c1.SetLogy()
			hist_Base = rt.TH1F("hB",DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Events",histLimsDict[var][0],histLimsDict[var][1],histLimsDict[var][2])
			hist_Half = rt.TH1F("hH",DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Events",histLimsDict[var][0],histLimsDict[var][1],histLimsDict[var][2])
			hist_Doub = rt.TH1F("hD",DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Events",histLimsDict[var][0],histLimsDict[var][1],histLimsDict[var][2])
			hist_Base.Sumw2()
			hist_Half.Sumw2()
			hist_Doub.Sumw2()
			treeBase.Draw(var+">>hB",preSelc,"hist SAME")
			treeHalf.Draw(var+">>hH",preSelc,"hist SAME")
			treeDoub.Draw(var+">>hD",preSelc,"hist SAME")
			hist_BaseError = hist_Base.Clone("hBE")
			for iBin in range(1, hist_BaseError.GetNbinsX()+1):
				hist_BaseError.SetBinError(iBin, hist_BaseError.GetBinContent(iBin)*error)
			#hist_Base = rt.gDirectory.Get("hB")
			#hist_Half = rt.gDirectory.Get("hH")
			#hist_Doub = rt.gDirectory.Get("hD")
			hist_Base.SetLineColor(rt.kBlack)
			hist_Half.SetLineColor(rt.kRed)
			hist_Doub.SetLineColor(rt.kBlue)
			hist_BaseError.SetFillColorAlpha(rt.kBlack,0.5)
			hist_BaseError.Draw("same e2")
			hist_Base.SetLineWidth(2)
			hist_Half.SetLineWidth(2)
			hist_Doub.SetLineWidth(2)
			hist_Base.SetTitle(DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Events")
			hist_Half.SetTitle(DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Ratio")
			hist_Doub.SetTitle(DSET+" "+idString.replace("-","=").replace("_",", ")+";"+var+";Ratio")
			c1.Update()
			leg = rt.TLegend(0.7,0.7,0.9,0.9)
			leg.AddEntry(hist_Doub, "Double")
			leg.AddEntry(hist_Base, "Baseline")
			leg.AddEntry(hist_Half, "Half")
			leg.Draw("same")
			c1.SaveAs(saveDir+var+"_"+DSET+"_"+idString.replace(".","p")+".png")
			for iBin in range(1, hist_Doub.GetNbinsX()+1):
				print("{} {} {} {} {} {} {} {} {} {} {}".format(dSet, mZprime, mDark, rInv, Alpha, iBin, hist_Base.GetBinLowEdge(iBin), hist_Base.GetBinContent(iBin),hist_Doub.GetBinContent(iBin), hist_Half.GetBinContent(iBin), error))
			iH = hist_Half.Integral()
			iB = hist_Base.Integral()
			iD = hist_Doub.Integral()
			c1.SetLogy(0)
			hist_Doub.Divide(hist_Base)
			hist_Doub.SetMarkerStyle(5)
			hist_Doub.SetMarkerSize(2)
			hist_Doub.SetMarkerColor(rt.kBlue)
			hist_Half.Divide(hist_Base)
			lineDoub = rt.TF1("lD","[0]")
			lineHalf = rt.TF1("lH","[0]")
			line_Low = rt.TF1("lLow",str(1-error),1500,3800)
			line_Hi = rt.TF1("lHi",str(1+error),1500,3800)
			line_Low.SetLineColor(rt.kBlack)
			line_Hi.SetLineColor(rt.kBlack)
			hist_Doub.Fit(lineDoub,"LQCN")
			hist_Half.Fit(lineHalf,"LCQN")
			fD = lineDoub.GetParameter(0)
			fH = lineHalf.GetParameter(0)
			hist_Doub.Draw("p e1")
			hist_Half.Draw("p e1 same")
			line_Low.Draw("same")
			line_Hi.Draw("same")
			leg = rt.TLegend(0.7,0.7,0.9,0.9)
			leg.AddEntry(hist_Doub, "Double","lpe")
			leg.AddEntry(hist_Half, "Half","lpe")
			leg.AddEntry(line_Low, "{}% Error".format(error*100),"l")
			leg.Draw("same")
			c1.SaveAs(saveDir+var+"_"+DSET+"_"+idString.replace(".","p")+"_ratio.png")
			#print("{} {} {} {} {} {} {} {} {} {} {} {}".format(dSet, mZprime, mDark, rInv, Alpha, iB, iD, fD, float(iB)/iD, iH, fH, float(iB)/iH))
			









