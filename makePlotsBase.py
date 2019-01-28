# script to plot various histograms in .root file
import ROOT as rt
import sys, os, argparse

rt.gROOT.SetBatch(True)
def makePng(LoH, name, doLeg = True):
	if len(LoH) == 0:
		return False
	c = rt.TCanvas("c1","c1",1200,900)
	stack = rt.THStack()
	for i, hist in enumerate(LoH):
		LoH[i].SetLineColor(i+1)
		if LoH[i].Integral() != 0:
			LoH[i].Scale(1/LoH[i].Integral())
		stack.Add(LoH[i])
	stack.Draw("nostack")
	if doLeg:
		c.BuildLegend(0.6,0.6,0.8,0.8)
	c.SaveAs(outPath+"/"+name+ext)

# arguments
# [0] - makePlots.py
# [-i] - path to .root input file
# [-o] - directory to save image files to
# [-pdf] - extension, defualts to '.png', '.pdf' is also an option
# []

parser = argparse.ArgumentParser(description="Plot histograms avaiable.")

parser.add_argument(
	'-i',
	dest='inFile',
	type=str,
	help='path to input file')

parser.add_argument(
	'-o',
	dest='outPath',
	type=str,
	help='path to output location')

parser.add_argument(
	'-pdf',
	dest='ext',
	action='store_const',
	default='.png',
	const='.pdf',
	help='Use to save as pdf, default png')

args = parser.parse_args()
inFile = args.inFile
outPath = args.outPath
ext = args.ext

print("-"*20)
print("Drawing Histograms from " + inFile)
print("- "*10)
print("Saving plots as " + ext + " files.")
print("- "*10)
print("Saving plots to " + outPath)
print("-"*20)


inTFile = rt.TFile.Open(inFile)

listOfKeys = inTFile.GetListOfKeys()

histVals = []
histMarkers = []

listOfHistLists = []
for hist in listOfKeys:
	hName = hist.GetName().split("_")
	# hist names are "hist_value_eventSeclection_subEventSelection"
	nDim = hName[0] # "hist" for 1d, "hist2d" for 2d
	if not hName[1] in histVals:
		histVals.append(hName[1])
	for mark in hName[2:]:
		if not mark in histMarkers:
			histMarkers.append(mark)

histListDict = {}
for Val in histVals:
	for Marker in histMarkers:
		histListDict[Val+"_"+Marker] = []

for hist in listOfKeys:
	if "2d" in hist.GetName():
		continue
	hName = hist.GetName().split("_")
	hVal = hName[1]
	hMarkers = hName[2:]
	for Marker in hMarkers:
		histListDict[hVal+"_"+Marker].append(hist)

for key, value in histListDict.items():
	if len(value) <= 1:
		continue
	hList = []
	for hKey in value:
		hList.append(rt.TH1F())
		inTFile.GetObject(hKey.GetName(),hList[-1])
		if hList[-1].GetEntries() == 0:
			hList.remove(hList[-1])
	makePng(hList, key)





