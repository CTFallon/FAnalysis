import ROOT as rt
import inspect
import sys
from array import array
from DataFormats.FWLite import Events, Handle

def trans_mass_Njet(jetsL, met, metPhi):
	#visible = rt.Math.LorentzVector<rt.Math.PxPyPzE4D>()
	jet0 = rt.TLorentzVector(jetsL[0].px(),jetsL[0].py(),jetsL[0].pz(),jetsL[0].energy())
	jet1 = rt.TLorentzVector(jetsL[1].px(),jetsL[1].py(),jetsL[1].pz(),jetsL[1].energy())
	#print(type(jet0), type(jet1))
	visible = jet0 + jet1
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	return rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))



# args, in order mZprime, mDark, Rinv, Alpha
dSet, mZprime, mDark, rInv, Alpha = sys.argv[1:]

print("dSet mZprime mDark rInv Alpha")
print("{} {} {} {} {}".format(dSet, mZprime, mDark, rInv, Alpha))

fileStart = "root://cmsxrootd.fnal.gov//store/user/cfallon/RFSV/"

rootOutFile = rt.TFile.Open("Skims/{}_mZprime-{}_mDark-{}_rinv-{}_alpha-{}.root".format(dSet,mZprime,mDark,rInv,Alpha), "recreate")

rootOutTree = rt.TTree("tree_"+dSet,"output of RFSV")

jet0Pt = array('f',[ 0 ])
jet0Eta = array('f',[ 0 ])
jet1Pt = array('f',[ 0 ])
jet1Eta = array('f',[ 0 ])

eventMET = array('f', [ 0 ])
eventRT = array('f', [ 0 ])
eventMT = array('f', [ 0 ])

rootOutTree.Branch("jet0Pt", jet0Pt, 'jet0Pt/F')
rootOutTree.Branch("jet0Eta", jet0Eta, 'jet0Eta/F')
rootOutTree.Branch("jet1Pt", jet1Pt, 'jet1Pt/F')
rootOutTree.Branch("jet1Eta", jet1Eta, 'jet1Eta/F')

rootOutTree.Branch("eventMET", eventMET, 'eventMET/F')
rootOutTree.Branch("eventRT", eventRT, 'eventRT/F')
rootOutTree.Branch("eventMT", eventMT, 'eventMT/F')

listOfFiles = [fileStart+"step1_GEN_{}_mZprime-{}_mDark-{}_rinv-{}_alpha-{}_n-500_part-{}.root".format(dSet,mZprime,mDark,rInv,Alpha,iPart) for iPart in range(1,21)]

events = Events ( listOfFiles )

handleJets = Handle("std::vector<reco::GenJet>")
handleMET = Handle("std::vector<reco::GenMET>")

labelJets = ("ak8GenJetsNoNu")
labelMET = ("genMetTrue")

rt.gROOT.SetBatch(True)
rt.gROOT.SetStyle("Plain")

for event in events:
	event.getByLabel(labelJets,handleJets)
	event.getByLabel(labelMET,handleMET)
	jets = handleJets.product()
	MET = handleMET.product()
	numJets = len(jets)
	if numJets < 2:
		continue
	
	jet0Pt[0] = jets[0].pt()
	jet0Eta[0] = jets[0].eta()
	jet1Pt[0] = jets[1].pt()
	jet1Eta[0] = jets[1].eta()
	eventMET[0] = MET[0].pt()
	eventMT[0] = trans_mass_Njet([jets[0],jets[1]],MET[0].pt(), MET[0].phi())
	eventRT[0] = eventMET[0]/eventMT[0]
	rootOutTree.Fill()
rootOutFile.Write()
rootOutFile.Close()
	
