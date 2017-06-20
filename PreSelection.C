#define PreSelection_cxx
#include "PreSelection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PreSelection::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PreSelection.C
//      root> PreSelection t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("GenElectrons",1); // vec<TLorentzVector>
   fChain->SetBranchStatus("GenElectrons_fromTau",1);//vector<bool>
   fChain->SetBranchStatus("GenElectrons_MT2Activity",1);// vector<double>
   fChain->SetBranchStatus("GenElectrons_RecoTrkAct",1);//vector<double>
   fChain->SetBranchStatus("GenElectrons_RecoTrkd3",1);//vector<double>
   fChain->SetBranchStatus("GenElectrons_RecoTrkIso",1);//vector<double>
   fChain->SetBranchStatus("GenHT",1); // float
   fChain->SetBranchStatus("GenJets",1); // vector<TLorentzVector>
   fChain->SetBranchStatus("GenJets_HTMask",1); // vector<bool>
   fChain->SetBranchStatus("GenJets_MHTMask",1); // vector<bool>
   fChain->SetBranchStatus("GenMET",1); // float
   fChain->SetBranchStatus("GenMETPhi",1); // float
   fChain->SetBranchStatus("GenMHT",1); // float
   fChain->SetBranchStatus("GenMHTPhi",1); // float
   fChain->SetBranchStatus("GenMuons",1); // vector<TLorentzVector>
   fChain->SetBranchStatus("GenMuons_fromTau",1); //vector<bool>
   fChain->SetBranchStatus("GenMuons_MT2Activity",1);//vector<double> 
   fChain->SetBranchStatus("GenMuons_RecoTrkAct",1);//vector<double>
   fChain->SetBranchStatus("GenMuons_RecoTrkd3",1);//vector<double>
   fChain->SetBranchStatus("GenMuons_RecoTrkIso",1);//vector<double>
   fChain->SetBranchStatus("GenParticles",1); //vector<TLorentzVector>
   fChain->SetBranchStatus("GenParticles_ParentId",1); //vector<int>
   fChain->SetBranchStatus("GenParticles_ParentIdx",1);// vector<int>
   fChain->SetBranchStatus("GenParticles_PdgId",1); //vector<int>
   fChain->SetBranchStatus("GenParticles_Status",1);//vector<int>
   fChain->SetBranchStatus("GenTaus",1); //vector<TLorentzVector>
   fChain->SetBranchStatus("GenTaus_had",1); //vector<bool>
   fChain->SetBranchStatus("GenTaus_LeadRecoTrkAct",1); //vector<double>
   fChain->SetBranchStatus("GenTaus_LeadRecoTrkd3",1);//vector<double>
   fChain->SetBranchStatus("GenTaus_LeadRecoTrkIso",1);//vector<double>
   fChain->SetBranchStatus("GenTaus_LeadTrk",1);//vector<TLorentzVector>
   fChain->SetBranchStatus("GenTaus_MT2Activity",1);//vector<double>
   fChain->SetBranchStatus("GenTaus_NNeutralHadrons",1);//vector<int>
   fChain->SetBranchStatus("GenTaus_NProngs",1);//vector<int>
   fChain->SetBranchStatus("GenTaus_Nu",1);//vector<TLorentzVector>
   fChain->SetBranchStatus("GenTops",1);//vector<TLorentzVector>
   fChain->SetBranchStatus("GenTopWeight",1);//float
    
   TCanvas * c1 = new TCanvas();
      std::cout << "Here" << std::endl;
   // Create Proper Histograms //////////////////////////////////////
   TH1F * h_GenHT = new TH1F();
   //////////////////////////////////////////////////////////////////

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      //Fill Histograms, be carfule of vector branches //////////////
      h_GenHT->Fill(GenHT);
      ///////////////////////////////////////////////////////////////
      }
   h_GenHT->Draw();
   c1->SaveAs("Test.pdf");
}
