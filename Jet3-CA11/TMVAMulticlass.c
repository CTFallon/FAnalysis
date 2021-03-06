#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAMultiClassGui.h"

using namespace TMVA;
void TMVAMulticlass( TString myMethodList = "" )
{
   // This loads the library
   TMVA::Tools::Instance();
   // to get access to the GUI and all tmva macros
   //
   //     TString tmva_dir(TString(gRootDir) + "/tmva");
   //     if(gSystem->Getenv("TMVASYS"))
   //        tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   //     gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   //     gROOT->ProcessLine(".L TMVAMultiClassGui.C");
   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 1;
   Use["BDTG"]            = 1;
   Use["DNN"]             = 0;
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 0;
   //---------------------------------------------------------------
   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlass" << std::endl;
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);
         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }
   // Create a new root output file.
   TString outfileName = "TMVAMulticlass.root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   dataloader->AddVariable( "var1", 'F' );
   dataloader->AddVariable( "var2", "Variable 2", "", 'F' );
   dataloader->AddVariable( "var3", "Variable 3", "units", 'F' );
   dataloader->AddVariable( "var4", "Variable 4", "units", 'F' );
   TFile *input(0);
   TString fname = "TMVATrainingInput.root";
   if (!gSystem->AccessPathName( fname )) {
      // first we try to find the file in the local directory
      std::cout << "--- TMVAMulticlass   : Accessing " << fname << std::endl;
      input = TFile::Open( fname );
   }
   else {
      std::cout << "Creating testdata...." << std::endl;
      TString createDataMacro = TString(gROOT->GetTutorialsDir()) + "/tmva/createData.C";
      gROOT->ProcessLine(TString::Format(".L %s",createDataMacro.Data()));
      gROOT->ProcessLine("create_MultipleBackground(2000)");
      std::cout << " created tmva_example_multiple_background.root for tests of the multiclass features"<<std::endl;
      input = TFile::Open( fname );
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   TTree *signalTree  = (TTree*)input->Get("jet12");
   TTree *background0 = (TTree*)input->Get("jet1");
   TTree *background1 = (TTree*)input->Get("jet13");
   TTree *background2 = (TTree*)input->Get("jet123");
   TTree *background3 = (TTree*)input->Get("jet1234");
   TTree *background4 = (TTree*)input->Get("jetOther");

   gROOT->cd( outfileName+TString(":/") );
   dataloader->AddTree    (signalTree,"jets12");
   dataloader->AddTree    (background0,"jets1");
   dataloader->AddTree    (background1,"jets13");
   dataloader->AddTree    (background2,"jets123");
   dataloader->AddTree    (background3,"jets1234");
   dataloader->AddTree    (background4,"jetsOther");

   dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );
   if (Use["BDTG"]) // gradient boosted decision trees
      factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   if (Use["MLP"]) // neural network
      factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
   if (Use["FDA_GA"]) // functional discriminant with GA minimizer
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   if (Use["PDEFoam"]) // PDE-Foam approach
      factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam", "!H:!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
   if (Use["DNN"]) {
       TString layoutString ("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
       TString training0 ("LearningRate=1e-1, Momentum=0.5, Repetitions=1, ConvergenceSteps=10,"
                          " BatchSize=256, TestRepetitions=10, Multithreading=True");
       TString training1 ("LearningRate=1e-2, Momentum=0.0, Repetitions=1, ConvergenceSteps=10,"
                          " BatchSize=256, TestRepetitions=7, Multithreading=True");
       TString trainingStrategyString ("TrainingStrategy=");
       trainingStrategyString += training0 + "|" + training1;
       TString nnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM:Architecture=STANDARD");
       nnOptions.Append (":"); nnOptions.Append (layoutString);
       nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
       factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN", nnOptions );
   }
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   // --------------------------------------------------------------
   // Save the output
   outputFile->Close();
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAMultiClassGui( outfileName );
}
int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAMulticlass(methodList);
   return 0;
}
