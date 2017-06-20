#include <iostream>
#include <TFile.h>
#include <TTree.h>

#include "PreSelection.h"
#include "PreSelection.C"

using namespace std;
 
void doLoop(){
  PreSelection analysis = PreSelection();
  analysis.Loop();
}
