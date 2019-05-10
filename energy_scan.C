

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLegend.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
// The event tree containing labelled MC data with sqrt_s and mcSignal features
const TString fname = "/home/besuser1/Tommaso/MC/data/finalData/feat2/inclmc12.root";

Double_t cutLow = 2.5, peak = 3.097, cutHigh = 3.5;    //cut limits in GeV
Double_t cutStep = 0.005;
Int_t nCutsL = (Int_t)((peak - cutLow) / cutStep);
Int_t nCutsR = (Int_t)((cutHigh - peak) / cutStep);
auto cutLeft  = new Double_t[nCutsL];
auto cutRight = new Double_t[nCutsR];

auto S = new Int_t[nCutsL * nCutsR];
auto B = new Int_t[nCutsL * nCutsR];
auto significance = new Double_t[nCutsL * nCutsR];



void energy_scan(){


   std::cout << " nCutsL: " << nCutsL << " nCutsR: " << nCutsR << std::endl;

   // Initialise the cuts and the other arrays
   int x = 0;
   for(int l = 0; l < nCutsL; l++){
     cutLeft[l] = cutLow;
     cutLow += cutStep;
   }
   for(int r = 0; r < nCutsR; r++){
     cutRight[r] = peak;
     peak += cutStep;;
   }
   for(int x = 0; x < nCutsL * nCutsR; x++){     
     S[x] = 0;
     B[x] = 0;
     significance[x] = 0;
   } 
   

   for(int l = 0; l < nCutsL; l++){
     for(int r = 0; r < nCutsR; r++){
     }
   }


   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   Float_t energy;
   reader -> AddVariable("sqrt_s", &energy);

   // Prepare input tree 
   TFile *input(0);
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Prepare the event tree
   Double_t denergy=energy;
   Int_t signal;
   TTree* theTree = (TTree*)input->Get("ntp1");
   theTree->SetBranchAddress( "sqrt_s",       &denergy    );
   theTree->SetBranchAddress( "mcSignal",     &signal     );


   // Loop over the events in the tree
   for(int ievt = 0; ievt < theTree->GetEntries(); ievt++){
     theTree->GetEntry(ievt);
     if(ievt % 10000 == 0) std::cout << "---Processing event " << ievt << std::endl;
     energy = denergy;
     x = 0;
     for(int l = 0; l < nCutsL; l++){  //for any possible cut interval
       for(int r = 0; r < nCutsR; r++){  
	 if(cutLeft[l] < energy && energy < cutRight[r]){  // if the event falls into the interval 
	   if(signal == 1) S[x]++;     /// if it is signal increment S
	   if(signal == 0) B[x]++;     /// if it is background increment B
	 }
	   x++;
       }
     }
   }


   // Calculate the significance of the intervals
   for(int y = 0; y < nCutsR * nCutsL; y++){  
     significance[y] = ((double)S[y]) / sqrt((double)(S[y]+B[y]));
   }
   
   // Find the best significance
   Double_t bestCutL, bestCutR, bestS = 0;
   Int_t L, R;
   x = 0;
   for(int l = 0; l < nCutsL; l++){  
     for(int r = 0; r < nCutsR; r++){
       if(significance[x] > bestS){
	 bestS = significance[x];
	 L = l;
	 R = r;
       }
	 x++;
     }
   }
   bestCutL = cutLeft[L];
   bestCutR = cutRight[R];

   std::cout << "The best separating interval is [" << bestCutL << ", " << bestCutR << "]" << std::endl;




}
