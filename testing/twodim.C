

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
// The event tree containing labelled MC data with sqrt_s, chisq and mcSignal features
const TString fname = "/home/besuser1/Tommaso/MC/data/finalData/feat4/inclmc12.root";


Double_t cutL = 3.065, cutR = 3.112, cutChi = 15.94;   //cuts on energy and cut on chisq
Int_t S = 0, B = 0;
Double_t significance;
Double_t SoverB;



void twodim(){

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   Float_t energy,chi;
   reader -> AddVariable("sqrt_s", &energy);
   reader -> AddVariable("chisq4C", &chi);

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
   Double_t denergy=energy, dchi=chi;
   Int_t signal;
   TTree* theTree = (TTree*)input->Get("ntp1");
   theTree->SetBranchAddress( "sqrt_s",       &denergy    );
   theTree->SetBranchAddress( "chisq4C",       &dchi    );
   theTree->SetBranchAddress( "mcSignal",     &signal     );

   Int_t totS = 0, totB = 0;

   // Loop over the events in the tree
   for(int ievt = 0; ievt < theTree->GetEntries(); ievt++){
     theTree->GetEntry(ievt);
     energy = denergy; chi = dchi;
     if(energy > cutL && energy < cutR && chi < cutChi){ // if the event is a positive
       if(signal == 1) S++;
       if(signal == 0) B++;
     }
     if(signal == 0) totB++;
     if(signal == 1) totS++;

   }
			 


   // Calculate the significance 
   significance = ((Double_t)S) / sqrt((Double_t)(S+B));
   SoverB = ((Double_t)S) / ((Double_t)B); 
   std::cout << "\nTotal signal events: " << totS << ", total background events: " << totB << std::endl;
   std::cout << "Cuts on invariant mass: interval I = [" << cutL << ", " << cutR << "] GeV" << std::endl;
   std::cout << "Cut on chi squared at C = " << cutChi << std::endl;

   std::cout << "The statistical significance after the bi-variate cut is S/sqrt(S+B) = " << significance<< std::endl;
   std::cout << "The signal to background ratio after the bi-variate cut is S/B = " << SoverB << "\n" << std::endl;




}
