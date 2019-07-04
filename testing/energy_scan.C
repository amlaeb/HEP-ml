

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
const TString fname = "./inclmc12.root";

Double_t cutLow = 2.5, peak = 3.097, cutHigh = 3.5;    //cut limits in GeV
Double_t cutStep = 0.005;
Int_t nCutsL = (Int_t)((peak - cutLow) / cutStep);
Int_t nCutsR = (Int_t)((cutHigh - peak) / cutStep);
auto cutLeft  = new Double_t[nCutsL];
auto cutRight = new Double_t[nCutsR];

auto S = new Int_t[nCutsL * nCutsR];
auto B = new Int_t[nCutsL * nCutsR];
auto sigEff = new Double_t[nCutsL * nCutsR];
auto bkgRej = new Double_t[nCutsL * nCutsR];
auto significance = new Double_t[nCutsL * nCutsR];




void energy_scan(){

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

   // some more statistics
   Int_t totS = 0, totB = 0;


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
     if(signal == 0) totB++;
     if(signal == 1) totS++;
   }

   // Calculate the significance of the intervals
   auto signi = new TGraph2D();   
   x = 0;
   for(int l = 0; l < nCutsL; l++){  //for any possible cut interval
     for(int r = 0; r < nCutsR; r++){ 
       significance[x] = ((double)S[x]) / sqrt((double)(S[x]+B[x]));
       sigEff[x] = (double)S[x] / (double)totS;
       bkgRej[x] = (double)(totB-B[x])/ (double)totB;
       signi -> SetPoint(x, cutLeft[l], cutRight[r], significance[x]);
       x++;
     }
   }
   
   Double_t sEff, bRej;
   // Find the best significance
   Double_t bestCutL, bestCutR, bestS = 0;
   Int_t L, R;
   Double_t SoverB;
   x = 0;
   for(int l = 0; l < nCutsL; l++){  
     for(int r = 0; r < nCutsR; r++){
       if(significance[x] > bestS){
	 bestS = significance[x];
	 sEff = sigEff[x];
	 bRej = bkgRej[x];
	 SoverB = ((Double_t)S[x]) / ((Double_t)B[x]);
	 L = l;
	 R = r;
       }
	 x++;
     }
   }
   bestCutL = cutLeft[L];
   bestCutR = cutRight[R];

   auto hs = new TH1F("signal", "signal", 100, bestCutL, bestCutR);
   auto hb = new TH1F("bakground", "background", 100, bestCutL, bestCutR);

   cout << "Filling histogram..." << endl;
   for(int ievt = 0; ievt < theTree->GetEntries(); ievt++){
     theTree->GetEntry(ievt);
     if(ievt % 10000 == 0) std::cout << "---Processing event " << ievt << std::endl;
     energy = denergy;
     if(energy < bestCutR && energy > bestCutL){
       if(signal == 0) hb -> Fill(energy);
       if(signal == 1) hs -> Fill(energy);
     }
   }
   std::cout << "\nTotal signal events: " << totS << ", total background events: " << totB << std::endl;
   std::cout << "The best separating interval is I = [" << bestCutL << ", " << bestCutR << "] GeV" << std::endl;
   std::cout << "The statistical significance after the cut is S/sqrt(S+B) = " << bestS << std::endl;
   std::cout << "The signal to background ratio after the cut is S/B = " << SoverB << "\n" << std::endl;
   cout << "The signal efficiency is " << sEff << endl;
   cout << "The background rejection is " << bRej << endl; 

   TApplication *app = new TApplication("app",0,NULL);
   TCanvas c1;
   c1.cd();
   signi -> SetLineWidth(2);
   signi -> Draw("TRI2");
   TCanvas c2;
   c2.cd();
   hs->SetLineColor(kRed);
   hs->SetFillColor(kRed);
   hs->SetFillStyle(3004);
   hb->SetLineColor(kBlue);
   hb->SetFillColor(kBlue);
   hb->SetFillStyle(3005);
   auto STACK = new THStack;
   STACK->Add(hs);
   STACK->Add(hb);
   STACK->Draw("NOSTACK");
   app -> Run(true);




}
