

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

// The event tree containing labelled MC data with chisq, energy and mcSignal features
const TString fname = "./inclmc12.root";

Double_t cutLow = 0, cutHigh = 200;     //cut limits
Double_t cutStep = 0.01;
Double_t sEff, bRej;
auto nCuts = (Int_t)((cutHigh - cutLow) / cutStep);
auto cut = new Double_t[nCuts];

auto S = new Int_t[nCuts];
auto B = new Int_t[nCuts];
auto significance = new Double_t[nCuts];

auto hs = new TH1F("signal", "signal", 100, 2.95, 3.2);
auto hb = new TH1F("bakground", "background", 100, 2.95, 3.2);


void chisq_scan(){

  for(int i = 0; i < nCuts; i++){
      S[i] = 0;
      B[i] = 0;
      cut[i] = cutLow;
      cutLow += cutStep;
    }

  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  Float_t chi, energy;
   reader -> AddVariable("chisq4C", &chi);
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
   Double_t dchi = chi, denergy = energy;
   Int_t signal;
   TTree* theTree = (TTree*)input->Get("ntp1");
   theTree->SetBranchAddress( "chisq4C",       &dchi    );
   theTree->SetBranchAddress( "sqrt_s",       &denergy   );
   theTree->SetBranchAddress( "mcSignal",     &signal     );

   Int_t totS = 0, totB = 0;


   // Loop over the events in the tree
   for(int ievt = 0; ievt < theTree->GetEntries(); ievt++){
     theTree->GetEntry(ievt);
     if(ievt % 10000 == 0) std::cout << "---Processing event " << ievt << std::endl;
     chi = dchi; energy = denergy;
     for(int i = 0; i < nCuts; i++){
       if(chi < cut[i]){    // if the event fall within the cut window
	 if(signal == 1){
	   S[i]++;    /// if it is signal classify it as such
	 }
	 if(signal == 0){
	   B[i]++;    /// otherwise as background
	 }
       }
     }
     if(signal == 0){
       totB++;
     }
     if(signal == 1){
       totS++;
     }
   }

  
   // Calculate the significance of each cut
   for(int y = 0; y < nCuts; y++){  
     significance[y] = ((double)S[y]) / sqrt((double)(S[y]+B[y]));
   }

   // Find the best significance
   Double_t bestCut, bestS = 0, SoverB;
   Int_t cutIndex;
   for(int i = 0; i < nCuts; i++){
     if(significance[i] > bestS){
       bestS = significance[i];
       sEff = (double)S[i] / (double)totS;
       bRej = 1 - (double)B[i] / (double)totB;
       SoverB = ((Double_t)S[i]) / ((Double_t)B[i]);
       cutIndex = i;
     }
   }
   bestCut = cut[cutIndex];

   cout << "Filling histogram..." << endl;
   for(int ievt = 0; ievt < theTree->GetEntries(); ievt++){
     theTree->GetEntry(ievt);
     chi = dchi; energy = denergy;
     if(chi < bestCut){
       if(signal == 0){
	 hb->Fill(energy);
       }
       if(signal == 1){
	 hs->Fill(energy);
       }
     }
   }
   std::cout << "\nTotal signal events: " << totS << ", total background events: " << totB << std::endl;
   std::cout << "The statistical significance before the cut is S/sqrt(S+B) = " << ((double)totS) / sqrt(((double)totS)+((double)totB)) << std::endl;
   std::cout << "The best cut is at C = " << bestCut << std::endl;
   std::cout << "The statistical significance after the cut is S/sqrt(S+B) = " << bestS << std::endl;
   std::cout << "The signal to background ratio after the cut is S/B = " << SoverB << "\n" << std::endl;
   std::cout << "The signal efficiency is " << sEff << std::endl;
   std::cout << "The background rejection is " << bRej << std::endl;


   TApplication *app = new TApplication("app",0,NULL);
   TCanvas c1;
   c1.cd();
   TGraph *signi = new TGraph(nCuts, cut, significance);
   signi -> SetLineWidth(2);
   signi -> SetLineColor(kBlack);
   signi -> Draw("ALP");
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
