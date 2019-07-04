/*
ROOT macro for classification of MC data. This includes a scan over a range of cuts for maximising the
efficiency of classification. This is possible since MC data has labels. Once the best cut value is found,
real data can be classified using the macro (still non-existing) RDclassification.C  (=RealDataClassification)

This macro also prepares the output of the MVAs to be used as inputs for still another training and 
classification step.

Classification macro to be used with algorithms trained on the 4-momenta feature space.
*/

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

//Change this depending on the feature you're using:
int featSpace = 3;
// Also change these according to the feature space and the trained algorithms
// Some strings are initialised here for ease of change, including the algorithm weight files and the 
// input and output file names
const TString weightsBDT = "/home/besuser1/Tommaso/MC/finalfinalRound/training/dataset/weights/TMVAfactory_BDT3.weights.xml";
const TString weightsANN = "/home/besuser1/Tommaso/MC/finalfinalRound/training/dataset/weights/TMVAfactory_ANN3.weights.xml";
const TString fname = "/home/besuser1/Tommaso/MC/data/finalData/allFeatNew/inclmc12.root";
const TString outFileName = "/home/besuser1/Tommaso/MC/finalfinalRound/classification/cl3.root";
// Cut array values
Double_t BDT_cut_min = -0.8, BDT_cut_step = 0.01, BDT_cut_max = 0.12;
Double_t ANN_cut_min = 0, ANN_cut_step = 0.01, ANN_cut_max = 0.92;
const Int_t BDT_n_cuts = (Int_t)((BDT_cut_max - BDT_cut_min) / BDT_cut_step), ANN_n_cuts = (Int_t)((ANN_cut_max - ANN_cut_min) / ANN_cut_step);
Double_t TMVA_BEST_CUT_ANN = 0.29, TMVA_TRUE_SIG = 0, TMVA_TRUE_BKG = 0;


// The function to calculate the center of mass energy given the relevant four momenta
Double_t comEnergy(Double_t pxp, Double_t pyp, Double_t pzp,  Double_t ep, Double_t pxm, Double_t pym, Double_t pzm,  Double_t em, Double_t pxg, Double_t pyg, Double_t pzg,  Double_t eg){
  
  Double_t pxt   = pxp + pxm + pxg;
  Double_t pyt   = pyp + pym + pyg;
  Double_t pzt   = pzp + pzm + pzg;
  Double_t et    = ep + em + eg;
  Double_t s     = et * et - ( pxt * pxt + pyt * pyt + pzt * pzt );
  Double_t sqrtS = sqrt(s);
  return sqrtS;
}


//The function to calculate the missing mass, assuming the e+ and e- are colliding head-on
Double_t missMass(Double_t pxp, Double_t pyp, Double_t pzp,  Double_t ep, Double_t pxm, Double_t pym, Double_t pzm,  Double_t em){

  Double_t mJPsi = 3.097;
  Double_t sumX = pxp + pxm;
  Double_t sumY = pyp + pym;
  Double_t sumZ = pzp + pzm;
  Double_t sumE = ep + em - mJPsi;
  Double_t sqMissMass = sumE * sumE - sumX * sumX - sumY * sumY - sumZ * sumZ;
  return sqrt(sqMissMass);

}









// The actual classification macro
void cl( TString myMethodList = ""){

  // Uninteresting code---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;
  
  
  
   // Neurl Networks
   int fANN          = 1; 
   // Boosted Decision Trees
   int fBDT             = 0; 

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;


   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
     
      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
	   std::cout << "Method \"" << regMethod
		     << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	     std::cout << it->first << " ";
	   }
	   std::cout << std::endl;
	   return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------


   // Variables to gather statistics about the classification. These will be used to
   // find the optimal cut and to plot the ROC curve and other stuff.
   //// Arrays containing the cut values. 
   Double_t *BDT_cuts = new Double_t[BDT_n_cuts];
   Double_t *ANN_cuts = new Double_t[ANN_n_cuts];
   for(int i = 0; i < BDT_n_cuts; i++){
     BDT_cuts[i] = BDT_cut_min;
     BDT_cut_min += BDT_cut_step;
   }
   for(int i = 0; i < ANN_n_cuts; i++){
     ANN_cuts[i] = ANN_cut_min;
     ANN_cut_min += ANN_cut_step;
   }
   //// Arrays containing the number of true/false positives/negatives for both classifiers
   //// and for each cut value. They are of course filled with zeroes.
   Int_t *truePosBDT = new Int_t[BDT_n_cuts];
   Int_t *trueNegBDT = new Int_t[BDT_n_cuts];
   Int_t *truePosANN = new Int_t[ANN_n_cuts];
   Int_t *trueNegANN = new Int_t[ANN_n_cuts];
   Int_t *falsePosBDT = new Int_t[BDT_n_cuts];
   Int_t *falseNegBDT = new Int_t[BDT_n_cuts];
   Int_t *falsePosANN = new Int_t[ANN_n_cuts];
   Int_t *falseNegANN = new Int_t[ANN_n_cuts];
   Int_t *totPosBDT = new Int_t[BDT_n_cuts];
   Int_t *totNegBDT = new Int_t[BDT_n_cuts];
   Int_t *totPosANN = new Int_t[BDT_n_cuts];
   Int_t *totNegANN = new Int_t[BDT_n_cuts];
 for(int i = 0; i < BDT_n_cuts; i++){
     truePosBDT[i] = 0;
     trueNegBDT[i] = 0;
     falsePosBDT[i] = 0;
     falseNegBDT[i] = 0;
     totPosBDT[i] = 0;
     totNegBDT[i] = 0;
   } 
   for(int i = 0; i < ANN_n_cuts; i++){
     truePosANN[i] = 0;
     trueNegANN[i] = 0;
     falsePosANN[i] = 0;
     falseNegANN[i] = 0;
     totPosANN[i] = 0;
     totNegANN[i] = 0;
   } 
   // Arrays to gather statistics for the ROC curve and for the cuts; no need to initialise.
   Double_t *sigEffBDT = new Double_t[BDT_n_cuts];
   Double_t *bkgRejBDT = new Double_t[BDT_n_cuts];
   Double_t *sigEffANN = new Double_t[ANN_n_cuts];
   Double_t *bkgRejANN = new Double_t[ANN_n_cuts];
   Double_t *significanceBDT = new Double_t[BDT_n_cuts];  //estimator for bdt
   Double_t *significanceANN = new Double_t[ANN_n_cuts];  //estimator for ann

   // Other statistic variables
   Double_t best_cut_BDT = 0;
   Double_t best_cut_ANN = 0;
   Int_t best_cut_index_BDT = 0;
   Int_t best_cut_index_ANN = 0;
   Int_t BACKGROUND_TOTAL = 0;
   Int_t SIGNAL_TOTAL = 0;







   // Create the MVA Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   Float_t Pp_px, Pp_py, Pp_pz, Pp_e, Pm_px, Pm_py, Pm_pz, Pm_e;
   Float_t gamma_px, gamma_py, gamma_pz, gamma_e, energy, chi;

   reader -> AddVariable("if4CPp_px",    &Pp_px);
   reader -> AddVariable("if4CPp_py",    &Pp_py);
   reader -> AddVariable("if4CPp_pz",    &Pp_pz);
   reader -> AddVariable("if4CPp_e",     &Pp_e);
   reader -> AddVariable("if4CPm_px",    &Pm_px);
   reader -> AddVariable("if4CPm_py",    &Pm_py);
   reader -> AddVariable("if4CPm_pz",    &Pm_pz);
   reader -> AddVariable("if4CPm_e",     &Pm_e);
   reader -> AddVariable("if4Cgamma_px", &gamma_px);
   reader -> AddVariable("if4Cgamma_py", &gamma_py);
   reader -> AddVariable("if4Cgamma_pz", &gamma_pz);
   reader -> AddVariable("if4Cgamma_e",  &gamma_e);
   if(featSpace == 2){
     reader -> AddVariable("sqrt_s",  &energy);
   }
   if(featSpace == 3){
     reader -> AddVariable("chisq4C",  &chi);
   };
   // Book method(s)

   if (fANN)   reader->BookMVA( "ANN", weightsANN );
   if (fBDT)      reader->BookMVA( "BDT",    weightsBDT );
  



   // Book output histograms and graphs for significance and S/B ratio
   UInt_t nbin = 800;    //in order to be able to rebin later
   TH1F *histANN(0);         //contains the full MVA output
   TH1F *histSigANN(0);      //contains the classified signal 
   TH1F *histBkgANN(0);      //contains the classified background
   TH1F *histOutSigANN(0);   //contains the MC signal in the MVA output
   TH1F *histOutBkgANN(0);//contains the MC background in the MVA output
   TH1F *histBDT(0);
   TH1F *histSigBDT(0);
   TH1F *histBkgBDT(0);
   TH1F *histOutSigBDT(0);
   TH1F *histOutBkgBDT(0);
   Double_t Emin = 2.8, Emax = 3.5;   //in GeV
   TH1F *totData = new TH1F("total_Data", "total_Data", nbin, Emin, Emax);
   TH1F *true_com_sig = new TH1F("true data signal", "true data signal", nbin, Emin, Emax);
   TH1F *true_com_bkg = new TH1F("true data background", "true data background", nbin, Emin, Emax);
   TH1F *after_cut_sig_BDT = new TH1F("signal after cut for BDT", "signal after cut for BDT", nbin, 3.0, 3.15);
   TH1F *after_cut_bkg_BDT = new TH1F("background after cut for BDT", "background after cut for BDT", nbin, 3.0, 3.15);
   TH1F *after_cut_sig_ANN = new TH1F("signal after cut for ANN", "signal after cut for ANN", nbin, 3.0, 3.15);
   TH1F *after_cut_bkg_ANN = new TH1F("background after cut for ANN", "background after cut for ANN", nbin, 3.0, 3.15);
   TH1F *after_cut_chi_sig_ANN = new TH1F("chi signal after cut ANN", "chi signal after cut ANN", nbin, 0, 200);
   TH1F *after_cut_chi_bkg_ANN = new TH1F("chi bkg after cut ANN", "chi bkg after cut ANN", nbin, 0, 200);
   TH1F *after_cut_chi_sig_BDT = new TH1F("chi signal after cut BDT", "chi signal after cut BDT", nbin, 0, 200);
   TH1F *after_cut_chi_bkg_BDT = new TH1F("chi bkg after cut BDT", "chi bkg after cut BDT", nbin, 0, 200);
   TH1F *missMassSigANN = new TH1F("missMassSigANN", "missMassSigANN", nbin, -.1, .5);
   TH1F *missMassBkgANN = new TH1F("missMassBkgANN", "missMassBkgANN", nbin, -.1, .5);
   TH1F *missMassSigBDT = new TH1F("missMassSigBDT", "missMassSigBDT", nbin, -.1, .5);
   TH1F *missMassBkgBDT = new TH1F("missMassBkgBDT", "missMassBkgBDT", nbin, -.1, .5);
   TGraph *SignificanceANN = new TGraph();
   TGraph *ratioANN = new TGraph();
   TGraph *SignificanceBDT = new TGraph();
   TGraph *ratioBDT = new TGraph();
   TGraph *rocANN = new TGraph();
   TGraph *rocBDT = new TGraph();

   //for each MVA there are 3 histograms: one for the classifier response, and one each for plotting the energy of SIG and BKG events 
   if (fANN)           histANN     = new TH1F( "ANN-output",           "ANN-output",           nbin, 0, 1 );

   if (fANN)      histOutSigANN   = new TH1F( "ANN-output-signal",        "ANN-output-signal",        nbin, 0, 1.0);
   if (fANN)      histOutBkgANN   = new TH1F( "ANN-output-background",        "ANN-output-background",        nbin, 0, 1.0);
   if (fBDT)       histOutSigBDT   = new TH1F( "BDT-output-signal",        "BDT-output-signal",        nbin, -0.4, 0.2);
   if (fBDT)       histOutBkgBDT   = new TH1F( "BDT-output-background",        "BDT-output-background",        nbin, -0.4, 0.2);
   if (fANN)        histSigANN   = new TH1F( "ANN-signal",    "ANN-signal",    nbin, Emin, Emax);
   if (fANN)        histBkgANN   = new TH1F( "ANN-background",    "ANN-background",    nbin, Emin, Emax );
   if (fBDT)           histBDT     = new TH1F( "BDT-output",           "BDT-output",           nbin, -0.4, 0.2 );
   if (fBDT)           histSigBDT  = new TH1F( "BDT-signal",       "BDT-signal",       nbin, Emin, Emax );
   if (fBDT)           histBkgBDT  = new TH1F( "BDT-background",       "BDT-background",       nbin, Emin,Emax );

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
   std::cout << "--- Select signal sample" << std::endl;
   Double_t dPp_px=Pp_px, dPp_py=Pp_py, dPp_pz=Pp_pz, dPp_e=Pp_e, dPm_px=Pm_px, dPm_py=Pm_py, dPm_pz=Pm_pz, dPm_e=Pm_e;
   Double_t dgamma_px=gamma_px, dgamma_py=gamma_py, dgamma_pz=gamma_pz, dgamma_e=gamma_e, denergy=energy, dchi=chi;
   Int_t signal;

   TTree* theTree = (TTree*)input->Get("ntp1");
   theTree->SetBranchAddress( "if4CPp_px",   &dPp_px    );
   theTree->SetBranchAddress( "if4CPp_py",   &dPp_py    );
   theTree->SetBranchAddress( "if4CPp_pz",   &dPp_pz    );
   theTree->SetBranchAddress( "if4CPp_e",    &dPp_e     );
   theTree->SetBranchAddress( "if4CPm_px",   &dPm_px    );
   theTree->SetBranchAddress( "if4CPm_py",   &dPm_py    );
   theTree->SetBranchAddress( "if4CPm_pz",   &dPm_pz    );
   theTree->SetBranchAddress( "if4CPm_e",    &dPm_e     );
   theTree->SetBranchAddress( "if4Cgamma_px", &dgamma_px );
   theTree->SetBranchAddress( "if4Cgamma_py", &dgamma_py );
   theTree->SetBranchAddress( "if4Cgamma_pz", &dgamma_pz );
   theTree->SetBranchAddress( "if4Cgamma_e",  &dgamma_e  );
   theTree->SetBranchAddress( "mcSignal",   &signal     );
   theTree->SetBranchAddress( "sqrt_s", &denergy);
   theTree->SetBranchAddress("chisq4C", &dchi);



   // Now the first event loop, to find the optimal cut value and to gather statistics on the cuts.
   // Histograms will be filled later on in another event loop.
   for(Long64_t ievt = 0; ievt < theTree->GetEntries(); ievt++){

     // Get the entry and assign the variables
     theTree->GetEntry(ievt);     
     Pp_px = dPp_px; Pp_py = dPp_py; Pp_pz = dPp_pz; Pp_e=dPp_e; Pm_px=dPm_px;
     Pm_py = dPm_py; Pm_pz = dPm_pz; Pm_e = dPm_e;
     gamma_px = dgamma_px; gamma_py=dgamma_py; gamma_pz = dgamma_pz;
     gamma_e = dgamma_e;
     if(featSpace == 2) energy = denergy; 
     if(featSpace == 3) chi = dchi;
     // Calculate sqrt(s) for the event and fill in the general histogram (SIG+BKG)
     Double_t comEn = denergy;
     totData -> Fill( comEn );
     if(signal == 0){
       BACKGROUND_TOTAL++;
       true_com_bkg -> Fill(comEn);
     }
     if(signal == 1){
       SIGNAL_TOTAL++;
       true_com_sig -> Fill(comEn);
     }

     // Statistics gathering
     if(fANN){
       double score = reader -> EvaluateMVA("ANN");   //get the MVA output for the current event
       //now loop on cut values and gather statistics
       histANN->Fill(score);
       for(int i = 0; i < ANN_n_cuts; i++){
	 if(score > ANN_cuts[i]){  // if this cut classifies the event as signal
	   totPosANN[i]++;         /// the event is a positive
	   if(signal == 1){        /// if the event is actually signal
	     truePosANN[i]++;      //// the event is a true positive
	   }
	   if(signal == 0){        /// if the event is actually background
	     falsePosANN[i]++;     //// the event is a false positive
	   }
	 }
	 else{                     // if this cut classifies the event as background
	   totNegANN[i]++;         /// the event is a negative
	   if(signal == 1){        /// if the event is actually signal
	     falseNegANN[i]++;     //// the event is a false negative
	   }
	   if(signal == 0){        /// if the event is actually background
	     trueNegANN[i]++;      //// the event is a true negative
	   }
	 }
       }
     }

     if(fBDT){
       double score = reader -> EvaluateMVA("BDT");   //get the MVA output for the current event
       histBDT -> Fill(score); // fill the MVA output histogram
       //now loop on cut values and gather statistics
       for(int i = 0; i < BDT_n_cuts; i++){
	 if(score > BDT_cuts[i]){  // if this cut classifies the event as signal
	   totPosBDT[i]++;         /// the event is a positive
	   if(signal == 1){        /// if the event is actually signal
	     truePosBDT[i]++;      //// the event is a true positive
	   }
	   if(signal == 0){        /// if the event is actually background
	     falsePosBDT[i]++;     //// the event is a false positive
	   }
	 }
	 else{                     // if this cut classifies the event as background
	   totNegBDT[i]++;         /// the event is a negative
	   if(signal == 1){        /// if the event is actually signal
	     falseNegBDT[i]++;     //// the event is a false negative
	   }
	   if(signal == 0){        /// if the event is actually background
	     trueNegBDT[i]++;      //// the event is a true negative
	   }
	 }
       }
     }

     // fill in the output trees
     if (ievt%10000 == 0) std::cout << "--- ... Processed event: " << ievt << std::endl;
   }


   // Now find the best cut by looking at the statistics
   for(int i = 0; i < BDT_n_cuts; i++){
     significanceBDT[i] = (double)truePosBDT[i] / sqrt((double)totPosBDT[i]);
     SignificanceBDT -> SetPoint(i, BDT_cuts[i], significanceBDT[i]);
     ratioBDT -> SetPoint(i, BDT_cuts[i], (double)truePosBDT[i] / (double)falsePosBDT[i]);
     bkgRejBDT[i] = (double)trueNegBDT[i] / (double)(falsePosBDT[i] + trueNegBDT[i]);
     sigEffBDT[i] = (double)truePosBDT[i] / (double)(truePosBDT[i] + falseNegBDT[i]);
     rocBDT -> SetPoint(i, sigEffBDT[i], bkgRejBDT[i]);
     if(significanceBDT[i] >= significanceBDT[best_cut_index_BDT]){ // find the maximum significance
       best_cut_index_BDT = i;
     }
   }

   best_cut_BDT = BDT_cuts[best_cut_index_BDT];
   for(int i = 0; i < ANN_n_cuts; i++){
     significanceANN[i] = (double)truePosANN[i] / sqrt((double)totPosANN[i]);
     SignificanceANN -> SetPoint(i, ANN_cuts[i], significanceANN[i]);
     ratioANN -> SetPoint(i, ANN_cuts[i], (double)truePosANN[i] / (double)falsePosANN[i]);
     bkgRejANN[i] = (double)trueNegANN[i] / (double)(falsePosANN[i] + trueNegANN[i]);
     sigEffANN[i] = (double)truePosANN[i] / (double)(truePosANN[i] + falseNegANN[i]);
     rocANN -> SetPoint(i, sigEffANN[i], bkgRejANN[i]);
     if(significanceANN[i] >= significanceANN[best_cut_index_ANN]){ // find the maximum significance
       best_cut_index_ANN = i;
     }
   }
   best_cut_ANN = ANN_cuts[best_cut_index_ANN];

   // Print some statistics 
   std::cout << "\nStatistics for the classification-----------------------------\n" << std::endl;
   std::cout << "Total signal events " << SIGNAL_TOTAL << std::endl;
   std::cout << "Total background events " << BACKGROUND_TOTAL << "\n" << std::endl;
   std::cout << "MVA method: BDT"<< std::endl;
   std::cout << "Best cut "<< best_cut_BDT << std::endl;
   std::cout << "S = "<< truePosBDT[best_cut_index_BDT] << std::endl;
   std::cout << "B = "<< falsePosBDT[best_cut_index_BDT] << std::endl;
   std::cout << "Statistical significance " << significanceBDT[best_cut_index_BDT] << std::endl;
   std::cout << "Background rejection "<< bkgRejBDT[best_cut_index_BDT] << std::endl;
   std::cout << "Signal efficiency " << sigEffBDT[best_cut_index_BDT] << "\n" << std::endl;
   std::cout << "MVA method: ANN"<< std::endl;
   std::cout << "Best cut "<< best_cut_ANN << std::endl;
   std::cout << "S = "<< truePosANN[best_cut_index_ANN] << std::endl;
   std::cout << "B = "<< falsePosANN[best_cut_index_ANN] << std::endl;
   std::cout << "Statistical significance " << significanceANN[best_cut_index_ANN] << std::endl;
   std::cout << "Background rejection "<< bkgRejANN[best_cut_index_ANN] << std::endl;
   std::cout << "Signal efficiency " << sigEffANN[best_cut_index_ANN] << "\n" << std::endl;




 std::cout << "\n\nNow filling sqrt(s) histograms " << std::endl;

   // Now another event loop with the best cut values to fill the sqrt(s) histograms

 for (Long64_t ievt=0; ievt< theTree->GetEntries();ievt++) {

     if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     // Get the entry and assign the variables
     theTree->GetEntry(ievt);     
     Pp_px = dPp_px; Pp_py = dPp_py; Pp_pz = dPp_pz; Pp_e=dPp_e; Pm_px=dPm_px;
     Pm_py = dPm_py; Pm_pz = dPm_pz; Pm_e = dPm_e;
     gamma_px = dgamma_px; gamma_py=dgamma_py; gamma_pz = dgamma_pz;
     gamma_e = dgamma_e;
     if(featSpace == 2) energy = denergy; 
     if(featSpace == 3) chi = dchi;

     // Calculate sqrt(s) for the event and fill in the general histogram (SIG+BKG)
     Double_t comEn = comEnergy(Pp_px, Pp_py, Pp_pz, Pp_e, Pm_px, Pm_py,Pm_pz, Pm_e, gamma_px, gamma_py, gamma_pz, gamma_e);
     // Calculate the missing mass
     Double_t missM = missMass(Pp_px, Pp_py, Pp_pz, Pp_e, Pm_px, Pm_py,Pm_pz, Pm_e);
     if(fBDT){
       Double_t score = reader -> EvaluateMVA("BDT");
       if(score > best_cut_BDT){
	 histSigBDT -> Fill(comEn);
	 if(signal==1){
	   after_cut_sig_BDT->Fill(comEn);
	   after_cut_chi_sig_BDT->Fill(dchi);
	   histOutSigBDT -> Fill(score); // fill the MVA output histogram
	   missMassSigBDT->Fill(missM);
	 }
	 if(signal==0){
	   after_cut_bkg_BDT->Fill(comEn);
	   after_cut_chi_bkg_BDT->Fill(dchi);
	   histOutBkgBDT -> Fill(score); // fill the MVA output histogram
	   missMassBkgBDT->Fill(missM);

	 }
       }
       else{
	 histBkgBDT -> Fill(comEn);
       }     
     }

     if(fANN){
       Double_t score = reader -> EvaluateMVA("ANN");
       if(score > best_cut_ANN){
	 histSigANN -> Fill(comEn);
	 if(signal==1){
	   after_cut_sig_ANN->Fill(comEn);
	   after_cut_chi_sig_ANN->Fill(dchi);
	   histOutSigANN -> Fill(score); // fill the MVA output histogram
	   missMassSigANN->Fill(missM);

	 }
	 if(signal==0){
	   after_cut_bkg_ANN->Fill(comEn);
	   after_cut_chi_bkg_ANN->Fill(dchi);
	   histOutBkgANN -> Fill(score); // fill the MVA output histogram
	   missMassBkgANN->Fill(missM);

	 }
       }
       else{
	 histBkgANN -> Fill(comEn);
       }     
     }


   }


 // Write the control histograms

 TFile *target  = new TFile( outFileName ,"RECREATE" );
 totData -> Write();
 true_com_sig ->Write();
 true_com_bkg ->Write();
 
 if (fANN){    
   histANN  ->Write();
   histSigANN  ->Write();
   histOutSigANN  ->Write();
   histOutBkgANN  ->Write();
   histBkgANN  ->Write();
   after_cut_sig_ANN ->Write();
   after_cut_bkg_ANN ->Write();
   after_cut_chi_sig_ANN ->Write();
   after_cut_chi_bkg_ANN ->Write();
   missMassSigANN->Write();
   missMassBkgANN->Write();
   SignificanceANN -> SetName("ANN significance");
   SignificanceANN -> Write();
   ratioANN -> SetName("ANN ratio");
   ratioANN -> Write();
   rocANN -> SetName("roc ANN");
   rocANN -> Write();
 }

 if (fBDT){
   histBDT    ->Write();
   histSigBDT    ->Write();
   histBkgBDT    ->Write();
   histOutSigBDT  ->Write();
   histOutBkgBDT  ->Write();
   after_cut_sig_BDT ->Write();
   after_cut_bkg_BDT ->Write();
   after_cut_chi_sig_BDT ->Write();
   after_cut_chi_bkg_BDT ->Write();
   missMassSigBDT->Write();
   missMassBkgBDT->Write();
   SignificanceBDT -> SetName("BDT significance");
   SignificanceBDT -> Write();
   ratioBDT -> SetName("BDT ratio");
   ratioBDT -> Write();
   rocBDT -> SetName("roc BDT");
   rocBDT -> Write();
 }



   std::cout << "\n\n\n--- Created root file: " + outFileName + " containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;

}


// I don't precisely know why this is here
int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   cl(methodList);
   return 0;
}
