#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCut.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>

//
// A steering Root script for the SHMS calorimeter calibration.
//


 //void dc_cal_hms(Int_t runNo, Int_t eventNo) {
 // TFile *f =  new TFile(Form("../../ROOTfiles/hms_coin_replay_production_%d_%d.root", runNo,eventNo));

  //TTree *t = (TTree*)f->Get("T");
  
 void cal_cal_shms(int numruns=1, int RunNumber=1791, int nstop=50000) {

  Double_t D_CALO_FP= 275. ;
  Double_t XMIN= -60.;
  Double_t XMAX=  60.;
  Double_t YMIN= -58.;
  Double_t YMAX=  58.;

  Double_t DELTA_MIN= -10;   //SHMS nominal acceptance
  Double_t DELTA_MAX=  22;

  Double_t BETA_MIN= 0.5;
  Double_t BETA_MAX= 1.5;

  Double_t HGCER_NPE_MIN= 3;
  Double_t NGCER_NPE_MIN= 1.;   //



  //  TTree* fTree;
  UInt_t fNentries;
  UInt_t fNstart;
  UInt_t fNstop;
 
  fNstart = 0;
  fNstop = nstop;
  // Declaration of leaves types

  // Track parameters.


  Double_t        P_tr_x;   //X FP
  Double_t        P_tr_xp;
  Double_t        P_tr_y;   //Y FP
  Double_t        P_tr_yp;


  Double_t        P_hgcer_npe;
  Double_t        P_ngcer_npe;
  Double_t        P_tr_beta;

  Double_t        pdelta; 
  Double_t        pcaleprnorm;
  Double_t        pcaletotnorm;
  Double_t        pcalepr;
  Double_t        pcaletrack;


  TH1F* hEcal;
  TH1F* hEpr;
  TH2F* hDPvsEcal;
  TH2F* hESHvsEPR;

  // Initialize the analysis clock
  clock_t t = clock();
 
  cout << "Plot run " << RunNumber << endl;

//Whole calorimeter fid. limits



  //Reset ROOT and connect tree file.

  gROOT->Reset();

  TChain *fTree = new TChain("T");
  for (int i=0;i<numruns;i++){
   TString fname = "ROOTfiles/shms_coin_replay_production_";
   fname += RunNumber;
   fname +="_";
   fname +=nstop;
   fname +=".root";
   fTree->Add(fname);
   cout << "THcPShowerCalib::Init: Root file name = " << fname << endl;
  }

  fNentries = fTree->GetEntries();
  cout << "THcPShowerCalib::Init: fNentries= " << fNentries << endl;

  fTree->SetBranchAddress("P.dc.x_fp", &P_tr_x);
  fTree->SetBranchAddress("P.dc.y_fp", &P_tr_y);
  fTree->SetBranchAddress("P.dc.xp_fp",&P_tr_xp);
  fTree->SetBranchAddress("P.dc.yp_fp",&P_tr_yp);
 
  fTree->SetBranchAddress("P.hgcer.npeSum", &P_hgcer_npe);
  fTree->SetBranchAddress("P.ngcer.npeSum", &P_ngcer_npe);

  fTree->SetBranchAddress("P.gtr.beta", &P_tr_beta);
  fTree->SetBranchAddress("P.gtr.dp", &pdelta); 
  fTree->SetBranchAddress("P.cal.eprtracknorm", &pcaleprnorm);
  fTree->SetBranchAddress("P.cal.etottracknorm", &pcaletotnorm);
  fTree->SetBranchAddress("P.cal.eprtrack", &pcalepr); 
  fTree->SetBranchAddress("P.cal.etot", &pcaletrack);    
     
 
 TCut pid_cut;
      pid_cut = "P.gtr.beta>=0.8&&P.gtr.beta<=1.3&&P.hgcer.npeSum<=0.&&P.gtr.dp>=-15&&P.gtr.dp<=15&&P.cal.etottracknorm>=0.";

  // Histogram declarations.
  gStyle->SetOptFit(2);

    Int_t cnt=0;
    hEcal = new TH1F("hEcal", "Edep/P calibrated", 200, 0., 2.);
    hEpr = new TH1F("hEpr", "Epr calibrated", 200, 0., 2.);
    hDPvsEcal = new TH2F("hDPvsEcal", "#DeltaP versus Edep/P ",
  		       400,0.,2., 440,DELTA_MIN-1.,DELTA_MAX+1.);
    hESHvsEPR = new TH2F("hESHvsEPR", "E_{SH} versus E_{PR}",
		       300,0.0,1.4, 300,0.,1.4);

    for (UInt_t ientry=TMath::Max(UInt_t(0),fNstart);
       ientry<TMath::Min(fNstop,fNentries); ientry++) {

     fTree->GetEntry(ientry);


     if (ientry%100000 == 0) cout << "   ReadTree: " << ientry << endl;

  // Request single electron track in calorimeter's fid. volume.
  //


   bool good_trk =   pdelta > DELTA_MIN &&
		    pdelta < DELTA_MAX &&
		    P_tr_x + P_tr_xp*D_CALO_FP > XMIN &&
		    P_tr_x + P_tr_xp*D_CALO_FP < XMAX &&
                    P_tr_y + P_tr_yp*D_CALO_FP > YMIN &&
		    P_tr_y + P_tr_yp*D_CALO_FP < YMAX ;
   if (!good_trk) continue;

   bool good_ngcer = P_ngcer_npe > NGCER_NPE_MIN ;
   if(!good_ngcer) continue;

   bool good_hgcer = P_hgcer_npe > HGCER_NPE_MIN ;
   if(!good_hgcer) continue;

   bool good_beta = P_tr_beta > BETA_MIN &&
                   P_tr_beta < BETA_MAX ;
   if(!good_beta) continue;
   hEcal->Fill(pcaletotnorm);   
   hEpr->Fill(pcaleprnorm);
   hDPvsEcal->Fill(pcaletotnorm,pdelta);
   hESHvsEPR->Fill(pcaleprnorm,(pcaletotnorm-pcaleprnorm));
   cnt++;
   }
     cout << cnt << " passed all cuts  " << endl;


  // Plot histograms
   
  TCanvas* Canvas =
    new TCanvas("Canvas", "SHMS Shower Counter performance", 1000, 667);
  Canvas->Divide(2,2);


  Canvas->cd(1);
  hESHvsEPR->Draw("colz");
  // Normalized energy deposition after calibration.

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(2);
  Canvas->cd(2);

  hEcal->Fit("gaus","","",0.9,1.1);
  hEcal->GetFunction("gaus")->SetLineColor(2);
  hEcal->GetFunction("gaus")->SetLineWidth(1);
  hEcal->GetFunction("gaus")->SetLineStyle(1);

  // SHMS delta(P) versus the calibrated energy deposition.

  Canvas->cd(3);
  hDPvsEcal->Draw("colz");

  //  Canvas->cd(4);
  //hEpr->Draw();

  // Save canvas in a pdf format.
  Canvas->Print(Form("%d_%d.pdf",RunNumber,nstop));

  // Calculate the analysis rate
  t = clock() - t;
  printf ("The analysis took %.1f seconds", ((float) t) / CLOCKS_PER_SEC);
}

