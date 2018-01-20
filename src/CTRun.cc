// Filename: CTRun.cc
// Description: Main CTRun class. 
// Author: Latiful Kabir < latif@jlab.org >
// Created: Fri Apr 10 14:50:00 2017 (-0400)
// URL: jlab.org/~latif

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TSystem.h"
#include "Constants.h"
#include "CTRun.h"

using namespace std;

ClassImp(CTRun)

CTRun::CTRun(int runNumber)
{
    fRunNumber = runNumber;
    fFirstRun = runNumber;
    fLastRun = runNumber;
    fRunCounter = 0;
    Init();
}


CTRun::CTRun(int firstRun, int lastRun)
{
    fRunNumber = -1;
    fFirstRun = firstRun;
    fLastRun = lastRun;
    fRunCounter = 0;
    Init();
}


CTRun::CTRun()
{
    
}
CTRun::~CTRun()
{
    fFile->Close();
}

void CTRun::Init()
{
    fChain = new TChain("ct_chain","CT Chain");

    Bool_t runNotFound = kTRUE;
    fRunExist = kFALSE;
    for(Int_t run = fFirstRun; run <= fLastRun; ++run)
    {
	fFileName = ROOT_FILE_PATH;
	fFileName += FILE_PREFIX;
	fFileName += run;
	fFileName +=  ".root";

	runNotFound = gSystem->AccessPathName(fFileName);
	if(runNotFound)
	{	    
	    cout << "Requested run " << run << " does NOT exist" <<endl;
	    cout << "File name with complete path: "<<fFileName<<endl;
	    continue;
	}
	fFileName +=  "/T";	
	fChain->Add(fFileName);
	++fRunCounter;
	cout << "Adding the run number:"<< run <<endl;
    }
    if(fRunCounter < 1)
	return;

    SetBranchAddressCT();
    fRunExist = kTRUE;
}

void CTRun::Print()
{
    if(!fRunExist)
	return;
    fChain->Print();
}

TChain* CTRun::GetChain()
{
    if(!fRunExist)
	return NULL;
    return fChain;
}


TTree* CTRun::GetTree()
{
    if(!fRunExist)
	return NULL;
    return (TTree*)fChain;
}


Int_t CTRun::GetRunNumber()
{
    return fRunNumber;
}

void CTRun::ActivateCTBranches()
{
    fChain->SetBranchStatus("*",0);
    fChain->SetBranchStatus("P.kin.primary.*",1);
    fChain->SetBranchStatus("H.kin.secondary.*",1);
    fChain->SetBranchStatus("P.gtr.*",1);
    fChain->SetBranchStatus("H.gtr.*",1);

    fChain->SetBranchStatus("P.dc.x*",1);
    fChain->SetBranchStatus("P.dc.y*",1);

    fChain->SetBranchStatus("H.dc.x*",1);
    fChain->SetBranchStatus("H.dc.y*",1);

    fChain->SetBranchStatus("P.hgcer.*",1);
    fChain->SetBranchStatus("H.cer.*",1);

    fChain->SetBranchStatus("P.hod.*",1);
    fChain->SetBranchStatus("H.hod.*",1);

    fChain->SetBranchStatus("P.cal.*",1);
    fChain->SetBranchStatus("H.cal.*",1);

    fChain->SetBranchStatus("T.coin.pTRIG*",1);
    fChain->SetBranchStatus("T.coin.hTRIG*",1);
}


// For the buffer fCTEventData
void CTRun::SetBranchAddressCT()
{
   //=============================== Set branch addresses.==================
    //------------ Kinematics ----------------------------
    //SHMS
    fChain->SetBranchAddress("P.kin.primary.W", &fP_kin_Primary_W);
    fChain->SetBranchAddress("P.kin.primary.Q2", &fP_kin_Primary_Q2);
    fChain->SetBranchAddress("P.kin.primary.nu", &fP_kin_Primary_nu);

    //HMS
    fChain->SetBranchAddress("H.kin.secondary.pmiss", &fH_kin_secondary_pmiss);
    fChain->SetBranchAddress("H.kin.secondary.pmiss_x", &fH_kin_secondary_pmiss_x);
    fChain->SetBranchAddress("H.kin.secondary.emiss", &fH_kin_secondary_emiss);

    //------------ Golden Track Quantities ----------------
    //SHMS 
    fChain->SetBranchAddress("P.gtr.p", &fP_gtr_p);
    fChain->SetBranchAddress("P.gtr.dp", &fP_gtr_dp);
    fChain->SetBranchAddress("P.gtr.beta", &fP_gtr_beta);            
    
    //HMS
    fChain->SetBranchAddress("H.gtr.p", &fH_gtr_p);
    fChain->SetBranchAddress("H.gtr.dp", &fH_gtr_dp);
    fChain->SetBranchAddress("H.gtr.beta", &fH_gtr_beta);

    // ---------------- Drift Chamber ---------------
    //SHMS
    fChain->SetBranchAddress("P.dc.x_fp", &fP_dc_x_fp);
    fChain->SetBranchAddress("P.dc.xp_fp", &fP_dc_xp_fp);
    fChain->SetBranchAddress("P.dc.y_fp", &fP_dc_y_fp);
    fChain->SetBranchAddress("P.dc.yp_fp", &fP_dc_yp_fp);    

    //HMS
    fChain->SetBranchAddress("H.dc.x_fp", &fH_dc_x_fp);
    fChain->SetBranchAddress("H.dc.xp_fp", &fH_dc_xp_fp);
    fChain->SetBranchAddress("H.dc.y_fp", &fH_dc_y_fp);
    fChain->SetBranchAddress("H.dc.yp_fp", &fH_dc_yp_fp);    

    //------------ Cerenkov Counter ----------    
    fChain->SetBranchAddress("P.hgcer.npeSum", &fP_hgcer_npeSum);
    fChain->SetBranchAddress("H.cer.npeSum", &fH_cer_npeSum);

    //------------ Hodoscope ----------
    //SHMS
    fChain->SetBranchAddress("P.hod.goodstarttime", &fP_hod_goodstarttime);
    fChain->SetBranchAddress("P.hod.starttime", &fP_hod_starttime);
    fChain->SetBranchAddress("P.hod.fpHitsTime", &fP_hod_fpHitsTime);

    //HMS
    fChain->SetBranchAddress("H.hod.goodstarttime", &fH_hod_goodstarttime);  
    fChain->SetBranchAddress("H.hod.starttime", &fH_hod_starttime);    
    fChain->SetBranchAddress("H.hod.fpHitsTime", &fH_hod_fpHitsTime);

    //------------ Calorimeter ----------    
    //SHMS
    fChain->SetBranchAddress("P.cal.eprtracknorm", &fP_cal_eprtracknorm);
    fChain->SetBranchAddress("P.cal.etottracknorm", &fP_cal_etottracknorm);
    
    //HMS
    fChain->SetBranchAddress("H.cal.eprtracknorm", &fH_cal_eprtracknorm);
    fChain->SetBranchAddress("H.cal.etottracknorm", &fH_cal_etottracknorm);

    
    //--------- Trigger ----------------
    //SHMS
    fChain->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw", &fT_coin_pTRIG1_ROC1_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw", &fT_coin_pTRIG4_ROC1_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw", &fT_coin_pTRIG1_ROC2_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw", &fT_coin_pTRIG4_ROC2_tdcTimeRaw);

    //HMS
    fChain->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTimeRaw", &fT_coin_hTRIG1_ROC1_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.hTRIG4_ROC1_tdcTimeRaw", &fT_coin_hTRIG4_ROC1_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.hTRIG1_ROC2_tdcTimeRaw", &fT_coin_hTRIG1_ROC2_tdcTimeRaw);
    fChain->SetBranchAddress("T.coin.hTRIG4_ROC2_tdcTimeRaw", &fT_coin_hTRIG4_ROC2_tdcTimeRaw);

    //fChain->SetBranchAddress("fEvtHdr.fEvtType", &EvtData);  // This branch has some problem at this moment, SetBranchAddress() does not work. It's a nested branch.
   
}
