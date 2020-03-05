/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 *
 * BBC Event Plane Builder
 * \author Ding Chen
 *
 * Updated to Calculate v2
 * \date Dec 23, 2019
 *
 * Tweak to analyze BES-II FXT 3GeV, 7.2GeV and more
 * \author Ding Chen
 * \date Feb 19, 2020
 */

// This is needed for calling standalone classes (not needed on RACF)
//#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <stdio.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
//EPD Utilizer
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"


// Define global constants
const Int_t daynumber     = 6;
const Int_t Ncentralities = 10;
const Int_t order         = 20;
const Int_t twoorder      = 2 * order;
const double d_mean  = 1.02087;//primary
const double d_sigma = 0.00444335;//primary
// BBC azimuthal distribution
Double_t BBC_GetPhi(Int_t e_w,Int_t iTile,Double_t Vx,Double_t Vy) {
    // Get phi of BBC tile
    Double_t Tx = 0.0, Ty = 0.0;
    TRandom *r = new TRandom();
    switch(iTile) {
        case 0: {Ty = 9.64;Tx = 0.0;}
            break;
        case 1: {Ty = 9.64/2.0;Tx = 9.64;}
            break;
        case 2: {Ty = -9.64/2.0;Tx = 9.64;}
            break;
        case 3: {Ty = -9.64;Tx = 0.0;}
            break;
        case 4: {Ty = -9.64/2.0;Tx = -9.64;}
            break;
        case 5: {Ty = 9.64/2.0;Tx = -9.64;}
            break;
        case 6: {Ty = 9.64*3.0/2.0;Tx = (r->Rndm() > 0.5)? 9.64:-9.64;}
            break;
        case 7: {Ty = 9.64*2.0;Tx = 0.0;}
            break;
        case 8: {Ty = 9.64;Tx = 9.64*2.0;}
            break;
        case 9: {Ty = 0.0;Tx = 9.64*2.0;}
            break;
        case 10: {Ty = -9.64;Tx = 9.64*2.0;}
            break;
        case 11: {Ty = -9.64*3.0/2.0;Tx = (r->Rndm() > 0.5)? 9.64:-9.64;}
            break;
        case 12: {Ty = -9.64*2.0;Tx = 0.0;}
            break;
        case 13: {Ty = -9.64;Tx = -9.64*2.0;}
            break;
        case 14: {Ty = 0.0;Tx = -9.64*2.0;}
            break;
        case 15: {Ty = 9.64;Tx = -9.64*2.0;}
            break;
    }
    delete r;
    Double_t bbc_phi = TMath::ATan2(Ty-Vy,Tx-Vx);
    //if(e_w == 1) bbc_phi = bbc_phi - TMath::Pi();
    if(e_w == 0) {
        if (bbc_phi > -0.001) bbc_phi = TMath::Pi() - bbc_phi;
        //if (bbc_phi > 0) bbc_phi = TMath::Pi() - bbc_phi;
        else                  bbc_phi= -TMath::Pi() - bbc_phi;
    }
    if(bbc_phi < 0.0            ) bbc_phi += 2.0*TMath::Pi();
    if(bbc_phi > 2.0*TMath::Pi()) bbc_phi -= 2.0*TMath::Pi();
    return bbc_phi;
}

//////////////////////////////// Main Function /////////////////////////////////
void PicoDstAnalyzer3(const Char_t *inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_BBCEP",
                      Int_t   inputp1 = 1
                    )
{
  // orders of Fourier expansion of momentum space distribution
  Int_t EpOrder = inputp1; // EpOrder = 1, 2, 3

  // Set particle masses
  Double_t Mass_Pion     = 0.13957061;
  Double_t Mass_Kaon     = 0.493677;
  Double_t Mass_Proton   = 0.938272;
  // Set rapidity/pseudorapidity range
  Int_t rapidityBins = 15; Double_t rapidityLow = -2.9/*-3.0*/; Double_t rapidityHigh = 0.1/*0.0*/;
  // Set transverse momentum range
  Int_t ptBins = 15; Double_t ptLow = 0.0; Double_t ptHigh = 3.0;
  // BBC ADC gain factors
  Double_t egain[16];
  Double_t emean[16];
  Double_t esum[16];
  Double_t emean_c[16];
  Double_t esum_c[16];
  Double_t wgain[16];
  Double_t wmean[16];
  Double_t wsum[16];
  Double_t wmean_c[16];
  Double_t wsum_c[16];

  // Initialize BBC ADC gain parameters
  for(Int_t i=0;i<16;i++) {
      egain[i] = 1.0; emean[i] = 0.0; esum[i] = 0.0; emean_c[i] = 0.0; esum_c[i] = 0.0;
      wgain[i] = 1.0; wmean[i] = 0.0; wsum[i] = 0.0; wmean_c[i] = 0.0; wsum_c[i] = 0.0;
  }

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);
  picoReader->SetStatus("StPicoEpdHit",1);
  picoReader->SetStatus("EpdHit",1);

  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !picoReader->chain() ) {
      std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  // EPD EP finder to get EPD event plane
  TString EpdEpOutputName = "EpdEpCorrectionHistograms_OUTPUT_";
  EpdEpOutputName += outFile;
  EpdEpOutputName += ".root";
  StEpdEpFinder *mEpFinder = new StEpdEpFinder(1,EpdEpOutputName,"/star/u/dchen/GitHub/fxtPicoAna/EpdEpCorrectionHistograms_INPUT.root");
  int format = 2;
  mEpFinder->SetEpdHitFormat(format);    // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
  mEpFinder->SetnMipThreshold(0.3);    // recommended by EPD group
  mEpFinder->SetMaxTileWeight(1.0);     // recommended by EPD group
  // mEpFinder->SetEtaWeights(2, TH2D EtaWeight);   // histogram is binned in |eta| and centrality
  // mEpFinder->SetRingWeights(2, double* RingWeights);    // RingWeights is a 1D array of 16 elements.
  TClonesArray * mEpdHits = new TClonesArray("StPicoEpdHit");
  unsigned int found;
  // Retrieve picoDst TChain*
  TChain *mPicoDst = picoReader->chain();
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  std::cout << "EpdHit Branch returned found= " << found << std::endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);

  outFile.Append(".picoDst.result.root");
  // outFile.Prepend(Form("%d_%d_%d_",cutID,variationID,versionID));
  TFile *outputFile = new TFile(outFile,"recreate");

  int ievtcut[5] = {0};
  int itrkcut[6] = {0};

  // Histograms

  // QA Plots

  TH1D *  h_phi_y = new TH1D("h_phi_y","Rapidity Distribution of K^{+}K^{-} Pairs",60,-3.0,3.0);
  h_phi_y->GetXaxis()->SetTitle("y");
  h_phi_y->GetYaxis()->SetTitle("# of K^{+}K^{-} pairs");

  TH1D *  hist_trackmult_proton = new TH1D("hist_trackmult_proton","hist_trackmult_proton",100,-0.5,99.5);
  hist_trackmult_proton->GetXaxis()->SetTitle("# of Protons");
  hist_trackmult_proton->GetYaxis()->SetTitle("# of Events");

  TH1D *  hist_trackmult_pionPlus = new TH1D("hist_trackmult_pionPlus","hist_trackmult_pionPlus",100,-0.5,99.5);
  hist_trackmult_pionPlus->GetXaxis()->SetTitle("# of #pi^{#plus}");
  hist_trackmult_pionPlus->GetYaxis()->SetTitle("# of Events");

  TH1D *  hist_trackmult_pionMinus = new TH1D("hist_trackmult_pionMinus","hist_trackmult_pionMinus",100,-0.5,99.5);
  hist_trackmult_pionMinus->GetXaxis()->SetTitle("# of #pi^{#minus}");
  hist_trackmult_pionMinus->GetYaxis()->SetTitle("# of Events");

  TH1D *  hist_trackmult_kaonPlus = new TH1D("hist_trackmult_kaonPlus","hist_trackmult_kaonPlus",100,-0.5,99.5);
  hist_trackmult_kaonPlus->GetXaxis()->SetTitle("# of K^{#plus}");
  hist_trackmult_kaonPlus->GetYaxis()->SetTitle("# of Events");

  TH1D *  hist_trackmult_kaonMinus = new TH1D("hist_trackmult_kaonMinus","hist_trackmult_kaonMinus",100,-0.5,99.5);
  hist_trackmult_kaonMinus->GetXaxis()->SetTitle("# of K^{#minus}");
  hist_trackmult_kaonMinus->GetYaxis()->SetTitle("# of Events");



  TH1D *  hist_ratio = new TH1D("hist_ratio","hist_ratio",100,0,2);
  hist_ratio->GetXaxis()->SetTitle("nHitsFit/nHitsPoss");
  hist_ratio->GetYaxis()->SetTitle("# of Tracks");

  TH1D *  hist_nHits = new TH1D("hist_nHits","hist_nHits",100,-0.5,99.5);
  hist_nHits->GetXaxis()->SetTitle("nHits");
  hist_nHits->GetYaxis()->SetTitle("# of Tracks");

  TH1D *  hist_ndEdx = new TH1D("hist_ndEdx","hist_ndEdx",100,-0.5,99.5);
  hist_ndEdx->GetXaxis()->SetTitle("nDedx");
  hist_ndEdx->GetYaxis()->SetTitle("# of Tracks");

  TH1D *  hist_DCA = new TH1D("hist_DCA","hist_DCA",100,0,10.0);
  hist_DCA->GetXaxis()->SetTitle("DCA [cm]");
  hist_DCA->GetYaxis()->SetTitle("# of Tracks");

  TH1D *  h_evt_vs_cut = new TH1D("h_evt_vs_cut","h_evt_vs_cut",11,-0.5,10.5);
  h_evt_vs_cut->GetXaxis()->SetTitle("Event-Level cuts");
  h_evt_vs_cut->GetYaxis()->SetTitle("# of Events");

  TH1D *  h_trk_vs_cut = new TH1D("h_trk_vs_cut","h_trk_vs_cut",11,-0.5,10.5);
  h_trk_vs_cut->GetXaxis()->SetTitle("Track level cuts");
  h_trk_vs_cut->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_runId = new TH1D("hist_runId","Event runId",4001,-0.5,4000.5);
  hist_runId->GetXaxis()->SetTitle("RunId");
  hist_runId->GetYaxis()->SetTitle("# of events");

  TH1D *hist_triggerID = new TH1D("hist_triggerID","Event TriggerId",801,-0.5,800.5);
  hist_triggerID->GetXaxis()->SetTitle("TriggerID");
  hist_triggerID->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Vz = new TH1D("hist_Vz","V_{Z} [cm]",6000,-300.0,300.0);
  hist_Vz->GetXaxis()->SetTitle("V_{Z} [cm]");
  hist_Vz->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Vz_cut = new TH1D("hist_Vz_cut","V_{Z} after cut [cm]",6000,-300.0,300.0);
  hist_Vz_cut->GetXaxis()->SetTitle("V_{Z} [cm]");
  hist_Vz_cut->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Vr = new TH1D("hist_Vr","V_{R} [cm]",500,0.0,20.0);
  hist_Vr->GetXaxis()->SetTitle("V_{R} [cm]");
  hist_Vr->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Vr_cut = new TH1D("hist_Vr_cut","V_{R} after cut [cm]",500,0.0,20.0);
  hist_Vr_cut->GetXaxis()->SetTitle("V_{R} [cm]");
  hist_Vr_cut->GetYaxis()->SetTitle("# of events");

  TH2D *hist_VyVx = new TH2D("hist_VyVx","V_{Y} [cm] vs. V_{X} [cm]",500,-5.0,5.0,500,-5.0,5.0);
  hist_VyVx->GetXaxis()->SetTitle("V_{X} [cm]");
  hist_VyVx->GetYaxis()->SetTitle("V_{Y} [cm]");

  TH2D *hist_VyVx_cut = new TH2D("hist_VyVx_cut","V_{Y} [cm] vs. V_{X} after cut [cm]",500,-5.0,5.0,500,-5.0,5.0);
  hist_VyVx_cut->GetXaxis()->SetTitle("V_{X} [cm]");
  hist_VyVx_cut->GetYaxis()->SetTitle("V_{Y} [cm]");

  TH1D *hist_realTrackMult = new TH1D("hist_realTrackMult","Actual track multiplicity",1001,-0.5,1000.5);
  hist_realTrackMult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult->GetXaxis()->SetTitle("# of events");

  TH2D *hist_realTrackMult_refmult = new TH2D("hist_realTrackMult_refmult","Actual track multiplicity vs. RefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  hist_realTrackMult_refmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_refmult->GetXaxis()->SetTitle("RefMult");

  TH2D *hist_realTrackMult_grefmult = new TH2D("hist_realTrackMult_grefmult","Actual track multiplicity vs. gRefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  hist_realTrackMult_grefmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_grefmult->GetXaxis()->SetTitle("gRefMult");

  TH2D *hist_realTrackMult_tofmult = new TH2D("hist_realTrackMult_tofmult","Actual track multiplicity vs. TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  hist_realTrackMult_tofmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_tofmult->GetXaxis()->SetTitle("TofMult");

  TH2D *hist_realTrackMult_trackmult = new TH2D("hist_realTrackMult_trackmult","Actual track multiplicity vs. TrackMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  hist_realTrackMult_trackmult->GetXaxis()->SetTitle("Actual TrackMult");
  hist_realTrackMult_trackmult->GetXaxis()->SetTitle("TrackMult");

  TH1D *hist_cent = new TH1D("hist_cent","Centrality",Ncentralities+1,-0.5,Ncentralities+0.5);
  hist_cent->GetXaxis()->SetTitle("Centrality bin");
  hist_cent->GetYaxis()->SetTitle("# of events");

  //Track histograms
  TH1D *hist_pt = new TH1D("hist_pt","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta = new TH1D("hist_eta","#eta",200,-3.0,0.5);
  hist_eta->GetXaxis()->SetTitle("#eta");
  hist_eta->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_phi = new TH1D("hist_phi","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx = new TH2D("hist_dEdx","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta = new TH2D("hist_beta","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass = new TH2D("hist_mass","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH1D *hist_pt_proton = new TH1D("hist_pt_proton","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt_proton->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_proton->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta_proton = new TH1D("hist_eta_proton","#eta",200,-3.0,0.5);
  hist_eta_proton->GetXaxis()->SetTitle("#eta");
  hist_eta_proton->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_y_proton = new TH1D("hist_y_proton","Rapidity y",200,-3.0,0.5);
  hist_y_proton->GetXaxis()->SetTitle("Rapidity y");
  hist_y_proton->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_rap_eta_proton = new TH2D("hist_rap_eta_proton","proton y versus #eta",250,-2.5,0,250,-2.5,0);
  hist_rap_eta_proton->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_proton->GetYaxis()->SetTitle("Rapidity y");

  TH2D *hist_pt_y_proton = new TH2D("hist_pt_y_proton","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_proton->GetXaxis()->SetTitle("y");
  hist_pt_y_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_eta_proton = new TH2D("hist_pt_eta_proton","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_eta_proton->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D *hist_phi_proton = new TH1D("hist_phi_proton","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi_proton->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_proton->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx_proton = new TH2D("hist_dEdx_proton","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_proton->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta_proton = new TH2D("hist_beta_proton","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_proton->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass_proton = new TH2D("hist_mass_proton","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_proton->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH1D *hist_pt_pionPlus = new TH1D("hist_pt_pionPlus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt_pionPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_pionPlus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta_pionPlus = new TH1D("hist_eta_pionPlus","#eta",200,-3.0,0.5);
  hist_eta_pionPlus->GetXaxis()->SetTitle("#eta");
  hist_eta_pionPlus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_y_pionPlus = new TH1D("hist_y_pionPlus","y",200,-3.0,0.5);
  hist_y_pionPlus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_pionPlus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_rap_eta_pionPlus = new TH2D("hist_rap_eta_pionPlus","pionPlus y versus #eta",250,-2.5,0,250,-2.5,0);
  hist_rap_eta_pionPlus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_pionPlus->GetYaxis()->SetTitle("Rapidity y");

  TH2D *hist_pt_y_pionPlus = new TH2D("hist_pt_y_pionPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_pionPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_pionPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_eta_pionPlus = new TH2D("hist_pt_eta_pionPlus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_eta_pionPlus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_pionPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D *hist_phi_pionPlus = new TH1D("hist_phi_pionPlus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi_pionPlus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_pionPlus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx_pionPlus = new TH2D("hist_dEdx_pionPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_pionPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta_pionPlus = new TH2D("hist_beta_pionPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_pionPlus->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass_pionPlus = new TH2D("hist_mass_pionPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_pionPlus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH1D *hist_pt_pionMinus = new TH1D("hist_pt_pionMinus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt_pionMinus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_pionMinus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta_pionMinus = new TH1D("hist_eta_pionMinus","#eta",200,-3.0,0.5);
  hist_eta_pionMinus->GetXaxis()->SetTitle("#eta");
  hist_eta_pionMinus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_y_pionMinus = new TH1D("hist_y_pionMinus","y",200,-3.0,0.5);
  hist_y_pionMinus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_pionMinus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_rap_eta_pionMinus = new TH2D("hist_rap_eta_pionMinus","pionMinus y versus #eta",250,-2.5,0,250,-2.5,0);
  hist_rap_eta_pionMinus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_pionMinus->GetYaxis()->SetTitle("Rapidity y");

  TH2D *hist_pt_y_pionMinus = new TH2D("hist_pt_y_pionMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_pionMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_pionMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_eta_pionMinus = new TH2D("hist_pt_eta_pionMinus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_eta_pionMinus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_pionMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D *hist_phi_pionMinus = new TH1D("hist_phi_pionMinus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi_pionMinus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_pionMinus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx_pionMinus = new TH2D("hist_dEdx_pionMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_pionMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta_pionMinus = new TH2D("hist_beta_pionMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_pionMinus->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass_pionMinus = new TH2D("hist_mass_pionMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_pionMinus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH1D *hist_pt_kaonPlus = new TH1D("hist_pt_kaonPlus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt_kaonPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_kaonPlus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta_kaonPlus = new TH1D("hist_eta_kaonPlus","#eta",200,-3.0,0.5);
  hist_eta_kaonPlus->GetXaxis()->SetTitle("#eta");
  hist_eta_kaonPlus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_y_kaonPlus = new TH1D("hist_y_kaonPlus","y",200,-3.0,0.5);
  hist_y_kaonPlus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_kaonPlus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_rap_eta_kaonPlus = new TH2D("hist_rap_eta_kaonPlus","kaonPlus y versus #eta",250,-2.5,0,250,-2.5,0);
  hist_rap_eta_kaonPlus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_kaonPlus->GetYaxis()->SetTitle("Rapidity y");

  TH2D *hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_kaonPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_eta_kaonPlus = new TH2D("hist_pt_eta_kaonPlus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_eta_kaonPlus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D *hist_phi_kaonPlus = new TH1D("hist_phi_kaonPlus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi_kaonPlus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_kaonPlus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx_kaonPlus = new TH2D("hist_dEdx_kaonPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_kaonPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta_kaonPlus = new TH2D("hist_beta_kaonPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_kaonPlus->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass_kaonPlus = new TH2D("hist_mass_kaonPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_kaonPlus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH1D *hist_pt_kaonMinus = new TH1D("hist_pt_kaonMinus","p_{T} [GeV/c]",1000,0.0,5.0);
  hist_pt_kaonMinus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_kaonMinus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_eta_kaonMinus = new TH1D("hist_eta_kaonMinus","#eta",200,-3.0,0.5);
  hist_eta_kaonMinus->GetXaxis()->SetTitle("#eta");
  hist_eta_kaonMinus->GetYaxis()->SetTitle("# of tracks");

  TH1D *hist_y_kaonMinus = new TH1D("hist_y_kaonMinus","y",200,-3.0,0.5);
  hist_y_kaonMinus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_kaonMinus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_rap_eta_kaonMinus = new TH2D("hist_rap_eta_kaonMinus","kaonMinus y versus #eta",250,-2.5,0,250,-2.5,0);
  hist_rap_eta_kaonMinus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_kaonMinus->GetYaxis()->SetTitle("Rapidity y");

  TH2D *hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_kaonMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_eta_kaonMinus = new TH2D("hist_pt_eta_kaonMinus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_eta_kaonMinus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D *hist_phi_kaonMinus = new TH1D("hist_phi_kaonMinus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_phi_kaonMinus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_kaonMinus->GetYaxis()->SetTitle("# of tracks");

  TH2D *hist_dEdx_kaonMinus = new TH2D("hist_dEdx_kaonMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_kaonMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *hist_beta_kaonMinus = new TH2D("hist_beta_kaonMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  hist_beta_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_kaonMinus->GetYaxis()->SetTitle("1/#beta");

  TH2D *hist_mass_kaonMinus = new TH2D("hist_mass_kaonMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  hist_mass_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_kaonMinus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");

  TH2D *hist_pt_y_Phi = new TH2D("hist_pt_y_Phi","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_Phi->GetXaxis()->SetTitle("y");
  hist_pt_y_Phi->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D * h_dip_angle = new TH1D("h_dip_angle","h_dip_angle",1000,-1,1.0);

  //////////////////////////////////////////////////////////////////////////////

  TH2D *GPC_request_m2_vs_dEdx_Pos = new TH2D("GPC_request_m2_vs_dEdx_Pos"
                                                ,"Positively charged track m^2 (GeV/c^2)^2 vs dE/dx (keV/cm)"
                                                ,1000,0.0,10.0
                                                ,1000,-0.6,4.0
                                                );
  GPC_request_m2_vs_dEdx_Pos->GetXaxis()->SetTitle("TPC dE/dx (keV/cm)");
  GPC_request_m2_vs_dEdx_Pos->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_dEdx_Neg = new TH2D("GPC_request_m2_vs_dEdx_Neg"
                                              ,"Negatively charged track m^2 (GeV/c^2)^2 vs dE/dx (keV/cm)"
                                              ,1000,0.0,10.0
                                              ,1000,-0.6,4.0
                                              );
  GPC_request_m2_vs_dEdx_Neg->GetXaxis()->SetTitle("TPC dE/dx (keV/cm)");
  GPC_request_m2_vs_dEdx_Neg->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_nSigmaProton = new TH2D("GPC_request_m2_vs_nSigmaProton"
                                                  ,"Positively charged track m^2 (GeV/c^2)^2 vs nSigmaProton"
                                                  ,1000,-50.0,50.0
                                                  ,1000,-0.6,4.0
                                                  );
  GPC_request_m2_vs_nSigmaProton->GetXaxis()->SetTitle("TPC dE/dx n#sigma_{p}");
  GPC_request_m2_vs_nSigmaProton->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_nSigmaPionPlus = new TH2D("GPC_request_m2_vs_nSigmaPionPlus"
                                                  ,"Positively charged track m^2 (GeV/c^2)^2 vs nSigmaPion"
                                                  ,1000,-50.0,50.0
                                                  ,1000,-0.6,4.0
                                                  );
  GPC_request_m2_vs_nSigmaPionPlus->GetXaxis()->SetTitle("TPC dE/dx n#sigma_{p}");
  GPC_request_m2_vs_nSigmaPionPlus->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_nSigmaPionMinus = new TH2D("GPC_request_m2_vs_nSigmaPionMinus"
                                                    ,"Negatively charged track m^2 (GeV/c^2)^2 vs nSigmaPion"
                                                    ,1000,-50.0,50.0
                                                    ,1000,-0.6,4.0
                                                    );
  GPC_request_m2_vs_nSigmaPionMinus->GetXaxis()->SetTitle("TPC dE/dx n#sigma_{p}");
  GPC_request_m2_vs_nSigmaPionMinus->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_nSigmaKaonPlus = new TH2D("GPC_request_m2_vs_nSigmaKaonPlus"
                                                    ,"Positively charged track m^2 (GeV/c^2)^2 vs nSigmaKaon"
                                                    ,1000,-50.0,50.0
                                                    ,1000,-0.6,4.0
                                                    );
  GPC_request_m2_vs_nSigmaKaonPlus->GetXaxis()->SetTitle("TPC dE/dx n#sigma_{p}");
  GPC_request_m2_vs_nSigmaKaonPlus->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  TH2D *GPC_request_m2_vs_nSigmaKaonMinus = new TH2D("GPC_request_m2_vs_nSigmaKaonMinus"
                                                     ,"Negatively charged track m^2 (GeV/c^2)^2 vs nSigmaKaon"
                                                     ,1000,-50.0,50.0
                                                     ,1000,-0.6,4.0
                                                     );
  GPC_request_m2_vs_nSigmaKaonMinus->GetXaxis()->SetTitle("TPC dE/dx n#sigma_{p}");
  GPC_request_m2_vs_nSigmaKaonMinus->GetYaxis()->SetTitle("TOF reconstructed m^{2} (GeV/c^2)^2");

  //////////////////////////////////////////////////////////////////////////////

  // END QA plots

  // Flow Histograms
  TProfile3D *profile3D_proton_v2 = new TProfile3D("profile3D_proton_v2","Proton v_{2}",Ncentralities,0.5,Ncentralities+0.5,ptBins,ptLow,ptHigh,rapidityBins,rapidityLow,rapidityHigh,"");
  profile3D_proton_v2->BuildOptions(-1,1,"");
  profile3D_proton_v2->GetXaxis()->SetTitle("Centrality bin");
  profile3D_proton_v2->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  profile3D_proton_v2->GetZaxis()->SetTitle("y");
  profile3D_proton_v2->Sumw2();

  TProfile3D *profile3D_pionPlus_v2 = new TProfile3D("profile3D_pionPlus_v2","#pi^{#plus} v_{2}",Ncentralities,0.5,Ncentralities+0.5,ptBins,ptLow,ptHigh,rapidityBins,rapidityLow,rapidityHigh,"");
  profile3D_pionPlus_v2->BuildOptions(-1,1,"");
  profile3D_pionPlus_v2->GetXaxis()->SetTitle("Centrality bin");
  profile3D_pionPlus_v2->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  profile3D_pionPlus_v2->GetZaxis()->SetTitle("y");
  profile3D_pionPlus_v2->Sumw2();

  TProfile3D *profile3D_pionMinus_v2 = new TProfile3D("profile3D_pionMinus_v2","#pi^{#minus} v_{2}",Ncentralities,0.5,Ncentralities+0.5,ptBins,ptLow,ptHigh,rapidityBins,rapidityLow,rapidityHigh,"");
  profile3D_pionMinus_v2->BuildOptions(-1,1,"");
  profile3D_pionMinus_v2->GetXaxis()->SetTitle("Centrality bin");
  profile3D_pionMinus_v2->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  profile3D_pionMinus_v2->GetZaxis()->SetTitle("y");
  profile3D_pionMinus_v2->Sumw2();

  TProfile3D *profile3D_kaonPlus_v2 = new TProfile3D("profile3D_kaonPlus_v2","K^{#plus} v_{2}",Ncentralities,0.5,Ncentralities+0.5,ptBins,ptLow,ptHigh,rapidityBins,rapidityLow,rapidityHigh,"");
  profile3D_kaonPlus_v2->BuildOptions(-1,1,"");
  profile3D_kaonPlus_v2->GetXaxis()->SetTitle("Centrality bin");
  profile3D_kaonPlus_v2->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  profile3D_kaonPlus_v2->GetZaxis()->SetTitle("y");
  profile3D_kaonPlus_v2->Sumw2();

  TProfile3D *profile3D_kaonMinus_v2 = new TProfile3D("profile3D_kaonMinus_v2","K^{#minus} v_{2}",Ncentralities,0.5,Ncentralities+0.5,ptBins,ptLow,ptHigh,rapidityBins,rapidityLow,rapidityHigh,"");
  profile3D_kaonMinus_v2->BuildOptions(-1,1,"");
  profile3D_kaonMinus_v2->GetXaxis()->SetTitle("Centrality bin");
  profile3D_kaonMinus_v2->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  profile3D_kaonMinus_v2->GetZaxis()->SetTitle("y");
  profile3D_kaonMinus_v2->Sumw2();

  // EP flow Plots
  TH2D *h2_proton_v2 = new TH2D("h2_proton_v2","Proton v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_proton_v2->GetYaxis()->SetTitle("v_{2}");
  h2_proton_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *h2_pionPlus_v2 = new TH2D("h2_pionPlus_v2","#pi^{#plus} v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_pionPlus_v2->GetYaxis()->SetTitle("v_{2}");
  h2_pionPlus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *h2_pionMinus_v2 = new TH2D("h2_pionMinus_v2","#pi^{#minus} v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_pionMinus_v2->GetYaxis()->SetTitle("v_{2}");
  h2_pionMinus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *h2_kaonPlus_v2 = new TH2D("h2_kaonPlus_v2","K^{#plus} v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_kaonPlus_v2->GetYaxis()->SetTitle("v_{2}");
  h2_kaonPlus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *h2_kaonMinus_v2 = new TH2D("h2_kaonMinus_v2","K^{#minus} v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_kaonMinus_v2->GetYaxis()->SetTitle("v_{2}");
  h2_kaonMinus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *h2_pions_v2 = new TH2D("h2_pions_v2","#pi v_{2}",ptBins,ptLow,ptHigh,1000,-1.0,1.0);
  h2_pions_v2->GetYaxis()->SetTitle("v_{2}");
  h2_pions_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  TProfile *profile_proton_v2 = new TProfile("profile_proton_v2","Proton v_{2}",ptBins,ptLow,ptHigh,"");
  profile_proton_v2->BuildOptions(-1,1,"");
  profile_proton_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_proton_v2->Sumw2();

  TProfile *profile_pionPlus_v2 = new TProfile("profile_pionPlus_v2","#pi^{#plus} v_{2}",ptBins,ptLow,ptHigh,"");
  profile_pionPlus_v2->BuildOptions(-1,1,"");
  profile_pionPlus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_pionPlus_v2->Sumw2();

  TProfile *profile_pionMinus_v2 = new TProfile("profile_pionMinus_v2","#pi^{#minus} v_{2}",ptBins,ptLow,ptHigh,"");
  profile_pionMinus_v2->BuildOptions(-1,1,"");
  profile_pionMinus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_pionMinus_v2->Sumw2();

  TProfile *profile_kaonPlus_v2 = new TProfile("profile_kaonPlus_v2","K^{#plus} v_{2}",ptBins,ptLow,ptHigh,"");
  profile_kaonPlus_v2->BuildOptions(-1,1,"");
  profile_kaonPlus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_kaonPlus_v2->Sumw2();

  TProfile *profile_kaonMinus_v2 = new TProfile("profile_kaonMinus_v2","K^{#minus} v_{2}",ptBins,ptLow,ptHigh,"");
  profile_kaonMinus_v2->BuildOptions(-1,1,"");
  profile_kaonMinus_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_kaonMinus_v2->Sumw2();


  TProfile *profile_pions_v2 = new TProfile("profile_pions_v2","#pi^{#plus} v_{2}",ptBins,ptLow,ptHigh,"");
  profile_pions_v2->BuildOptions(-1,1,"");
  profile_pions_v2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  profile_pions_v2->Sumw2();

  // Phi meson plots
  TH1D *  h_K_DCA_r      = new TH1D("h_K_DCA_r","h_K_DCA_r",400,0.0,4.0);
  TH1D *  h_K_obj_DCA_r  = new TH1D("h_K_obj_DCA_r","h_K_obj_DCA_r",400,0.0,4.0);
  TH1D *  h_K_diff_DCA_r = new TH1D("h_K_diff_DCA_r","h_K_diff_DCA_r",800,-4.0,4.0);

  TH1D *  h_prim_inv_m_PHI    = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",1000,0.9,1.1);

  TH2D  *h2_phi_v2_vs_invM = new TH2D("h2_phi_v2_vs_invM","h2_phi_v2_vs_invM",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_bin2 = new TH2D("h2_phi_v2_vs_invM_bin2","h2_phi_v2_vs_invM_bin2, -1.5 <= y < -1.0",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_bin2->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_bin2->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_bin3 = new TH2D("h2_phi_v2_vs_invM_bin3","h2_phi_v2_vs_invM_bin3, -1.0 <= y < -0.5",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_bin3->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_bin3->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_bin4 = new TH2D("h2_phi_v2_vs_invM_bin4","h2_phi_v2_vs_invM_bin4, -0.5 <= y <= 0.0",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_bin4->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_bin4->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin1 = new TH2D("h2_phi_v2_vs_invM_pTbin1","h2_phi_v2_vs_invM_pTbin1, 0 <= p_{T} < 0.5 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin1->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin1->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin2 = new TH2D("h2_phi_v2_vs_invM_pTbin2","h2_phi_v2_vs_invM_pTbin2, 0.5 <= p_{T} < 1.0 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin2->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin2->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin3 = new TH2D("h2_phi_v2_vs_invM_pTbin3","h2_phi_v2_vs_invM_pTbin3, 1.0 <= p_{T} < 1.5 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin3->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin3->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin4 = new TH2D("h2_phi_v2_vs_invM_pTbin4","h2_phi_v2_vs_invM_pTbin4, 1.5 <= p_{T} < 2.0 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin4->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin4->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin5 = new TH2D("h2_phi_v2_vs_invM_pTbin5","h2_phi_v2_vs_invM_pTbin5, 2.0 <= p_{T} < 2.5 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin5->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin5->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin6 = new TH2D("h2_phi_v2_vs_invM_pTbin6","h2_phi_v2_vs_invM_pTbin6, 2.5 <= p_{T} < 3.0 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin6->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin6->GetYaxis()->SetTitle("v_{2}");

  TH2D  *h2_phi_v2_vs_invM_pTbin7 = new TH2D("h2_phi_v2_vs_invM_pTbin7","h2_phi_v2_vs_invM_pTbin7, 3.0 <= p_{T} < 3.5 (GeV/c)",100,0.9,1.1,1000,-1.0,1.0);
  h2_phi_v2_vs_invM_pTbin7->GetXaxis()->SetTitle("K^{+}K^{-} Invariant Mass(GeV/c^{2})");
  h2_phi_v2_vs_invM_pTbin7->GetYaxis()->SetTitle("v_{2}");

  TProfile  *TP_phi_v2_vs_invM      = nullptr;
  TProfile  *TP_phi_v2_vs_invM_bin2 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_bin3 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_bin4 = nullptr;

  TProfile  *TP_phi_v2_vs_invM_pTbin1 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin2 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin3 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin4 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin5 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin6 = nullptr;
  TProfile  *TP_phi_v2_vs_invM_pTbin7 = nullptr;

  // END Phi meson plots

  // END Flow Histograms
  Char_t name[100], description[200];

  // EPD EPs
  TH2D *hist2_Epd_east_Qy_Qx_raw = new TH2D("hist2_Epd_east_Qy_Qx_raw","EPD east Qy vs Qx",600,-3.0,3.0,600,-3.0,3.0);
  hist2_Epd_east_Qy_Qx_raw->GetXaxis()->SetTitle("Q_x^{EPD east}_{2} ");
  hist2_Epd_east_Qy_Qx_raw->GetYaxis()->SetTitle("Q_y^{EPD east}_{2} ");

  TH2D *hist2_Epd_east_Qy_Qx_Weighted = new TH2D("hist2_Epd_east_Qy_Qx_Weighted","EPD east Qy vs Qx (Weighted)",600,-3.0,3.0,600,-3.0,3.0);
  hist2_Epd_east_Qy_Qx_Weighted->GetXaxis()->SetTitle("Q_x^{EPD east}_{2} ");
  hist2_Epd_east_Qy_Qx_Weighted->GetYaxis()->SetTitle("Q_y^{EPD east}_{2} ");

  TH2D *hist2_Epd_west_Qy_Qx_raw = new TH2D("hist2_Epd_west_Qy_Qx_raw","EPD west Qy vs Qx",600,-3.0,3.0,600,-3.0,3.0);
  hist2_Epd_west_Qy_Qx_raw->GetXaxis()->SetTitle("Q_x^{EPD west}_{2} ");
  hist2_Epd_west_Qy_Qx_raw->GetYaxis()->SetTitle("Q_y^{EPD west}_{2} ");

  TH2D *hist2_Epd_west_Qy_Qx_Weighted = new TH2D("hist2_Epd_west_Qy_Qx_Weighted","EPD west Qy vs Qx (Weighted)",600,-3.0,3.0,600,-3.0,3.0);
  hist2_Epd_west_Qy_Qx_Weighted->GetXaxis()->SetTitle("Q_x^{EPD west}_{2} ");
  hist2_Epd_west_Qy_Qx_Weighted->GetYaxis()->SetTitle("Q_y^{EPD west}_{2} ");

  TH1D *hist_Epd_east_psi_raw = new TH1D("hist_Epd_east_psi_raw","EPD east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_east_psi_raw->GetXaxis()->SetTitle("#psi^{EPD east}_{2} [Radian]");
  hist_Epd_east_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Epd_east_psi_Weighted = new TH1D("hist_Epd_east_psi_Weighted","EPD east EP (Weighted)",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_east_psi_Weighted->GetXaxis()->SetTitle("#psi^{EPD east}_{2} [Radian]");
  hist_Epd_east_psi_Weighted->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Epd_east_psi_Shifted = new TH1D("hist_Epd_east_psi_Shifted","EPD east EP (Weighted & Shifted)",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_east_psi_Shifted->GetXaxis()->SetTitle("#psi^{EPD east}_{2} [Radian]");
  hist_Epd_east_psi_Shifted->GetYaxis()->SetTitle("# of events");


  TH1D *hist_Epd_west_psi_raw = new TH1D("hist_Epd_west_psi_raw","EPD west EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_west_psi_raw->GetXaxis()->SetTitle("#psi^{EPD west}_{2} [Radian]");
  hist_Epd_west_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_Epd_west_psi_Weighted = new TH1D("hist_Epd_west_psi_Weighted","EPD west EP (Weighted)",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_west_psi_Weighted->GetXaxis()->SetTitle("#psi^{EPD west}_{2} [Radian]");
  hist_Epd_west_psi_Weighted->GetYaxis()->SetTitle("# of events");


  TH1D *hist_Epd_west_psi_Shifted = new TH1D("hist_Epd_west_psi_Shifted","EPD west EP (Weighted & Shifted)",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_Epd_west_psi_Shifted->GetXaxis()->SetTitle("#psi^{EPD west}_{2} [Radian]");
  hist_Epd_west_psi_Shifted->GetYaxis()->SetTitle("# of events");

  // TPC EPs
  TH1D *hist_tpc_east_psi_raw = new TH1D("hist_tpc_east_psi_raw","TPC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_east_psi_raw->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  hist_tpc_east_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_tpc_east_psi_recentered = new TH1D("hist_tpc_east_psi_recentered","TPC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_east_psi_recentered->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  hist_tpc_east_psi_recentered->GetYaxis()->SetTitle("# of events");

  TH1D *hist_tpc_east_psi_flattened = new TH1D("hist_tpc_east_psi_flattened","TPC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_east_psi_flattened->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  hist_tpc_east_psi_flattened->GetYaxis()->SetTitle("# of events");

  TH1D *hist_tpc_west_psi_raw = new TH1D("hist_tpc_west_psi_raw","TPC west EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_west_psi_raw->GetXaxis()->SetTitle("#psi^{TPC west}_{2} [Radian]");
  hist_tpc_west_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_tpc_west_psi_recentered = new TH1D("hist_tpc_west_psi_recentered","TPC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_west_psi_recentered->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  hist_tpc_west_psi_recentered->GetYaxis()->SetTitle("# of events");

  TH1D *hist_tpc_west_psi_flattened = new TH1D("hist_tpc_west_psi_flattened","TPC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_tpc_west_psi_flattened->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  hist_tpc_west_psi_flattened->GetYaxis()->SetTitle("# of events");

  TProfile3D *profile3D_tpc_east_Qx_Qy = new TProfile3D("profile3D_tpc_east_Qx_Qy","profile3D_tpc_east_Qx_Qy",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,2,0.5,2.5,"");
  profile3D_tpc_east_Qx_Qy->BuildOptions(0.0,0.0,"");

  TProfile3D *profile3D_tpc_east_psiShift = new TProfile3D("profile3D_tpc_east_psiShift","profile3D_tpc_east_psiShift",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,twoorder,0.5,twoorder+0.5,"");
  profile3D_tpc_east_psiShift->BuildOptions(-1.0,1.0,"");

  TProfile3D *profile3D_tpc_west_Qx_Qy = new TProfile3D("profile3D_tpc_west_Qx_Qy","profile3D_tpc_west_Qx_Qy",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,2,0.5,2.5,"");
  profile3D_tpc_west_Qx_Qy->BuildOptions(0.0,0.0,"");

  TProfile3D *profile3D_tpc_west_psiShift = new TProfile3D("profile3D_tpc_west_psiShift","profile3D_tpc_west_psiShift",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,twoorder,0.5,twoorder+0.5,"");
  profile3D_tpc_west_psiShift->BuildOptions(-1.0,1.0,"");

  // END TPC EPs

  // BBC EPs
  TH1D *bbc_east_adc_profile = new TH1D("bbc_east_adc_profile","bbc_east_adc_profile",16,0.5,16.5);
  TH1D *bbc_east_gain_corrected_adc_profile = new TH1D("bbc_east_gain_corrected_adc_profile","bbc_east_gain_corrected_adc_profile",16,0.5,16.5);

  TH1D *hist_bbc_east_psi_raw = new TH1D("hist_bbc_east_psi_raw","BBC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_east_psi_raw->GetXaxis()->SetTitle("#psi^{BBC east}_{2} [Radian]");
  hist_bbc_east_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_bbc_east_psi_recentered = new TH1D("hist_bbc_east_psi_recentered","BBC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_east_psi_recentered->GetXaxis()->SetTitle("#psi^{BBC east}_{2} [Radian]");
  hist_bbc_east_psi_recentered->GetYaxis()->SetTitle("# of events");

  TH1D *hist_bbc_east_psi_flattened = new TH1D("hist_bbc_east_psi_flattened","BBC east EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_east_psi_flattened->GetXaxis()->SetTitle("#psi^{BBC east}_{2} [Radian]");
  hist_bbc_east_psi_flattened->GetYaxis()->SetTitle("# of events");

  TProfile3D *profile3D_bbc_east_Qx_Qy = new TProfile3D("profile3D_bbc_east_Qx_Qy","profile3D_bbc_east_Qx_Qy",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,2,0.5,2.5,"");
  profile3D_bbc_east_Qx_Qy->BuildOptions(0.0,0.0,"");

  TProfile3D *profile3D_bbc_east_psiShift = new TProfile3D("profile3D_bbc_east_psiShift","profile3D_bbc_east_psiShift",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,twoorder,0.5,twoorder+0.5,"");
  profile3D_bbc_east_psiShift->BuildOptions(-1.0,1.0,"");

  TH1D *bbc_west_adc_profile = new TH1D("bbc_west_adc_profile","bbc_west_adc_profile",16,0.5,16.5);
  TH1D *bbc_west_gain_corrected_adc_profile = new TH1D("bbc_west_gain_corrected_adc_profile","bbc_west_gain_corrected_adc_profile",16,0.5,16.5);

  TH1D *hist_bbc_west_psi_raw = new TH1D("hist_bbc_west_psi_raw","BBC west EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_west_psi_raw->GetXaxis()->SetTitle("#psi^{BBC west}_{2} [Radian]");
  hist_bbc_west_psi_raw->GetYaxis()->SetTitle("# of events");

  TH1D *hist_bbc_west_psi_recentered = new TH1D("hist_bbc_west_psi_recentered","BBC west EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_west_psi_recentered->GetXaxis()->SetTitle("#psi^{BBC west}_{2} [Radian]");
  hist_bbc_west_psi_recentered->GetYaxis()->SetTitle("# of events");

  TH1D *hist_bbc_west_psi_flattened = new TH1D("hist_bbc_west_psi_flattened","BBC west EP",500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  hist_bbc_west_psi_flattened->GetXaxis()->SetTitle("#psi^{BBC west}_{2} [Radian]");
  hist_bbc_west_psi_flattened->GetYaxis()->SetTitle("# of events");

  TProfile3D *profile3D_bbc_west_Qx_Qy = new TProfile3D("profile3D_bbc_west_Qx_Qy","profile3D_bbc_west_Qx_Qy",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,2,0.5,2.5,"");
  profile3D_bbc_west_Qx_Qy->BuildOptions(0.0,0.0,"");

  TProfile3D *profile3D_bbc_west_psiShift = new TProfile3D("profile3D_bbc_west_psiShift","profile3D_bbc_west_psiShift",daynumber,0.5,daynumber+0.5,Ncentralities,0.5,Ncentralities+0.5,twoorder,0.5,twoorder+0.5,"");
  profile3D_bbc_west_psiShift->BuildOptions(-1.0,1.0,"");

  // EP resolutions plots
  Double_t xbin[8] = {0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0};

  TProfile *profile_correlation_tpc_east_epd_east = new TProfile("profile_correlation_tpc_east_epd_east","<cos(#psi^{TPC east}_{1} #minus #psi^{epd east}_{1})>",Ncentralities,0.5,Ncentralities+0.5,-1.0,1.0,"");
  profile_correlation_tpc_east_epd_east->GetXaxis()->SetTitle("Centrality (%)");
  profile_correlation_tpc_east_epd_east->GetYaxis()->SetTitle("Correlation");

  TProfile *profile_correlation_tpc_west_epd_east = new TProfile("profile_correlation_tpc_west_epd_east","<cos(#psi^{TPC west}_{1} #minus #psi^{epd east}_{1})>",Ncentralities,0.5,Ncentralities+0.5,-1.0,1.0,"");
  profile_correlation_tpc_west_epd_east->GetXaxis()->SetTitle("Centrality (%)");
  profile_correlation_tpc_west_epd_east->GetYaxis()->SetTitle("Correlation");

  TProfile *profile_correlation_tpc_east_tpc_west = new TProfile("profile_correlation_tpc_east_tpc_west","<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>",Ncentralities,0.5,Ncentralities+0.5,-1.0,1.0,"");
  profile_correlation_tpc_east_tpc_west->GetXaxis()->SetTitle("Centrality (%)");
  profile_correlation_tpc_east_tpc_west->GetYaxis()->SetTitle("Correlation");

  TProfile *profile_correlation_tpc_east_bbc_east = new TProfile("profile_correlation_tpc_east_bbc_east","<cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>",Ncentralities,0.5,Ncentralities+0.5,-1.0,1.0,"");
  profile_correlation_tpc_east_bbc_east->GetXaxis()->SetTitle("Centrality (%)");
  profile_correlation_tpc_east_bbc_east->GetYaxis()->SetTitle("Correlation");

  TProfile *profile_correlation_tpc_west_bbc_east = new TProfile("profile_correlation_tpc_west_bbc_east","<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})>",Ncentralities,0.5,Ncentralities+0.5,-1.0,1.0,"");
  profile_correlation_tpc_west_bbc_east->GetXaxis()->SetTitle("Centrality (%)");
  profile_correlation_tpc_west_bbc_east->GetYaxis()->SetTitle("Correlation");

  // END EP resolutions plots
  // EP correlations in 2D
  TH2D *correlation2D_epd_east_tpc_east = new TH2D("correlation2D_epd_east_tpc_east","#psi^{EPD east}_{1} vs. #psi^{TPC east}_{1}",50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  correlation2D_epd_east_tpc_east->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  correlation2D_epd_east_tpc_east->GetYaxis()->SetTitle("#psi^{EPD east}_{2} [Radian]");

  TH2D *correlation2D_epd_east_tpc_west = new TH2D("correlation2D_epd_east_tpc_west","#psi^{EPD east}_{1} vs. #psi^{TPC west}_{1}",50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  correlation2D_epd_east_tpc_west->GetXaxis()->SetTitle("#psi^{TPC west}_{2} [Radian]");
  correlation2D_epd_east_tpc_west->GetYaxis()->SetTitle("#psi^{EPD east}_{2} [Radian]");

  TH2D *correlation2D_tpc_sub = new TH2D("correlation2D_tpc_sub","#psi^{TPC east}_{1} vs. #psi^{TPC west}_{1}",50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  correlation2D_tpc_sub->GetXaxis()->SetTitle("#psi^{TPC west}_{2} [Radian]");
  correlation2D_tpc_sub->GetYaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");

  TH2D *correlation2D_bbc_east_tpc_east = new TH2D("correlation2D_bbc_east_tpc_east","#psi^{BBC east}_{1} vs. #psi^{TPC east}_{1}",50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  correlation2D_bbc_east_tpc_east->GetXaxis()->SetTitle("#psi^{TPC east}_{2} [Radian]");
  correlation2D_bbc_east_tpc_east->GetYaxis()->SetTitle("#psi^{BBC east}_{2} [Radian]");

  TH2D *correlation2D_bbc_east_tpc_west = new TH2D("correlation2D_bbc_east_tpc_west","#psi^{BBC east}_{1} vs. #psi^{TPC west}_{1}",50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  correlation2D_bbc_east_tpc_west->GetXaxis()->SetTitle("#psi^{TPC west}_{2} [Radian]");
  correlation2D_bbc_east_tpc_west->GetYaxis()->SetTitle("#psi^{BBC east}_{2} [Radian]");

  // END EP correlations in 2D


  // END Histograms

  // Read input files
  TFile *eventPlanes_input   = new TFile("/star/u/dchen/GitHub/fxtPicoAna/TpcEpCorrectionHistograms_INPUT.root","read");

  // Read resolution for systematic analysis
  TFile *t_resolution_iniput = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_event_plane_resolution/event_plane_resolution_primary.root","read");

  TH1D *h_resolution_input = (TH1D*)t_resolution_iniput->Get("h_resolution");
  Double_t d_resolution_input[Ncentralities]={0};
  for(int i = 0; i<Ncentralities;i++)
  {
    d_resolution_input[i]=h_resolution_input->GetBinContent(i+1);
    std::cout << Form("d_resolution_input[%d] = ", i) << d_resolution_input[i] << std::endl;
  }

  TProfile3D *profile3D_tpc_east_Qx_Qy_input = nullptr;
  TProfile3D *profile3D_tpc_east_psiShift_input = nullptr;
  TProfile3D *profile3D_tpc_west_Qx_Qy_input = nullptr;
  TProfile3D *profile3D_tpc_west_psiShift_input = nullptr;

  TH1D       *bbc_east_adc_profile_input = nullptr;
  TProfile3D *profile3D_bbc_east_Qx_Qy_input = nullptr;
  TProfile3D *profile3D_bbc_east_psiShift_input = nullptr;

  TH1D       *bbc_west_adc_profile_input = nullptr;
  TProfile3D *profile3D_bbc_west_Qx_Qy_input = nullptr;
  TProfile3D *profile3D_bbc_west_psiShift_input = nullptr;

  // EP Resolution
  TProfile *profile_correlation_tpc_east_tpc_west_input = nullptr;
  TProfile *profile_correlation_tpc_east_bbc_east_input = nullptr;
  TProfile *profile_correlation_tpc_west_bbc_east_input = nullptr;

  Double_t d_correlation_tpc_east_tpc_west_input[Ncentralities]={0};
  Double_t d_correlation_bbc_east_tpc_east_input[Ncentralities]={0};
  Double_t d_correlation_bbc_east_tpc_west_input[Ncentralities]={0};
  Double_t d_resolution[Ncentralities]={0};

  // Read particle efficiency table
  TH3D *ProtonEfficiencyTable;
  TH3D *PionPlusEfficiencyTable;
  TH3D *PionMinusEfficiencyTable;
  TH3D *KaonMinusEfficiencyTable;
  TH3D *KaonPlusEfficiencyTable;

  if(!eventPlanes_input->IsOpen()) std::cout<<"No EP file!"<<std::endl;
  if( eventPlanes_input->IsOpen())
  {
    std::cout<<"EP file loaded successfully!"<<std::endl;
    // Read TPC EP histograms
    profile3D_tpc_east_Qx_Qy_input    = (TProfile3D*)eventPlanes_input->Get("profile3D_tpc_east_Qx_Qy");
    profile3D_tpc_east_psiShift_input = (TProfile3D*)eventPlanes_input->Get("profile3D_tpc_east_psiShift");

    profile3D_tpc_west_Qx_Qy_input    = (TProfile3D*)eventPlanes_input->Get("profile3D_tpc_west_Qx_Qy");
    profile3D_tpc_west_psiShift_input = (TProfile3D*)eventPlanes_input->Get("profile3D_tpc_west_psiShift");

    // Read BBC EP histograms
    bbc_east_adc_profile_input        = (TH1D*)eventPlanes_input->Get("bbc_east_adc_profile");
    profile3D_bbc_east_Qx_Qy_input    = (TProfile3D*)eventPlanes_input->Get("profile3D_bbc_east_Qx_Qy");
    profile3D_bbc_east_psiShift_input = (TProfile3D*)eventPlanes_input->Get("profile3D_bbc_east_psiShift");

    bbc_west_adc_profile_input        = (TH1D*)eventPlanes_input->Get("bbc_west_adc_profile");
    profile3D_bbc_west_Qx_Qy_input    = (TProfile3D*)eventPlanes_input->Get("profile3D_bbc_west_Qx_Qy");
    profile3D_bbc_west_psiShift_input = (TProfile3D*)eventPlanes_input->Get("profile3D_bbc_west_psiShift");

    // Read EP correlation histograms
    profile_correlation_tpc_east_tpc_west_input = (TProfile*)eventPlanes_input->Get("profile_correlation_tpc_east_tpc_west");
    profile_correlation_tpc_east_bbc_east_input = (TProfile*)eventPlanes_input->Get("profile_correlation_tpc_east_bbc_east");
    profile_correlation_tpc_west_bbc_east_input = (TProfile*)eventPlanes_input->Get("profile_correlation_tpc_west_bbc_east");


    // Compute BBC ADC gain parameters
    if(bbc_east_adc_profile_input && bbc_east_adc_profile_input->GetEntries() > 0) {
        Double_t ADCtotal_inner = 0.0, ADCtotal_outer1 = 0.0, ADCtotal_outer2 = 0.0;
        for(Int_t i=1;i<=16;i++) {
            Double_t value1 = bbc_east_adc_profile_input->GetBinContent(i);
            if( i <= 6 ) ADCtotal_inner += value1;
            if( i == 7 || i == 10 || i == 12 || i == 15 ) ADCtotal_outer1 += value1;
            if( i == 8 || i ==  9 || i == 11 || i == 13 || i == 14 || i == 16 ) ADCtotal_outer2 += value1;
        }
        for(Int_t i=1;i<=16;i++) {
            Double_t value1 = bbc_east_adc_profile_input->GetBinContent(i);
            if( i <=  6 ) egain[i-1] = (value1 > 0.0)? ADCtotal_inner/(6.0*value1):0.0;
            if( i ==  7 || i == 12 ) egain[i-1] = (value1 > 0.0)? ADCtotal_outer1/(3.0*value1):0.0;
            if( i == 10 || i == 15 ) egain[i-1] = (value1 > 0.0)? ADCtotal_outer1/(6.0*value1):0.0;
            if( i ==  8 || i == 9 || i == 11 || i == 13 || i == 14 || i == 16 ) egain[i-1] = (value1 > 0.0)? ADCtotal_outer2/(6.0*value1):0.0;
        }
    }
    if(bbc_west_adc_profile_input && bbc_west_adc_profile_input->GetEntries() > 0) {
        Double_t ADCtotal_inner = 0.0, ADCtotal_outer1 = 0.0, ADCtotal_outer2 = 0.0;
        for(Int_t i=1;i<=16;i++) {
            Double_t value1 = bbc_west_adc_profile_input->GetBinContent(i);
            if( i <= 6 ) ADCtotal_inner += value1;
            if( i == 7 || i == 10 || i == 12 || i == 15 ) ADCtotal_outer1 += value1;
            if( i == 8 || i ==  9 || i == 11 || i == 13 || i == 14 || i == 16 ) ADCtotal_outer2 += value1;
        }
        for(Int_t i=1;i<=16;i++) {
            Double_t value1 = bbc_west_adc_profile_input->GetBinContent(i);
            if( i <=  6 ) wgain[i-1] = (value1 > 0.0)? ADCtotal_inner/(6.0*value1):0.0;
            if( i ==  7 || i == 12 ) wgain[i-1] = (value1 > 0.0)? ADCtotal_outer1/(3.0*value1):0.0;
            if( i == 10 || i == 15 ) wgain[i-1] = (value1 > 0.0)? ADCtotal_outer1/(6.0*value1):0.0;
            if( i ==  8 || i == 9 || i == 11 || i == 13 || i == 14 || i == 16 ) wgain[i-1] = (value1 > 0.0)? ADCtotal_outer2/(6.0*value1):0.0;
        }
    }
    // Read particle efficiency table

    TFile *PiPlusEffTableFile = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/PiPlusEffTable.root","read");
    if( !PiPlusEffTableFile->IsOpen() ) std::cout<<"No Pi+ efficiency table input!"<<std::endl;
    if(  PiPlusEffTableFile->IsOpen() ) {
        std::cout<<"Pi+ efficiency table loaded successfully!"<<std::endl;
        PionPlusEfficiencyTable = (TH3D*)PiPlusEffTableFile->Get("PiPlusEffTable");
    }
    TFile *PiMinusEffTableFile = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/PiMinusEffTable.root","read");
    if( !PiMinusEffTableFile->IsOpen() ) std::cout<<"No Pi- efficiency table input!"<<std::endl;
    if(  PiMinusEffTableFile->IsOpen() ) {
        std::cout<<"Pi- efficiency table loaded successfully!"<<std::endl;
		    PionMinusEfficiencyTable = (TH3D*)PiMinusEffTableFile->Get("PiMinusEffTable");
    }
    TFile *KMinusEffTableFile = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/KMinusEffTable.root","read");
    if( !KMinusEffTableFile->IsOpen() ) std::cout<<"No K- efficiency table input!"<<std::endl;
    if(  KMinusEffTableFile->IsOpen() ) {
        std::cout<<"K- efficiency table loaded successfully!"<<std::endl;
        KaonMinusEfficiencyTable = (TH3D*)KMinusEffTableFile->Get("KMinusEffTable");
    }
    TFile *KPlusEffTableFile = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/KPlusEffTable.root","read");
    if( !KPlusEffTableFile->IsOpen() ) std::cout<<"No K+ efficiency table input!"<<std::endl;
    if(  KPlusEffTableFile->IsOpen() ) {
        std::cout<<"K+ efficiency table loaded successfully!"<<std::endl;
        KaonPlusEfficiencyTable = (TH3D*)KPlusEffTableFile->Get("KPlusEffTable");
    }
	  TFile *ProtonEffTableFile = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/ProtonEffTable.root","read");
    if( !ProtonEffTableFile->IsOpen() ) std::cout<<"No Proton efficiency table input!"<<std::endl;
    if(  ProtonEffTableFile->IsOpen() ) {
        std::cout<<"Proton efficiency table loaded successfully!"<<std::endl;
		    ProtonEfficiencyTable = (TH3D*)ProtonEffTableFile->Get("ProtonEffTable");
    }

    // Compute EP resolutions
    if(profile_correlation_tpc_east_tpc_west_input&&profile_correlation_tpc_east_bbc_east_input&&profile_correlation_tpc_west_bbc_east_input)
    {
      for(Int_t i=0;i<Ncentralities;i++)
      {
        d_correlation_tpc_east_tpc_west_input[i] = profile_correlation_tpc_east_tpc_west_input->GetBinContent(i+1);
        d_correlation_bbc_east_tpc_east_input[i] = profile_correlation_tpc_east_bbc_east_input->GetBinContent(i+1);
        d_correlation_bbc_east_tpc_west_input[i] = profile_correlation_tpc_west_bbc_east_input->GetBinContent(i+1);
        Double_t value = d_correlation_bbc_east_tpc_east_input[i]*d_correlation_bbc_east_tpc_west_input[i]/d_correlation_tpc_east_tpc_west_input[i];
        if(value > 0) d_resolution[i] = TMath::Sqrt(value);
        std::cout<< Form("resolution %d = ",i+1) << d_resolution[i]<<std::endl;
      }
    }
std::cout<<ProtonEfficiencyTable<<std::endl;
std::cout<<PionPlusEfficiencyTable<<std::endl;
std::cout<<PionMinusEfficiencyTable<<std::endl;
std::cout<<KaonMinusEfficiencyTable<<std::endl;
std::cout<<KaonPlusEfficiencyTable<<std::endl;

  }
  TH1D *BBC_ADC_EAST[16], *BBC_ADC_WEST[16];
  // BBC ADC channels
  for(Int_t i=0;i<16;i++) {
      sprintf(name,"BBC_ADC_EAST%d",i+1);
      BBC_ADC_EAST[i] = new TH1D(name,name,4500,0.0,4500.0);
      BBC_ADC_EAST[i]->GetXaxis()->SetTitle("ADC");
      BBC_ADC_EAST[i]->GetYaxis()->SetTitle("# of events");

      sprintf(name,"BBC_ADC_WEST%d",i+1);
      BBC_ADC_WEST[i] = new TH1D(name,name,4500,0.0,4500.0);
      BBC_ADC_WEST[i]->GetXaxis()->SetTitle("ADC");
      BBC_ADC_WEST[i]->GetYaxis()->SetTitle("# of events");
  }
  // std::cout << "test 1 "<<std::endl;
  //////////////////////////// Event Loop //////////////////////////////////////
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
  {
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1)
    << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..."
        << std::endl;
        break;
    }
    // std::cout << "test 2 "<<std::endl;

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();


    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
        break;
    }

    ievtcut[0] += 1;

    const Float_t B = event->bField();
    double d_MagField = event->bField();

    Int_t nTracks = dst->numberOfTracks();
    Int_t runId = event->runId();
    // std::cout << "runId = " << runId << std::endl;

    Double_t Day = (Double_t)runId - 19151028.0;//- 20160023.0;
    hist_runId->Fill(Day);

    Double_t primaryVertex_X = (Double_t)event->primaryVertex().X();
    Double_t primaryVertex_Y = (Double_t)event->primaryVertex().Y();
    Double_t primaryVertex_Z = (Double_t)event->primaryVertex().Z();
    Double_t primaryVertex_perp = (Double_t)event->primaryVertex().Perp();
    hist_Vz->Fill(primaryVertex_Z);
    hist_Vr->Fill(primaryVertex_perp);
    hist_VyVx->Fill(primaryVertex_X,primaryVertex_Y);

    double  d_zvtx  = -9999.0;
    double  d_xvtx  = -9999.0;
    double  d_yvtx  = -9999.0;

    std::vector <unsigned int> triggerIDs;

    //============================ Trigger Selection ===========================
    triggerIDs.clear();
    triggerIDs       = event->triggerIds();
    bool b_bad_trig = true;

    // loop for the trigger ids and see if any == 1
    for(unsigned int i=0; i < triggerIDs.size(); i++)
      {
        Double_t d_trigger = (Double_t)triggerIDs[i] - 620050.0;
        hist_triggerID->Fill(d_trigger);
        if(triggerIDs[i] == 630802) b_bad_trig = false; // hlt_fixedTargetGood 7.2GeV
      }

    //=========================== End Trigger Slection =========================
    TVector3 pVtx     = event->primaryVertex();
    // Primary Vertex

    //=========================== Z-VTX Selection ==============================
    TVector3 v3D_vtx  = event->primaryVertex();
    d_zvtx = pVtx.z();
    d_xvtx = pVtx.x();
    d_yvtx = pVtx.y();
    // bool b_bad_zvtx  = ((d_zvtx < 210.0) || (d_zvtx > 212.0));
    // Insert systematic check cuts
    // bool b_bad_zvtx   = (cutID == 1) ? TMath::Abs(d_zvtx - 211.0)>(0.8 + 0.04*variationID) : ((d_zvtx < 210.0) || (d_zvtx > 212.0));
    bool b_bad_zvtx   =  ((d_zvtx < 199.0) || (d_zvtx > 202.0)); //FXT_26p5_2018
    bool b_bad_xvtx   =  ((d_xvtx < -1.0) || (d_xvtx > 1.0)); //FXT_26p5_2018
    bool b_bad_yvtx   =  ((d_yvtx < -3.0) || (d_yvtx > -0.5)); //FXT_26p5_2018
    bool b_bad_rvtx   =  primaryVertex_perp > 3.0;

    // bool b_bad_zvtx   =  ((d_zvtx < 210.0) || (d_zvtx > 212.0)); //FXT 4.5GeV 2016


    //======================== END Z-VTX Selection =============================

    bool b_bad_evt  = b_bad_zvtx || b_bad_trig || b_bad_xvtx || b_bad_yvtx || b_bad_rvtx;
    if(b_bad_evt) continue;

    hist_Vz_cut->Fill(primaryVertex_Z);
    hist_Vr_cut->Fill(primaryVertex_perp);
    hist_VyVx_cut->Fill(primaryVertex_X,primaryVertex_Y);
    // Bad Event Cut
    ievtcut[1] += 1;

    Int_t refMult = event->refMult();
    Int_t grefMult = event->grefMult();

    /***************** Centrality Track Loop to determine centrality ******************/
    Int_t centrality = 0;

    bool b_cent_01  = false;
    bool b_cent_02  = false;
    bool b_cent_03  = false;
    bool b_cent_04  = false;
    bool b_cent_05  = false;
    bool b_cent_06  = false;
    bool b_cent_07  = false;
    bool b_cent_08  = false;
    bool b_cent_09  = false;
    bool b_cent_10  = false;

    int nGoodTracks = 0;
    int nTrkvsCuts  = 0;
    // std::cout <<" test 3" << std::endl;
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

        // Retrieve i-th pico track
        StPicoTrack *picoTrack = dst->track(iTrk);

        if(!picoTrack) continue;

        hist_ratio->Fill(((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()));
        hist_nHits->Fill((double)picoTrack->nHitsFit());
        hist_ndEdx->Fill(picoTrack->nHitsDedx());
        hist_DCA->Fill(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z));
        //============== === Track Level Cuts ==============
        bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);

        bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
        //bool b_not_enough_hits = false;//((double)picoTrack->nHitsFit()) < 15;
        bool b_bad_track    = b_bad_dEdx || b_bad_tracking /*|| b_not_enough_hits*/;
        if(b_bad_track) continue;

        //============== END Track Level Cuts ==============
        //============= Primary Track Cut ===================
        if(!picoTrack->isPrimary()) continue;

        StPicoBTofPidTraits *trait = NULL;
        if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
        double d_tofBeta0        = -999;
        if(trait) d_tofBeta0 = trait->btofBeta();

        double d_px0      = picoTrack->gMom().x();
        double d_py0      = picoTrack->gMom().y();
        double d_pz0      = picoTrack->gMom().z();
        double d_pT0      = picoTrack->gPt();
        double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
        double mass2      = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

        nGoodTracks++;
        if(d_tofBeta0 == -999) continue;

        nTrkvsCuts++;
    }
    ////////////////////// END Centrality Track Loop for 7p2 ///////////////////////////

    bool b_pileup   = (nGoodTracks >  270);
    b_cent_01       = (nGoodTracks >= 200);// 0 - 10%
    b_cent_02       = (nGoodTracks >= 150 && nGoodTracks < 200);// 10 - 20%
    b_cent_03       = (nGoodTracks >= 124  && nGoodTracks < 150);// 20 - 30%
    b_cent_04       = (nGoodTracks >= 100  && nGoodTracks < 124);// 30 - 40%
    b_cent_05       = (nGoodTracks >= 72  && nGoodTracks < 100);// 40 - 50%
    b_cent_06       = (nGoodTracks >= 50  && nGoodTracks < 72);// 50 - 60%
    b_cent_07       = (nGoodTracks >= 40  && nGoodTracks < 50);// 60 - 70%
    b_cent_08       = (nGoodTracks >= 30  && nGoodTracks < 40);// 70 - 80%
    b_cent_09       = (nGoodTracks >= 20  && nGoodTracks < 30);// 80 - 90%
    b_cent_10       = (nGoodTracks >= 10  && nGoodTracks < 20);// >90%
    bool b_low_mult = (nGoodTracks <= 10);

    if(b_cent_01) centrality = 1;
    if(b_cent_02) centrality = 2;
    if(b_cent_03) centrality = 3;
    if(b_cent_04) centrality = 4;
    if(b_cent_05) centrality = 5;
    if(b_cent_06) centrality = 6;
    if(b_cent_07) centrality = 7;
    if(b_cent_08) centrality = 8;
    if(b_cent_09) centrality = 9;
    if(b_cent_10) centrality = 10;

    hist_realTrackMult->Fill(nGoodTracks);
    hist_cent->Fill(centrality);
    hist_realTrackMult_refmult->Fill(nGoodTracks,refMult);
    hist_realTrackMult_grefmult->Fill(nGoodTracks,grefMult);
    // hist_realTrackMult_tofmult->Fill(nGoodTracks,tofMult);
    if(b_pileup) continue;
    // pile-up cut
    ievtcut[2] += 1;

    //EPD EP result
    StEpdEpInfo result = mEpFinder->Results(mEpdHits,pVtx,0);  // and now you have all the EP info you could ever want :-)
    Double_t EastRawQx = (Double_t) result.EastRawQ(EpOrder).X();
    Double_t EastRawQy = (Double_t) result.EastRawQ(EpOrder).Y();
    Double_t WestRawQx = (Double_t) result.WestRawQ(EpOrder).X();
    Double_t WestRawQy = (Double_t) result.WestRawQ(EpOrder).Y();

    Double_t EastWeightedQx = (Double_t) result.EastPhiWeightedQ(EpOrder).X();
    Double_t EastWeightedQy = (Double_t) result.EastPhiWeightedQ(EpOrder).Y();
    Double_t WestWeightedQx = (Double_t) result.WestPhiWeightedQ(EpOrder).X();
    Double_t WestWeightedQy = (Double_t) result.WestPhiWeightedQ(EpOrder).Y();


    if(EastRawQx!=0 || EastRawQy!=0 )
    {
      hist2_Epd_east_Qy_Qx_raw->Fill(EastRawQx,EastRawQy);
      hist2_Epd_east_Qy_Qx_Weighted->Fill(EastWeightedQx,EastWeightedQy);

      hist_Epd_east_psi_raw->Fill(result.EastRawPsi(EpOrder));
      hist_Epd_east_psi_Weighted->Fill(result.EastPhiWeightedPsi(EpOrder));
      hist_Epd_east_psi_Shifted->Fill(result.EastPhiWeightedAndShiftedPsi(EpOrder));

    }
    if(WestRawQx!=0 || WestRawQy!=0 )
    {
      hist2_Epd_west_Qy_Qx_raw->Fill(WestRawQx,WestRawQy);
      hist2_Epd_west_Qy_Qx_Weighted->Fill(WestWeightedQx,WestWeightedQy);

      hist_Epd_west_psi_raw->Fill(result.WestRawPsi(EpOrder));
      hist_Epd_west_psi_Weighted->Fill(result.WestPhiWeightedPsi(EpOrder));
      hist_Epd_west_psi_Shifted->Fill(result.WestPhiWeightedAndShiftedPsi(EpOrder));

    }


    // Define event plane parameters
    Int_t N_tpc_east = 0, N_tpc_west = 0, N_thirdEP = 0;
    Double_t tpc_east_Qx = 0.0, tpc_east_Qy = 0.0, tpc_east_Qweight = 0.0;
    Double_t tpc_west_Qx = 0.0, tpc_west_Qy = 0.0, tpc_west_Qweight = 0.0;
    Double_t tpc_full_Qx = 0.0, tpc_full_Qy = 0.0;
    Double_t tpc_east_plane1 = -999.0, tpc_east_plane2 = -999.0, tpc_east_plane3 = -999.0;
    Double_t tpc_west_plane1 = -999.0, tpc_west_plane2 = -999.0, tpc_west_plane3 = -999.0;
    Double_t tpc_full_plane1 = -999.0, tpc_full_plane2 = -999.0, tpc_full_plane3 = -999.0;
    Double_t thirdEP_Qx = 0.0, thirdEP_Qy = 0.0, thirdEP_Qweight = 0.0;
    Double_t thirdEP_plane1 = -999.0, thirdEP_plane2 = -999.0, thirdEP_plane3 = -999.0;


    /////////// Primary Tracks loop to get TPC EP parameters //////////////
    for(Int_t itr=0;itr<nTracks;itr++)
    {
      // Get Track pointer
      StPicoTrack *picoTrack = dst->track(itr);
      itrkcut[0] += 1;

      if(!picoTrack) continue;
      // PicoTrack Cut
      itrkcut[1] += 1;

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut
      itrkcut[2] += 1;


      // bool b_bad_dEdx     = (cutID == 2) ? (picoTrack->nHitsDedx() <= 5*variationID):(picoTrack->nHitsDedx() <= 0);
      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);

      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      // bool b_bad_DCA      = (cutID == 3) ? (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= variationID):(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);
      // if(cutID == 3 && variationID == 0) b_bad_DCA = false;

      bool b_bad_DCA      = (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);

      bool b_bad_track    = b_bad_dEdx || b_bad_tracking || b_bad_DCA;

      if(b_bad_track) continue;
      itrkcut[3] += 1;

      Double_t pt = picoTrack->pPt();
      hist_pt->Fill(pt);
      Double_t pz = picoTrack->pMom().Z();
      Double_t eta = picoTrack->pMom().Eta();
      hist_eta->Fill(eta);

      // Get PID parameters
      Int_t charge = picoTrack->charge();
      Double_t trackP = picoTrack->pPtot();
      Double_t mass2 = 0.0;

      // Check if TOF info available
      StPicoBTofPidTraits *trait        = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      double d_tofBeta                  = -999;
      if(trait) d_tofBeta               = trait->btofBeta();
      if(d_tofBeta != -999) mass2 = trackP * trackP *( ( 1.0 / ( d_tofBeta*d_tofBeta ) ) - 1.0 );
      Double_t phi = picoTrack->pMom().Phi();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();
      hist_phi->Fill(phi);
      // Define flow weight
      Double_t w0 = 0.0, w1 = 0.0;

      // PID for TPC Event plane
      // Protons
      if(
        TMath::Abs(picoTrack->nSigmaProton()) <  2.0
        && (
            d_tofBeta != -999.0
            && mass2 > 0.5
            && mass2 < 1.5
        )
        && charge > 0
        && pt > 0.4
        && pt < 2.0
      )
      {
        // Get particle track rapidity
        Double_t energy_Proton = TMath::Sqrt(trackP*trackP + Mass_Proton*Mass_Proton);
        Double_t rap_Proton    = 0.5*TMath::Log( (energy_Proton + pz) / (energy_Proton - pz) );

        Double_t d_mT_Proton   = TMath::Sqrt(pt*pt + Mass_Proton*Mass_Proton);
        // Get eff corr
        Double_t efficiency = 0.0;//1.0
        Int_t i_ybin_Proton = ProtonEfficiencyTable->GetYaxis()->FindBin(rap_Proton);
        Int_t i_zbin_Proton = ProtonEfficiencyTable->GetZaxis()->FindBin(d_mT_Proton-Mass_Proton);

        efficiency = ProtonEfficiencyTable ->GetBinContent(centrality,i_ybin_Proton,i_zbin_Proton);
        efficiency = (efficiency >= 0.001 && efficiency <= 1.0)? 1.0 / efficiency : 0.0;
        efficiency = 1.0;
        // Get phi weight
        Double_t phi_weight = 1.0;

        if(efficiency > 0.0)
        {
          w0 = pt; // pt as weight to calculate v2. // rap_Proton + 1.52;//1.4144;
          w1 = efficiency * phi_weight;
          //w1 = phi_weight;
        }
      }

      // Pions
      if( //TMath::Abs( track->nSigmaPion() ) < 2.0
         picoTrack->nSigmaPion() > -2.0 && picoTrack->nSigmaPion() < 2.0
         && ( d_tofBeta != -999.0
               && ( (  picoTrack->nSigmaPion() > 2.0 && picoTrack->nSigmaPion() < 6.0
                      && ( mass2 < -0.005 || mass2 > 0.005 )
                    ) || (picoTrack->nSigmaPion() > -4.0 && picoTrack->nSigmaPion() <= 2.0) )
               && mass2 > -0.1
               && mass2 < 0.15
              )
         && pt > 0.2
         && pt < 1.6

         )
         {
           // Get particle track rapidity
           Double_t energy_Pion = TMath::Sqrt(trackP*trackP + Mass_Pion*Mass_Pion);
           Double_t rap_Pion = 0.5*TMath::Log( (energy_Pion + pz) / (energy_Pion - pz) );
           if(charge > 0)
           {
             // Get eff corr
             // Double_t efficiency = 1.0;

             Double_t d_mT_Pion   = TMath::Sqrt(pt*pt + Mass_Pion*Mass_Pion);
             // Get eff corr
             Double_t efficiency = 0.0;//1.0
             Int_t i_ybin_Pion = PionPlusEfficiencyTable->GetYaxis()->FindBin(rap_Pion);
             Int_t i_zbin_Pion = PionPlusEfficiencyTable->GetZaxis()->FindBin(d_mT_Pion-Mass_Pion);

             efficiency = PionPlusEfficiencyTable ->GetBinContent(centrality,i_ybin_Pion,i_zbin_Pion);
             efficiency = (efficiency >= 0.001 && efficiency <= 1.0)? 1.0 / efficiency : 0.0;
             efficiency = 1.0;

             if(efficiency > 0.0)
             {
               w0 = pt; // pt as weight to calculate v2. // rap_Pion + 1.52;//1.4144;
               w1 = efficiency /* -1.0*/;
             }
           }
           if(charge < 0)
           {
             // Get eff corr
             // Double_t efficiency = 1.0;

             Double_t d_mT_Pion   = TMath::Sqrt(pt*pt + Mass_Pion*Mass_Pion);
             // Get eff corr
             Double_t efficiency = 0.0;//1.0
             Int_t i_ybin_Pion = PionMinusEfficiencyTable->GetYaxis()->FindBin(rap_Pion);
             Int_t i_zbin_Pion = PionMinusEfficiencyTable->GetZaxis()->FindBin(d_mT_Pion-Mass_Pion);

             efficiency = PionMinusEfficiencyTable ->GetBinContent(centrality,i_ybin_Pion,i_zbin_Pion);
             efficiency = (efficiency >= 0.001 && efficiency <= 1.0)? 1.0 / efficiency : 0.0;
             efficiency = 1.0;

             if(efficiency > 0.0)
             {
               w0 = pt; // pt as weight to calculate v2. // rap_Pion + 1.52;//1.4144;
               w1 = efficiency /* -1.0*/;
             }
           }
         }

      // Kaons
      // Remove Kaons for EP to get rid of auto correlation
      // END PID
      Bool_t IsEast = kFALSE, IsWest = kFALSE;

      Double_t randomNumber = gRandom->Uniform(1);
      //if( w0 < -0.05 ) IsEast = kTRUE;
      //if( w0 >  0.05 ) IsWest = kTRUE;
      if( randomNumber < 0.5 ) IsEast = kTRUE;
      if( randomNumber >= 0.5 ) IsWest = kTRUE;
      // East
      if( IsEast )
      {
          N_tpc_east++;
          tpc_east_Qx += w0*w1 * TMath::Cos(EpOrder * phi);
          tpc_east_Qy += w0*w1 * TMath::Sin(EpOrder * phi);
          tpc_east_Qweight += TMath::Abs(w0*w1);
      }
      // West
      if( IsWest )
      {
          N_tpc_west++;
          tpc_west_Qx += w0*w1 * TMath::Cos(EpOrder * phi);
          tpc_west_Qy += w0*w1 * TMath::Sin(EpOrder * phi);
          tpc_west_Qweight += TMath::Abs(w0*w1);
      }

    }
    /////// END Primary Tracks loop to get event plane parameters //////////////

    // Compute TPC east event planes
    if(N_tpc_east >= 5 && tpc_east_Qweight > 0.0)
    {
      tpc_east_Qx /= tpc_east_Qweight;
      tpc_east_Qy /= tpc_east_Qweight;
      if(tpc_east_Qx || tpc_east_Qy)
      {
        tpc_east_plane1 = (1. / EpOrder) * TMath::ATan2(tpc_east_Qy,tpc_east_Qx);
        if(tpc_east_plane1 < 0.0                             ) tpc_east_plane1 += (1. / EpOrder) * 2.0*TMath::Pi();
        if(tpc_east_plane1 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_east_plane1 -= (1. / EpOrder) * 2.0*TMath::Pi();
        hist_tpc_east_psi_raw->Fill(tpc_east_plane1);
        // Recenter reaction plane vector
        profile3D_tpc_east_Qx_Qy->Fill(Day,centrality,1,tpc_east_Qx);
        profile3D_tpc_east_Qx_Qy->Fill(Day,centrality,2,tpc_east_Qy);
        if(eventPlanes_input->IsOpen() && profile3D_tpc_east_Qx_Qy_input && profile3D_tpc_east_Qx_Qy_input->GetEntries() > 0)
        {
          tpc_east_Qx -= profile3D_tpc_east_Qx_Qy_input->GetBinContent(Day,centrality,1);
          tpc_east_Qy -= profile3D_tpc_east_Qx_Qy_input->GetBinContent(Day,centrality,2);
        }
        if(tpc_east_Qx || tpc_east_Qy)
        {
          tpc_east_plane2 = (1. / EpOrder) * TMath::ATan2(tpc_east_Qy,tpc_east_Qx);
          if(tpc_east_plane2 < 0.0                             ) tpc_east_plane2 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(tpc_east_plane2 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_east_plane2 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_tpc_east_psi_recentered->Fill(tpc_east_plane2);
          // Shift recenterd event plane to flat
          Double_t reaction_plane_new = tpc_east_plane2;
          for(Int_t k=0;k<order;k++)
          {
            profile3D_tpc_east_psiShift->Fill(Day,centrality,1+2*k,TMath::Cos(EpOrder * reaction_plane_new*(k+1)));
            profile3D_tpc_east_psiShift->Fill(Day,centrality,2+2*k,TMath::Sin(EpOrder * reaction_plane_new*(k+1)));
          }
          Double_t psi_mean[twoorder];
          for(Int_t i=0;i<twoorder;i++)
          {
            psi_mean[i] = 0.0;
          }
          if(eventPlanes_input->IsOpen() && profile3D_tpc_east_psiShift_input && profile3D_tpc_east_psiShift_input->GetEntries() > 0)
          {
            for(Int_t k=0;k<order;k++)
            {
              psi_mean[0+2*k] = profile3D_tpc_east_psiShift_input->GetBinContent(Day,centrality,1+2*k);
              psi_mean[1+2*k] = profile3D_tpc_east_psiShift_input->GetBinContent(Day,centrality,2+2*k);
            }

          }
          for(Int_t k=0;k<order;k++)
          {
            reaction_plane_new += (1. / EpOrder) * ( -2.0*psi_mean[1+2*k] * TMath::Cos( EpOrder * tpc_east_plane2*(k+1) )
                                                    + 2.0*psi_mean[0+2*k] * TMath::Sin( EpOrder * tpc_east_plane2*(k+1) ) ) / (k+1);
          }
          tpc_east_plane3 = reaction_plane_new;
          if(tpc_east_plane3 < 0.0                             ) tpc_east_plane3 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(tpc_east_plane3 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_east_plane3 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_tpc_east_psi_flattened->Fill(tpc_east_plane3);
        }
      }
    }
    // END Compute TPC east event planes

    // Compute TPC west event planes
    if(N_tpc_west >= 5 && tpc_west_Qweight > 0.0)
    {
      tpc_west_Qx /= tpc_west_Qweight;
      tpc_west_Qy /= tpc_west_Qweight;
      if(tpc_west_Qx || tpc_west_Qy)
      {
        tpc_west_plane1 = (1. / EpOrder) * TMath::ATan2(tpc_west_Qy,tpc_west_Qx);
        if(tpc_west_plane1 < 0.0                             ) tpc_west_plane1 += (1. / EpOrder) * 2.0*TMath::Pi();
        if(tpc_west_plane1 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_west_plane1 -= (1. / EpOrder) * 2.0*TMath::Pi();
        hist_tpc_west_psi_raw->Fill(tpc_west_plane1);
        // Recenter reaction plane vector
        profile3D_tpc_west_Qx_Qy->Fill(Day,centrality,1,tpc_west_Qx);
        profile3D_tpc_west_Qx_Qy->Fill(Day,centrality,2,tpc_west_Qy);
        if(eventPlanes_input->IsOpen() && profile3D_tpc_west_Qx_Qy_input && profile3D_tpc_west_Qx_Qy_input->GetEntries() > 0)
        {
          tpc_west_Qx -= profile3D_tpc_west_Qx_Qy_input->GetBinContent(Day,centrality,1);
          tpc_west_Qy -= profile3D_tpc_west_Qx_Qy_input->GetBinContent(Day,centrality,2);
        }
        if(tpc_west_Qx || tpc_west_Qy)
        {
          tpc_west_plane2 = (1. / EpOrder) * TMath::ATan2(tpc_west_Qy,tpc_west_Qx);
          if(tpc_west_plane2 < 0.0                             ) tpc_west_plane2 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(tpc_west_plane2 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_west_plane2 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_tpc_west_psi_recentered->Fill(tpc_west_plane2);
          // Shift recenterd event plane to flat
          Double_t reaction_plane_new = tpc_west_plane2;
          for(Int_t k=0;k<order;k++)
          {
            profile3D_tpc_west_psiShift->Fill(Day,centrality,1+2*k,TMath::Cos(EpOrder * reaction_plane_new*(k+1)));
            profile3D_tpc_west_psiShift->Fill(Day,centrality,2+2*k,TMath::Sin(EpOrder * reaction_plane_new*(k+1)));
          }
          Double_t psi_mean[twoorder];
          for(Int_t i=0;i<twoorder;i++)
          {
              psi_mean[i] = 0.0;
          }
          if(eventPlanes_input->IsOpen() && profile3D_tpc_west_psiShift_input && profile3D_tpc_west_psiShift_input->GetEntries() > 0)
          {
            for(Int_t k=0;k<order;k++)
            {
              psi_mean[0+2*k] = profile3D_tpc_west_psiShift_input->GetBinContent(Day,centrality,1+2*k);
              psi_mean[1+2*k] = profile3D_tpc_west_psiShift_input->GetBinContent(Day,centrality,2+2*k);
            }
          }
          for(Int_t k=0;k<order;k++)
          {
            reaction_plane_new += (1. / EpOrder) * ( -2.0*psi_mean[1+2*k] * TMath::Cos( EpOrder * tpc_west_plane2*(k+1) )
                                                    + 2.0*psi_mean[0+2*k] * TMath::Sin( EpOrder * tpc_west_plane2*(k+1) ) ) / (k+1);
          }
          tpc_west_plane3 = reaction_plane_new;
          if(tpc_west_plane3 < 0.0                             ) tpc_west_plane3 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(tpc_west_plane3 > (1. / EpOrder) * 2.0*TMath::Pi()) tpc_west_plane3 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_tpc_west_psi_flattened->Fill(tpc_west_plane3);
        }
      }
    }

    // END Compute TPC west event planes

    // Compute bbc event planes
    Int_t N_bbc_east = 0, N_bbc_west = 0;
    Double_t bbc_east_Qx = 0.0, bbc_east_Qy = 0.0, bbc_east_Qweight = 0.0;
    Double_t bbc_west_Qx = 0.0, bbc_west_Qy = 0.0, bbc_west_Qweight = 0.0;
    Double_t bbc_full_Qx = 0.0, bbc_full_Qy = 0.0;
    Double_t bbc_east_plane1 = -999.0, bbc_east_plane2 = -999.0, bbc_east_plane3 = -999.0;
    Double_t bbc_west_plane1 = -999.0, bbc_west_plane2 = -999.0, bbc_west_plane3 = -999.0;
    Double_t bbc_full_plane1 = -999.0, bbc_full_plane2 = -999.0, bbc_full_plane3 = -999.0;

    // Get BBC psi
    Double_t nHitE[16], nHitW[16];
    for(Int_t i=0;i<16;i++)
    {
        nHitE[i] = 0.0;
        nHitW[i] = 0.0;
    }
    Double_t etotal = 0.0, wtotal = 0.0;
    for(Int_t i=0;i<16;i++)
    {
        Double_t eadc = event->bbcAdcEast(i); // BBC east
        BBC_ADC_EAST[i]->Fill(eadc);
        Double_t wadc = event->bbcAdcWest(i); // BBC west
        BBC_ADC_WEST[i]->Fill(wadc);
        if(eadc > 30 && eadc < 4000)
        {
            nHitE[i] = (egain[i] > 0.0)? eadc * egain[i] : eadc;
            etotal += nHitE[i];

            emean[i] += eadc;
            esum[i]++;

            emean_c[i] += nHitE[i];
            esum_c[i]++;
        }
        if(wadc > 30)
        {
            nHitW[i] = (wgain[i] > 0.0)? wadc * wgain[i] : wadc;
            wtotal += nHitW[i];

            wmean[i] += wadc;
            wsum[i]++;

            wmean_c[i] += nHitW[i];
            wsum_c[i]++;
        }
    }

    for(Int_t i=0;i<16;i++)
    {
        nHitE[i] *= (etotal > 0.0)? 1.0 / etotal : 0.0;
        nHitW[i] *= (wtotal > 0.0)? 1.0 / wtotal : 0.0;
    }
    for(Int_t iTile=0;iTile<16;iTile++)
    {
      // BBC east
      if(nHitE[iTile] > 0.0)
      {
          N_bbc_east++;
          bbc_east_Qx += TMath::Cos( EpOrder * BBC_GetPhi(0,iTile,primaryVertex_X,primaryVertex_Y) ) * nHitE[iTile];
          bbc_east_Qy += TMath::Sin( EpOrder * BBC_GetPhi(0,iTile,primaryVertex_X,primaryVertex_Y) ) * nHitE[iTile];
          bbc_east_Qweight += nHitE[iTile];
      }
      // BBC west
      if(nHitW[iTile] > 0.0)
      {
          N_bbc_west++;
          bbc_west_Qx += TMath::Cos( EpOrder * BBC_GetPhi(1,iTile,primaryVertex_X,primaryVertex_Y) ) * nHitW[iTile];
          bbc_west_Qy += TMath::Sin( EpOrder * BBC_GetPhi(1,iTile,primaryVertex_X,primaryVertex_Y) ) * nHitW[iTile];
          bbc_west_Qweight += nHitW[iTile];
      }
    }
    // std::cout << "test 4 "<<std::endl;

    // BBC east
    if(N_bbc_east >= 5 && bbc_east_Qweight > 0.0)
    {
      bbc_east_Qx /= bbc_east_Qweight;
      bbc_east_Qy /= bbc_east_Qweight;
      if(bbc_east_Qx || bbc_east_Qy)
      {
        bbc_east_plane1 = (1. / EpOrder) * TMath::ATan2(bbc_east_Qy,bbc_east_Qx);
        if(bbc_east_plane1 < 0.0                             ) bbc_east_plane1 += (1. / EpOrder) * 2.0*TMath::Pi();
        if(bbc_east_plane1 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_east_plane1 -= (1. / EpOrder) * 2.0*TMath::Pi();
        hist_bbc_east_psi_raw->Fill(bbc_east_plane1);

        // Recenter event plane vector
        profile3D_bbc_east_Qx_Qy->Fill(Day,centrality,1,bbc_east_Qx);
        profile3D_bbc_east_Qx_Qy->Fill(Day,centrality,2,bbc_east_Qy);
        if(eventPlanes_input->IsOpen() && profile3D_bbc_east_Qx_Qy_input && profile3D_bbc_east_Qx_Qy_input->GetEntries() > 0)
        {
            bbc_east_Qx -= profile3D_bbc_east_Qx_Qy_input->GetBinContent(Day,centrality,1);
            bbc_east_Qy -= profile3D_bbc_east_Qx_Qy_input->GetBinContent(Day,centrality,2);
        }
        if(bbc_east_Qx || bbc_east_Qy)
        {
          bbc_east_plane2 = (1. / EpOrder) * TMath::ATan2(bbc_east_Qy,bbc_east_Qx);
          if(bbc_east_plane2 < 0.0                             ) bbc_east_plane2 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(bbc_east_plane2 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_east_plane2 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_bbc_east_psi_recentered->Fill(bbc_east_plane2);

          // Shift recentered event plane to flat
          Double_t reaction_plane_new = bbc_east_plane2;
          for(Int_t k=0;k<order;k++)
          {
              profile3D_bbc_east_psiShift->Fill(Day,centrality,1+2*k,TMath::Cos(EpOrder * reaction_plane_new*(k+1)));
              profile3D_bbc_east_psiShift->Fill(Day,centrality,2+2*k,TMath::Sin(EpOrder * reaction_plane_new*(k+1)));
          }
          Double_t psi_mean[twoorder] = {0.0};
          if(eventPlanes_input->IsOpen() && profile3D_bbc_east_psiShift_input && profile3D_bbc_east_psiShift_input->GetEntries() > 0)
          {
            for(Int_t k=0;k<order;k++)
            {
              psi_mean[0+2*k] = profile3D_bbc_east_psiShift_input->GetBinContent(Day,centrality,1+2*k);
              psi_mean[1+2*k] = profile3D_bbc_east_psiShift_input->GetBinContent(Day,centrality,2+2*k);
            }
          }
          for(Int_t k=0;k<order;k++)
          {
            reaction_plane_new += (1. / EpOrder) * ( -2.0*psi_mean[1+2*k] * TMath::Cos( EpOrder * bbc_east_plane2*(k+1) )
                                                    + 2.0*psi_mean[0+2*k] * TMath::Sin( EpOrder * bbc_east_plane2*(k+1) ) ) / (k+1);

          }
          bbc_east_plane3 = reaction_plane_new;
          if(bbc_east_plane3 < 0.0                             ) bbc_east_plane3 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(bbc_east_plane3 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_east_plane3 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_bbc_east_psi_flattened->Fill(bbc_east_plane3);
        }
      }
    }
    // BBC west
    if(N_bbc_west >= 2 && bbc_west_Qweight > 0.0)
    {
      bbc_west_Qx /= bbc_west_Qweight;
      bbc_west_Qy /= bbc_west_Qweight;
      if(bbc_west_Qx || bbc_west_Qy)
      {
        bbc_west_plane1 = (1. / EpOrder) * TMath::ATan2(bbc_west_Qy,bbc_west_Qx);
        if(bbc_west_plane1 < 0.0                             ) bbc_west_plane1 += (1. / EpOrder) * 2.0*TMath::Pi();
        if(bbc_west_plane1 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_west_plane1 -= (1. / EpOrder) * 2.0*TMath::Pi();
        hist_bbc_west_psi_raw->Fill(bbc_west_plane1);

        // Recenter event plane vector
        profile3D_bbc_west_Qx_Qy->Fill(Day,centrality,1,bbc_west_Qx);
        profile3D_bbc_west_Qx_Qy->Fill(Day,centrality,2,bbc_west_Qy);
        if(eventPlanes_input->IsOpen() && profile3D_bbc_west_Qx_Qy_input && profile3D_bbc_west_Qx_Qy_input->GetEntries() > 0)
        {
            bbc_west_Qx -= profile3D_bbc_west_Qx_Qy_input->GetBinContent(Day,centrality,1);
            bbc_west_Qy -= profile3D_bbc_west_Qx_Qy_input->GetBinContent(Day,centrality,2);
        }
        if(bbc_west_Qx || bbc_west_Qy)
        {
          bbc_west_plane2 = (1. / EpOrder) * TMath::ATan2(bbc_west_Qy,bbc_west_Qx);
          if(bbc_west_plane2 < 0.0                             ) bbc_west_plane2 += (1. / EpOrder) * 2.0*TMath::Pi();
          if(bbc_west_plane2 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_west_plane2 -= (1. / EpOrder) * 2.0*TMath::Pi();
          hist_bbc_west_psi_recentered->Fill(bbc_west_plane2);
        }
        // Shift recentered event plane to flat
        Double_t reaction_plane_new = bbc_west_plane2;
        for(Int_t k=0;k<order;k++)
        {
            profile3D_bbc_west_psiShift->Fill(Day,centrality,1+2*k,TMath::Cos(EpOrder * reaction_plane_new*(k+1)));
            profile3D_bbc_west_psiShift->Fill(Day,centrality,2+2*k,TMath::Sin(EpOrder * reaction_plane_new*(k+1)));
        }
        Double_t psi_mean[twoorder] = {0.0};
        if(eventPlanes_input->IsOpen() && profile3D_bbc_west_psiShift_input && profile3D_bbc_west_psiShift_input->GetEntries() > 0)
        {
          for(Int_t k=0;k<order;k++)
          {
            psi_mean[0+2*k] = profile3D_bbc_west_psiShift_input->GetBinContent(Day,centrality,1+2*k);
            psi_mean[1+2*k] = profile3D_bbc_west_psiShift_input->GetBinContent(Day,centrality,2+2*k);
          }
        }
        for(Int_t k=0;k<order;k++)
        {
          reaction_plane_new += (1. / EpOrder) * ( -2.0*psi_mean[1+2*k] * TMath::Cos( EpOrder * bbc_west_plane2*(k+1) )
                                                  + 2.0*psi_mean[0+2*k] * TMath::Sin( EpOrder * bbc_west_plane2*(k+1) ) ) / (k+1);
        }
        bbc_west_plane3 = reaction_plane_new;
        if(bbc_west_plane3 < 0.0                             ) bbc_west_plane3 += (1. / EpOrder) * 2.0*TMath::Pi();
        if(bbc_west_plane3 > (1. / EpOrder) * 2.0*TMath::Pi()) bbc_west_plane3 -= (1. / EpOrder) * 2.0*TMath::Pi();
        hist_bbc_west_psi_flattened->Fill(bbc_west_plane3);
      }

    }

    // Accumulate Event Plane correlations
    // Double_t centTMP = Ncentralities-centrality+1;
    if( tpc_east_plane3 >= 0.0 && tpc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() )
    {
      if( tpc_west_plane3 >= 0.0 && tpc_west_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() )
      {
        profile_correlation_tpc_east_tpc_west->Fill(centrality,TMath::Cos(EpOrder * (tpc_east_plane3 - tpc_west_plane3)));
        correlation2D_tpc_sub->Fill(tpc_west_plane3,tpc_east_plane3);
      }
      if( bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() )
      {
        profile_correlation_tpc_east_bbc_east->Fill(centrality,TMath::Cos(EpOrder * (tpc_east_plane3 - bbc_east_plane3 - TMath::Pi())));
        correlation2D_bbc_east_tpc_east->Fill(tpc_east_plane3,bbc_east_plane3);
      }
      if( result.EastPhiWeightedAndShiftedPsi(EpOrder) >= 0.0 && result.EastPhiWeightedAndShiftedPsi(EpOrder) <= (1. / EpOrder) * 2.0*TMath::Pi() )
      {
        profile_correlation_tpc_east_epd_east->Fill(centrality,TMath::Cos(EpOrder * (tpc_east_plane3 - result.EastPhiWeightedAndShiftedPsi(EpOrder) )));
        correlation2D_epd_east_tpc_east->Fill(tpc_east_plane3,result.EastPhiWeightedAndShiftedPsi(EpOrder));
      }
    }
    if( tpc_west_plane3 >= 0.0 && tpc_west_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() )
    {
      if( bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() )
      {
        profile_correlation_tpc_west_bbc_east->Fill(centrality,TMath::Cos(EpOrder * (tpc_west_plane3 - bbc_east_plane3 - TMath::Pi())));
        correlation2D_bbc_east_tpc_west->Fill(tpc_west_plane3,bbc_east_plane3);
      }
      if( result.EastPhiWeightedAndShiftedPsi(EpOrder) >= 0.0 && result.EastPhiWeightedAndShiftedPsi(EpOrder) <= (1. / EpOrder) * 2.0*TMath::Pi() )
      {
        profile_correlation_tpc_west_epd_east->Fill(centrality,TMath::Cos(EpOrder * (tpc_west_plane3 - result.EastPhiWeightedAndShiftedPsi(EpOrder) )));
        correlation2D_epd_east_tpc_west->Fill(tpc_west_plane3,result.EastPhiWeightedAndShiftedPsi(EpOrder));
      }
    }

    // END Accumulate Event Plane correlations
    // std::cout << "test 5 "<<std::endl;

    // Get Three EP resolutions
    Double_t d_resolution2 = TMath::Cos(EpOrder * (tpc_west_plane3 - bbc_east_plane3 - TMath::Pi()))*TMath::Cos(EpOrder * (tpc_east_plane3 - bbc_east_plane3 - TMath::Pi()))
                            /TMath::Cos(EpOrder * (tpc_east_plane3 - tpc_west_plane3));
    Double_t d_resolution_test = 0;
    if(d_resolution2>0) d_resolution_test = TMath::Sqrt(d_resolution2);



    // Define PID parameters
    Int_t Nprotons = 0, NpionPlus = 0, NpionMinus = 0, NkaonPlus = 0, NkaonMinus = 0;
    //////////////////// Primary Tracks loop to get flow ///////////////////////
    for(Int_t itr=0;itr<nTracks;itr++)
    {
      // Get Track pointer
      StPicoTrack *picoTrack = dst->track(itr);

      if(!picoTrack) continue;
      // PicoTrack Cut

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      // bool b_bad_dEdx     = (cutID == 2)?  (picoTrack->nHitsDedx() <= 5*variationID):(picoTrack->nHitsDedx() <= 0);
      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);

      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      // bool b_bad_DCA      = (cutID == 3) ? (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= variationID):(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);
      // if(cutID == 3 && variationID == 0) b_bad_DCA = false;
      bool b_bad_DCA      = (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);

      bool b_bad_track    = b_bad_dEdx || b_bad_tracking || b_bad_DCA;

      if(b_bad_track) continue;

      // Get track-wise parameter histograms
      Double_t pt = picoTrack->pPt();
      Double_t pz = picoTrack->pMom().Z();
      Double_t eta = picoTrack->pMom().Eta();
      Double_t phi = picoTrack->pMom().Phi();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();

      // Get PID parameters
      Int_t charge = picoTrack->charge();
      Double_t Beta = -999.0;
      Double_t trackP = picoTrack->pPtot();
      Double_t mass2 = 0.0;
      hist_dEdx->Fill(charge*trackP,picoTrack->dEdx());

      // Check if TOF info available
      if( picoTrack->isTofTrack() )
      {
          // Retrieve corresponding trait
          StPicoBTofPidTraits *trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
          if( trait ) {
              Beta = trait->btofBeta();
          }
      }

      if(Beta != -999.0)
      {
        mass2 = trackP*trackP * ( ( 1.0 / ( Beta*Beta ) ) - 1.0 );
        hist_beta->Fill(charge*trackP,1.0/Beta); //cout<<charge*trackP<<endl;
        hist_mass->Fill(charge*trackP,mass2);
        if(charge > 0)
        {
            GPC_request_m2_vs_nSigmaProton->Fill(picoTrack->nSigmaProton(),mass2);
            GPC_request_m2_vs_nSigmaPionPlus->Fill(picoTrack->nSigmaPion(),mass2);
            GPC_request_m2_vs_nSigmaKaonPlus->Fill(picoTrack->nSigmaKaon(),mass2);

            GPC_request_m2_vs_dEdx_Pos->Fill(picoTrack->dEdx(),mass2);
        }
        if(charge < 0)
        {
            GPC_request_m2_vs_nSigmaPionMinus->Fill(picoTrack->nSigmaPion(),mass2);
            GPC_request_m2_vs_nSigmaKaonMinus->Fill(picoTrack->nSigmaKaon(),mass2);

            GPC_request_m2_vs_dEdx_Neg->Fill(picoTrack->dEdx(),mass2);
        }
      }
      // 2nd particle identifications for flow calculation (1st PID is for TPC EP calculation)
      // Protons
      if(
        TMath::Abs(picoTrack->nSigmaProton()) <  2.0
        && (
            Beta != -999.0
            && mass2 > 0.5
            && mass2 < 1.5
        )
        && charge > 0
        && pt > 0.4
        && pt < 2.0
      )
      {
        // Count proton tracks number
        Nprotons++;
        // Get particle track rapidity
        Double_t energy_Proton = TMath::Sqrt(trackP*trackP + Mass_Proton*Mass_Proton);
        Double_t rap_Proton = 0.5*TMath::Log( (energy_Proton + pz) / (energy_Proton - pz) );

        // Fill histograms
        hist_pt_proton->Fill(pt);
        hist_eta_proton->Fill(eta);
        hist_y_proton->Fill(rap_Proton);
        hist_rap_eta_proton->Fill(eta,rap_Proton);
        hist_pt_y_proton->Fill(rap_Proton,pt,1/*efficiency*/);
        hist_pt_eta_proton->Fill(eta,pt,1/*efficiency*/);
        hist_phi_proton->Fill(phi);
        hist_dEdx_proton->Fill(charge*trackP,picoTrack->dEdx());
        if(Beta != -999.0)
        {
            hist_beta_proton->Fill(charge*trackP,1.0/Beta);
            hist_mass_proton->Fill(charge*trackP,mass2);
        }

        // Accumulate flows
        if( /*chooseBBCeastEP &&*/ bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() /*&& d_resolution[centrality-1]>0*//*&& efficiency > 0.0 */) {
            profile3D_proton_v2->Fill(centrality,pt,rap_Proton/*track->eta()*/,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi()))/*/d_resolution[centrality-1]*/,/*efficiency*/1.0);
            h2_proton_v2       ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
        }
        // if( /*chooseTPCEP &&*/ reaction_plane3_ex_west[itr] >= 0.0 && reaction_plane3_ex_west[itr] <= 2.0*TMath::Pi() && efficiency > 0.0 ) {
        //     profile3D_proton_v2_tpc->Fill(centrality,pt,rap_Proton/*track->eta()*/,TMath::Cos(1.0*phi - reaction_plane3_ex_west[itr]),/*efficiency*/1.0);
        // }

      }

      // Pions
      if(
        picoTrack->nSigmaPion() > -2.0 && picoTrack->nSigmaPion() < 2.0
        && ( Beta != -999.0
            && ( (  picoTrack->nSigmaPion() > 2.0 && picoTrack->nSigmaPion() < 6.0
                  && ( mass2 < -0.005 || mass2 > 0.005 )
                  ) || (picoTrack->nSigmaPion() > -4.0 && picoTrack->nSigmaPion() <= 2.0) )
            && mass2 > -0.1
            && mass2 < 0.15
           )
            //&& ( ( trackP < 0.62 && TMath::Abs(mass2) > 0.005 ) || trackP > 0.62 )
        && pt > 0.2
        && pt < 1.6
      )
      {
        // Get particle track rapidity
        Double_t energy_Pion = TMath::Sqrt(trackP*trackP + Mass_Pion*Mass_Pion);
        Double_t rap_Pion = 0.5*TMath::Log( (energy_Pion + pz) / (energy_Pion - pz) );
        // pionPlus
        if(charge > 0)
        {
          // Get eff corr
          Double_t efficiency = 1.0;

          // Count pionPlus tracks number
          NpionPlus++;

          // Fill histograms
          hist_pt_pionPlus->Fill(pt);
          hist_eta_pionPlus->Fill(eta);
          hist_y_pionPlus->Fill(rap_Pion);
          hist_rap_eta_pionPlus->Fill(eta,rap_Pion);
          hist_pt_y_pionPlus->Fill(rap_Pion,pt,efficiency);
          hist_pt_eta_pionPlus->Fill(eta,pt,efficiency);
          hist_phi_pionPlus->Fill(phi);
          hist_dEdx_pionPlus->Fill(charge*trackP,picoTrack->dEdx());
          if(Beta != -999.0)
          {
              hist_beta_pionPlus->Fill(charge*trackP,1.0/Beta);
              hist_mass_pionPlus->Fill(charge*trackP,mass2);
          }
          // Accumulate flows
          if( /*chooseBBCeastEP &&*/ bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() && efficiency > 0.0  /*&& d_resolution[centrality-1]>0*/) {
              profile3D_pionPlus_v2->Fill(centrality,pt,rap_Pion/*track->eta()*/,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi()))/*/d_resolution[centrality-1]*/,efficiency/*1.0*/);
              h2_pionPlus_v2       ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
              h2_pions_v2          ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
          }

        }
        // pionMinus
        if(charge < 0)
        {
          // Get eff corr
          Double_t efficiency = 1.0;
          // Count pionMinus tracks number
          NpionMinus++;
          // Fill histograms
          hist_pt_pionMinus->Fill(pt);
          hist_eta_pionMinus->Fill(eta);
          hist_y_pionMinus->Fill(rap_Pion);
          hist_rap_eta_pionMinus->Fill(eta,rap_Pion);
          hist_pt_y_pionMinus->Fill(rap_Pion,pt,efficiency);
          hist_pt_eta_pionMinus->Fill(eta,pt,efficiency);
          hist_phi_pionMinus->Fill(phi);
          hist_dEdx_pionMinus->Fill(charge*trackP,picoTrack->dEdx());
          if(Beta != -999.0) {
              hist_beta_pionMinus->Fill(charge*trackP,1.0/Beta);
              hist_mass_pionMinus->Fill(charge*trackP,mass2);
          }
          // Accumulate flows
          if( /*chooseBBCeastEP &&*/ bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() && efficiency > 0.0 /*&& d_resolution[centrality-1]>0*/) {
              profile3D_pionMinus_v2->Fill(centrality,pt,rap_Pion/*track->eta()*/,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi()))/*/d_resolution[centrality-1]*/,efficiency/*1.0*/);
              h2_pionMinus_v2       ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
              h2_pions_v2           ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
          }

        }
      }

      // Kaons
      if( picoTrack->nSigmaKaon() > -2.0 && picoTrack->nSigmaKaon() < 2.0
         && ( Beta != -999.0
             && mass2 > 0.19//0.17
             && mass2 < 0.3//0.33
            )
         && pt > 0.2
         // && pt < 1.6
      )
      {
        // Get particle track rapidity
        Double_t energy_Kaon = TMath::Sqrt(trackP*trackP + Mass_Kaon*Mass_Kaon);
        Double_t rap_Kaon = 0.5*TMath::Log( (energy_Kaon + pz) / (energy_Kaon - pz) );
        // kaonPlus
        if(charge > 0)
        {
          // Get eff corr
          Double_t efficiency = 1.0;

          // Count kaonPlus tracks number
          NkaonPlus++;
          // Fill histograms
          hist_pt_kaonPlus->Fill(pt);
          hist_eta_kaonPlus->Fill(eta);
          hist_y_kaonPlus->Fill(rap_Kaon);
          hist_rap_eta_kaonPlus->Fill(eta,rap_Kaon);
          // hist_pt_y_kaonPlus->Fill(rap_Kaon,pt,efficiency);
          hist_pt_eta_kaonPlus->Fill(eta,pt,efficiency);
          hist_phi_kaonPlus->Fill(phi);
          hist_dEdx_kaonPlus->Fill(charge*trackP,picoTrack->dEdx());

          if(Beta != -999.0) {
              hist_beta_kaonPlus->Fill(charge*trackP,1.0/Beta);
              hist_mass_kaonPlus->Fill(charge*trackP,mass2);
          }
          // Accumulate flows
          if( /*chooseBBCeastEP &&*/ bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() && efficiency > 0.0 /*&& d_resolution[centrality-1]>0*/) {
              profile3D_kaonPlus_v2->Fill(centrality,pt,rap_Kaon/*track->eta()*/,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi()))/*/d_resolution[centrality-1]*/,efficiency/*1.0*/);
              h2_kaonPlus_v2       ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
          }

        }
        // kaonMinus
        if(charge < 0)
        {
          // Get eff corr
          Double_t efficiency = 1.0;
          NkaonMinus++;
          // fill histograms
          hist_pt_kaonMinus->Fill(pt);
          hist_eta_kaonMinus->Fill(eta);
          hist_y_kaonMinus->Fill(rap_Kaon);
          hist_rap_eta_kaonMinus->Fill(eta,rap_Kaon);
          // hist_pt_y_kaonMinus->Fill(rap_Kaon,pt,efficiency);
          hist_pt_eta_kaonMinus->Fill(eta,pt,efficiency);
          hist_phi_kaonMinus->Fill(phi);
          hist_dEdx_kaonMinus->Fill(charge*trackP,picoTrack->dEdx());
          if(Beta != -999.0) {
              hist_beta_kaonMinus->Fill(charge*trackP,1.0/Beta);
              hist_mass_kaonMinus->Fill(charge*trackP,mass2);
          }
          // Accumulate flows
          if( /*chooseBBCeastEP &&*/ bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() && efficiency > 0.0 /*&& d_resolution[centrality-1]>0*/) {
              profile3D_kaonMinus_v2->Fill(centrality,pt,rap_Kaon/*track->eta()*/,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi()))/*/d_resolution[centrality-1]*/,efficiency/*1.0*/);
              h2_kaonMinus_v2       ->Fill(pt,TMath::Cos(EpOrder * (phi - bbc_east_plane3 - TMath::Pi())),1);
          }

        }

      }

    }
    //////////////// END Primary Tracks loop to get flow ///////////////////////
    hist_trackmult_proton->Fill(Nprotons);
    hist_trackmult_pionPlus->Fill(NpionPlus);
    hist_trackmult_pionMinus->Fill(NpionMinus);
    hist_trackmult_kaonPlus->Fill(NkaonPlus);
    hist_trackmult_kaonMinus->Fill(NkaonMinus);

    // Phi meson flow analysis
    //========================= Track Cut Settings ===============================================
    int i_Nhits_min = 9;  //minimum number of hits

    //Mon Jul  3 09:17:54 EDT 2017 from phi paper
    double d_pT_min = 0.1;//0.05;
    double d_pT_max = 10.0;//0.05;
    double d_mom_min = 0.1;//0.05;
    double d_mom_max = 10.0;// 1.0;//0.05;
    double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut // Mon Jul  3 09:16:46 EDT 2017

    //--- Not being used Fri Jun 16 10:48:30 EDT 2017
    // double d_PRO_daughter_DCA = 0.3; // trk must be greater than this dca
    // double d_PI_daughter_DCA  = 1.0; // trk must be greater than this dca

    // Lambda: 1.0, K0s 1.2
    //  double d_cut_dca_daughters    = 2.0; //distance between daughters, must be less than this
    // double d_cut_dca_daughters_lam    = 1.0;//2.0; //distance between daughters, must be less than this
    // double d_cut_dca_daughters_k0s    = 1.2; //distance between daughters, must be less than this


    // Lambda: 1.5 Anti Lambda: 2.0 K0s: 1.0
    //double d_cut_dca_mother       = 5.0; // must be less than this
    // double d_cut_dca_mother_lam       = 1.5;//5.0; // must be less than this
    // double d_cut_dca_mother_k0s       = 1.0; // must be less than this

    // Lambda: 3.0, K0s: 2.0
    //  double d_cut_mother_decay_length = 3.0; // must be greater than this
    // double d_cut_mother_decay_length_lam = 3.0; // must be greater than this
    // double d_cut_mother_decay_length_k0s = 2.0; // must be greater than this
    // double d_cut_mother_decay_length_RHO = 0.5; // must be LESS    than this
    // double d_cut_mother_decay_length_PHI = 0.5; // must be LESS    than this
    //======================== END Track Settings ================================================

    //======================= Primary Track Loop to kaon tracks =================================================
    vector<StPicoTrack *> v_pri_tracks;
    vector<StPicoTrack *> v_pri_tracks_pl;
    vector<StPicoTrack *> v_pri_tracks_mi;

    int index = 0;

    double d_PI_m    = 0.13957018;
    double d_PRO_m   = 0.9382720813;
    double d_K_m     = 0.493677;

    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);
      // Retrieve i-th pico track

      if(picoTrack == NULL)       continue;

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      // To get beta from Btof to calculate mass2
      StPicoBTofPidTraits *trait = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );

      //nHits minimum cut
      unsigned short nHits = picoTrack->nHits();
      if(nHits < i_Nhits_min)     continue;

      bool b_bad_dEdx      = false;
      bool b_bad_tracking  = false;
      bool b_bad_ToF       = false;

      // b_bad_dEdx     = (cutID == 20)?  (picoTrack->nHitsDedx() <= 5*variationID):(picoTrack->nHitsDedx() <= 0);
      b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);

      if(trait) b_bad_ToF       = (trait->btof() <= 0.0);
      bool b_bad_track          = b_bad_dEdx || b_bad_tracking;
      if(b_bad_track)              continue;
      // Bad Track Cut

      double d_TPCnSigmaElectron = fabs(picoTrack->nSigmaElectron());
      double d_TPCnSigmaPion   = fabs(picoTrack->nSigmaPion());
      double d_TPCnSigmaProton = fabs(picoTrack->nSigmaProton());
      double d_TPCnSigmaKaon   = fabs(picoTrack->nSigmaKaon());

      double tofBeta       = -999;
      if(trait) tofBeta    = trait->btofBeta();

      double d_tofBeta0    = -999;
      if(trait) d_tofBeta0 = trait->btofBeta();

      // test new sets of PID cuts
      bool b_PI  = false; //fabs(d_TPCnSigmaPion)   < d_SigmaCutLevel;
      bool b_PRO = false; //fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
      bool b_K   = false; //fabs(d_TPCnSigmaKaon)   < d_SigmaCutLevel;
      bool b_E   = false; //fabs(d_TPCnSigmaElectron)< d_SigmaCutLevel;

      double d_charge     = picoTrack->charge();
      double d_px0        = picoTrack->pMom().x();
      double d_py0        = picoTrack->pMom().y();
      double d_pz0        = picoTrack->pMom().z();
      double d_pT0        = picoTrack->pPt();
      double d_mom0       = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

      double mass2        = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
      double d_E_PRO      = sqrt(d_PRO_m*d_PRO_m + d_mom0*d_mom0);
      double d_E_PI       = sqrt(d_PI_m*d_PI_m + d_mom0*d_mom0);
      double d_E_K        = sqrt(d_K_m*d_K_m + d_mom0*d_mom0);

      double d_y_PRO      = 0.5*TMath::Log((d_E_PRO + d_pz0)/(d_E_PRO - d_pz0));
      double d_y_PI       = 0.5*TMath::Log((d_E_PI + d_pz0)/(d_E_PI - d_pz0));
      double d_y_K        = 0.5*TMath::Log((d_E_K + d_pz0)/(d_E_K - d_pz0));

      if(d_charge > 0.0)
        {
          if( b_E
            && (fabs(d_TPCnSigmaElectron) < 3.0)
            && (mass2 < 0.005)
            && (mass2 > -0.008)
            && (d_pT0 < 0.3)){
              b_E = true; b_PI = false; b_PRO = false; b_K = false;
            // h_E_plus_pT -> Fill(d_pT0);
            // h_E_plus_y  -> Fill(d_y_K);
            // h2_E_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_E_plus_mT_Diff -> Fill((d_mT_E - d_E_m));
          }
          // kaon PID
                  if( picoTrack->nSigmaKaon() > -2.0 && picoTrack->nSigmaKaon() < 2.0
                     && ( d_tofBeta0 != -999.0
                         && mass2 > 0.19//0.17
                         && mass2 < 0.3//0.33
                        )
                     && d_pT0 > 0.2
                     && d_pT0 < 1.6
                  ){
                       b_E = false; b_PI = false; b_PRO = false; b_K = true;
                    // h_K_plus_pT -> Fill(d_pT0);
                    // h_K_plus_y  -> Fill(d_y_K);
                    // h2_K_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
                    // h_K_plus_mT_Diff -> Fill((d_mT_K - d_K_m));
                  }

          if(
            picoTrack->nSigmaPion() > -2.0 && picoTrack->nSigmaPion() < 2.0
            && ( d_tofBeta0 != -999.0
                && ( (  picoTrack->nSigmaPion() > 2.0 && picoTrack->nSigmaPion() < 6.0
                      && ( mass2 < -0.005 || mass2 > 0.005 )
                      ) || (picoTrack->nSigmaPion() > -4.0 && picoTrack->nSigmaPion() <= 2.0) )
                && mass2 > -0.1
                && mass2 < 0.15
               )
                //&& ( ( trackP < 0.62 && TMath::Abs(mass2) > 0.005 ) || trackP > 0.62 )
            && d_pT0 > 0.2
            && d_pT0 < 1.6
          ){
            b_E = false; b_PI = true; b_PRO = false; b_K = false;
            // h_PI_plus_pT -> Fill(d_pT0);
            // h_PI_plus_y  -> Fill(d_y_PI);
            // h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            // h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
          }
          if(
            TMath::Abs(picoTrack->nSigmaProton()) <  2.0
            && (
                d_tofBeta0 != -999.0
                && mass2 > 0.5
                && mass2 < 1.5
            )
            && d_pT0 > 0.4
            && d_pT0 < 2.0
          )
            {
              b_E = false; b_PI = false; b_PRO = true; b_K = false;
            // h_PRO_plus_pT -> Fill(d_pT0);
            // h_PRO_plus_y  -> Fill(d_y_PRO);
            // h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
            // h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
          }
        }
      else
        {
          if( b_E
            && (fabs(d_TPCnSigmaElectron) < 3.0)
            && (mass2 < 0.005)
            && (mass2 > -0.008)
            && (d_pT0 < 0.3)){
              b_E = true; b_PI = false; b_PRO = false; b_K = false;
            // h_E_minus_pT -> Fill(d_pT0);
            // h_E_minus_y  -> Fill(d_y_K);
            // h2_E_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_E_minus_mT_Diff -> Fill((d_mT_E - d_E_m));
          }

          // kaon PID
                  if( picoTrack->nSigmaKaon() > -2.0 && picoTrack->nSigmaKaon() < 2.0
                     && ( d_tofBeta0 != -999.0
                         && mass2 > 0.19//0.17
                         && mass2 < 0.3//0.33
                        )
                     && d_pT0 > 0.2
                     && d_pT0 < 1.6
                  ){
                       b_E = false; b_PI = false; b_PRO = false; b_K = true;
                    // h_K_plus_pT -> Fill(d_pT0);
                    // h_K_plus_y  -> Fill(d_y_K);
                    // h2_K_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
                    // h_K_plus_mT_Diff -> Fill((d_mT_K - d_K_m));
                  }

          if(
            picoTrack->nSigmaPion() > -2.0 && picoTrack->nSigmaPion() < 2.0
            && ( d_tofBeta0 != -999.0
                && ( (  picoTrack->nSigmaPion() > 2.0 && picoTrack->nSigmaPion() < 6.0
                      && ( mass2 < -0.005 || mass2 > 0.005 )
                      ) || (picoTrack->nSigmaPion() > -4.0 && picoTrack->nSigmaPion() <= 2.0) )
                && mass2 > -0.1
                && mass2 < 0.15
               )
                //&& ( ( trackP < 0.62 && TMath::Abs(mass2) > 0.005 ) || trackP > 0.62 )
            && d_pT0 > 0.2
            && d_pT0 < 1.6
          ){
            b_E = false; b_PI = true; b_PRO = false; b_K = false;
            // h_PI_minus_pT -> Fill(d_pT0);
            // h_PI_minus_y  -> Fill(d_y_PI);
            // h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            // h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
          }
          if(
            TMath::Abs(picoTrack->nSigmaProton()) <  2.0
            && (
                d_tofBeta0 != -999.0
                && mass2 > 0.5
                && mass2 < 1.5
            )
            && d_pT0 > 0.4
            && d_pT0 < 2.0
          )
            {
              b_E = false; b_PI = false; b_PRO = true; b_K = false;
            // h_PRO_minus_pT -> Fill(d_pT0);
            // h_PRO_minus_y  -> Fill(d_y_PRO);
            // h2_PRO_minus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
            // h_PRO_minus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
          }
        }

        if(d_charge > 0.0)
          {
            if(b_PI||b_PRO||b_E)    continue;
            if(!b_K)                continue;
          }
        else
          {
            if(b_E||b_PI||b_PRO)    continue;
            if(!b_K)                continue;
          }
        // Kaon Cut

        if( (d_pT0<d_pT_min)   || (d_pT0 > d_pT_max))   continue;
        if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max)) continue;
        //pT Min, Mom Min Cut

        bool b_bad_DCA       = false;
        StPicoPhysicalHelix trackhelix = picoTrack->helix(B);
        double helixpathl              = trackhelix.pathLength(v3D_vtx, false);
        TVector3 v3D_dca               = trackhelix.at(helixpathl)-v3D_vtx;
        double d_helix_DCA_r           = v3D_dca.Mag();
        // double d_DCA_r_cut             = (cutID == 21) ? variationID : 4.0;
        double d_DCA_r_cut             = 3.0;
        TVector3 v3D_obj_DCA           = picoTrack->gDCA(pVtx);
        double d_obj_DCA               = v3D_obj_DCA.Mag();

        h_K_DCA_r      -> Fill(d_helix_DCA_r);
        h_K_obj_DCA_r  -> Fill(d_obj_DCA);
        h_K_diff_DCA_r -> Fill(d_helix_DCA_r-d_obj_DCA);

        if(d_helix_DCA_r > d_DCA_r_cut) b_bad_DCA == true;
        // if(cutID == 21 && variationID == 0) b_bad_DCA == false;
        if(b_bad_DCA)              continue;
        //Kaon DCA Cut

        Double_t efficiency = 1.0 ;
        if(d_charge > 0.0)
        {
          hist_pt_y_kaonPlus->Fill(d_y_K,d_pT0,efficiency);
          v_pri_tracks_pl.push_back(picoTrack);// check the kaon QA plots when put into the TTree
        }
        else if(d_charge < 0.0)
        {
          hist_pt_y_kaonMinus->Fill(d_y_K,d_pT0,efficiency);
          v_pri_tracks_mi.push_back(picoTrack);
        }
        v_pri_tracks.push_back(picoTrack);

        index++;

    }
    //=================== END Primary Track Loop =================================================

    //======== Invariant Mass Nested Primary Track Loop to get Phi meson flow ====================

    vector<StPicoTrack *> v_pri_tracks0 = v_pri_tracks_pl;
    vector<StPicoTrack *> v_pri_tracks1 = v_pri_tracks_mi;

    for(unsigned int i = 0; i < v_pri_tracks0.size();i++)
    {
      StPicoTrack * picoTrack0 = v_pri_tracks0[i];
      if(!picoTrack0) continue;

      // double d_TPCnSigmaPion0   = fabs(picoTrack0->nSigmaPion());
      // double d_TPCnSigmaProton0 = fabs(picoTrack0->nSigmaProton());
      double d_TPCnSigmaKaon0   = fabs(picoTrack0->nSigmaKaon());

      // bool b_PI0  = fabs(d_TPCnSigmaPion0)   < d_SigmaCutLevel;
      // bool b_PRO0 = fabs(d_TPCnSigmaProton0) < d_SigmaCutLevel;
      // bool b_K0   = fabs(d_TPCnSigmaKaon0)   < d_SigmaCutLevel;

      // if( b_PI0
      //    && (d_TPCnSigmaPion0 < d_TPCnSigmaProton0)
      //    && (d_TPCnSigmaPion0 < d_TPCnSigmaKaon0) )
      //   { b_PI0 = true; b_PRO0 = false;}// b_K0 = false;}
      //
      // if( b_PRO0
      //    && (d_TPCnSigmaProton0 < d_TPCnSigmaPion0)
      //    && (d_TPCnSigmaProton0 < d_TPCnSigmaKaon0) )
      //   { b_PRO0 = true; b_PI0 = false;}// b_K0 = false;}
      //
      // if( b_K0
      //    && (d_TPCnSigmaKaon0 < d_TPCnSigmaProton0)
      //    && (d_TPCnSigmaKaon0 < d_TPCnSigmaPion0) )
      //   { b_K0 = true; b_PRO0 = false; b_PI0 = false;}

      double d_M0       = -9999.0;
      double d_px0      = picoTrack0->pMom().x();
      double d_py0      = picoTrack0->pMom().y();
      double d_pz0      = picoTrack0->pMom().z();
      double d_pT0      = picoTrack0->pPt();
      double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
      double d_charge0  = picoTrack0->charge();
      double d_mc0      = d_mom0/d_charge0;
      double d_eta0     = picoTrack0->pMom().Eta();
      double d_phi0     = picoTrack0->pMom().Phi();
      TVector3 v3D_obj_p0 = picoTrack0->pMom();

      /*
      StPicoBTofPidTraits *trait = NULL;
      if(picoTrack0->isTofTrack()) trait = dst->btofPidTraits(picoTrack0->bTofPidTraitsIndex());
      double d_tofBeta0 = -999;
      if(trait) d_tofBeta0 = trait->btofBeta();

      StPicoPhysicalHelix    trackhelix0 = picoTrack0->helix(B);

      if(b_PRO0)
        {
          d_M0 = d_PRO_m;
        }
      else if(b_PI0)
        {
          d_M0 = d_PI_m;
        }
      else if(b_K0)
        {
          d_M0 = d_K_m;
        }*/
      d_M0 = d_K_m;//test

      for(unsigned int j = 0; j < v_pri_tracks1.size(); j++)
      {
        StPicoTrack * picoTrack1 = v_pri_tracks1[j];

        if(!picoTrack1 || (picoTrack0->id() == picoTrack1->id())) continue;
        /*
        double d_TPCnSigmaPion1   = fabs(picoTrack1->nSigmaPion());
        double d_TPCnSigmaProton1 = fabs(picoTrack1->nSigmaProton());
        */
        double d_TPCnSigmaKaon1   = fabs(picoTrack1->nSigmaKaon());

        /*
        bool b_PI1  = fabs(d_TPCnSigmaPion1) < d_SigmaCutLevel;
        bool b_PRO1 = fabs(d_TPCnSigmaProton1) < d_SigmaCutLevel;
        bool b_K1   = fabs(d_TPCnSigmaKaon1) < d_SigmaCutLevel;
        if( b_PI1
           && (d_TPCnSigmaPion1 < d_TPCnSigmaProton1)
           && (d_TPCnSigmaPion1 < d_TPCnSigmaKaon1) )
          { b_PI1 = true; b_PRO1 = false;}// b_K1 = false;}

        if( b_PRO1
           && (d_TPCnSigmaProton1 < d_TPCnSigmaPion1)
           && (d_TPCnSigmaProton1 < d_TPCnSigmaKaon1) )
          { b_PRO1 = true; b_PI1 = false;}// b_K1 = false;}

        if( b_K1
           && (d_TPCnSigmaKaon1 < d_TPCnSigmaProton1)
           && (d_TPCnSigmaKaon1 < d_TPCnSigmaPion1) )
          { b_K1 = true; b_PRO1 = false; b_PI1 = false;}

        double d_M1 = -9999.0;

        if(b_PRO1) d_M1 = d_PRO_m;
        else if(b_PI1)  d_M1 = d_PI_m;
        else if(b_K1)   d_M1 = d_K_m;
        */
        double d_M1 = d_K_m;

        double d_px1      = picoTrack1->pMom().x();
        double d_py1      = picoTrack1->pMom().y();
        double d_pz1      = picoTrack1->pMom().z();
        double d_pT1      = picoTrack1->pPt();
        double d_mom1     = sqrt(d_pT1*d_pT1 + d_pz1*d_pz1);
        double d_charge1    = picoTrack1->charge();
        double d_mc1        = d_mom1/d_charge1;
        double d_eta1       = picoTrack1->pMom().Eta();
        double d_phi1       = picoTrack1->pMom().Phi();
        TVector3 v3D_obj_p1 = picoTrack1->pMom();
        /*
        StPicoBTofPidTraits *trait = NULL;
        if(picoTrack1->isTofTrack()) trait = dst->btofPidTraits(picoTrack1->bTofPidTraitsIndex());
        double d_tofBeta1 = -999;
        if(trait) d_tofBeta1 = trait->btofBeta();
        */
        if(d_charge0 == d_charge1) continue;
        bool b_PHI    = true;//test

        /*
        bool b_PHI    = b_K0 && b_K1;
        bool b_RHO    = b_PI0 && b_PI1;
        bool b_K0S    = b_RHO;


        bool b_LAMBDA = (b_PRO0 && b_PI1)||(b_PRO1 && b_PI0);
        bool b_V0 = b_PHI || b_RHO || b_LAMBDA;
        if(!b_V0) continue;

        StPicoPhysicalHelix    trackhelix1 = picoTrack1->helix(B);

        pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);

        TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, d_MagField*kilogauss);
        TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, d_MagField*kilogauss);

        TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
        TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);

        double d_dca_products =(v3D_x_daughter0-v3D_x_daughter1).Mag();

        if(d_dca_products > d_cut_dca_daughters_lam) b_LAMBDA = false;
        if(d_dca_products > d_cut_dca_daughters_k0s) b_RHO    = false;

        TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)*0.5;
        TVector3 v3D_xvec_decayl = v3D_x_mother - v3D_vtx;
        TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;

        double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);

        double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );

        if(d_dca_mother > d_cut_dca_mother_lam) b_LAMBDA = false;
        if(d_dca_mother > d_cut_dca_mother_k0s) b_RHO    = false;

        double d_mother_decay_length =  v3D_xvec_decayl.Mag();

        if(d_mother_decay_length < d_cut_mother_decay_length_lam) b_LAMBDA = false;
        if(d_mother_decay_length < d_cut_mother_decay_length_k0s) b_K0S    = false;
        if(d_mother_decay_length > d_cut_mother_decay_length_RHO) b_RHO    = false;
        */
        double d_E0 = sqrt(v3D_obj_p0.Mag2()+d_M0*d_M0);
        double d_E1 = sqrt(v3D_obj_p1.Mag2()+d_M1*d_M1);

        double d_y0 = 0.5*TMath::Log((d_E0 + d_pz0)/(d_E0 - d_pz0));
        double d_y1 = 0.5*TMath::Log((d_E1 + d_pz1)/(d_E1 - d_pz1));

        double d_mT0        = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
        double d_mT1        = sqrt(d_pT1*d_pT1 + d_M1*d_M1);
        /*
        // efficiency correction
        Int_t i_ybin0 = KaonPlusEfficiencyTable->GetYaxis()->FindBin(d_y0);
        Int_t i_zbin0 = KaonPlusEfficiencyTable->GetZaxis()->FindBin(d_mT0-d_M0);

        Int_t i_ybin1 = KaonMinusEfficiencyTable->GetYaxis()->FindBin(d_y1);
        Int_t i_zbin1 = KaonMinusEfficiencyTable->GetZaxis()->FindBin(d_mT1-d_M1);
        */
        double d_eff_corr0 = 1;
        double d_eff_corr1 = 1;

        /*
        double d_eff_corr0 = KaonPlusEfficiencyTable ->GetBinContent(centrality,i_ybin0,i_zbin0);
        double d_eff_corr1 = KaonMinusEfficiencyTable->GetBinContent(centrality,i_ybin1,i_zbin1);

        // efficiency corrections
        d_eff_corr0 = (d_eff_corr0 <= 0.01 || d_eff_corr0 >= 1) ? 1 : 1 / d_eff_corr0;
        d_eff_corr1 = (d_eff_corr1 <= 0.01 || d_eff_corr1 >= 1) ? 1 : 1 / d_eff_corr1;

        */
        double d_inv_m = sqrt(d_M0*d_M0
                              +d_M1*d_M1
                              +2.0*d_E0*d_E1
                              -2.0*(v3D_obj_p0.Dot(v3D_obj_p1)) );

        // if(d_mother_decay_length > d_cut_mother_decay_length_PHI) b_PHI    = false;
        // Decay Length Cut

        double d_dip_angle = TMath::ACos((d_pT0*d_pT1+d_pz0*d_pz1) / (d_mom0*d_mom1) );
        h_dip_angle->Fill(d_dip_angle);
        // if(d_dip_angle < 0.04) b_PHI = false;
        // Dip Angle Cut

        // double d_mother_m = -9999;

        double d_pT_phi = sqrt(d_px0*d_px0 + d_py0*d_py0 +d_px1*d_px1 +d_py1+d_py1 + 2.*d_px0*d_px1 + 2.*d_py0*d_py1);
    	  double m_phi = 1.019455;
    	  double d_mT_phi = sqrt(d_pT_phi*d_pT_phi + m_phi*m_phi );

    	  double d_phi_pz = d_pz0+d_pz1;
    	  double d_phi_E  = d_E0+d_E1;
    	  double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;

        // double d_phi_azimuth = v3D_p_mother.Phi();
        /*
        if(d_phi_azimuth < 0.0            ) d_phi_azimuth += 2.0*TMath::Pi();
        if(d_phi_azimuth > 2.0*TMath::Pi()) d_phi_azimuth -= 2.0*TMath::Pi();
        */
        if(b_PHI) h_prim_inv_m_PHI    -> Fill(d_inv_m);
        double d_v2_raw_phi = -999.0;
        /*
        if(  bbc_east_plane3 >= 0.0 && bbc_east_plane3 <= (1. / EpOrder) * 2.0*TMath::Pi() ) {
           d_v2_raw_phi = TMath::Cos(EpOrder * (d_phi_azimuth - bbc_east_plane3 - TMath::Pi()));
        }
        */
        if(/*b_cent_01||b_cent_02||*/b_cent_07) continue;
        if(d_v2_raw_phi == -999.0) continue;
        // std::cout<< d_v2_raw_phi <<std::endl;
        // if(b_cent_03) d_v2_raw_phi /= 0.398419;
        // if(b_cent_04) d_v2_raw_phi /= 0.454645;
        // if(b_cent_05) d_v2_raw_phi /= 0.498949;
        // if(b_cent_06) d_v2_raw_phi /= 0.596082;

        // Add EP resolution later

        // if(b_cent_01) d_v2_raw_phi /= d_resolution[0];
        // if(b_cent_02) d_v2_raw_phi /= d_resolution[1];
        // if(b_cent_03) d_v2_raw_phi /= d_resolution[2];
        // if(b_cent_04) d_v2_raw_phi /= d_resolution[3];
        // if(b_cent_05) d_v2_raw_phi /= d_resolution[4];
        // if(b_cent_06) d_v2_raw_phi /= d_resolution[5];

        // std::cout<< d_v2_raw_phi <<std::endl;

        h2_phi_v2_vs_invM                                        -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_phi_y>=-1.5 && d_phi_y<-1.0) h2_phi_v2_vs_invM_bin2 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_phi_y>=-1.0 && d_phi_y<-0.5) h2_phi_v2_vs_invM_bin3 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_phi_y>=-0.5 && d_phi_y<=0.0) h2_phi_v2_vs_invM_bin4 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);

        // fill in pT bin
        if(d_pT_phi>=0 && d_pT_phi<0.5) h2_phi_v2_vs_invM_pTbin1 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=0.5 && d_pT_phi<1.0) h2_phi_v2_vs_invM_pTbin2 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=1.0 && d_pT_phi<1.5) h2_phi_v2_vs_invM_pTbin3 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=1.5 && d_pT_phi<2.0) h2_phi_v2_vs_invM_pTbin4 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=2.0 && d_pT_phi<2.5) h2_phi_v2_vs_invM_pTbin5 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=2.5 && d_pT_phi<3.0) h2_phi_v2_vs_invM_pTbin6 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);
        if(d_pT_phi>=3.0 && d_pT_phi<3.5) h2_phi_v2_vs_invM_pTbin7 -> Fill(d_inv_m,d_v2_raw_phi,d_eff_corr0*d_eff_corr1);


        if(d_inv_m > (d_mean + 3*d_sigma) || d_inv_m < (d_mean - 3*d_sigma)) continue;
        h_phi_y->Fill(d_phi_y);
        hist_pt_y_Phi->Fill(d_phi_y,d_pT_phi);



      }
    }
    //=================== END Invariant Mass Nested Primary Track Loop ===========================

    // END Phi meson flow analysis


  }
  //////////////////////// END Event Loop //////////////////////////////////////

  // Compute BBC ADCs mean
  for(Int_t i=0;i<16;i++) {
      if(esum[i] > 0.0) {
          emean[i] /= esum[i];
          bbc_east_adc_profile->SetBinContent(i+1,emean[i]);
      }
      if(esum_c[i] > 0.0) {
          emean_c[i] /= esum_c[i];
          bbc_east_gain_corrected_adc_profile->SetBinContent(i+1,emean_c[i]);
      }
      if(wsum[i] > 0.0) {
          wmean[i] /= wsum[i];
          bbc_west_adc_profile->SetBinContent(i+1,wmean[i]);
      }
      if(wsum_c[i] > 0.0) {
          wmean_c[i] /= wsum_c[i];
          bbc_west_gain_corrected_adc_profile->SetBinContent(i+1,wmean_c[i]);
      }
  }

  h_evt_vs_cut->SetBinContent(1,ievtcut[0]);
  h_evt_vs_cut->SetBinContent(2,ievtcut[1]);
  h_evt_vs_cut->SetBinContent(3,ievtcut[2]);
  // h_evt_vs_cut->SetBinContent(4,ievtcut[3]);
  // h_evt_vs_cut->SetBinContent(5,ievtcut[4]);

  h_trk_vs_cut->SetBinContent(1,itrkcut[0]);
  h_trk_vs_cut->SetBinContent(2,itrkcut[1]);
  h_trk_vs_cut->SetBinContent(3,itrkcut[2]);
  h_trk_vs_cut->SetBinContent(4,itrkcut[3]);

  h_evt_vs_cut->GetXaxis()->SetBinLabel(1,"No Cut");
  h_evt_vs_cut->GetXaxis()->SetBinLabel(2,"ZVtx & trigger");
  h_evt_vs_cut->GetXaxis()->SetBinLabel(3,"final");
  // h_evt_vs_cut->GetXaxis()->SetBinLabel(4,"Primary track");
  // h_evt_vs_cut->GetXaxis()->SetBinLabel(5,"bad track");

  h_trk_vs_cut->GetXaxis()->SetBinLabel(1,"No Cut");
  h_trk_vs_cut->GetXaxis()->SetBinLabel(2,"picoTrack cut");
  h_trk_vs_cut->GetXaxis()->SetBinLabel(3,"Primary Track cut");
  h_trk_vs_cut->GetXaxis()->SetBinLabel(4,"bad track cut");

  outputFile->cd();
  h_phi_y->Write();
  hist_mass->Write();
  hist_dEdx->Write();
  hist_trackmult_proton->Write();
  hist_trackmult_pionPlus->Write();
  hist_trackmult_pionMinus->Write();
  hist_trackmult_kaonPlus->Write();
  hist_trackmult_kaonMinus->Write();

  hist_ratio->Write();
  hist_nHits->Write();
  hist_ndEdx->Write();
  hist_DCA->Write();

  h_evt_vs_cut->Write();
  h_trk_vs_cut->Write();

  hist_runId->Write();
  hist_triggerID->Write();
  hist_Vz->Write();
  hist_Vz_cut->Write();
  hist_Vr->Write();
  hist_Vr_cut->Write();
  hist_VyVx->Write();
  hist_VyVx_cut->Write();

  hist_realTrackMult->Write();
  hist_cent->Write();
  hist_realTrackMult_refmult->Write();
  hist_realTrackMult_grefmult->Write();
  // hist_realTrackMult_tofmult->Write();
  hist_cent->Write();
  hist_pt->Write();
  hist_eta->Write();
  hist_phi->Write();

  hist_pt_proton->Write();
  hist_eta_proton->Write();
  hist_y_proton->Write();
  hist_rap_eta_proton->Write();
  hist_pt_y_proton->Write();
  hist_pt_eta_proton->Write();
  hist_phi_proton->Write();
  hist_dEdx_proton->Write();
  hist_beta_proton->Write();
  hist_mass_proton->Write();

  hist_pt_pionPlus->Write();
  hist_eta_pionPlus->Write();
  hist_y_pionPlus->Write();
  hist_rap_eta_pionPlus->Write();
  hist_pt_y_pionPlus->Write();
  hist_pt_eta_pionPlus->Write();
  hist_phi_pionPlus->Write();
  hist_dEdx_pionPlus->Write();
  hist_beta_pionPlus->Write();
  hist_mass_pionPlus->Write();

  hist_pt_pionMinus->Write();
  hist_eta_pionMinus->Write();
  hist_y_pionMinus->Write();
  hist_rap_eta_pionMinus->Write();
  hist_pt_y_pionMinus->Write();
  hist_pt_eta_pionMinus->Write();
  hist_phi_pionMinus->Write();
  hist_dEdx_pionMinus->Write();
  hist_beta_pionMinus->Write();
  hist_mass_pionMinus->Write();

  hist_pt_kaonPlus->Write();
  hist_eta_kaonPlus->Write();
  hist_y_kaonPlus->Write();
  hist_rap_eta_kaonPlus->Write();
  hist_pt_y_kaonPlus->Write();
  hist_pt_eta_kaonPlus->Write();
  hist_phi_kaonPlus->Write();
  hist_dEdx_kaonPlus->Write();
  hist_beta_kaonPlus->Write();
  hist_mass_kaonPlus->Write();

  hist_pt_kaonMinus->Write();
  hist_eta_kaonMinus->Write();
  hist_y_kaonMinus->Write();
  hist_rap_eta_kaonMinus->Write();
  hist_pt_y_kaonMinus->Write();
  hist_pt_eta_kaonMinus->Write();
  hist_phi_kaonMinus->Write();
  hist_dEdx_kaonMinus->Write();
  hist_beta_kaonMinus->Write();
  hist_mass_kaonMinus->Write();

  hist2_Epd_east_Qy_Qx_raw->Write();
  hist2_Epd_west_Qy_Qx_raw->Write();
  hist2_Epd_east_Qy_Qx_Weighted->Write();
  hist2_Epd_west_Qy_Qx_Weighted->Write();

  hist_Epd_east_psi_raw->Write();
  hist_Epd_west_psi_raw->Write();
  hist_Epd_east_psi_Weighted->Write();
  hist_Epd_west_psi_Weighted->Write();
  hist_Epd_east_psi_Shifted->Write();
  hist_Epd_west_psi_Shifted->Write();

  profile_correlation_tpc_east_epd_east->Write();
  correlation2D_epd_east_tpc_east->Write();

  profile_correlation_tpc_west_epd_east->Write();
  correlation2D_epd_east_tpc_west->Write();

  profile_correlation_tpc_east_tpc_west->Write();
  correlation2D_tpc_sub->Write();

  profile_correlation_tpc_east_bbc_east->Write();
  correlation2D_bbc_east_tpc_east->Write();

  profile_correlation_tpc_west_bbc_east->Write();
  correlation2D_bbc_east_tpc_west->Write();

  hist_tpc_east_psi_raw->Write();
  hist_tpc_east_psi_recentered->Write();
  hist_tpc_east_psi_flattened->Write();

  hist_tpc_west_psi_raw->Write();
  hist_tpc_west_psi_recentered->Write();
  hist_tpc_west_psi_flattened->Write();

  hist_bbc_east_psi_raw->Write();
  hist_bbc_east_psi_recentered->Write();
  hist_bbc_east_psi_flattened->Write();

  profile3D_tpc_east_Qx_Qy->Write();
  profile3D_tpc_east_psiShift->Write();

  profile3D_tpc_west_Qx_Qy->Write();
  profile3D_tpc_west_psiShift->Write();

  profile3D_bbc_east_Qx_Qy->Write();
  profile3D_bbc_east_psiShift->Write();

  bbc_east_adc_profile->Write();
  bbc_east_gain_corrected_adc_profile->Write();

  hist_bbc_west_psi_raw->Write();
  hist_bbc_west_psi_recentered->Write();
  hist_bbc_west_psi_flattened->Write();
  profile3D_bbc_west_Qx_Qy->Write();
  profile3D_bbc_west_psiShift->Write();
  bbc_west_adc_profile->Write();
  bbc_west_gain_corrected_adc_profile->Write();

  h2_proton_v2->Write();
  h2_pionPlus_v2->Write();
  h2_pionMinus_v2->Write();
  h2_kaonPlus_v2->Write();
  h2_kaonMinus_v2->Write();
  h2_pions_v2->Write();

  profile_proton_v2    = h2_proton_v2->ProfileX();
  profile_pionPlus_v2  = h2_pionPlus_v2->ProfileX();
  profile_pionMinus_v2 = h2_pionMinus_v2->ProfileX();
  profile_kaonPlus_v2  = h2_kaonPlus_v2->ProfileX();
  profile_kaonMinus_v2 = h2_kaonMinus_v2->ProfileX();
  profile_pions_v2     = h2_pions_v2->ProfileX();

  profile_proton_v2->Write();
  profile_pionPlus_v2->Write();
  profile_pionMinus_v2->Write();
  profile_kaonPlus_v2->Write();
  profile_kaonMinus_v2->Write();
  profile_pions_v2->Write();

  profile3D_proton_v2->Write();
  profile3D_pionPlus_v2->Write();
  profile3D_pionMinus_v2->Write();
  profile3D_kaonPlus_v2->Write();
  profile3D_kaonMinus_v2->Write();

  h_K_DCA_r      ->Write();
  h_K_obj_DCA_r  ->Write();
  h_K_diff_DCA_r ->Write();

  hist_pt_y_Phi->Write();
  h_dip_angle->Write();
  h_prim_inv_m_PHI ->Write();
  h2_phi_v2_vs_invM->Write();
  h2_phi_v2_vs_invM_bin2->Write();
  h2_phi_v2_vs_invM_bin3->Write();
  h2_phi_v2_vs_invM_bin4->Write();

  TP_phi_v2_vs_invM      = h2_phi_v2_vs_invM->ProfileX();
  TP_phi_v2_vs_invM_bin2 = h2_phi_v2_vs_invM_bin2->ProfileX();
  TP_phi_v2_vs_invM_bin3 = h2_phi_v2_vs_invM_bin3->ProfileX();
  TP_phi_v2_vs_invM_bin4 = h2_phi_v2_vs_invM_bin4->ProfileX();

  TP_phi_v2_vs_invM     ->Write();
  TP_phi_v2_vs_invM_bin2->Write();
  TP_phi_v2_vs_invM_bin3->Write();
  TP_phi_v2_vs_invM_bin4->Write();

  h2_phi_v2_vs_invM_pTbin1->Write();
  h2_phi_v2_vs_invM_pTbin2->Write();
  h2_phi_v2_vs_invM_pTbin3->Write();
  h2_phi_v2_vs_invM_pTbin4->Write();
  h2_phi_v2_vs_invM_pTbin5->Write();
  h2_phi_v2_vs_invM_pTbin6->Write();
  h2_phi_v2_vs_invM_pTbin7->Write();

  TP_phi_v2_vs_invM_pTbin1 = h2_phi_v2_vs_invM_pTbin1->ProfileX();
  TP_phi_v2_vs_invM_pTbin2 = h2_phi_v2_vs_invM_pTbin2->ProfileX();
  TP_phi_v2_vs_invM_pTbin3 = h2_phi_v2_vs_invM_pTbin3->ProfileX();
  TP_phi_v2_vs_invM_pTbin4 = h2_phi_v2_vs_invM_pTbin4->ProfileX();
  TP_phi_v2_vs_invM_pTbin5 = h2_phi_v2_vs_invM_pTbin5->ProfileX();
  TP_phi_v2_vs_invM_pTbin6 = h2_phi_v2_vs_invM_pTbin6->ProfileX();
  TP_phi_v2_vs_invM_pTbin7 = h2_phi_v2_vs_invM_pTbin7->ProfileX();

  TP_phi_v2_vs_invM_pTbin1->Write();
  TP_phi_v2_vs_invM_pTbin2->Write();
  TP_phi_v2_vs_invM_pTbin3->Write();
  TP_phi_v2_vs_invM_pTbin4->Write();
  TP_phi_v2_vs_invM_pTbin5->Write();
  TP_phi_v2_vs_invM_pTbin6->Write();
  TP_phi_v2_vs_invM_pTbin7->Write();

  mEpFinder->Finish();

}
//////////////////////////// END Main Function /////////////////////////////////
