//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 13:53:10 2014 by ROOT version 5.34/10
// from TTree T/event
// found on file: ../MG2Dalone_200fC_shaping63_gasstarted_140818_11H23_025_signal.root
//////////////////////////////////////////////////////////

#ifndef Tsignal_R_h
#define Tsignal_R_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <map>
#include "tomography.h"
#include "detector.h"

using std::vector;
using std::map;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Tsignal_R {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   map<Tomography::det_type,unsigned short> det_N;
   Int_t           Nevent;
   Double_t        evttime;
   Float_t         (*StripAmpl_MG)[MG_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_MG_ped)[MG_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_MG_corr)[MG_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_MGv2)[MGv2_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_MGv2_ped)[MGv2_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_MGv2_corr)[MGv2_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_CM)[CM_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_CM_ped)[CM_Detector::Nchannel][Tomography::Nsample];
   Float_t         (*StripAmpl_CM_corr)[CM_Detector::Nchannel][Tomography::Nsample];

   // List of branches
   TBranch        *b_Nevent;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_StripAmpl_MG;   //!
   TBranch        *b_StripAmpl_MG_ped;   //!
   TBranch        *b_StripAmpl_MG_corr;   //!
   TBranch        *b_StripAmpl_MGv2;   //!
   TBranch        *b_StripAmpl_MGv2_ped;   //!
   TBranch        *b_StripAmpl_MGv2_corr;   //!
   TBranch        *b_StripAmpl_CM;   //!
   TBranch        *b_StripAmpl_CM_ped;   //!
   TBranch        *b_StripAmpl_CM_corr;   //!

   Tsignal_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   Tsignal_R();
   virtual ~Tsignal_R();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   vector<vector<double> > get_ampl(Tomography::det_type type_, unsigned short det_n_);

};

#endif
