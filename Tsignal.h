//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 13:53:10 2014 by ROOT version 5.34/10
// from TTree T/event
// found on file: ../MG2Dalone_200fC_shaping63_gasstarted_140818_11H23_025_signal.root
//////////////////////////////////////////////////////////

#ifndef Tsignal_h
#define Tsignal_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

using std::vector;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Tsignal {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Nevent;
   Int_t           TsampleNum[32];
   int CM_n;
   int MG_n;
   Float_t         (*StripAmpl_MG)[61][32];
   Float_t         (*StripAmpl_MG_ped)[61][32];
   Float_t         (*StripAmpl_MG_corr)[61][32];
   Float_t         (*StripAmpl_CM)[64][32];
   Float_t         (*StripAmpl_CM_ped)[64][32];
   Float_t         (*StripAmpl_CM_corr)[64][32];

   // List of branches
   TBranch        *b_Nevent;   //!
   TBranch        *b_TsampleNum;   //!
   TBranch        *b_StripAmpl_MG;   //!
   TBranch        *b_StripAmpl_MG_ped;   //!
   TBranch        *b_StripAmpl_MG_corr;   //!
   TBranch        *b_StripAmpl_CM;   //!
   TBranch        *b_StripAmpl_CM_ped;   //!
   TBranch        *b_StripAmpl_CM_corr;   //!

   Tsignal(TTree *tree, int CMN, int MGN);
   Tsignal();
   virtual ~Tsignal();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, int CMN, int MGN);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   vector<vector<double> > get_mg_ampl(int mg_n);
   vector<vector<double> > get_cm_ampl(int cm_n);

};

#endif
