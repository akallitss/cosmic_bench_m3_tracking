#define Tsignal_cxx
#include "Tsignal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>

using std::vector;

#include <iostream>
using std::cout;
using std::endl;

Tsignal::Tsignal(){
   
}

Tsignal::Tsignal(TTree *tree, int CMN_, int MGN_) : fChain(0) 
{
   Init(tree, CMN_, MGN_);
}

Tsignal::~Tsignal()
{
   if(MGN>0){
      delete StripAmpl_MG;
      delete StripAmpl_MG_ped;
      delete StripAmpl_MG_corr;
   }
   if(CMN>0){
      delete StripAmpl_CM;
      delete StripAmpl_CM_ped;
      delete StripAmpl_CM_corr;
   }
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tsignal::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tsignal::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Tsignal::Init(TTree *tree, int CMN_, int MGN_)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   MGN = MGN_;
   CMN = CMN_;
   if(MGN>0){
      StripAmpl_MG = new Float_t[MGN][61][Tomography::Nsample];
      StripAmpl_MG_ped = new Float_t[MGN][61][Tomography::Nsample];
      StripAmpl_MG_corr = new Float_t[MGN][61][Tomography::Nsample];
   }
   if(CMN>0){
      StripAmpl_CM = new Float_t[CMN][64][Tomography::Nsample];
      StripAmpl_CM_ped = new Float_t[CMN][64][Tomography::Nsample];
      StripAmpl_CM_corr = new Float_t[CMN][64][Tomography::Nsample];
   }

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("evttime",&evttime, &b_evttime);
   fChain->SetBranchAddress("TsampleNum", TsampleNum, &b_TsampleNum);
   if(MGN>0){
	   fChain->SetBranchAddress("StripAmpl_MG", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_MG_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_MG_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
	if(CMN>0){
	   fChain->SetBranchAddress("StripAmpl_CM", StripAmpl_MG, &b_StripAmpl_MG);
	   fChain->SetBranchAddress("StripAmpl_CM_ped", StripAmpl_MG_ped, &b_StripAmpl_MG_ped);
	   fChain->SetBranchAddress("StripAmpl_CM_corr", StripAmpl_MG_corr, &b_StripAmpl_MG_corr);
	}
   Notify();
}

Bool_t Tsignal::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tsignal::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tsignal::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

vector<vector<double> > Tsignal::get_mg_ampl(int mg_n){
	vector<vector<double> > return_array(61,vector<double>(Tomography::Nsample,0));
	for(int i=0;i<61;i++){
		for(int j=0;j<Tomography::Nsample;j++){
			return_array[i][j] = StripAmpl_MG_corr[mg_n][i][j];
		}
	}
	return return_array;
}
vector<vector<double> > Tsignal::get_cm_ampl(int cm_n){
	vector<vector<double> > return_array(64,vector<double>(Tomography::Nsample,0));
	for(int i=0;i<64;i++){
		for(int j=0;j<Tomography::Nsample;j++){
			return_array[i][j] = StripAmpl_MG_corr[cm_n][i][j];
		}
	}
	return return_array;
}
