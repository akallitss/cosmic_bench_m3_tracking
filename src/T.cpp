#define T_cxx
#include "T.h"
#include <iostream>

using std::cout;
using std::endl;

T::T(){
   
}

T::T(TTree *tree, int CM_n, int MG_n)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/* if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run6_Pb_analyse.root");
      if (!f) {
         f = new TFile("run6_Pb_analyse.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }*/
   Init(tree, CM_n, MG_n);
}

T::~T()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t T::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t T::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void T::Init(TTree *tree, int CM_n, int MG_n)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evn", &evn, &b_evn);
   evttime = 0;
   //fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   if(CM_n>0){
      CM_NClus = new int[CM_n];
      fChain->SetBranchAddress("CM_NClus", CM_NClus, &b_CM_NClus);
      CM_Spark = new int[CM_n];
      for(int i = 0;i<CM_n;i++){
         CM_Spark[i] = 0;
      }
      //fChain->SetBranchAddress("CM_Spark", CM_Spark, &b_CM_Spark);
      CM_ClusAmpl = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusAmpl", CM_ClusAmpl, &b_CM_ClusAmpl);
      CM_ClusSize = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusSize", CM_ClusSize, &b_CM_ClusSize);
      CM_ClusPos = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusPos", CM_ClusPos, &b_CM_ClusPos);
      CM_ClusMaxStripAmpl = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusMaxStripAmpl", CM_ClusMaxStripAmpl, &b_CM_ClusMaxStripAmpl);
      CM_ClusMaxStrip = new Int_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusMaxStrip", CM_ClusMaxStrip, &b_CM_ClusMaxStrip);
      CM_ClusMaxSample = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusMaxSample", CM_ClusMaxSample, &b_CM_ClusMaxSample);
      CM_ClusTOT = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusTOT", CM_ClusTOT, &b_CM_ClusTOT);
      CM_ClusT = new Double_t[CM_n][600];
      fChain->SetBranchAddress("CM_ClusT", CM_ClusT, &b_CM_ClusT);
      CM_StripMaxAmpl = new Double_t[CM_n][32];
      fChain->SetBranchAddress("CM_StripMaxAmpl", CM_StripMaxAmpl, &b_CM_StripMaxAmpl);
   }
   if(MG_n>0){
      MG_NClus = new int[MG_n];
      fChain->SetBranchAddress("MG_NClus", MG_NClus, &b_MG_NClus);
      MG_Spark = new int[MG_n];
      for(int i = 0;i<MG_n;i++){
         MG_Spark[i] = 0;
      }
      //fChain->SetBranchAddress("MG_Spark", MG_Spark, &b_MG_Spark);
      MG_ClusAmpl = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusAmpl", MG_ClusAmpl, &b_MG_ClusAmpl);
      MG_ClusSize = new Double_t[MG_n][300];
      MG_ClusPos = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusPos", MG_ClusPos, &b_MG_ClusPos);
      fChain->SetBranchAddress("MG_ClusSize", MG_ClusSize, &b_MG_ClusSize);
      MG_ClusMaxStripAmpl = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusMaxStripAmpl", MG_ClusMaxStripAmpl, &b_MG_ClusMaxStripAmpl);
      MG_ClusMaxSample = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusMaxSample", MG_ClusMaxSample, &b_MG_ClusMaxSample);
      MG_ClusTOT = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusTOT", MG_ClusTOT, &b_MG_ClusTOT);
      MG_ClusT = new Double_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusT", MG_ClusT, &b_MG_ClusT);
      MG_ClusMaxStrip = new Int_t[MG_n][300];
      fChain->SetBranchAddress("MG_ClusMaxStrip", MG_ClusMaxStrip, &b_MG_ClusMaxStrip);
      MG_StripMaxAmpl = new Double_t[MG_n][61];
      fChain->SetBranchAddress("MG_StripMaxAmpl", MG_StripMaxAmpl, &b_MG_StripMaxAmpl);
   }
   Notify();
}

Bool_t T::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void T::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}