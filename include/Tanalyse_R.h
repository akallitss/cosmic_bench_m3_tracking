//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 18 13:34:43 2014 by ROOT version 5.18/00b
// from TTree T/event
// found on file: run6_Pb_analyse.root
//////////////////////////////////////////////////////////

#ifndef Tanalyse_R_h
#define Tanalyse_R_h

#include <TChain.h>
#include <TFile.h>

#include "tomography.h"
#include "detector.h"
#include <map>

using std::map;

class CM_Detector;
class MG_Detector;
class MGv2_Detector;

class Tanalyse_R{
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   map<Tomography::det_type,unsigned short> det_N;
   Int_t           evn;
   Double_t        evttime;
   Int_t           *CM_NClus;
   Int_t           *CM_Spark;
   Double_t        (*CM_ClusAmpl)[600];
   Double_t        (*CM_ClusSize)[600];
   Double_t        (*CM_ClusPos)[600];
   Double_t        (*CM_ClusMaxStripAmpl)[600];
   Int_t           (*CM_ClusMaxStrip)[600];
   Double_t        (*CM_ClusMaxSample)[600];
   Double_t        (*CM_ClusTOT)[600];
   Double_t        (*CM_ClusT)[600];
   Double_t        (*CM_StripMaxAmpl)[CM_Detector::Nchannel/2];
   Int_t           *MG_NClus;
   Int_t           *MG_Spark;
   Double_t        (*MG_ClusAmpl)[300];
   Double_t        (*MG_ClusSize)[300];
   Double_t        (*MG_ClusPos)[300];
   Double_t        (*MG_ClusMaxStripAmpl)[300];
   Double_t        (*MG_ClusMaxSample)[300];
   Double_t        (*MG_ClusTOT)[300];
   Double_t        (*MG_ClusT)[300];
   Int_t           (*MG_ClusMaxStrip)[300];
   Double_t        (*MG_StripMaxAmpl)[MG_Detector::Nchannel];
   Int_t           *MGv2_NClus;
   Int_t           *MGv2_Spark;
   Double_t        (*MGv2_ClusAmpl)[300];
   Double_t        (*MGv2_ClusSize)[300];
   Double_t        (*MGv2_ClusPos)[300];
   Double_t        (*MGv2_ClusMaxStripAmpl)[300];
   Double_t        (*MGv2_ClusMaxSample)[300];
   Double_t        (*MGv2_ClusTOT)[300];
   Double_t        (*MGv2_ClusT)[300];
   Int_t           (*MGv2_ClusMaxStrip)[300];
   Double_t        (*MGv2_StripMaxAmpl)[MGv2_Detector::Nchannel];

   // List of branches
   TBranch        *b_evn;   //!
   TBranch        *b_evttime;
   TBranch        *b_CM_NClus;   //!
   TBranch        *b_CM_Spark;
   TBranch        *b_CM_ClusAmpl;   //!
   TBranch        *b_CM_ClusSize;   //!
   TBranch        *b_CM_ClusPos;   //!
   TBranch        *b_CM_ClusMaxStripAmpl;   //!
   TBranch        *b_CM_ClusMaxStrip;   //!
   TBranch        *b_CM_ClusMaxSample;   //!
   TBranch        *b_CM_ClusTOT;   //!
   TBranch        *b_CM_ClusT;   //!
   TBranch        *b_CM_StripMaxAmpl;   //!
   TBranch        *b_MG_NClus;   //!
   TBranch        *b_MG_Spark;
   TBranch        *b_MG_ClusAmpl;   //!
   TBranch        *b_MG_ClusSize;   //!
   TBranch        *b_MG_ClusPos;   //!
   TBranch        *b_MG_ClusMaxStripAmpl;   //!
   TBranch        *b_MG_ClusMaxSample;   //!
   TBranch        *b_MG_ClusTOT;   //!
   TBranch        *b_MG_ClusT;   //!
   TBranch        *b_MG_ClusMaxStrip;   //!
   TBranch        *b_MG_StripMaxAmpl;   //!
   TBranch        *b_MGv2_NClus;   //!
   TBranch        *b_MGv2_Spark;
   TBranch        *b_MGv2_ClusAmpl;   //!
   TBranch        *b_MGv2_ClusSize;   //!
   TBranch        *b_MGv2_ClusPos;   //!
   TBranch        *b_MGv2_ClusMaxStripAmpl;   //!
   TBranch        *b_MGv2_ClusMaxSample;   //!
   TBranch        *b_MGv2_ClusTOT;   //!
   TBranch        *b_MGv2_ClusT;   //!
   TBranch        *b_MGv2_ClusMaxStrip;   //!
   TBranch        *b_MGv2_StripMaxAmpl;   //!

   Tanalyse_R();
   Tanalyse_R(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual ~Tanalyse_R();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, map<Tomography::det_type,unsigned short> det_N_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

