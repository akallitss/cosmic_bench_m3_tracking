//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 13:53:10 2014 by ROOT version 5.34/10
// from TTree T/event
// found on file: ../MG2Dalone_200fC_shaping63_gasstarted_140818_11H23_025_signal.root
//////////////////////////////////////////////////////////

#ifndef Tsignal_W_h
#define Tsignal_W_h

//#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <map>
#include "tomography.h"

using std::vector;
using std::map;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Tsignal_W {
public :
   TTree          *T;
   TFile          *saveFile;

   // Declaration of leaf types
   int CMN;
   int MGN;
   Int_t           evn;
   Double_t        evttime;
   Float_t         (*StripAmpl_MG)[61][Tomography::Nsample];
   Float_t         (*StripAmpl_MG_ped)[61][Tomography::Nsample];
   Float_t         (*StripAmpl_MG_corr)[61][Tomography::Nsample];
   Float_t         (*StripAmpl_CM)[64][Tomography::Nsample];
   Float_t         (*StripAmpl_CM_ped)[64][Tomography::Nsample];
   Float_t         (*StripAmpl_CM_corr)[64][Tomography::Nsample];


   Tsignal_W(string saveFileName, int CM_n, int MG_n);
   ~Tsignal_W();
   void Init();
   TTree * getTree() const;
   void Write();
   void CloseFile();
   void Reset_raw();
   void Reset_ped();
   void Reset_corr();
   void fillTree_raw(int evn_, double evttime_, vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl);
   void fillTree_ped(vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl);
   void fillTree_corr(vector<vector<vector<float> > > mg_ampl, vector<vector<vector<float> > > cm_ampl);
   map<Tomography::det_type,vector<vector<vector<float> > > > read_raw(long entry);
   map<Tomography::det_type,vector<vector<vector<float> > > > read_ped(long entry);
   map<Tomography::det_type,vector<vector<vector<float> > > > read_corr(long entry);
};

#endif
