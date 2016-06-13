//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 15 13:53:10 2014 by ROOT version 5.34/10
// from TTree T/event
// found on file: ../MG2Dalone_200fC_shaping63_gasstarted_140818_11H23_025_signal.root
//////////////////////////////////////////////////////////

#ifndef Tsignal_W_h
#define Tsignal_W_h

#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <map>

#include "tomography.h"
#include "detector.h"

using std::vector;
using std::map;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Tsignal_W {
public :
   TTree          *T;
   TFile          *saveFile;

   // Declaration of leaf types
   map<Tomography::det_type,unsigned short> det_N;
   Int_t           Nevent;
   Double_t        evttime;
   /*
   Float_t         *StripAmpl_MG;
   Float_t         *StripAmpl_MG_ped;
   Float_t         *StripAmpl_MG_corr;
   Float_t         *StripAmpl_MGv2;
   Float_t         *StripAmpl_MGv2_ped;
   Float_t         *StripAmpl_MGv2_corr;
   Float_t         *StripAmpl_CM;
   Float_t         *StripAmpl_CM_ped;
   Float_t         *StripAmpl_CM_corr;
   */
   map<Tomography::det_type,Float_t*> StripAmpl;
   map<Tomography::det_type,Float_t*> StripAmpl_ped;
   map<Tomography::det_type,Float_t*> StripAmpl_corr;


   Tsignal_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_);
   ~Tsignal_W();
   void Init();
   TTree * getTree() const;
   void Write();
   void CloseFile();
   void Reset_raw();
   void Reset_ped();
   void Reset_corr();
   void disable_data_branches();
   void enable_all_branches();
   void enable_raw_branches();
   void enable_ped_branches();
   void enable_corr_branches();
   template<typename A> void fillTree_raw(int evn_, double evttime_, map<Tomography::det_type,vector<vector<vector<A> > > > ampl);
   template<typename A> void fillTree_ped(map<Tomography::det_type,vector<vector<vector<A> > > > ampl);
   template<typename A> void fillTree_corr(map<Tomography::det_type,vector<vector<vector<A> > > > ampl);
   template<typename A> map<Tomography::det_type,vector<vector<vector<A> > > > read_raw(long entry);
   template<typename A> map<Tomography::det_type,vector<vector<vector<A> > > > read_ped(long entry);
   template<typename A> map<Tomography::det_type,vector<vector<vector<A> > > > read_corr(long entry);
};

#endif
