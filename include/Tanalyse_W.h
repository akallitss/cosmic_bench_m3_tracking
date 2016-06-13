#ifndef Tanalyse_W_h
#define Tanalyse_W_h

#include <TTree.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <map>

#include "tomography.h"
#include "detector.h"

using std::string;
using std::vector;

class Event;

class Tanalyse_W{
   public:
      TTree          *T;
      TFile          *saveFile;

      // Declaration of leaf types
      map<Tomography::det_type,unsigned short> det_N;

      Int_t           evn;
      Double_t        evttime;
      /*
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
      */
      map<Tomography::det_type,Int_t*> NClus;
      map<Tomography::det_type,Int_t*> Spark;
      map<Tomography::det_type,Double_t*> ClusAmpl;
      map<Tomography::det_type,Double_t*> ClusSize;
      map<Tomography::det_type,Double_t*> ClusPos;
      map<Tomography::det_type,Double_t*> ClusMaxStripAmpl;
      map<Tomography::det_type,Int_t*> ClusMaxStrip;
      map<Tomography::det_type,Double_t*> ClusMaxSample;
      map<Tomography::det_type,Double_t*> ClusTOT;
      map<Tomography::det_type,Double_t*> ClusT;
      map<Tomography::det_type,Double_t*> StripMaxAmpl;

      //Tanalyse(string saveFileName);
      Tanalyse_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_);
      ~Tanalyse_W();
      void Init();
      TTree * getTree() const;
      void Write();
      void CloseFile();
      void fillTree(int evn_, double evttime_, map<Tomography::det_type,vector<Event*> > events);
};

#endif

