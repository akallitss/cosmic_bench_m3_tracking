#ifndef Tanalyse_h
#define Tanalyse_h

//#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include "event.h"
#include <string>
#include <map>

using std::string;
using std::map;

class MG_Event;
class CM_Event;

class Tanalyse{
   public:
      TTree          *T;
      TFile          *saveFile;

      // Declaration of leaf types
      int MGN;
      int CMN;

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
      Double_t        (*CM_StripMaxAmpl)[32];
      Int_t           *MG_NClus;
      Int_t           *MG_Spark;
      Double_t        (*MG_ClusAmpl)[300];
      Double_t        (*MG_ClusSize)[300];
      Double_t        (*MG_ClusPos)[300];
      Double_t        (*MG_ClusMaxStripAmpl)[300];
      Double_t        (*MG_ClusMaxSample)[300];
      Double_t        (*MG_ClusTOT)[300];
      Double_t        (*MG_ClusT)[300];
      Double_t        (*MG_StripMaxAmpl)[61];

      Tanalyse(string saveFileName);
      Tanalyse(string saveFileName, int CM_n, int MG_n);
      ~Tanalyse();
      void Init();
      TTree * getTree() const;
      void Write();
      void CloseFile();
      void fillTree(int evn_, double evttime_, map<int,MG_Event> mg_events, map<int,CM_Event> cm_events);
};

#endif

