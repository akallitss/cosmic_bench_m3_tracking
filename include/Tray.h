#ifndef Trays_h
#define Trays_h

#include <TTree.h>
#include <TFile.h>
#include "ray.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

class Tray{
   public:
      TTree          *T;
      TFile          *saveFile;

      // Declaration of leaf types

      Int_t           evn;
      Double_t        evttime;
      Double_t        Z_Up;
      Double_t        Z_Down;
      Int_t           rayN;
      vector<double>  X_Up;
      vector<double>  Y_Up;
      vector<double>  X_Down;
      vector<double>  Y_Down;
      vector<double>  Chi2X;
      vector<double>  Chi2Y;

      Tray(string saveFileName);
      ~Tray();
      void Init();
      TTree * getTree() const;
      void Write();
      void CloseFile();
      void fillTree(int evn_, double evttime_, vector<Ray> rays, double Z_Up_, double Z_Down_);
};

#endif