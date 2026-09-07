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
      vector<double>  NClusX;
      vector<double>  NClusY;
      // Per-layer measured cluster positions (mm, aligned) for each ray, NaN when the
      // layer has no cluster in that coordinate.  Layers 0..3 map to fixed station z's
      // (read from the run config on the analysis side).  Added to enable a per-plane
      // telescope self-resolution / DUT-pointing study (does not change the fit).
      vector<double>  mX0, mX1, mX2, mX3;
      vector<double>  mY0, mY1, mY2, mY3;
      // z (mm, aligned) of the X/Y plane actually used at each layer, NaN when absent.
      vector<double>  zX0, zX1, zX2, zX3;
      vector<double>  zY0, zY1, zY2, zY3;

      Tray(string saveFileName);
      ~Tray();
      void Init();
      TTree * getTree() const;
      void Write();
      void CloseFile();
      void fillTree(int evn_, double evttime_, vector<Ray> rays, double Z_Up_, double Z_Down_);
};

#endif