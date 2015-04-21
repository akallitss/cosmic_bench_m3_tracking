#define Tray_cxx
#include "Tray.h"
#include <iostream>

using std::cout;
using std::endl;

Tray::Tray(string saveFileName)
{
   saveFile = new TFile(saveFileName.c_str(),"RECREATE");
   T = new TTree("T","event");
   T->SetMaxTreeSize(10000000000000LL);
   Init();
}

Tray::~Tray()
{
   //delete T;
   delete saveFile;
}

void Tray::Init()
{
   T->Branch("evn", &evn, "evn/I");
   T->Branch("evttime", &evttime, "evttime/D");
   T->Branch("rayN",&rayN,"rayN/I");
   T->Branch("Z_Up",&Z_Up,"Z_Up/D");
   T->Branch("Z_Down",&Z_Down,"Z_Down/D");
   T->Branch("X_Up", "vector<double>", &X_Up);
   T->Branch("Y_Up", "vector<double>", &Y_Up);
   T->Branch("X_Down", "vector<double>", &X_Down);
   T->Branch("Y_Down", "vector<double>", &Y_Down);
   T->Branch("Chi2X", "vector<double>", &Chi2X);
   T->Branch("Chi2Y", "vector<double>", &Chi2Y);
}
TTree * Tray::getTree() const{
   return T->CloneTree();
}
void Tray::Write(){
   saveFile->cd();
   T->Write();
}
void Tray::CloseFile(){
   saveFile->Close();
}
void Tray::fillTree(int evn_, double evttime_, vector<Ray> rays, double Z_Up_, double Z_Down_){
   evn = evn_;
   evttime = evttime_;
   Z_Up = Z_Up_;
   Z_Down = Z_Down_;
   rayN = rays.size();
   X_Up = vector<double>(rayN,0);
   X_Down = vector<double>(rayN,0);
   Y_Up = vector<double>(rayN,0);
   Y_Down = vector<double>(rayN,0);
   Chi2X = vector<double>(rayN,-1);
   Chi2Y = vector<double>(rayN,-1);
   for(int i=0;i<rayN;i++){
      Chi2X[i] = rays[i].get_chiSquare_X();
      Chi2Y[i] = rays[i].get_chiSquare_Y();
      X_Up[i] = rays[i].eval_X(Z_Up);
      Y_Up[i] = rays[i].eval_Y(Z_Up);
      X_Down[i] = rays[i].eval_X(Z_Down);
      Y_Down[i] = rays[i].eval_Y(Z_Down);
   }
   T->Fill();
}