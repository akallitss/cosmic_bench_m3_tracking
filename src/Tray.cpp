#define Tray_cxx
#include "Tray.h"
#include "cluster.h"
#include <iostream>
#include <cmath>

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
   T->Branch("NClusX", "vector<double>", &NClusX);
   T->Branch("NClusY", "vector<double>", &NClusY);
   T->Branch("mX0", "vector<double>", &mX0);
   T->Branch("mX1", "vector<double>", &mX1);
   T->Branch("mX2", "vector<double>", &mX2);
   T->Branch("mX3", "vector<double>", &mX3);
   T->Branch("mY0", "vector<double>", &mY0);
   T->Branch("mY1", "vector<double>", &mY1);
   T->Branch("mY2", "vector<double>", &mY2);
   T->Branch("mY3", "vector<double>", &mY3);
   T->Branch("zX0", "vector<double>", &zX0);
   T->Branch("zX1", "vector<double>", &zX1);
   T->Branch("zX2", "vector<double>", &zX2);
   T->Branch("zX3", "vector<double>", &zX3);
   T->Branch("zY0", "vector<double>", &zY0);
   T->Branch("zY1", "vector<double>", &zY1);
   T->Branch("zY2", "vector<double>", &zY2);
   T->Branch("zY3", "vector<double>", &zY3);
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
   NClusX = vector<double>(rayN,0);
   NClusY = vector<double>(rayN,0);
   const double NaN = std::numeric_limits<double>::quiet_NaN();
   vector<vector<double>*> mX = {&mX0,&mX1,&mX2,&mX3};
   vector<vector<double>*> mY = {&mY0,&mY1,&mY2,&mY3};
   vector<vector<double>*> zX = {&zX0,&zX1,&zX2,&zX3};
   vector<vector<double>*> zY = {&zY0,&zY1,&zY2,&zY3};
   for(int L=0;L<4;L++){ mX[L]->assign(rayN,NaN); mY[L]->assign(rayN,NaN);
                         zX[L]->assign(rayN,NaN); zY[L]->assign(rayN,NaN); }
   for(int i=0;i<rayN;i++){
      Chi2X[i] = rays[i].get_chiSquare_X();
      Chi2Y[i] = rays[i].get_chiSquare_Y();
      NClusX[i] = rays[i].get_clus_x_n();
      NClusY[i] = rays[i].get_clus_y_n();
      X_Up[i] = rays[i].eval_X(Z_Up);
      Y_Up[i] = rays[i].eval_Y(Z_Up);
      X_Down[i] = rays[i].eval_X(Z_Down);
      Y_Down[i] = rays[i].eval_Y(Z_Down);
      // Per-layer measured positions used in this ray's fit (get_clus returns clones).
      vector<Cluster*> clus = rays[i].get_clus();
      for(vector<Cluster*>::iterator it=clus.begin();it!=clus.end();++it){
         int L = (*it)->get_layer();
         if(L<0 || L>3) continue;
         if((*it)->get_is_X()){ (*mX[L])[i] = (*it)->get_pos_mm(); (*zX[L])[i] = (*it)->get_z(); }
         else                 { (*mY[L])[i] = (*it)->get_pos_mm(); (*zY[L])[i] = (*it)->get_z(); }
         delete *it;  // get_clus() hands back Clones; free them
      }
   }
   T->Fill();
}