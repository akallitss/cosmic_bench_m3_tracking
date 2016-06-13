#define Tanalyse_W_cxx
#include "Tanalyse_W.h"
#include <iostream>
#include "event.h"
#include "cluster.h"

using std::cout;
using std::endl;
/*
Tanalyse_W::Tanalyse_W(string saveFileName){
   saveFile = new TFile(saveFileName.c_str(),"RECREATE");
}
*/
Tanalyse_W::Tanalyse_W(string saveFileName, map<Tomography::det_type,unsigned short> det_N_)
{
   saveFile = new TFile(saveFileName.c_str(),"RECREATE");
   saveFile->SetCompressionLevel(5);
   T = new TTree("T","event");
   T->SetMaxTreeSize(10000000000000LL);
   det_N = det_N_;
   Init();
}

Tanalyse_W::~Tanalyse_W()
{
   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      delete NClus[type_it->first];
      delete Spark[type_it->first];
      delete ClusAmpl[type_it->first];
      delete ClusSize[type_it->first];
      delete ClusPos[type_it->first];
      delete ClusTOT[type_it->first];
      delete ClusT[type_it->first];
      delete ClusMaxStrip[type_it->first];
      delete ClusMaxSample[type_it->first];
      delete ClusMaxStripAmpl[type_it->first];
      delete StripMaxAmpl[type_it->first];
   }
   //delete T;
   delete saveFile;
}

void Tanalyse_W::Init()
{
   T->Branch("evn", &evn, "evn/I");
   T->Branch("evttime", &evttime, "evttime/D");
   //Tout->Branch("finetstp", &finetstp, "finetstp/I");

   for(map<Tomography::det_type,unsigned short>::iterator type_it=det_N.begin();type_it!=det_N.end();++type_it){
      int current_channel = Tomography::Static_Detector[type_it->first]->get_Nchannel();
      int current_MaxNClus = Tomography::Static_Detector[type_it->first]->get_MaxNClus();
      NClus[type_it->first] = new int[type_it->second];
      Spark[type_it->first] = new int[type_it->second];
      ClusAmpl[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusSize[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusPos[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusMaxStripAmpl[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusMaxStrip[type_it->first] = new Int_t[type_it->second*current_MaxNClus];
      ClusMaxSample[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusTOT[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      ClusT[type_it->first] = new Double_t[type_it->second*current_MaxNClus];
      StripMaxAmpl[type_it->first] = new Double_t[type_it->second*current_channel];
      for(unsigned int i=0;i<(type_it->second);i++){
         Spark[type_it->first][i] = 0;
         for(int j=0;j<current_channel;j++){
             reinterpret_cast<Double_t(*)[type_it->second]>(StripMaxAmpl[type_it->first])[i][j] = 0;
         }
      }

      string current_name = Tomography::Static_Detector[type_it->first]->Name();

      char leef_NClus[100];
      sprintf(leef_NClus,"%s_NClus[%d]/I",(current_name).c_str(),type_it->second);
      T->Branch((current_name+"_NClus").c_str(), NClus[type_it->first], leef_NClus);

      char leef_Spark[100];
      sprintf(leef_Spark,"%s_Spark[%d]/I",(current_name).c_str(),type_it->second);
      T->Branch((current_name+"_Spark").c_str(), Spark[type_it->first], leef_Spark);

      char leef_ClusAmpl[100];
      sprintf(leef_ClusAmpl,"%s_ClusAmpl[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusAmpl").c_str(), ClusAmpl[type_it->first], leef_ClusAmpl);

      char leef_ClusSize[100];
      sprintf(leef_ClusSize,"%s_ClusSize[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusSize").c_str(), ClusSize[type_it->first], leef_ClusSize);

      char leef_ClusPos[100];
      sprintf(leef_ClusPos,"%s_ClusPos[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusPos").c_str(), ClusPos[type_it->first], leef_ClusPos);

      char leef_ClusMaxStripAmpl[100];
      sprintf(leef_ClusMaxStripAmpl,"%s_ClusMaxStripAmpl[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusMaxStripAmpl").c_str(), ClusMaxStripAmpl[type_it->first], leef_ClusMaxStripAmpl);

      char leef_ClusMaxSample[100];
      sprintf(leef_ClusMaxSample,"%s_ClusMaxSample[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusMaxSample").c_str(), ClusMaxSample[type_it->first], leef_ClusMaxSample);

      char leef_ClusTOT[100];
      sprintf(leef_ClusTOT,"%s_ClusTOT[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusTOT").c_str(), ClusTOT[type_it->first], leef_ClusTOT);

      char leef_ClusT[100];
      sprintf(leef_ClusT,"%s_ClusT[%d][%d]/D",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusT").c_str(), ClusT[type_it->first], leef_ClusT);

      char leef_StripMaxAmpl[100];
      sprintf(leef_StripMaxAmpl,"%s_StripMaxAmpl[%d][%d]/D",(current_name).c_str(),type_it->second,current_channel);
      T->Branch((current_name+"_StripMaxAmpl").c_str(), StripMaxAmpl[type_it->first], leef_StripMaxAmpl);

      char leef_ClusMaxStrip[100];
      sprintf(leef_ClusMaxStrip,"%s_ClusMaxStrip[%d][%d]/I",(current_name).c_str(),type_it->second,current_MaxNClus);
      T->Branch((current_name+"_ClusMaxStrip").c_str(), ClusMaxStrip[type_it->first], leef_ClusMaxStrip);

   }
   
}
TTree * Tanalyse_W::getTree() const{
   return T->CloneTree();
}
void Tanalyse_W::Write(){
   saveFile->cd();
   T->Write();
}
void Tanalyse_W::CloseFile(){
   saveFile->Close();
}
void Tanalyse_W::fillTree(int evn_, double evttime_, map<Tomography::det_type,vector<Event*> > events){
   evn = evn_;
   evttime = evttime_;

   for(map<Tomography::det_type,vector<Event*> >::iterator type_it=events.begin();type_it!=events.end();++type_it){
      if(det_N.count(type_it->first)==0){
         cout << "detector type mismatch : " << type_it->first << endl;
         return;
      }
      if((type_it->second).size() != det_N[type_it->first]){
         cout << "problem in event number" << endl;
         return;
      }
      int current_MaxNClus = Tomography::Static_Detector[type_it->first]->get_MaxNClus();
      for(vector<Event*>::iterator event_it=(type_it->second).begin();event_it!=(type_it->second).end();++event_it){
         int i = (*event_it)->get_n_in_tree();
         NClus[type_it->first][i] = (*event_it)->get_NClus();
         vector<Cluster*> current_clusters = (*event_it)->get_clusters();
         if(current_clusters.size()!=static_cast<unsigned int>(NClus[type_it->first][i])){
            cout << "problem in cluster number" << endl;
            return;
         }
         for(int j=0;j<NClus[type_it->first][i];j++){
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusAmpl[type_it->first])[i][j] = current_clusters[j]->get_ampl();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusT[type_it->first])[i][j] = current_clusters[j]->get_t();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusTOT[type_it->first])[i][j] = current_clusters[j]->get_TOT();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusPos[type_it->first])[i][j] = current_clusters[j]->get_pos();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusSize[type_it->first])[i][j] = current_clusters[j]->get_size();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusMaxSample[type_it->first])[i][j] = current_clusters[j]->get_maxSample();
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusMaxStripAmpl[type_it->first])[i][j] = current_clusters[j]->get_maxStripAmpl();
            reinterpret_cast<Int_t(*)[current_MaxNClus]>(ClusMaxStrip[type_it->first])[i][j] = current_clusters[j]->get_maxStrip();
         }
         for(int j=NClus[type_it->first][i];j<current_MaxNClus;j++){
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusAmpl[type_it->first])[i][j] = 0;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusT[type_it->first])[i][j] = 0;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusTOT[type_it->first])[i][j] = 0;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusPos[type_it->first])[i][j] = 0;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusSize[type_it->first])[i][j] = -1;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusMaxSample[type_it->first])[i][j] = 0;
            reinterpret_cast<Double_t(*)[current_MaxNClus]>(ClusMaxStripAmpl[type_it->first])[i][j] = 0;
            reinterpret_cast<Int_t(*)[current_MaxNClus]>(ClusMaxStrip[type_it->first])[i][j] = 0;
         }
         for(vector<Cluster*>::iterator clus_it=current_clusters.begin();clus_it!=current_clusters.end();++clus_it){
            delete *clus_it;
         }
      }
   }
   T->Fill();
}
