#define CAEN_comm_cpp
#include "CAEN_comm.h"

#include <utility>
#include <iostream>
#include <iomanip>

#include <cstring>

using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::showpos;
using std::noshowpos;

#include <CAENHVWrapper.h>

string CAEN_Ch::get_Param_name(Param param){
	switch(param){
		case V0Set :
			return "V0Set";
		case I0Set :
			return "I0Set";
		case V1Set :
			return "V1Set";
		case I1Set :
			return "I1Set";
		case RUp :
			return "RUp";
		case RDown :
			return "RDown";
		case Trip :
			return "Trip";
		case SVMax :
			return "SVMax";
		case VMon :
			return "VMon";
		case IMon :
			return "IMon";
		case Status :
			return "Status";
		case Pw :
			return "Pw";
		case POn :
			return "POn";
		case PDwn :
			return "PDwn";
		case TripInt :
			return "TripInt";
		case TripExt :
			return "TripExt";
		default :
			return "";
	}
}
CAEN_Ch::Param CAEN_Ch::get_Param(string name){
	if(name == "V0Set"){
		return V0Set;
	}
	else if(name == "I0Set"){
		return I0Set;
	}
	else if(name == "V1Set"){
		return V1Set;
	}
	else if(name == "I1Set"){
		return I1Set;
	}
	else if(name == "RUp"){
		return RUp;
	}
	else if(name == "RDown"){
		return RDown;
	}
	else if(name == "Trip"){
		return Trip;
	}
	else if(name == "SVMax"){
		return SVMax;
	}
	else if(name == "VMon"){
		return VMon;
	}
	else if(name == "IMon"){
		return IMon;
	}
	else if(name == "Status"){
		return Status;
	}
	else if(name == "Pw"){
		return Pw;
	}
	else if(name == "POn"){
		return POn;
	}
	else if(name == "PDwn"){
		return PDwn;
	}
	else if(name == "TripInt"){
		return TripInt;
	}
	else if(name == "TripExt"){
		return TripExt;
	}
	else{
		cout << "unknown param" << endl;
		return Status;
	}
}

bool CAEN_Ch::Param_is_float(Param param){
	return (param < 100);
}
bool CAEN_Ch::Param_is_voltage(Param param){
	return (param == V0Set || param == V1Set || param == VMon || param == SVMax);
}

Board_spec::Board_spec(){
	index = -1;
	chan_nb = 0;
	is_pos = false;
}
Board_spec::Board_spec(int index_,int chan_nb_,bool is_pos_){
	if(index_<0 || chan_nb_<0){
		index = -1;
		chan_nb = 0;
		is_pos = false;
	}
	else{
		index = index_;
		chan_nb = chan_nb_;
		is_pos = is_pos_;
	}
}
Board_spec::~Board_spec(){

}
int Board_spec::get_index() const{
	return index;
}
int Board_spec::get_chan_nb() const{
	return chan_nb;
}
bool Board_spec::get_is_pos() const{
	return is_pos;
}

CAEN_Comm::CAEN_Comm(){
	handle = -1;
	boards.clear();
}
CAEN_Comm::CAEN_Comm(string IP, string username, string passwd){
	handle = -1;
	boards.clear();
	CAENHVRESULT ret = CAENHV_InitSystem(SY4527, LINKTYPE_TCPIP, const_cast<char*>(IP.c_str()), username.c_str(), passwd.c_str(), &handle);
	if( ret != CAENHV_OK ){
		cout << "Could Not Connect" << endl;
		handle = -1;
		return;
	}
	unsigned short SlotNb;
	unsigned short * ChanNbList;
	char * ModelListRaw;
	char * DescriptionListRaw;
	unsigned short * SerNumList;
	unsigned char * FMWRelMinList;
	unsigned char * FMWRelMaxList;
	ret = CAENHV_GetCrateMap(handle,&SlotNb,&ChanNbList,&ModelListRaw,&DescriptionListRaw,&SerNumList,&FMWRelMinList,&FMWRelMaxList);
	if( ret != CAENHV_OK ){
		cout << "Could Not Fetch CrateMap" << endl;
		handle = -1;
		return;
	}
	for(int i=0;i<SlotNb;i++){
		if(ChanNbList[i]>0){
			string current_description = DescriptionListRaw;
			boards[i] = Board_spec(i,ChanNbList[i],(current_description.find("Pos") != string::npos));
		}
		DescriptionListRaw += strlen(DescriptionListRaw) + 1;
	}
}
CAEN_Comm::~CAEN_Comm(){
	CAENHV_DeinitSystem(handle);
	handle = -1;
	boards.clear();
}

map<int,map<int,string> > CAEN_Comm::get_Ch_name(){
	map<int,map<int,string> > return_map;
	CAENHVRESULT ret;
	for(map<int,Board_spec>::iterator board_it = boards.begin();board_it!=boards.end();++board_it){
		int current_chan_nb = (board_it->second).get_chan_nb();
		unsigned short * listChanNb = new unsigned short[current_chan_nb];
		for(int i=0;i<current_chan_nb;i++){
			listChanNb[i] = i;
		}
		char (*listChanName)[MAX_CH_NAME] = new char[current_chan_nb][MAX_CH_NAME];
		ret = CAENHV_GetChName(handle,(board_it->second).get_index(),current_chan_nb,listChanNb,listChanName);
		if(ret != CAENHV_OK){
			cout << "Could not fetch board " << (board_it->second).get_index() << " channels name" << endl;
			delete listChanName;
		}
		else{
			for(int i=0;i<current_chan_nb;i++){
				return_map[(board_it->second).get_index()][i] = listChanName[i];
			}
		}
		delete listChanNb;
	}
	return return_map;
}

map<int,map<int,CAEN_Ch::param_value> > CAEN_Comm::get_Ch_param(map<int,vector<int> > Ch_list, CAEN_Ch::Param param){
	map<int,map<int,CAEN_Ch::param_value> > return_map;
	CAENHVRESULT ret;
	for(map<int,vector<int> >::iterator ch_list_it = Ch_list.begin();ch_list_it!=Ch_list.end();++ch_list_it){
		if(boards.count(ch_list_it->first)==0){
			cout << "Board " << ch_list_it->first << " unknown" << endl;
			continue;
		}
		unsigned short * listChanNb = new unsigned short[(ch_list_it->second).size()];
		bool wrong_ch = false;
		for(unsigned int i=0;i<(ch_list_it->second).size();i++){
			listChanNb[i] = (ch_list_it->second)[i];
			if(listChanNb[i] > (boards[ch_list_it->first].get_chan_nb() - 1)){
				wrong_ch = true;
				break;
			}
		}
		if(wrong_ch){
			cout << "Invalid channel index" << endl;
			continue;
		}
		CAEN_Ch::param_value (*listChanParam) = new CAEN_Ch::param_value[(ch_list_it->second).size()];
		ret = CAENHV_GetChParam(handle,ch_list_it->first,CAEN_Ch::get_Param_name(param).c_str(),(ch_list_it->second).size(),listChanNb,listChanParam);
		if(ret != CAENHV_OK){
			cout << "Could not fetch param value" << endl;
			delete listChanParam;
			continue;
		}
		for(unsigned int i=0;i<(ch_list_it->second).size();i++){
			return_map[ch_list_it->first][listChanNb[i]] = listChanParam[i];
			if(CAEN_Ch::Param_is_voltage(param) && !(boards[ch_list_it->first].get_is_pos())){
				return_map[ch_list_it->first][listChanNb[i]].real = -return_map[ch_list_it->first][listChanNb[i]].real;
			}

		}
		delete listChanParam;
		delete listChanNb;
	}
	return return_map;
}
