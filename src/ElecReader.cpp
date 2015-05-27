#ifndef elecreader_cpp
#define elecreader_cpp
#include "ElecReader.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include "dataline.h"
#include "TMath.h"
#include <stdint.h>

using std::ostringstream;
using std::setw;
using std::setfill;
using std::cout;
using std::endl;

using TMath::Ldexp;

RawData::RawData(){

}
RawData::~RawData(){

}
RawData::RawData(const RawData& other){
	
}
RawData& RawData::operator=(const RawData& other){
	return *this;
}

FeuData::FeuData(): RawData(){
	Nevent = -1;
	evttime = 0;
	file = NULL;
	current_index = -1;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = 0;
			}
		}
	}
}
FeuData::~FeuData(){
	if(file != NULL){
		file->close();
		delete file;
	}
}
FeuData::FeuData(const FeuData& other): RawData(other){
	Nevent = other.Nevent;
	evttime = other.evttime;
	file = other.file;
	current_index = other.current_index;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
}
FeuData& FeuData::operator=(const FeuData& other){
	RawData::operator=(other);
	Nevent = other.Nevent;
	evttime = other.evttime;
	file = other.file;
	current_index = other.current_index;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
	return *this;
}

FeminosData::FeminosData(): RawData(){
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = 0;
			}
		}
	}
}
FeminosData::~FeminosData(){

}
FeminosData::FeminosData(const FeminosData& other): RawData(other){
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
}
FeminosData& FeminosData::operator=(const FeminosData& other){
	RawData::operator=(other);
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
	return *this;
}

ElecReader::ElecReader(){
	first_index = -1;
	last_index = -1;
	base_name.clear();
}
ElecReader::~ElecReader(){
	base_name.clear();
}
ElecReader::ElecReader(string base_name_,int first_index_,int last_index_){
	first_index = first_index_;
	last_index = last_index_;
	base_name = base_name_;
}
ElecReader::ElecReader(const ElecReader& other){
	first_index = other.first_index;
	last_index = other.last_index;
	base_name = other.base_name;
}
ElecReader& ElecReader::operator=(const ElecReader& other){
	first_index = other.first_index;
	last_index = other.last_index;
	base_name = other.base_name;
	return *this;
}

DreamElecReader::DreamElecReader(): ElecReader(){
	feu_data.clear();
	feu_id_to_n.clear();
}
DreamElecReader::~DreamElecReader(){
	feu_data.clear();
}
DreamElecReader::DreamElecReader(string base_name_,map<int,int> feu_id_to_n_,int first_index_,int last_index_): ElecReader(base_name_,first_index_,last_index_){
	feu_id_to_n = feu_id_to_n_;
	for(map<int,int>::iterator feu_it=feu_id_to_n.begin();feu_it!=feu_id_to_n.end();++feu_it){
		feu_data[feu_it->first].Nevent = 1;
		feu_data[feu_it->first].evttime = 0;
		feu_data[feu_it->first].current_index = first_index;
		ostringstream current_name;
		current_name << base_name << setw(3) << setfill('0') << feu_data[feu_it->first].current_index << "_" << setw(2) << setfill('0') << feu_it->second << "." << Tomography::DreamExt;
		feu_data[feu_it->first].file = new ifstream(current_name.str().c_str(),ifstream::binary);
		if((feu_data[feu_it->first].file)->is_open()) cout << "\n" << current_name.str() << " loaded !" << endl;
		for(int i=0;i<Tomography::Nasic_FEU;i++){
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::Nsample;k++){
					feu_data[feu_it->first].data[i][j][k] = 0;
				}
			}
		}
	}
}
DreamElecReader::DreamElecReader(const DreamElecReader& other): ElecReader(other){
	feu_data = other.feu_data;
	feu_id_to_n = other.feu_id_to_n;
}
DreamElecReader& DreamElecReader::operator=(const DreamElecReader& other){
	ElecReader::operator=(other);
	feu_data = other.feu_data;
	feu_id_to_n = other.feu_id_to_n;
	return *this;
}
void DreamElecReader::read_next_event(){
	long current_Nevent = -1;
	bool desync = false;
	map<int,bool> is_sync;
	for(map<int,FeuData>::iterator evn_it = feu_data.begin();evn_it!=feu_data.end();++evn_it){
		is_sync[evn_it->first] = false;
	}
	do{
		for(map<int,FeuData>::iterator evn_it = feu_data.begin();evn_it!=feu_data.end();++evn_it){
			if(is_sync[evn_it->first]) continue;
			read_next_event_file(evn_it->first);
		}
		current_Nevent = ((feu_data.begin())->second).Nevent;
		desync = false;
		for(map<int,FeuData>::iterator evn_it = feu_data.begin();evn_it!=feu_data.end();++evn_it){
			if((evn_it->second).Nevent != current_Nevent) desync = true;
			if((evn_it->second).Nevent > current_Nevent){
				current_Nevent = (evn_it->second).Nevent;
			}
		}
		if(desync){
			cout << "Desync detected, current event number in files are : " << endl;
			for(map<int,FeuData>::iterator evn_it = feu_data.begin();evn_it!=feu_data.end();++evn_it){
				is_sync[evn_it->first] = ((evn_it->second).Nevent == current_Nevent);
				cout << "    " << evn_it->first << " | " << (evn_it->second).Nevent << endl;
			}
		}
	} while(desync);
}
void DreamElecReader::read_next_event_file(int feu_id){
	if(!(feu_data.count(feu_id)>0)){
		cout << "Feu ID not found for read (" << feu_id << ")" << endl;
		return;
	}
	int isample=-1; int isample_prev=-2;
	int isample_nb=0;
	int asic_nb = 0;
	int ichannel=0;
	int asicN=0;
	int FeuN=0;
	int FeuHeaderLine=0;
	int DataHeaderLine=0;
	bool zs_mode = false;
	bool got_channel_id=false;
	bool event_complete = false;
	bool has_bug = false;
	int current_event = -1;
	int current_event_old = -1;
	double current_evttime = 0;
	reset_data(feu_id);
	DataLineDream current_data;
	while(feu_data[feu_id].current_index <= last_index){
		if((feu_data[feu_id].file)->eof()){
			(feu_data[feu_id].file)->close();
			feu_data[feu_id].current_index++;
			ostringstream current_name;
			current_name << base_name << setw(3) << setfill('0') << feu_data[feu_id].current_index << "_" << setw(2) << setfill('0') << feu_id_to_n[feu_id] << "." << Tomography::DreamExt;
			feu_data[feu_id].file->open(current_name.str().c_str(),ifstream::binary);
			if((feu_data[feu_id].file)->is_open()) cout << current_name.str() << " loaded !" << endl;
			else{
				cout << "can't load : " << current_name.str() << endl;
				reset_data(feu_id);
				return;
			}
		}
		(feu_data[feu_id].file)->read((char*)&current_data,sizeof(current_data));
		current_data.ntohs_();
		while((feu_data[feu_id].file)->good()){
			if(FeuHeaderLine<8 && current_data.is_Feu_header()){
				if(FeuHeaderLine==0){
					asic_nb = 0;
					zs_mode = current_data.get_zs_mode();
					FeuN = current_data.get_Feu_ID();
					if(FeuN != feu_id) cout << "problem in FeuN to FeuID mapping" << endl;
				}
				else if(FeuHeaderLine==1){
					current_event = current_data.get_data();
				}
				else if(FeuHeaderLine==2){
					current_evttime = current_data.get_data();
				}
				else if(FeuHeaderLine==3){
					isample_prev = isample;
					isample = current_data.get_sample_ID();
					if(isample!=isample_prev+1){
						cout << "problem in sample index : " << isample << "; " << isample_prev << endl;
						has_bug = true;
						break;
					}
				}
				else if(FeuHeaderLine==4){
					current_event += static_cast<int>(current_data.get_data()) << 12;
				}
				else{
					current_evttime += Ldexp(static_cast<double>(current_data.get_data()),12*(FeuHeaderLine-4)); // data*2^(12*(FeuHeaderLine-4))
				}
				FeuHeaderLine++;
			}
			else if(FeuHeaderLine>7 && current_data.is_Feu_header()){
				cout << "problem in Feu header" << endl;
				has_bug = true;
				break;
			}
			else if(FeuHeaderLine>3){
				if(DataHeaderLine<4 && current_data.is_data_header()){
					asicN = current_data.get_dream_ID();
					DataHeaderLine++;
				}
				else if(DataHeaderLine>3 && current_data.is_data_header()){
					cout << "problem in data header" << endl;
					has_bug = true;
					break;
				}
				else if(DataHeaderLine>3){
					if(current_data.is_data() && !zs_mode){
						feu_data[feu_id].data[asicN][ichannel][isample] = current_data.get_data();
						ichannel++;
					}
					else if(current_data.is_data_zs() && zs_mode){
						if(!got_channel_id){
							ichannel = current_data.get_channel_ID();
							got_channel_id = true;
						}
						else{
							feu_data[feu_id].data[asicN][ichannel][isample] = current_data.get_data();
							got_channel_id = false;
						}
					}
					else if(current_data.is_data_trailer()){
						if(ichannel!=64 && !zs_mode){
							cout << "problem in channel number" << endl;
							has_bug = true;
							break;
						}
						if(got_channel_id){
							cout << "problem in ZS data" << endl;
							has_bug = true;
							break;
						}
						if(DataHeaderLine!=4 && DataHeaderLine!=1){
							cout << "problem in data header" << endl;
							has_bug = true;
							break;
						}
						ichannel=0;
						asicN=0;
						DataHeaderLine=0;
						asic_nb++;
					}
				}
				else if(current_data.is_final_trailer()){
					if(ichannel!=0){
						cout << "problem in channel number" << endl;
						has_bug = true;
						break;
					}
					if(got_channel_id){
						cout << "problem in ZS data" << endl;
						has_bug = true;
						break;
					}
					if(FeuHeaderLine!=4 && FeuHeaderLine!=8){
						cout << "problem in Feu Header" << endl;
						has_bug = true;
						break;
					}
					if(asic_nb!=Tomography::Nasic_FEU){
						cout << "Problem in number of asic (" << asic_nb << ")" << endl;
						has_bug = true;
						break;
					}
					if(isample==0) current_event_old = current_event;
					else if(current_event_old != current_event && !has_bug){
						cout << "problem in event ID (" << current_event_old << ";" << current_event << ")[sample=" << isample << "]" << endl;
						//has_bug = true;
						//break; //comment these 2 lines until bugfix by irakli :)
					}
					isample_nb++;
					zs_mode = false;
					FeuHeaderLine=0;
					(feu_data[feu_id].file)->ignore(sizeof(current_data));
					if(current_data.is_EOE()){
						if(isample_nb!=Tomography::Nsample) cout << "Reached EOE with less than " << Tomography::Nsample << " samples (" << isample_nb << ")" << endl;
						isample=-1; isample_prev=-2;
						event_complete = true;
						break;
					}
				}
			}
			(feu_data[feu_id].file)->read((char*)&current_data,sizeof(current_data));
			current_data.ntohs_();	
		}
		if(event_complete || has_bug) break;
	}
	if(event_complete){
		feu_data[feu_id].Nevent = current_event;
		feu_data[feu_id].evttime = current_evttime;
	}
	if(has_bug){
		feu_data[feu_id].Nevent = -1;
		feu_data[feu_id].evttime = 0;

	}
}
double DreamElecReader::get_data(int asic_n,int channel_n,int sample_n){
	if(!(feu_data.count(asic_n/Tomography::Nasic_FEU)>0)){
		cout << "Unknown Feu " << asic_n/Tomography::Nasic_FEU << endl;
		return 0;
	}
	return feu_data[asic_n/Tomography::Nasic_FEU].data[asic_n%Tomography::Nasic_FEU][channel_n][sample_n];
}
long DreamElecReader::get_event_n(){
	long current_Nevent = ((feu_data.begin())->second).Nevent;
	bool desync = false;
	for(map<int,FeuData>::iterator event_it = feu_data.begin();event_it!=feu_data.end();++event_it){
		if((event_it->second).Nevent != current_Nevent) desync = true;
		if((event_it->second).Nevent > current_Nevent) current_Nevent = (event_it->second).Nevent;
	}
	if(desync) cout << "Warning ! Event ID desync detected ! Returned event ID might be inaccurate" << endl;
	return current_Nevent;
}
double DreamElecReader::get_evttime(){
	return ((feu_data.begin())->second).evttime;
}
void DreamElecReader::reset_data(){
	for(map<int, FeuData>::iterator data_it=feu_data.begin();data_it!=feu_data.end();++data_it){
		(data_it->second).Nevent = -1;
		(data_it->second).evttime = 0;
		for(int i=0;i<Tomography::Nasic_FEU;i++){
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::Nsample;k++){
					(data_it->second).data[i][j][k] = 0;
				}
			}
		}
	}
}
void DreamElecReader::reset_data(int feu_id){
	if(!(feu_data.count(feu_id)>0)){
		cout << "Feu ID not found for reset (" << feu_id << ")" << endl;
		return;
	}
	feu_data[feu_id].Nevent = -1;
	feu_data[feu_id].evttime = 0;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			for(int k=0;k<Tomography::Nsample;k++){
				feu_data[feu_id].data[i][j][k] = 0;
			}
		}
	}
}
bool DreamElecReader::is_end(){
	for(map<int, FeuData>::iterator data_it=feu_data.begin();data_it!=feu_data.end();++data_it){
		if((data_it->second).current_index<last_index) return false;
		if(!((data_it->second).file->eof())) return false;
	}
	return true;
}

FeminosElecReader::FeminosElecReader(): ElecReader(){
	Nevent = -1;
	evttime = 0;
	file = NULL;
	current_index = -1;
	feminos_data.clear();
}
FeminosElecReader::~FeminosElecReader(){
	if(file != NULL){
		file->close();
		delete file;
	}
	feminos_data.clear();
}
FeminosElecReader::FeminosElecReader(string base_name_,vector<int> fem_id,int first_index_,int last_index_): ElecReader(base_name_,first_index_,last_index_){
	Nevent = 1;
	evttime = 0;
	current_index = first_index;
	ostringstream current_name;
	current_name << base_name << setw(3) << setfill('0') << current_index << "." << Tomography::FeminosExt;
	file = new ifstream(current_name.str().c_str(),ifstream::binary);
	if(file->is_open()){
		cout << current_name.str() << " loaded !" << endl;
		file->ignore(26); //You can read run UID here
		file->ignore(2); //Beginning of frame (28 bytes from the beginning)
	}
	for(vector<int>::iterator fem_it=fem_id.begin();fem_it!=fem_id.end();++fem_it){
		for(int i=0;i<Tomography::Nasic_Feminos;i++){
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::Nsample;k++){
					feminos_data[*fem_it].data[i][j][k] = 0;
				}
			}
		}
	}
}
FeminosElecReader::FeminosElecReader(const FeminosElecReader& other): ElecReader(other){
	Nevent = other.Nevent;
	evttime = other.evttime;
	file = other.file;
	current_index = other.current_index;
	feminos_data = other.feminos_data;
}
FeminosElecReader& FeminosElecReader::operator=(const FeminosElecReader& other){
	ElecReader::operator=(other);
	Nevent = other.Nevent;
	evttime = other.evttime;
	file = other.file;
	current_index = other.current_index;
	feminos_data = other.feminos_data;
	return *this;
}
void FeminosElecReader::read_next_event(){
	int card=0;
	int chip=0;
	int channel=0;
	int itime = 0;
	bool inEvent = false;
	bool inFrame = false;
	int event_started = 0;
	DataLineFeminos current_data;
	double evttime_tmp;
	uint32_t Nevent_tmp;
	file->read((char*)&current_data,sizeof(current_data));
	bool event_complete = false;
	bool has_bug = false;
	while(current_index <= last_index){
		if(file->eof()){
			file->close();
			current_index++;
			ostringstream current_name;
			current_name << base_name << setw(3) << setfill('0') << current_index << "." << Tomography::FeminosExt;
			file->open(current_name.str().c_str(),ifstream::binary);
			if(file->is_open()){
				cout << current_name.str() << " loaded !" << endl;
				file->ignore(26); //You can read run UID here
				file->ignore(2); //Beginning of frame (28 bytes from the beginning)
			}
			else{
				cout << "can't load : " << current_name.str() << endl;
				reset_data();
				return;
			}
		}
		while(file->good()){
			if(inEvent){
				if(inFrame){
					if(current_data.is_event_start()){
						event_started++;
						uint64_t evttime_int = 0;
						file->read((char*)&evttime_int,3*sizeof(DataLineFeminos));
						evttime_tmp = evttime_int >> 16;
						//file->ignore(3*sizeof(current_data)); //contain timestamp
						file->read((char*)&Nevent_tmp,sizeof(Nevent_tmp));
					}
					else if(current_data.is_end_of_event()){
						event_started--;
						file->ignore(sizeof(current_data)); //contain eventsize
					}
					else if(current_data.is_end_of_frame()){
						inFrame = false;
					}
					else if(current_data.is_info()){
						card = current_data.get_card_ID();
						chip = current_data.get_chip_ID();
						channel = current_data.get_channel_ID();
						itime = 0;
					}
					else if(current_data.is_time()){
						itime= current_data.get_time();
					}
					else if(current_data.is_data()){
						feminos_data[card].data[chip][channel][itime] = current_data.get_data();
						itime++;
					}

				}
				else if(current_data.is_end_of_built_event()){
					if(event_started != 0){
						cout << "problem in fem number" << endl;
						has_bug = true;
						break;
					}
					inEvent = false;
					event_complete = true;
					break;
				}
				else if(current_data.is_frame_start()){
					inFrame = true;
				}
			}
			else if(current_data.is_built_event_start()){
				inEvent = true;
				reset_data();
				chip=0;
				channel=0;
				itime = 0;
				event_started = 0;
			}
			file->read((char*)&current_data,sizeof(current_data));
		}
		if(event_complete){
			Nevent = Nevent_tmp;
			evttime = evttime_tmp;
			break;
		}
		if(has_bug){
			Nevent = -1;
			evttime = 0;
			break;
		}
	}
}
double FeminosElecReader::get_data(int asic_n,int channel_n,int sample_n){
	if(!(feminos_data.count(asic_n/Tomography::Nasic_Feminos)>0)){
		cout << "Unknown Fem " << asic_n/Tomography::Nasic_Feminos << endl;
		return 0;
	}
	return feminos_data[asic_n/Tomography::Nasic_Feminos].data[asic_n%Tomography::Nasic_Feminos][channel_n][sample_n];
}
long FeminosElecReader::get_event_n(){
	return Nevent;
}
double FeminosElecReader::get_evttime(){
	return evttime;
}
void FeminosElecReader::reset_data(){
	Nevent = 1;
	evttime = 0;
	for(map<int, FeminosData>::iterator data_it=feminos_data.begin();data_it!=feminos_data.end();++data_it){
		for(int i=0;i<Tomography::Nasic_Feminos;i++){
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::Nsample;k++){
					(data_it->second).data[i][j][k] = 0;
				}
			}
		}
	}
}
bool FeminosElecReader::is_end(){
	if(current_index<last_index) return false;
	if(!(file->eof())) return false;
	return true;
}

#endif