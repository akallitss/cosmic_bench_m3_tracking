#ifndef elecreader_cpp
#define elecreader_cpp
#include "ElecReader.h"
#include "dataline.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdint.h>
#include "TMath.h"

//#include <sys/msg.h>
#include <unistd.h>

using std::ostringstream;
using std::setw;
using std::setfill;
using std::cout;
using std::endl;

using TMath::Ldexp;

#include <set>
using std::set;
#include <boost/filesystem.hpp>
using boost::filesystem::is_directory;
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::filesystem::directory_entry;
using boost::filesystem::directory_iterator;
using boost::filesystem::is_regular_file;
#include <regex>
using std::regex;
using std::smatch;
using std::regex_match;
#include <time.h>

RawData::RawData(){

}
RawData::~RawData(){

}
RawData::RawData(const RawData& /*other*/){
	
}
RawData& RawData::operator=(const RawData& /*other*/){
	return *this;
}

FeuData::FeuData(): RawData(){
	Nevent = -1;
	evttime = 0;
	file = NULL;
	current_index = -1;
	dream_mask = 0;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
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
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			delete[] data[i][j];
		}
	}
}
FeuData::FeuData(const FeuData& other): RawData(other){
	Nevent = other.Nevent;
	evttime = other.evttime;
	file = other.file;
	current_index = other.current_index;
	dream_mask = other.dream_mask;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
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
	dream_mask = other.dream_mask;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
	return *this;
}

FeminosData::FeminosData(): RawData(){
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				data[i][j][k] = 0;
			}
		}
	}
}
FeminosData::~FeminosData(){
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			delete[] data[i][j];
		}
	}
}
FeminosData::FeminosData(const FeminosData& other): RawData(other){
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
}
FeminosData& FeminosData::operator=(const FeminosData& other){
	RawData::operator=(other);
	for(int i=0;i<Tomography::Nasic_Feminos;i++){
		for(int j=0;j<Tomography::Nchannel;j++){
			data[i][j] = new double[Tomography::get_instance()->get_Nsample()];
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				data[i][j][k] = other.data[i][j][k];
			}
		}
	}
	return *this;
}
status_message::status_message(){
	mtype = -1;
}
status_message::~status_message(){
	mtype = -1;
}
data_message::data_message(){
	mtype = -1;
	memset(data,0,max_length*sizeof(uint16_t));
}
data_message::~data_message(){
	mtype = -1;
}

ElecReader::ElecReader(){
	first_index = -1;
	last_index = -1;
	base_name = "";
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

LiveElecReader::LiveElecReader(): ElecReader(){
	//queue_id = -1;
	current_message = NULL;
	message_index = 0;
	current_event_data = NULL;
}
LiveElecReader::LiveElecReader(vector<int> used_asics_, string pipe_name): ElecReader(){
	used_asics = used_asics_;
	for(vector<int>::iterator asics_it=used_asics.begin();asics_it!=used_asics.end();++asics_it){
		//data[*asics_it] = vector<vector<double> >(Tomography::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
		if(dream_mask.count((*asics_it)/Tomography::Nasic_FEU)==0) dream_mask[(*asics_it)/Tomography::Nasic_FEU] = 0;
		dream_mask[(*asics_it)/Tomography::Nasic_FEU] |= 0x1 << (*asics_it)%Tomography::Nasic_FEU;
	}
	//queue_id = msgget(ftok("/proc/self/exe",getpid()), 0666 | IPC_CREAT);
	live_reader_task = new Read_Live_Task(pipe_name/*queue_id*/);
	live_reader_thread = new Reader_Thread(live_reader_task);
	live_reader_thread->start();
	//cout << "listening to queue : " << queue_id << endl;
	current_message = new data_message();
	message_index = data_message::max_length;
	current_event_data = NULL;
}
LiveElecReader::~LiveElecReader(){
	live_reader_thread->stop();
	delete live_reader_thread;
	delete live_reader_task;
	current_event_data = NULL;
	data.clear();
	dream_mask.clear();
	//if(queue_id>0) msgctl(queue_id, IPC_RMID, NULL);
	//queue_id = 0;
	delete current_message;
	message_index = 0;
	used_asics.clear();
}
LiveElecReader::LiveElecReader(const LiveElecReader& other): ElecReader(other){
	//queue_id = other.queue_id;
	live_reader_task = other.live_reader_task;
	live_reader_thread = other.live_reader_thread;
	data = other.data;
	current_event_data = other.current_event_data;
	dream_mask = other.dream_mask;
	current_message = other.current_message;
	message_index = other.message_index;
	used_asics = other.used_asics;
	cout << "Warning, the live elec reader should not be copied !" << endl;
}
LiveElecReader& LiveElecReader::operator=(const LiveElecReader& other){
	ElecReader::operator=(other);;
	//queue_id = other.queue_id;
	live_reader_task = other.live_reader_task;
	live_reader_thread = other.live_reader_thread;
	data = other.data;
	current_event_data = other.current_event_data;
	dream_mask = other.dream_mask;
	current_message = other.current_message;
	message_index = other.message_index;
	used_asics = other.used_asics;
	cout << "Warning, the live elec reader should not be copied !" << endl;
	return *this;
}
void LiveElecReader::build_data(long ev_id, double ev_time){
	data[ev_id].Nevent = ev_id;
	data[ev_id].evttime = ev_time;
	for(vector<int>::iterator asics_it=used_asics.begin();asics_it!=used_asics.end();++asics_it){
		data[ev_id].strip_data[*asics_it] = vector<vector<double> >(Tomography::Nchannel,vector<double>(Tomography::get_instance()->get_Nsample(),0));
	}
}
DataLineDream LiveElecReader::get_next_word(){
	if(message_index>=data_message::max_length){	
		delete current_message;
		//cout << "waiting data" << endl;
		current_message = live_reader_task->wait_new_data();
		message_index = 0;
	}
	//cout << "getting " << message_index << " line of message" << endl;
	DataLineDream current_line((current_message->data)[message_index]);
	current_line.ntohs_();
	message_index++;
	return current_line;
}
void LiveElecReader::get_next_message(){
	delete current_message;
	current_message = live_reader_task->wait_new_data();
	message_index = 0;
}

void LiveElecReader::read_next_event(){
	//cout << "reading event" << endl;
	int isample=-1; //int isample_prev=-2;
	//int isample_nb=0;
	int asic_nb = 0;
	//unsigned int feu_nb = 0;
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
	//int current_event_old = -1;
	double current_evttime = 0;
	if(current_event_data) data.erase(current_event_data->Nevent);
	current_event_data = NULL;
	DataLineDream current_data;
	current_data = get_next_word();
	while(Tomography::get_instance()->get_can_continue() && (live_reader_thread->is_working() || live_reader_task->has_new_data()) && !(event_complete || has_bug)){
		if(FeuHeaderLine<8 && current_data.is_Feu_header()){
			if(FeuHeaderLine==0){
				asic_nb = 0;
				zs_mode = current_data.get_zs_mode();
				FeuN = current_data.get_Feu_ID();
				if(dream_mask.count(FeuN)==0){
					cout << "unknown FEU " << FeuN << endl;
					has_bug = true;
					break;
				}
				//if(FeuN != feu_id) cout << "problem in FeuN to FeuID mapping" << endl;
			}
			else if(FeuHeaderLine==1){
				current_event = current_data.get_data();
			}
			else if(FeuHeaderLine==2){
				current_evttime = current_data.get_data();
			}
			else if(FeuHeaderLine==3){
				//if(feu_nb==0 && asic_nb==0) isample_prev = current_data.get_sample_ID();;
				isample = current_data.get_sample_ID();
				/*
				if(isample!=isample_prev){
					cout << "problem in sample index : " << isample << "; " << isample_prev << endl;
					has_bug = true;
					break;
				}
				*/
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
			/*
			if(feu_nb==0 && asic_nb==0 && isample==0) current_event_old = current_event;
			if(current_event != current_event_old){
				cout << "problem in event id : " << current_event << "; " << current_event_old << endl;
				has_bug = true;
				break;
			}
			*/
			if(data.count(current_event) == 0) build_data(current_event,current_evttime);
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
					//cout << "getting data for asic : " << asicN+(8*FeuN) << " channel : " << ichannel << " sample : " << isample << endl;
					data[current_event].strip_data[asicN+(8*FeuN)][ichannel][isample] = current_data.get_data();
					ichannel++;
				}
				else if(current_data.is_data_zs() && zs_mode){
					if(!got_channel_id){
						ichannel = current_data.get_channel_ID();
						got_channel_id = true;
					}
					else{
						data[current_event].strip_data[asicN+(8*FeuN)][ichannel][isample] = current_data.get_data();
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
					asic_nb |= 0x1 << asicN;
					asicN=0;
					DataHeaderLine=0;
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
				if(asic_nb!=dream_mask[FeuN]){
					cout << "Problem in number of asic (" << asic_nb << ")" << endl;
					has_bug = true;
					break;
				}
				//else feu_nb++;
				//isample_nb++;
				(data[current_event].data_retrieved)++;
				zs_mode = false;
				FeuHeaderLine=0;
				FeuN = 0;
				//cout << "trailer reached" << endl;
				if(current_data.is_EOE()){
					//if(((data[current_event].data_retrieved) % Tomography::get_instance()->get_Nsample()) != 0) cout << "Reached EOE with less than " << Tomography::get_instance()->get_Nsample() << " samples (" << data[current_event].data_retrieved << ")" << endl;
					//(data[current_event].data_retrieved)++;
					isample=-1; //isample_prev=-2;
					//cout << "EOE reached for " << data[current_event].data_retrieved/Tomography::get_instance()->get_Nsample() << " FEUs out of " << dream_mask.size() << endl;
					if((data[current_event].data_retrieved)==((dream_mask.size())*(Tomography::get_instance()->get_Nsample()))){
						event_complete = true;
						break;
					}
				}
				current_data = get_next_word();
				current_data = DataLineDream();
				/*
				if(feu_nb==dream_mask.size()){
					feu_nb=0;
				}
				*/
			}
		}
		current_data = get_next_word();
	}
	if(event_complete){
		//cout << "complete event read !" << endl;
		data[current_event].Nevent = current_event;
		data[current_event].evttime = current_evttime;
		current_event_data = &(data[current_event]);
	}
	if(has_bug){
		cout << "bugged event read !" << endl;
		cout << message_index << " | " << current_event << " | " << FeuN << " | " << asicN << " | " << ichannel << " | " << isample << endl;
		get_next_message();
		data[current_event].Nevent = -1;
		data[current_event].evttime = 0;
		//current_event_data = &(data[current_event]);
	}
}
double LiveElecReader::get_data(int asic_n,int channel_n,int sample_n){
	if(current_event_data) return current_event_data->strip_data[asic_n][channel_n][sample_n];
	else return -1;
}
long LiveElecReader::get_event_n(){
	if(current_event_data) return current_event_data->Nevent;
	else return -1;
}
double LiveElecReader::get_evttime(){
	if(current_event_data) return current_event_data->evttime;
	else return -1;
}
bool LiveElecReader::is_end() const{
	return !((live_reader_task->has_new_data()) || (live_reader_thread->is_working()));
}

DreamElecReader::DreamElecReader(): ElecReader(){
	feu_data.clear();
	feu_id_to_n.clear();
}
DreamElecReader::~DreamElecReader(){
	feu_data.clear();
}
DreamElecReader::DreamElecReader(string base_name_,vector<FeuInfo> feu_info,int first_index_,int last_index_): ElecReader(base_name_,first_index_,last_index_){
	for(vector<FeuInfo>::iterator feu_it=feu_info.begin();feu_it!=feu_info.end();++feu_it){
		feu_id_to_n[feu_it->id] = feu_it->n;
		feu_data[feu_it->id].Nevent = 1;
		feu_data[feu_it->id].evttime = 0;
		feu_data[feu_it->id].current_index = first_index;
		feu_data[feu_it->id].file = new ifstream();
		this->open_file(feu_it->id);
		feu_data[feu_it->id].dream_mask = 0;
		for(int i=0;i<Tomography::Nasic_FEU;i++){
			feu_data[feu_it->id].dream_mask |= (((feu_it->dream_mask)[i]) ? 0x1 : 0x0) << i;
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
					feu_data[feu_it->id].data[i][j][k] = 0;
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
			this->read_next_event_file(evn_it->first);
			while((evn_it->second).Nevent < 0 && !is_end_feu(evn_it->first)){
				this->seek_next_EOE(evn_it->first);
				this->read_next_event_file(evn_it->first);
			}
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
	///reset_data(feu_id);
	DataLineDream current_data;
	while(feu_data[feu_id].current_index <= last_index){
		this->check_file(feu_id);
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
						asic_nb |= 0x1 << asicN;
						asicN=0;
						DataHeaderLine=0;
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
					if(asic_nb!=feu_data[feu_id].dream_mask){
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
						if(isample_nb!=Tomography::get_instance()->get_Nsample()) cout << "Reached EOE with less than " << Tomography::get_instance()->get_Nsample() << " samples (" << isample_nb << ")" << endl;
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
void DreamElecReader::seek_next_EOE(int feu_id){
	if(!(feu_data.count(feu_id)>0)){
		cout << "Feu ID not found for read (" << feu_id << ")" << endl;
		return;
	}
	bool eoe_reached = false;
	int line_skipped = 0;
	reset_data(feu_id);
	DataLineDream current_data;
	while(feu_data[feu_id].current_index <= last_index && !eoe_reached){
		this->check_file(feu_id);
		/*
		(feu_data[feu_id].file)->read((char*)&current_data,sizeof(current_data));
		current_data.ntohs_();
		*/
		do {
			(feu_data[feu_id].file)->read((char*)&current_data,sizeof(current_data));
			current_data.ntohs_();
			line_skipped++;
			if(current_data.is_EOE()) eoe_reached = true;
		} while((feu_data[feu_id].file)->good() && !eoe_reached);
	}
	int unmasked_feu = 0;
	for(int i=0;i<Tomography::Nasic_FEU;i++){
		unmasked_feu += ((feu_data[feu_id].dream_mask >> i) & 0x1);
	}
	cout << "    skipped " << line_skipped << " lines (approx. " << line_skipped/(74+10) << " packets) to realign with dream packet" << endl;
}
void DreamElecReader::check_file(int feu_id){
	if((feu_data[feu_id].file)->eof()){
		if(feu_data[feu_id].current_index == last_index){
			cout << "end of data for FEU : " << feu_id_to_n[feu_id] << endl;
			reset_data(feu_id);
			return;
		}
		(feu_data[feu_id].file)->close();
		feu_data[feu_id].current_index++;
		this->open_file(feu_id);
	}
}
void DreamElecReader::open_file(int feu_id){
	while(!((feu_data[feu_id].file)->is_open()) && feu_data[feu_id].current_index <= last_index){
		ostringstream current_name;
		current_name << base_name << setw(3) << setfill('0') << feu_data[feu_id].current_index << "_" << setw(2) << setfill('0') << feu_id_to_n[feu_id] << "." << Tomography::DreamExt;
		feu_data[feu_id].file->open(current_name.str().c_str(),ifstream::binary);
		if((feu_data[feu_id].file)->is_open()) cout << "\n" << current_name.str() << " loaded !" << endl;
		else{
			cout << "\ncan't load : " << current_name.str() << endl;
			feu_data[feu_id].current_index++;
			reset_data(feu_id);
		}
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
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
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
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				feu_data[feu_id].data[i][j][k] = 0;
			}
		}
	}
}
bool DreamElecReader::is_end() const{
	for(map<int, FeuData>::const_iterator data_it=feu_data.begin();data_it!=feu_data.end();++data_it){
		if(is_end_feu(data_it->first)) return true;
	}
	return false;
}
bool DreamElecReader::is_end_feu(int feu_id) const{
	if(((feu_data.find(feu_id))->second).current_index<last_index) return false;
	if(!((((feu_data.find(feu_id))->second).file)->is_open())) return true;
	if(!((((feu_data.find(feu_id))->second).file)->eof())) return false;
	return true;
}

DreamElecWattoReader::DreamElecWattoReader(): DreamElecReader(){
	evttime_offset = 0;
}
DreamElecWattoReader::~DreamElecWattoReader(){
	timestamp_to_filename.clear();
}
DreamElecWattoReader::DreamElecWattoReader(string directory,vector<FeuInfo> feu_info): DreamElecReader(){
	evttime_offset = 0;
	set<int> feu_n_list;
	map<pair<string,int>,set<int> > feu_list_by_file;
	for(vector<FeuInfo>::iterator feu_it=feu_info.begin();feu_it!=feu_info.end();++feu_it){
		feu_n_list.insert(feu_it->n);
	}

	path directory_p(directory);

	if(!exists(directory_p)){
		cout << "directory : " << directory << " does not exist" << endl;
		return;
	}
	if(!is_directory(directory_p)){
		cout << "directory : " << directory << " is not a directory" << endl;
		return;
	}
	for(directory_iterator file_it = directory_iterator(directory_p);file_it!=directory_iterator();++file_it){
		path current_path = *file_it;
		if(current_path.extension() != ".fdf") continue;
		regex time_regexp(".*([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{2})H([0-9]{2})_([0-9]{3})_([0-9]{2}).*");
		smatch string_match;
		regex_match(current_path.string(),string_match,time_regexp);
		if(string_match.size()!=8){
			cout << "problem parsing filename : " << current_path << endl;
			continue;
		}
		struct tm current_tm;
		strptime(("20"+string_match[1].str()+string_match[2].str()+string_match[3].str()+string_match[4].str()+string_match[5].str()+"00").c_str(),"%Y%m%d%H%M%S",&current_tm);
		unsigned int current_timestamp = mktime(&current_tm);
		string base_file_name = current_path.string().substr(0,current_path.string().size()-10);
		int file_index = atoi(string_match[6].str().c_str());
		int current_feu_n = atoi(string_match[7].str().c_str());
		if(feu_n_list.count(current_feu_n) == 0){
			cout << "file : " << current_path.string() << " is not associated with any known FEU" << endl;
			continue;
		}
		if(timestamp_to_filename.count(current_timestamp) == 0){
			timestamp_to_filename[current_timestamp] = pair<string,int>(base_file_name,file_index);
		}
		else if(timestamp_to_filename[current_timestamp].second < file_index){
			timestamp_to_filename[current_timestamp].second = file_index;
		}
		if(feu_list_by_file.count(pair<string,int>(base_file_name,file_index))>0){
			feu_list_by_file[pair<string,int>(base_file_name,file_index)].insert(current_feu_n);
		}
		else{
			feu_list_by_file[pair<string,int>(base_file_name,file_index)] = set<int>();
			feu_list_by_file[pair<string,int>(base_file_name,file_index)].insert(current_feu_n);
		}
	}

	reading_status = timestamp_to_filename.begin();
	while(reading_status != timestamp_to_filename.end()){
		bool is_good = true;
		for(int i=0;i<((reading_status->second).second);i++){
			if(feu_list_by_file[reading_status->second] != feu_n_list){
				cout << "file missing for some FEU" << endl;
				is_good = false;
				if(i>0){
					cout << "    stopping " << (reading_status->second).first << " at index " << i-1 << endl;
					(reading_status->second).second = i-1;
					++reading_status;
				}
				else{
					cout << "    dropping : " << (reading_status->second).first << endl;
					reading_status = timestamp_to_filename.erase(reading_status);
				}
				break;
			}
		}
		if(is_good){
			cout << "run : " << (reading_status->second).first << " will be processed from 0 to " << ((reading_status->second).second) << endl;
			++reading_status;
		}
	}

	reading_status = timestamp_to_filename.begin();
	first_index = 0;
	last_index = (reading_status->second).second;
	
	for(vector<FeuInfo>::iterator feu_it=feu_info.begin();feu_it!=feu_info.end();++feu_it){
		feu_id_to_n[feu_it->id] = feu_it->n;
		feu_data[feu_it->id].Nevent = 1;
		feu_data[feu_it->id].evttime = 0;
		feu_data[feu_it->id].current_index = 0;
		feu_data[feu_it->id].file = new ifstream();
		this->open_file(feu_it->id);
		feu_data[feu_it->id].dream_mask = 0;
		for(int i=0;i<Tomography::Nasic_FEU;i++){
			feu_data[feu_it->id].dream_mask |= (((feu_it->dream_mask)[i]) ? 0x1 : 0x0) << i;
			for(int j=0;j<Tomography::Nchannel;j++){
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
					feu_data[feu_it->id].data[i][j][k] = 0;
				}
			}
		}
	}
}
DreamElecWattoReader::DreamElecWattoReader(const DreamElecWattoReader& other): DreamElecReader(other){
	timestamp_to_filename = other.timestamp_to_filename;
	reading_status = other.reading_status;
	evttime_offset = other.evttime_offset;

}
DreamElecWattoReader& DreamElecWattoReader::operator=(const DreamElecWattoReader& other){
	ElecReader::operator=(other);
	timestamp_to_filename = other.timestamp_to_filename;
	reading_status = other.reading_status;
	evttime_offset = other.evttime_offset;
	return *this;
}
void DreamElecWattoReader::check_file(int feu_id){
	if((feu_data[feu_id].file)->eof()){
		if(feu_data[feu_id].current_index == last_index){
			if(reading_status==timestamp_to_filename.end()){
				cout << "should not be there, there is no data left to process !" << endl;
				reset_data();
				return;
			}
			cout << "end of data for FEU : " << feu_id_to_n[feu_id] << " for run " << (reading_status->second).first << endl;
			change_run();
			return;
		}
		(feu_data[feu_id].file)->close();
		feu_data[feu_id].current_index++;
		this->open_file(feu_id);
	}
}
void DreamElecWattoReader::open_file(int feu_id){
	if(reading_status==timestamp_to_filename.end()) return;
	while(!((feu_data[feu_id].file)->is_open()) && feu_data[feu_id].current_index <= last_index){
		ostringstream current_name;
		current_name << (reading_status->second).first << setw(3) << setfill('0') << feu_data[feu_id].current_index << "_" << setw(2) << setfill('0') << feu_id_to_n[feu_id] << "." << Tomography::DreamExt;
		feu_data[feu_id].file->open(current_name.str().c_str(),ifstream::binary);
		if((feu_data[feu_id].file)->is_open()) cout << "\n" << current_name.str() << " loaded !" << endl;
		else{
			cout << "\ncan't load : " << current_name.str() << endl;
			feu_data[feu_id].current_index++;
			reset_data(feu_id);
		}
	}
	if(!((feu_data[feu_id].file)->is_open())) change_run();
}
double DreamElecWattoReader::get_evttime(){
	return (((feu_data.begin())->second).evttime + evttime_offset);
}
void DreamElecWattoReader::change_run(){
	reset_data();
	if(reading_status==timestamp_to_filename.end()){
		last_index = -1;
		return;
	}
	++reading_status;
	if(reading_status==timestamp_to_filename.end()){
		last_index = -1;
		return;
	}
	last_index = (reading_status->second).second;
	evttime_offset = ((reading_status->first) - ((timestamp_to_filename.begin())->first))*(Tomography::get_instance()->get_clock_rate());

	for(map<int,FeuData>::iterator feu_it=feu_data.begin();feu_it!=feu_data.end();++feu_it){
		(feu_it->second).current_index = 0;
		((feu_it->second).file)->close();
		this->open_file(feu_it->first);
	}

}
bool DreamElecWattoReader::is_end() const{
	map<unsigned long,pair<string,int> >::const_iterator last_run = timestamp_to_filename.end();
	if(reading_status==last_run) return true;
	--last_run;
	if(reading_status!=last_run) return false;
	for(map<int, FeuData>::const_iterator data_it=feu_data.begin();data_it!=feu_data.end();++data_it){
		if(is_end_feu(data_it->first)) return true;
	}
	return false;
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
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
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
	double evttime_tmp = 0;
	uint32_t Nevent_tmp;
	file->read((char*)&current_data,sizeof(current_data));
	bool event_complete = false;
	bool has_bug = false;
	while(current_index <= last_index){
		if(file->eof()){
			if(current_index == last_index){
				cout << "end of data" << endl;
				reset_data();
				return;
			}
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
				for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
					(data_it->second).data[i][j][k] = 0;
				}
			}
		}
	}
}
bool FeminosElecReader::is_end() const{
	if(current_index<last_index) return false;
	if(!(file->eof())) return false;
	return true;
}

#endif
