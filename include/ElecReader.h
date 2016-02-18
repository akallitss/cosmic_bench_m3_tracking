#ifndef elecreader_h
#define elecreader_h
#include "tomography.h"
#include <string>
#include <map>
#include <fstream>
#include <stdint.h>
#include "MT_tomography.h"
#include "task/read_live_task.h"
#include "dataline.h"

using std::string;
using std::map;
using std::ifstream;


class RawData{
	public:
		RawData();
		virtual ~RawData() = 0;
		RawData(const RawData& other);
		RawData& operator=(const RawData& other);
};
class FeuData: public RawData{
	public:
		FeuData();
		~FeuData();
		FeuData(const FeuData& other);
		FeuData& operator=(const FeuData& other);
		double * data[Tomography::Nasic_FEU][Tomography::Nchannel];
		long Nevent;
		double evttime;
		ifstream * file;
		int current_index;
		int dream_mask;
};
struct FeuInfo{
	int id;
	int n;
	bool dream_mask[Tomography::Nasic_FEU];
};

class status_message{
	public:
		status_message();
		~status_message();
		long mtype; //==1
		short status;
};
class data_message{
	public:
		static const unsigned int max_length = 1114;
		data_message();
		~data_message();
		long mtype; //==2
		uint16_t data[1114]; // 1 uint16_t == 1 16-bit word, if 0 dreams are masked and every channel is above ZS thr the packet is 1114-word long which is the largest possible packet
};

class FeminosData: public RawData{
	public:
		FeminosData();
		~FeminosData();
		FeminosData(const FeminosData& other);
		FeminosData& operator=(const FeminosData& other);
		double * data[Tomography::Nasic_Feminos][Tomography::Nchannel];
};

class ElecReader{
	public:
		ElecReader();
		virtual ~ElecReader();
		ElecReader(string base_name_,int first_index_,int last_index_);
		ElecReader(const ElecReader& other);
		ElecReader& operator=(const ElecReader& other);
		virtual void read_next_event() = 0;
		virtual double get_data(int asic_n,int channel_n,int sample_n) = 0;
		virtual long get_event_n() = 0;
		virtual double get_evttime() = 0;
		virtual bool is_end() const = 0;
	protected:
		int first_index;
		int last_index;
		string base_name;
};

class LiveElecReader: public ElecReader{
	public:
		LiveElecReader();
		~LiveElecReader();
		LiveElecReader(vector<int> used_asics);
		LiveElecReader(const LiveElecReader& other);
		LiveElecReader& operator=(const LiveElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n) const;
		long get_event_n() const;
		double get_evttime() const;
		bool is_end() const;
	protected:
		void reset_data();
		DataLineDream get_next_word();
		data_message * current_message;
		unsigned int message_index;
		int queue_id;
		Read_Live_Task * live_reader_task;
		Reader_Thread * live_reader_thread;
		map<int,vector<vector<double> > > data;
		map<int,int> dream_mask;
		long Nevent;
		double evttime;
};

class DreamElecReader: public ElecReader{
	public:
		DreamElecReader();
		~DreamElecReader();
		DreamElecReader(string base_name_,vector<FeuInfo> feu_info,int first_index_,int last_index_);
		DreamElecReader(const DreamElecReader& other);
		DreamElecReader& operator=(const DreamElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n); // asic_n = 8*feu_n + asic_n_in_feu
		long get_event_n();
		double get_evttime();
		bool is_end() const;
		bool is_end_feu(int feu_id) const;
	protected:
		void reset_data();
		void reset_data(int feu_id);
		void read_next_event_file(int feu_id);
		void seek_next_EOE(int feu_id);
		void check_file(int feu_id);
		void open_file(int feu_id);
		map<int,FeuData> feu_data;
		map<int,int> feu_id_to_n;
};

class FeminosElecReader: public ElecReader{
	public:
		FeminosElecReader();
		~FeminosElecReader();
		FeminosElecReader(string base_name_,vector<int> fem_id,int first_index_,int last_index_);
		FeminosElecReader(const FeminosElecReader& other);
		FeminosElecReader& operator=(const FeminosElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n); // asic_n = 4*feminos_n + asic_n_in_feminos
		long get_event_n();
		double get_evttime();
		bool is_end() const;
	protected:
		void reset_data();
		map<int,FeminosData> feminos_data;
		long Nevent;
		double evttime;
		ifstream * file;
		int current_index;
};

#endif