#ifndef elecreader_h
#define elecreader_h
#include "tomography.h"
#include <string>
#include <map>
#include <utility>
#include <fstream>
#include <stdint.h>
#include "MT_tomography.h"
#include "task/read_live_task.h"
#include "dataline.h"

using std::string;
using std::map;
using std::pair;
using std::ifstream;


//abstract helper class to store raw information about the data output of the electronics
class RawData{
	public:
		RawData();
		virtual ~RawData() = 0;
		RawData(const RawData& other);
		RawData& operator=(const RawData& other);
};
//implementation of the RawData for the FEU/Dream electronics
class FeuData: public RawData{
	public:
		FeuData();
		~FeuData();
		FeuData(const FeuData& other);
		FeuData& operator=(const FeuData& other);
		double * data[Tomography::Nasic_FEU][Tomography::Nchannel]; //possibly buggy
		long Nevent;
		double evttime;
		ifstream * file;
		int current_index;
		int dream_mask;
};
//store the FEU informations
struct FeuInfo{
	int id;
	int n;
	bool dream_mask[Tomography::Nasic_FEU];
};

//message queue compliant class to pass status information
class status_message{
	public:
		status_message();
		~status_message();
		long mtype; //==1
		short status;
};
//message queue compliant class to pass data
class data_message{
	public:
		static const unsigned int max_length = 1114;
		data_message();
		~data_message();
		long mtype; //==2
		uint16_t data[1114]; // 1 uint16_t == 1 16-bit word, if 0 dreams are masked and every channel is above ZS thr the packet is 1114-word long which is the largest possible packet
};
//implementation of the RawData for the Feminos electronics
class FeminosData: public RawData{
	public:
		FeminosData();
		~FeminosData();
		FeminosData(const FeminosData& other);
		FeminosData& operator=(const FeminosData& other);
		double * data[Tomography::Nasic_Feminos][Tomography::Nchannel]; //possibly buggy
};

//abstract class used to read electronics specific data packets
class ElecReader{
	public:
		ElecReader();
		virtual ~ElecReader();
		//constructor giving the base electronics file name and the index limits (because the data is splitted in many files of usually 1GB)
		ElecReader(string base_name_,int first_index_,int last_index_);
		//copy constructor
		ElecReader(const ElecReader& other);
		//copy assignment
		ElecReader& operator=(const ElecReader& other);
		//read the file until a complete event has been decoded
		virtual void read_next_event() = 0;
		//retrieve data of current event from given asic, channel and sample
		virtual double get_data(int asic_n,int channel_n,int sample_n) = 0;
		//retrieve event id of current event
		virtual long get_event_n() = 0;
		//retrieve event timestamp of current event
		virtual double get_evttime() = 0;
		//check if there is more data to read
		virtual bool is_end() const = 0;
	protected:
		int first_index;
		int last_index;
		string base_name;
};

//struct to store all of an event raw data (before mapping)
struct asic_event{
	public:
		long Nevent;
		double evttime;
		map<int,vector<vector<double> > > strip_data;
		unsigned long data_retrieved;
		asic_event(): Nevent(-1), evttime(0), strip_data(), data_retrieved(0){}
};

//implementation of the ElecReader for online event reading (message queue) using the -q option of FeuDataReader
class LiveElecReader: public ElecReader{
	public:
		LiveElecReader();
		~LiveElecReader();
		LiveElecReader(vector<int> used_asics_, string pipe_name);
		LiveElecReader(const LiveElecReader& other);
		LiveElecReader& operator=(const LiveElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n);
		long get_event_n();
		double get_evttime();
		bool is_end() const;
	protected:
		void build_data(long ev_id, double ev_time);
		DataLineDream get_next_word();
		void get_next_message();
		data_message * current_message;
		unsigned int message_index;
		vector<int> used_asics;
		//int queue_id;
		Read_Live_Task * live_reader_task;
		Reader_Thread * live_reader_thread;
		map<long, asic_event > data;
		asic_event * current_event_data;
		map<int,int> dream_mask;
};

//implementation of the ElecReader for offline reading of FeuUdpControl made .fdf files
class DreamElecReader: public ElecReader{
	public:
		DreamElecReader();
		virtual ~DreamElecReader();
		DreamElecReader(string base_name_,vector<FeuInfo> feu_info,int first_index_,int last_index_);
		DreamElecReader(const DreamElecReader& other);
		DreamElecReader& operator=(const DreamElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n); // asic_n = 8*feu_n + asic_n_in_feu
		virtual long get_event_n();
		virtual double get_evttime();
		virtual bool is_end() const;
		bool is_end_feu(int feu_id) const;
	protected:
		void reset_data();
		void reset_data(int feu_id);
		void read_next_event_file(int feu_id);
		void seek_next_EOE(int feu_id);
		virtual void check_file(int feu_id);
		virtual void open_file(int feu_id);
		map<int,FeuData> feu_data;
		map<int,int> feu_id_to_n;
};

//implementation of the ElecReader for offline reading of FeuUdpControl made .fdf files, synchronizing multiple small run from a directory in a single one with consistant timestamps
class DreamElecWattoReader: public DreamElecReader{
	public:
		DreamElecWattoReader();
		~DreamElecWattoReader();
		DreamElecWattoReader(string directory,vector<FeuInfo> feu_info);
		DreamElecWattoReader(const DreamElecWattoReader& other);
		DreamElecWattoReader& operator=(const DreamElecWattoReader& other);
		double get_evttime();
		bool is_end() const;
	protected:
		void check_file(int feu_id);
		void open_file(int feu_id);
		void change_run();
		map<unsigned long,pair<string,int> > timestamp_to_filename;
		map<unsigned long,pair<string,int> >::iterator reading_status;
		double evttime_offset;
};

//implementation of the ElecReader for offline reading of mclient made files
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
