#ifndef elecreader_h
#define elecreader_h
#include "tomography.h"
#include <string>
#include <map>
#include <fstream>

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
		double data[Tomography::Nasic_FEU][Tomography::Nchannel][Tomography::Nsample];
		long Nevent;
		double evttime;
		ifstream * file;
		int current_index;
};

class FeminosData: public RawData{
	public:
		FeminosData();
		~FeminosData();
		FeminosData(const FeminosData& other);
		FeminosData& operator=(const FeminosData& other);
		double data[Tomography::Nasic_Feminos][Tomography::Nchannel][Tomography::Nsample];
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
		virtual bool is_end() = 0;
	protected:
		int first_index;
		int last_index;
		string base_name;
};

class DreamElecReader: public ElecReader{
	public:
		DreamElecReader();
		~DreamElecReader();
		DreamElecReader(string base_name_,map<int,int> feu_id_to_n_,int first_index_,int last_index_);
		DreamElecReader(const DreamElecReader& other);
		DreamElecReader& operator=(const DreamElecReader& other);
		void read_next_event();
		double get_data(int asic_n,int channel_n,int sample_n); // asic_n = 8*feu_n + asic_n_in_feu
		long get_event_n();
		double get_evttime();
		bool is_end();
		bool is_end_feu(int feu_id);
	protected:
		void reset_data();
		void reset_data(int feu_id);
		void read_next_event_file(int feu_id);
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
		bool is_end();
	protected:
		void reset_data();
		map<int,FeminosData> feminos_data;
		long Nevent;
		double evttime;
		ifstream * file;
		int current_index;
};

#endif