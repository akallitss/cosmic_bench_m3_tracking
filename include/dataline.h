#ifndef dataline_h
#define dataline_h
#include <stdint.h>

class DataLineDream{
	public:
		DataLineDream();
		DataLineDream(uint16_t input);
		~DataLineDream();
		void ntohs_();
		bool is_final_trailer() const;
		bool is_data_trailer() const;
		bool is_first_line() const;
		bool is_data() const;
		bool is_data_zs() const;
		bool is_channel_ID() const;
		bool is_Feu_header() const;
		bool is_data_header() const;
		bool is_EOE() const;
		bool get_zs_mode() const;
		int get_Feu_ID() const;
		int get_sample_ID() const;
		int get_channel_ID() const;
		int get_dream_ID() const;
		short get_data() const;
		//unsigned short int data;
		uint16_t data;
};

class DataLineFeminos{
	public:
		DataLineFeminos();
		~DataLineFeminos();
		bool is_info() const;
		int get_card_ID() const;
		int get_chip_ID() const;
		int get_channel_ID() const;
		bool is_time() const;
		int get_time() const;
		bool is_data() const;
		short get_data() const;
		bool is_end_of_frame() const;
		bool is_end_of_event() const;
		bool is_event_start() const;
		bool is_frame_start() const;
		bool is_built_event_start() const;
		bool is_end_of_built_event() const;
		//unsigned short int data;
		uint16_t data;
};

#endif