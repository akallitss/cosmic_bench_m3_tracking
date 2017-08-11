#ifndef dataline_h
#define dataline_h
#include <stdint.h>

//helper class to read DREAM raw data
//see Dream packet of dream user manual for more details
class DataLineDream{
	public:
		DataLineDream();
		DataLineDream(uint16_t input);
		~DataLineDream();
		//change endianness
		void ntohs_();
		//check if data corresponds to a packet final trailer
		bool is_final_trailer() const;
		//check if data corresponds to a data trailer
		bool is_data_trailer() const;
		//check if data corresponds to the first line of a packet
		bool is_first_line() const;
		//check if data corresponds to signal amplitude
		bool is_data() const;
		//check if data corresponds to a signal amplitude in Zero Suppress Mode
		bool is_data_zs() const;
		//check if data corresponds to a channel ID
		bool is_channel_ID() const;
		//check if data corresponds to the FEU header of the packet
		bool is_Feu_header() const;
		//check if data corresponds to the data header
		bool is_data_header() const;
		//check if data corresponds to the End Of Event word
		bool is_EOE() const;
		//retrive the Zero Suppress Mode if is_data_zs()
		bool get_zs_mode() const;
		//retrieve the FEU id if is_Feu_header()
		int get_Feu_ID() const;
		//retrieve the sample id if is_Feu_header()
		int get_sample_ID() const;
		//retrieve the sample id if is_data_zs()
		int get_channel_ID() const;
		//retrieve the sample id if is_data_header()
		int get_dream_ID() const;
		//retrieve signal amplitude if is_data()
		short get_data() const;
		//unsigned short int data;
		uint16_t data;
};

//helper class to read FEMINOS raw data
//see event data frame of feminos user manual for more details
class DataLineFeminos{
	public:
		DataLineFeminos();
		~DataLineFeminos();
		//check if data corresponds to chip/channel id info
		bool is_info() const;
		//retrieve feminos id if is_frame_start()
		int get_card_ID() const;
		//retrieve AGET/AFTER chip id if is_info()
		int get_chip_ID() const;
		//retrieve channel id if is_info()
		int get_channel_ID() const;
		//check if data corresponds to sample id
		bool is_time() const;
		//retrieve first over threshold sample id if is_time() in zero suppress mode
		int get_time() const;
		//check if data corresponds to signal amplitude
		bool is_data() const;
		//retrieve signal amplitude if is_data()
		short get_data() const;
		//check if data corresponds to the End Of Frame word
		bool is_end_of_frame() const;
		//check if data corresponds to the end of event
		bool is_end_of_event() const;
		//check if data corresponds to the start of event
		bool is_event_start() const;
		//check if data corresponds to the Start Of Frame word
		bool is_frame_start() const;
//check if data corresponds to the Start Of Built Event word
		bool is_built_event_start() const;
		//check if data corresponds to the End Of Built Event word
		bool is_end_of_built_event() const;
		//unsigned short int data;
		uint16_t data;
};

#endif