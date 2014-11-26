#ifndef dataline_h
#define dataline_h

class DataLineDream{
	public:
		DataLineDream();
		~DataLineDream();
		void ntohs_();
		bool is_final_trailer() const;
		bool is_first_line() const;
		bool is_data() const;
		bool is_channel_ID() const;
		int get_channel_ID() const;
		int get_dream_ID() const;
		float get_data() const;
		unsigned short int data;
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
		float get_data() const;
		bool is_end_of_frame() const;
		bool is_end_of_event() const;
		bool is_event_start() const;
		bool is_frame_start() const;
		bool is_built_event_start() const;
		bool is_end_of_built_event() const;
		unsigned short int data;
};

#endif