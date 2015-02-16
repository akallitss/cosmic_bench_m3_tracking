#ifndef header_h
#define header_h

class HeaderC{
	public:
		HeaderC();
		~HeaderC();
		void ntohs_();
		int get_sampleIndex() const;
		bool check_type() const;
		bool get_ped_mode() const;
		bool get_cms_mode() const;
		bool get_zs_mode() const;
		unsigned short int param;
		unsigned short int eventid;
		unsigned short int timeStamp;
		unsigned short int sampleIndex;
};

#endif