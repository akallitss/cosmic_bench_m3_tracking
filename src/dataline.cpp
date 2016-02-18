#define dataline_cpp
#include "dataline.h"
#include <arpa/inet.h>

DataLineDream::DataLineDream(){
	data = 0;
}
DataLineDream::DataLineDream(uint16_t input): data(input){
	
}
DataLineDream::~DataLineDream(){

}
void DataLineDream::ntohs_(){
	data = ntohs(data);
}
bool DataLineDream::is_final_trailer() const{
	return (((data) & 0x7000)>>12)==7;
}
bool DataLineDream::is_data_trailer() const{
	return (((data) & 0x6000)>>13)==2;
}
bool DataLineDream::is_first_line() const{
	return (((data) & 0x7000)>>12)==3;
}
bool DataLineDream::is_Feu_header() const{
	return (((data) & 0x7000)>>12)==6;
}
bool DataLineDream::is_data_header() const{
	return (((data) & 0x6000)>>13)==1;
}
bool DataLineDream::is_data() const{
	return (((data) & 0x7000)>>12)==0;
}
bool DataLineDream::is_data_zs() const{
	return (((data) & 0x6000)>>13)==0;
}
bool DataLineDream::is_channel_ID() const{
	return (((data) & 0x7000)>>12)==1;
}
bool DataLineDream::is_EOE() const{
	return (((data) & 0x7800)>>11)==15;
	//return (((data) & 0x0800)>>11)==1;
}
bool DataLineDream::get_zs_mode() const{
	return (((data) & 0x400)>>10);
}
int DataLineDream::get_sample_ID() const{
	return (((data) & 0xFF8)>>3);
}
int DataLineDream::get_Feu_ID() const{
	return (((data) & 0xFF));
}
int DataLineDream::get_dream_ID() const{
	return (((data) & 0xE00)>>9);
}
int DataLineDream::get_channel_ID() const{
	return (((data) & 0x3F));
}
short DataLineDream::get_data() const{
	return (((data) & 0xFFF));
}

DataLineFeminos::DataLineFeminos(){
	data = 0;
}
DataLineFeminos::~DataLineFeminos(){

}
bool DataLineFeminos::is_info() const{
	return ((data & 0xC000)>>14)==3;
}
int DataLineFeminos::get_card_ID() const{
	return (data & 0x3E00)>>9;
}
int DataLineFeminos::get_chip_ID() const{
	return (data & 0x180)>>7;
}
int DataLineFeminos::get_channel_ID() const{
	return (data & 0x7F);
}
bool DataLineFeminos::is_time() const{
	return ((data & 0xFE00)>>9)==7;
}
int DataLineFeminos::get_time() const{
	return (data & 0x1FF);
}
bool DataLineFeminos::is_data() const{
	return ((data & 0x3000)>>12)==3;
}
short DataLineFeminos::get_data() const{
	return (data & 0xFFF);
}
bool DataLineFeminos::is_end_of_frame() const{
	return data==15;
}
bool DataLineFeminos::is_end_of_event() const{
	return ((data & 0xFFF0)>>4)==14;
}
bool DataLineFeminos::is_event_start() const{
	return ((data & 0xFFF0)>>4)==15;
}
bool DataLineFeminos::is_frame_start() const{
	return ((data & 0xFE00)>>9)==4;
}
bool DataLineFeminos::is_built_event_start() const{
	return data == 9;
}
bool DataLineFeminos::is_end_of_built_event() const{
	return data == 8;
}