#define dataline_cpp
#include "dataline.h"
#include "mygblink.h"
#include <arpa/inet.h>

DataLineDream::DataLineDream(){
	data = 0;
}
DataLineDream::~DataLineDream(){

}
void DataLineDream::ntohs_(){
	data = ntohs(data);
}
bool DataLineDream::is_final_trailer() const{
	return GET_TYPE(data)==7;
}
bool DataLineDream::is_first_line() const{
	return GET_TYPE(data)==3;
}
bool DataLineDream::is_final_header() const{
	return GET_TYPE(data)==2;
}
bool DataLineDream::is_data() const{
	return GET_TYPE(data)==0;
}
bool DataLineDream::is_channel_ID() const{
	return GET_TYPE(data)==1;
}
int DataLineDream::get_dream_ID() const{
	return GET_DREAMID(data);
}
int DataLineDream::get_channel_ID() const{
	return GET_CHANNELID(data);
}
float DataLineDream::get_data() const{
	return GET_DATA(data);
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
float DataLineFeminos::get_data() const{
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