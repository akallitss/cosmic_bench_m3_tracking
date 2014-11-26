#define header_cpp
#include "header.h"
#include "mygblink.h"
#include <arpa/inet.h>


HeaderC::HeaderC(){
	param = 0;
	eventid = 0;
	timeStamp = 0;
	sampleIndex = 0;
}
HeaderC::~HeaderC(){

}
void HeaderC::ntohs_(){
	param = ntohs(param);
	eventid = ntohs(eventid);
	timeStamp = ntohs(timeStamp);
	sampleIndex = ntohs(sampleIndex);
}
int HeaderC::get_sampleIndex() const{
	return GET_SAMPLEINDEX(sampleIndex);
}
bool HeaderC::get_ped_mode() const{
	return GET_PED(param);
}
bool HeaderC::get_cms_mode() const{
	return GET_CMS(param);
}
bool HeaderC::get_zs_mode() const{
	return GET_ZS(param);
}
bool HeaderC::check_type() const{
	return !(GET_TYPE(param)!=6 && GET_TYPE(eventid)!=6 && GET_TYPE(timeStamp)!=6 && GET_TYPE(sampleIndex)!=6);
}