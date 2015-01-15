#include "tomography.h"

ostream& Tomography::operator<<(ostream& os, const det_type& det){
	switch(det){
		case CM : os << "CM"; break;
		case CM_Demux : os << "CM_Demux"; break;
		case MG : os << "MG"; break;
		default : os << "unknown det";
	}
	return os;
}
ostream& Tomography::operator<<(ostream& os, const strip_type& strip){
	switch(strip){
		case Wide : os << "Wide"; break;
		case Thin : os << "Thin"; break;
		case Demux : os << "Demux"; break;
		default : os << "unknown strip type";
	}
	return os;
}
ostream& Tomography::operator<<(ostream& os, const elec_type& elec){
	switch(elec){
		case Dream : os << "Dream"; break;
		case Feminos : os << "Feminos"; break;
		default : os << "unknown electronic type";
	}
	return os;
}
Tomography::elec_type Tomography::str_to_elec(string str){
	elec_type return_value;
	if(str == "dream") return_value = Dream;
	else if(str == "feminos") return_value = Feminos;
	return return_value;
}