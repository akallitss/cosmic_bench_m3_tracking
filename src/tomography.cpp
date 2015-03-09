#include "tomography.h"
#include <map>

using std::map;

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

template<typename T,typename R>
ostream& operator<<(ostream& os, const map<T,R>& map_){
	if(map_.size()<1){
		os << "[]";
		return os;
	}
	os << "[ ";
	typename map<T,R>::const_iterator it=map_.begin();
	typename map<T,R>::const_iterator jt = map_.end();
	jt--;
	while(it!=jt){
		os << it->first << " -> " << it->second << " ; ";
		++it;
	}
	os << it->first << " -> " << it->second << " ]";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os,const vector<T>& vec_){
	if(vec_.empty()){
		os << "{}";
		return os;
	}
	os << "{ ";
	typename vector<T>::const_iterator it = vec_.begin();
	typename vector<T>::const_iterator jt = vec_.end();
	jt--;
	while(it!=jt){
		os << *it << " ; ";
		++it;
	}
	os << *it << " }";
	return os;
}

//template ostream& operator<<(ostream& os, const map<int,int>& map_);
template ostream& operator<<(ostream& os, const map<int,double>& map_);
template ostream& operator<<(ostream& os, const map<double,int>& map_);
template ostream& operator<<(ostream& os, const map<int,int>& map_);
template ostream& operator<<(ostream& os, const map<bool,map<int,int> >& map_);

Tomography::elec_type Tomography::str_to_elec(string str){
	elec_type return_value = unknown_elec;
	if(str == "dream") return_value = Dream;
	else if(str == "feminos") return_value = Feminos;
	return return_value;
}