#define ped_task_cpp

#include "task/ped_task.h"

#include "detector.h"

/*
Ped_Task::Ped_Task(map<Tomography::det_type,vector<vector<float> > > ped_): Typed_Task<raw_data>(){
	ped = ped_;
	next_task = NULL;
}
*/

template<typename A>
Ped_Task::Ped_Task(map<Tomography::det_type,vector<vector<A> > > ped_, Typed_Task<ped_data> * next_task_): Typed_Task<raw_data>(){
	for(typename map<Tomography::det_type,vector<vector<A> > >::iterator type_it=ped_.begin();type_it!=ped_.end();++type_it){
		ped[type_it->first] = vector<vector<double> >((type_it->second).size(),vector<double>((type_it->second).front().size(),0));
		vector<vector<double> >::iterator final_det_it = ped[type_it->first].begin();
		for(typename vector<vector<A> >::iterator det_it=(type_it->second).begin();det_it!=(type_it->second).end();++det_it){
			vector<double>::iterator final_strip_it = final_det_it->begin();
			for(typename vector<A>::iterator strip_it=det_it->begin();strip_it!=det_it->end();++strip_it){
				*final_strip_it = *strip_it;
				++final_strip_it;
			}
			++final_det_it;
		}
	}
	next_task = next_task_;
}
template Ped_Task::Ped_Task(map<Tomography::det_type,vector<vector<float> > > ped_, Typed_Task<ped_data> * next_task_);

template<>
Ped_Task::Ped_Task(map<Tomography::det_type,vector<vector<double> > > ped_, Typed_Task<ped_data> * next_task_): Typed_Task<raw_data>(){
	ped = ped_;
	next_task = next_task_;
}
Ped_Task::~Ped_Task(){
	delete next_task;
}
bool Ped_Task::do_task(){
	raw_data * current_data = get_next_data();
	if(current_data->Nevent < 0){
		delete current_data;
		return false;
	}
	ped_data * sub_data = new ped_data();
	sub_data->Nevent = current_data->Nevent;
	sub_data->evttime = current_data->evttime;
	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator type_it=(current_data->strip_data).begin();type_it!=(current_data->strip_data).end();++type_it){
		map<Tomography::det_type,vector<vector<double> > >::const_iterator ped_it = ped.find(type_it->first);
		vector<vector<double> >::const_iterator ped_det_it = (ped_it->second).begin();
		(sub_data->strip_data)[type_it->first] = vector<vector<vector<double> > >((type_it->second).size(),vector<vector<double> >(Tomography::Static_Detector[type_it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
		vector<vector<vector<double> > >::iterator sub_det_it = (sub_data->strip_data)[type_it->first].begin();
		for(vector<vector<vector<float> > >::iterator det_it=(type_it->second).begin();det_it!=(type_it->second).end();++det_it){
			vector<double>::const_iterator ped_strip_it = ped_det_it->begin();
			vector<vector<double> >::iterator sub_strip_it = sub_det_it->begin();
			for(vector<vector<float> >::iterator strip_it=det_it->begin();strip_it!=det_it->end();++strip_it){
				vector<double>::iterator sub_sample_it = sub_strip_it->begin();
				for(vector<float>::iterator sample_it=strip_it->begin();sample_it!=strip_it->end();++sample_it){
					*sub_sample_it = (*sample_it) - (*ped_strip_it);
					//*sample_it -= *ped_strip_it;
					++sub_sample_it;
				}
				++ped_strip_it;
				++sub_strip_it;
			}
			++ped_det_it;
			++sub_det_it;
		}
	}
	next_task->push_next_data(sub_data);
	delete current_data;
	return true;
}
void Ped_Task::update_task_list() const{
	add_task(next_task);
}
bool Ped_Task::is_queueable() const{
	return true;
}
/*
Corr_Task::Corr_Task(): Typed_Task<ped_data>(){
	next_task = NULL;
}
*/
Corr_Task::Corr_Task(Typed_Task<corr_data> * next_task_): Typed_Task<ped_data>(){
	next_task= next_task_;
}
Corr_Task::~Corr_Task(){
	delete next_task;
}
bool Corr_Task::do_task(){
	ped_data * current_data = get_next_data();
	if(current_data->Nevent < 0){
		delete current_data;
		return false;
	}
	corr_data * sub_data = new corr_data();
	sub_data->Nevent = current_data->Nevent;
	sub_data->evttime = current_data->evttime;
	for(map<Tomography::det_type,vector<vector<vector<double> > > >::iterator it = (current_data->strip_data).begin();it!=(current_data->strip_data).end();++it){
		(sub_data->strip_data)[it->first] = vector<vector<vector<double> > >((it->second).size(),vector<vector<double> >(Tomography::Static_Detector[it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
		vector<vector<vector<double> > >::iterator sub_det_it = (sub_data->strip_data)[it->first].begin();
		for(vector<vector<vector<double> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				for(int det_div=0;det_div<(Tomography::Static_Detector[it->first]->get_CMN_div());det_div++){
					int strip_nb = (Tomography::Static_Detector[it->first]->get_Nchannel())/(Tomography::Static_Detector[it->first]->get_CMN_div()) + (Tomography::Static_Detector[it->first]->get_Nchannel())%(Tomography::Static_Detector[it->first]->get_CMN_div());
					int strip_offset = det_div*strip_nb;
					vector<double> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						current_sample[j] = (*jt)[j+strip_offset][k];
					}
					sort(current_sample.begin(),current_sample.end());
					double median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						(*sub_det_it)[j+strip_offset][k] = (*jt)[j+strip_offset][k] - median;
						//(*jt)[j+strip_offset][k] -= median;
					}
				}
			}
			++sub_det_it;
		}
	}
	next_task->push_next_data(sub_data);
	delete current_data;
	return true;
}
void Corr_Task::update_task_list() const{
	add_task(next_task);
}
bool Corr_Task::is_queueable() const{
	return true;
}
/*
Ped_Corr_Task::Ped_Corr_Task(map<Tomography::det_type,vector<vector<float> > > ped_): Typed_Task<raw_data>(){
	ped = ped_;
	next_task = NULL;
}
*/
template<typename A>
Ped_Corr_Task::Ped_Corr_Task(map<Tomography::det_type,vector<vector<A> > > ped_, Typed_Task<corr_data> * next_task_): Typed_Task<raw_data>(){
	for(typename map<Tomography::det_type,vector<vector<A> > >::iterator type_it=ped_.begin();type_it!=ped_.end();++type_it){
		ped[type_it->first] = vector<vector<double> >((type_it->second).size(),vector<double>((type_it->second).front().size(),0));
		vector<vector<double> >::iterator final_det_it = ped[type_it->first].begin();
		for(typename vector<vector<A> >::iterator det_it=(type_it->second).begin();det_it!=(type_it->second).end();++det_it){
			vector<double>::iterator final_strip_it = final_det_it->begin();
			for(typename vector<A>::iterator strip_it=det_it->begin();strip_it!=det_it->end();++strip_it){
				*final_strip_it = *strip_it;
				++final_strip_it;
			}
			++final_det_it;
		}
	}
	next_task = next_task_;
}
template Ped_Corr_Task::Ped_Corr_Task(map<Tomography::det_type,vector<vector<float> > > ped_, Typed_Task<corr_data> * next_task_);
template Ped_Corr_Task::Ped_Corr_Task(map<Tomography::det_type,vector<vector<double> > > ped_, Typed_Task<corr_data> * next_task_);

/*
template<>
Ped_Corr_Task::Ped_Corr_Task(map<Tomography::det_type,vector<vector<double> > > ped_, Typed_Task<corr_data> * next_task_): Typed_Task<raw_data>(){
	ped = ped_;
	next_task = next_task_;
}
*/
Ped_Corr_Task::~Ped_Corr_Task(){
	delete next_task;
}
bool Ped_Corr_Task::do_task(){
	raw_data * current_data = get_next_data();
	if(current_data->Nevent < 0){
		delete current_data;
		return false;
	}
	corr_data * sub_data = new corr_data();
	sub_data->Nevent = current_data->Nevent;
	sub_data->evttime = current_data->evttime;

	for(map<Tomography::det_type,vector<vector<vector<float> > > >::iterator it = (current_data->strip_data).begin();it!=(current_data->strip_data).end();++it){
		if((it->second).empty()) continue;
		(sub_data->strip_data)[it->first] = vector<vector<vector<double> > >((it->second).size(),vector<vector<double> >(Tomography::Static_Detector[it->first]->get_Nchannel(),vector<double>(Tomography::get_instance()->get_Nsample(),0)));
		map<Tomography::det_type,vector<vector<double> > >::const_iterator ped_it = ped.find(it->first);
		vector<vector<double> >::const_iterator ped_jt = (ped_it->second).begin();
		vector<vector<vector<double> > >::iterator out_it = (sub_data->strip_data)[it->first].begin();
		for(vector<vector<vector<float> > >::iterator jt = (it->second).begin();jt!=(it->second).end();++jt){
			for(int k=0;k<Tomography::get_instance()->get_Nsample();k++){
				for(int det_div=0;det_div<(Tomography::Static_Detector[it->first]->get_CMN_div());det_div++){
					int strip_nb = (Tomography::Static_Detector[it->first]->get_Nchannel())/(Tomography::Static_Detector[it->first]->get_CMN_div()) + (Tomography::Static_Detector[it->first]->get_Nchannel())%(Tomography::Static_Detector[it->first]->get_CMN_div());
					int strip_offset = det_div*strip_nb;
					vector<float> current_sample(strip_nb,0);
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						current_sample[j] = (*jt)[j+strip_offset][k] - (*ped_jt)[j+strip_offset];
					}
					sort(current_sample.begin(),current_sample.end());
					float median = current_sample[strip_nb/2];
					for(int j=0;(j<strip_nb && (j+strip_offset)<(Tomography::Static_Detector[it->first]->get_Nchannel()));j++){
						(*out_it)[j+strip_offset][k] = (*jt)[j+strip_offset][k] - median - (*ped_jt)[j+strip_offset];
					}
				}
			}
			++ped_jt;
			++out_it;
		}
	}
	next_task->push_next_data(sub_data);
	delete current_data;
	return true;
}
void Ped_Corr_Task::update_task_list() const{
	add_task(next_task);
}
bool Ped_Corr_Task::is_queueable() const{
	return true;
}
