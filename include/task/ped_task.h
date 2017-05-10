#ifndef ped_task_h
#define ped_task_h

#include "MT_tomography.h"
#include "tomography.h"

#include <vector>
#include <map>

using std::vector;
using std::map;

class Ped_Task: public Typed_Task<raw_data>{
	public:
		//Ped_Task(map<Tomography::det_type,vector<vector<double> > > ped_);
		template<typename A> Ped_Task(map<Tomography::det_type,vector<vector<A> > > ped_, Typed_Task<ped_data> * next_task_);
		~Ped_Task();
		bool do_task();
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		map<Tomography::det_type,vector<vector<double> > > ped;
		Typed_Task<ped_data> * next_task;
};
class Corr_Task: public Typed_Task<ped_data>{
	public:
		//Corr_Task();
		Corr_Task(Typed_Task<corr_data> * next_task_);
		~Corr_Task();
		bool do_task();
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		Typed_Task<corr_data> * next_task;
};
class Ped_Corr_Task: public Typed_Task<raw_data>{
	public:
		//Ped_Corr_Task(map<Tomography::det_type,vector<vector<double> > > ped_);
		template<typename A> Ped_Corr_Task(map<Tomography::det_type,vector<vector<A> > > ped_, Typed_Task<corr_data> * next_task_);
		~Ped_Corr_Task();
		bool do_task();
		void update_task_list() const;
		bool is_queueable() const;
	protected:
		map<Tomography::det_type,vector<vector<double> > > ped;
		Typed_Task<corr_data> * next_task;
};

#endif
