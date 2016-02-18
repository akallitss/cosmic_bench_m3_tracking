#ifndef write_rays_task_h
#define write_rays_task_h

#include "MT_tomography.h"

class Tray;

class Write_Rays_Task: public Output_Task<ray_data>{
	public:
		Write_Rays_Task(Tray * writer_, double z_up_, double z_down_);
		//Write_Rays_Task(Tray * writer_, double z_up_, double z_down_, Typed_Task<ray_data> * next_task_);
		~Write_Rays_Task();
		bool do_task();
		bool can_exec() const;
		void update_task_list() const;
	protected:
		Tray * writer;
		double z_up;
		double z_down;
};
/*
class Write_Raypairs_Task: public IO_Task{
	public:
		Write_Raypairs_Task();
		Write_Raypairs_Task(Task * next_task_);
		~Write_Raypairs_Task();
		bool do_task();
		bool can_exec();
		void update_task_list();
};
*/

#endif