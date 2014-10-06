#ifndef detector_h
#define detector_h
#include <string>
#include <vector>
#include "event.h"

using std::string;
using std::vector;

class Event;
class CM_Event;
class CM_Demux_Event;
class MG_Event;
class CosmicBenchEvent;
class Ray_2D;
class Ray;
class RayPair;
class Cluster;
class CM_Cluster;
class CM_Demux_Cluster;
class MG_Cluster;
class Detector;
class CM_Detector;
class MG_Detector;
class CosmicBench;
class T;

class Detector{
	public:
		//getters
		double get_z() const;
		bool get_is_X() const;
		bool get_is_up() const;
		bool get_is_ref() const;
		virtual string get_type() const = 0;
		double get_offset() const;
		bool get_direction() const;
		//seters
		void set_ClusTOTCut_Min(double cut);
		void set_ClusMaxSampleCut_Min(double cut);
		void set_ClusMaxSampleCut_Max(double cut);
		virtual ~Detector();
	protected:
		Detector();	
		Detector(const Detector& other);
		Detector& operator=(const Detector& other);
		Detector(double z_, bool is_X_, bool is_up_, bool is_ref_, double offset_, bool direction_);	
		double z; //altitude inside cosmic bench
		bool is_X;//coordinate measured by the detector
		bool is_up;//bloc (up|down) which the detector is part of
		bool is_ref; //used to test detectors
		double offset; //used for alignement
		bool direction; // direction of the axis
		//Detector dependent Cuts
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;

};

class CM_Detector: public Detector{
	friend class CM_Cluster;
	friend class CM_Demux_Cluster;
	public:	
		CM_Detector();
		CM_Detector(const CM_Detector& other);
		CM_Detector& operator=(const CM_Detector& other);
		CM_Detector(double z_, bool is_X_, bool is_up_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_);
		~CM_Detector();
		//CosMulti general charac
		static const double thinStripPitch; // distance between the middle of 2 adjacent thin strips
		static const double wideStripPitch; // distance between the middle of 2 adjacent wide strips
		//Cut setters
		void set_ClusMaxStripAmplCut_Min_Wide(double cut);
		void set_ClusSizeCut_Max_Wide(double cut);
		int get_cm_n_in_tree() const;
		bool get_use_thin_strip() const;
		string get_type() const;
	protected:
		int cm_n_in_tree;
		bool use_thin_strip;
		//Detector dependent Cuts
		double ClusMaxStripAmplCut_Min_Wide;
		double ClusSizeCut_Max_Wide;
};

class MG_Detector: public Detector{
	friend class MG_Cluster;
	public:
		MG_Detector();
		MG_Detector(const MG_Detector& other);
		MG_Detector& operator=(const MG_Detector& other);
		MG_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_);
		~MG_Detector();
		//MultiGen general charac
		static const double StripPitch; // distance between the middle of 2 adjacent strips
		//Cut setters
		void set_ClusSizeCut_Min(double cut);
		int get_mg_n_in_tree() const;
		string get_type() const;
	protected:
		int mg_n_in_tree;
		//Detector dependant cuts
		double ClusSizeCut_Min;
};

class CosmicBench{
	public:
		CosmicBench();
		~CosmicBench();
		//void add_MM(Detector * det);
		int get_CM_N() const;
		int get_MG_N() const;
		Detector * get_detector(unsigned int i) const;
	protected:
		vector<Detector*> detectors;
		int CM_N;
		int MG_N;
};

bool operator==(Detector const &det1, Detector const &det2);
bool operator!=(Detector const &det1, Detector const &det2);

#endif