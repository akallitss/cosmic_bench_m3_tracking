#ifndef event_h
#define event_h
#include <string>
#include <vector>
#include <map>
#include "cluster.h"
#include "ray.h"
#include "T.h"
#include "detector.h"
#include "analyse.h"

using std::string;
using std::vector;
using std::map;

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

class Event{
	public:
		int get_evn() const;
		string get_type() const;
		int get_n_in_tree() const;
		bool get_is_ref() const;
		double get_z() const;
		int get_NClus() const;
		virtual ~Event();
		virtual void MultiCluster() = 0;
	protected:
		Event();
		Event(const Event& other);
		Event& operator=(const Event& other);
		Event(T * treeObject, bool use_srf_,int entry = -1);
		int evn;
		double z;
		string type;
		int n_in_tree;
		bool has_spark;
		bool is_ref;
		int NClus;
		bool use_srf;
		vector<vector<double> > strip_ampl;
		struct StripInfo {
			double MaxAmpl;
			int MaxSample;
			int TOT;
			double Time;
		};
};

class CM_Event: public Event{
	friend class CM_Demux_Event;
	public:
		CM_Event();
		CM_Event(const CM_Event& other);
		CM_Event& operator=(const CM_Event& other);
		CM_Event(T * treeObject,CM_Detector * det, bool use_srf_, int entry = -1);
		CM_Event(CM_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_);
		~CM_Event();
		void MultiCluster();
		vector<CM_Cluster> get_clusters() const;
	protected:
		CM_Detector detector;
		bool use_thin_strip;
		vector<CM_Cluster> clusters;
};

class CM_Demux_Event: public Event{
	friend class Analyse;
	public:
		CM_Demux_Event();
		CM_Demux_Event(const CM_Demux_Event& other);
		CM_Demux_Event& operator=(const CM_Demux_Event& other);
		CM_Demux_Event(const CM_Event& rawEvent);
		vector<CM_Demux_Cluster> get_clusters() const;
		~CM_Demux_Event();
		void MultiCluster();
	protected:
		vector<CM_Demux_Cluster> clusters;
};

class MG_Event: public Event{
	friend class Analyse;
	public:
		MG_Event();
		MG_Event(const MG_Event& other);
		MG_Event& operator=(const MG_Event& other);
		MG_Event(T * treeObject,MG_Detector * det, bool use_srf_, int entry = -1);
		MG_Event(MG_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_);
		~MG_Event();
		void MultiCluster();
		void HoughCluster();
		vector<MG_Cluster> get_clusters() const;
	protected:
		MG_Detector detector;
		vector<MG_Cluster> clusters;
};

//Group events objects of a same event
class CosmicBenchEvent{
	friend class Analyse;
	public:
		CosmicBenchEvent();
		CosmicBenchEvent(const CosmicBenchEvent& other);
		CosmicBenchEvent& operator=(const CosmicBenchEvent& other);
		CosmicBenchEvent(CosmicBench * detectors, T * treeObject, bool use_srf_, int entry = -1);
		~CosmicBenchEvent();
		void createPairs();
		RayPair get_rayPair(unsigned int i) const;
		unsigned int get_rayPairs_N() const;
		unsigned int get_event_N() const;
		unsigned int get_clus_N() const;
		unsigned int get_clus_N_by_det(Detector * det) const;
		void Demux_CM();
		vector<Ray> get_absorption_rays();
		static vector<map<double,int> > combinaisons(map<double,int> sizes);
		//void MultiCluster();
	protected:
		int evn;
		vector<Event*> events;
		//After processing
		vector<RayPair> rayPairs;
};

#endif