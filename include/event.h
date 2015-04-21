#ifndef event_h
#define event_h
#include <string>
#include <vector>
#include <map>
#include "cluster.h"
#include "ray.h"
#include "Tanalyse_R.h"
#include "detector.h"
#include "analyse.h"
#include <TH1D.h>
#include <TCanvas.h>
#include "tomography.h"

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
class Tanalyse_R;

class Event{
	public:
		int get_evn() const;
		Tomography::det_type get_type() const;
		int get_n_in_tree() const;
		bool get_is_ref() const;
		double get_z() const;
		int get_NClus() const;
		bool get_is_X() const;
		virtual ~Event();
		virtual void MultiCluster() = 0;
		virtual void do_cuts() = 0;
		virtual void set_strip_ampl(vector<vector<double> > strip_ampl_) = 0;
	protected:
		Event();
		Event(const Event& other);
		Event& operator=(const Event& other);
		Event(Tanalyse_R * treeObject, bool use_srf_,long entry = -1);
		int evn;
		double z;
		bool is_X;
		Tomography::det_type type;
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
			bool signal_sample[Tomography::Nsample];
		};
};

class CM_Event: public Event{
	friend class CM_Demux_Event;
	public:
		CM_Event();
		CM_Event(const CM_Event& other);
		CM_Event& operator=(const CM_Event& other);
		CM_Event(Tanalyse_R * treeObject,CM_Detector * det, bool use_srf_, long entry = -1);
		CM_Event(CM_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_, int evn_);
		~CM_Event();
		void MultiCluster();
		void do_cuts();
		vector<CM_Cluster> get_clusters() const;
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
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
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~CM_Demux_Event();
		void MultiCluster();
		void do_cuts();
	protected:
		vector<CM_Demux_Cluster> clusters;
};

class MG_Event: public Event{
	friend class Analyse;
	public:
		MG_Event();
		MG_Event(const MG_Event& other);
		MG_Event& operator=(const MG_Event& other);
		MG_Event(Tanalyse_R * treeObject,MG_Detector * det, bool use_srf_, long entry = -1);
		MG_Event(MG_Detector detector_, vector<vector<double> > strip_ampl_, bool use_srf_, int evn_);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~MG_Event();
		void MultiCluster();
		void HoughCluster();
		void do_cuts();
		vector<MG_Cluster> get_clusters() const;
		TH1D * get_ampl_hist() const;
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
		CosmicBenchEvent(CosmicBench * detectors, Tanalyse_R * treeObject, bool use_srf_, long entry = -1);
		CosmicBenchEvent(CosmicBench * detectors, vector<Event*> events_);
		~CosmicBenchEvent();
		void createPairs();
		void EventDisplay(TCanvas * c1 = 0);
		RayPair get_rayPair(unsigned int i) const;
		unsigned int get_rayPairs_N() const;
		unsigned int get_event_N() const;
		unsigned int get_clus_N() const;
		unsigned int get_clus_N_by_det(Detector * det) const;
		void Demux_CM();
		void do_cuts();
		vector<Ray> get_absorption_rays(double chiSquare_threshold = Tomography::chisquare_threshold);
		vector<Ray> get_hough_rays(double chiSquare_threshold = Tomography::chisquare_threshold);
		template<typename T>
		static vector<map<T,int> > combinaisons(map<T,int> sizes, bool allow_drop = false);
		//void MultiCluster();
	protected:
		int evn;
		vector<Event*> events;
		//After processing
		vector<RayPair> rayPairs;
};

#endif