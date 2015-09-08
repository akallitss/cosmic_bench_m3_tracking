#ifndef event_h
#define event_h
#include <string>
#include <vector>
#include <map>

#include "ray.h"
#include "tomography.h"

using std::string;
using std::vector;
using std::map;

class Cluster;
class Detector;
class CM_Detector;
class MG_Detector;
class MGv2_Detector;
class CosmicBench;
class Tanalyse_R;
class Analyse;

class TH1D;
class TCanvas;

class Event{
	friend class Analyse;
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
		virtual void HoughCluster() = 0;
		void do_cuts();
		virtual void set_strip_ampl(vector<vector<double> > strip_ampl_) = 0;
		vector<Cluster*> get_clusters() const;
		Detector * get_det() const;
		double get_StripPitch() const;
		virtual TH1D * get_ampl_hist() const = 0;
		virtual Event * Clone() const = 0;
	protected:
		struct StripInfo {
			double MaxAmpl;
			int MaxSample;
			int TOT;
			double Time;
			bool signal_sample[Tomography::Nsample];
		};
		Event(int evn_ = -1);
		Event(const Event& other);
		Event& operator=(const Event& other);
		Event(Tanalyse_R * treeObject,const Detector * const det,long entry = -1);
		Event(const Detector * const detector_,int evn_);
		int evn;
		Tomography::det_type type;
		bool has_spark;
		vector<vector<double> > strip_ampl;
		vector<Cluster*> clusters;
		Detector * detector;
};

class CM_Event: public Event{
	friend class CM_Demux_Event;
	public:
		CM_Event();
		CM_Event(const CM_Event& other);
		CM_Event& operator=(const CM_Event& other);
		CM_Event(Tanalyse_R * treeObject,const CM_Detector * const det, long entry = -1);
		CM_Event(const CM_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_);
		~CM_Event();
		void MultiCluster();
		void HoughCluster();
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		TH1D * get_ampl_hist() const;
		Event * Clone() const;
};

class CM_Demux_Event: public Event{
	public:
		CM_Demux_Event();
		CM_Demux_Event(const CM_Demux_Event& other);
		CM_Demux_Event& operator=(const CM_Demux_Event& other);
		CM_Demux_Event(const CM_Event& rawEvent);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~CM_Demux_Event();
		void MultiCluster();
		void HoughCluster();
		TH1D * get_ampl_hist() const;
		Event * Clone() const;
};

class MG_Event: public Event{
	public:
		MG_Event();
		MG_Event(const MG_Event& other);
		MG_Event& operator=(const MG_Event& other);
		MG_Event(Tanalyse_R * treeObject,const MG_Detector * const det, long entry = -1, bool use_srf_ = false);
		MG_Event(const MG_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, bool use_srf_ = false);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~MG_Event();
		void MultiCluster();
		void HoughCluster();
		TH1D * get_ampl_hist() const;
		Event * Clone() const;
	protected:
		bool use_srf;
};

class MGv2_Event: public Event{
	public:
		MGv2_Event();
		MGv2_Event(const MGv2_Event& other);
		MGv2_Event& operator=(const MGv2_Event& other);
		MGv2_Event(Tanalyse_R * treeObject,const MGv2_Detector * const det, long entry = -1, bool use_srf_ = false);
		MGv2_Event(const MGv2_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, bool use_srf_ = false);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~MGv2_Event();
		void MultiCluster();
		void HoughCluster();
		TH1D * get_ampl_hist() const;
		Event * Clone() const;
	protected:
		bool use_srf;
};

//Group events objects of a same event
class CosmicBenchEvent{
	friend class Analyse;
	friend class Carac;
	public:
		CosmicBenchEvent();
		CosmicBenchEvent(const CosmicBenchEvent& other);
		CosmicBenchEvent& operator=(const CosmicBenchEvent& other);
		CosmicBenchEvent(const CosmicBench * const detectors, Tanalyse_R * treeObject, long entry = -1);
		CosmicBenchEvent(const CosmicBench * const detectors, vector<Event*> events_);
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