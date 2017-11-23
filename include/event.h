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

//abstract class to store a whole event data in a single detector
class Event{
	friend class Analyse;
	public:
		//retrieve the event id
		int get_evn() const;
		//retrieve the event timestamp
		double get_evttime() const;
		//retrieve the type of the detector measuring this event
		Tomography::det_type get_type() const;
		//retrieve the detector id measuring this event
		int get_n_in_tree() const;
		//check if the detector measuring this event is included in tracking
		bool get_is_ref() const;
		//retrieve the altitude of the detector measuring this event
		double get_z() const;
		//retrieve the number of reconstructed cluster inside the detector for this event
		int get_NClus() const;
		//check if the detector measuring this event gives information about the X or Y coordinate
		bool get_is_X() const;
		virtual ~Event();
		//do the clustering using the strip amplitudes
		virtual void MultiCluster() = 0;
		//do the clustering using the strip amplitudes with the convolution method
		virtual void ConvCluster() = 0;
		//do the clustering using the strip amplitudes with the Hough method
		virtual void HoughCluster(int hole_nb) = 0;
		//delete clusters which do not pass the detector specific cuts
		void do_cuts();
		//set the strip amplitudes to do the clustering
		virtual void set_strip_ampl(vector<vector<double> > strip_ampl_) = 0;
		//retrieve the built clusters
		vector<Cluster*> get_clusters() const;
		//retrieve the associated detector
		Detector * get_det() const;
		//retrieve the strip pitch of the detector measuring this event
		double get_StripPitch() const;
		//build the histogram of maximum amplitude by strip
		virtual TH1D * get_ampl_hist() const = 0;
		//build the histogram of Time Over Threshold by strip
		virtual TH1D * get_TOT_hist() const = 0;
		//make a copy of this event
		virtual Event * Clone() const = 0;
	protected:
		//helper struct to store strip informations
		struct StripInfo {
			double MaxAmpl;
			int MaxSample;
			int TOT;
			double Time;
			bool signal_sample[Tomography::Max_Nsample];
		};
		Event(int evn_ = -1);
		//copy constructor
		Event(const Event& other);
		//copy assignment
		Event& operator=(const Event& other);
		//build event for the given detector by using the treeObject data at the given entry
		Event(Tanalyse_R * treeObject,const Detector * const det,long entry);
		//build event for the given detector by using the treeObject data at the current loaded entry
		Event(const Tanalyse_R * const treeObject,const Detector * const det);
		//build event for the given detector and giving explicitly the event id and timestamp
		Event(const Detector * const detector_,int evn_, double evttime_);
		int evn;
		double evttime;
		Tomography::det_type type;
		bool has_spark;
		vector<vector<double> > strip_ampl;
		vector<Cluster*> clusters;
		Detector * detector;
};

//event implementation for a single full connector dummy detector
class dummy_Event: public Event{
	public:
		dummy_Event();
		dummy_Event(const dummy_Event& other);
		dummy_Event& operator=(const dummy_Event& other);
		~dummy_Event();
		void MultiCluster();
		void ConvCluster();
		void HoughCluster(int hole_nb);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		TH1D * get_ampl_hist() const;
		TH1D * get_TOT_hist() const;
		Event * Clone() const;
};

//implementation for CosMulti detector to store still multiplexed informations
class CM_Event: public Event{
	friend class CM_Demux_Event;
	public:
		CM_Event();
		CM_Event(const CM_Event& other);
		CM_Event& operator=(const CM_Event& other);
		CM_Event(Tanalyse_R * treeObject,const CM_Detector * const det, long entry);
		CM_Event(const Tanalyse_R * const treeObject,const CM_Detector * const det);
		CM_Event(const CM_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_);
		~CM_Event();
		void MultiCluster();
		void ConvCluster();
		void HoughCluster(int hole_nb);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		TH1D * get_ampl_hist() const;
		TH1D * get_TOT_hist() const;
		Event * Clone() const;
};

//implementation for CosMulti detector to store demultiplexed informations
class CM_Demux_Event: public Event{
	public:
		CM_Demux_Event();
		CM_Demux_Event(const CM_Demux_Event& other);
		CM_Demux_Event& operator=(const CM_Demux_Event& other);
		CM_Demux_Event(const CM_Event& rawEvent);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~CM_Demux_Event();
		void MultiCluster();
		void ConvCluster();
		void HoughCluster(int hole_nb);
		TH1D * get_ampl_hist() const;
		TH1D * get_TOT_hist() const;
		Event * Clone() const;
};

//implementation to store MultiGen V1 informations (50x50cm^2; 1024 strips)
class MG_Event: public Event{
	public:
		MG_Event();
		MG_Event(const MG_Event& other);
		MG_Event& operator=(const MG_Event& other);
		MG_Event(Tanalyse_R * treeObject,const MG_Detector * const det, long entry);
		MG_Event(const Tanalyse_R * const treeObject,const MG_Detector * const det);
		MG_Event(const MG_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~MG_Event();
		void MultiCluster();
		void ConvCluster();
		void HoughCluster(int hole_nb);
		TH1D * get_ampl_hist() const;
		TH1D * get_TOT_hist() const;
		Event * Clone() const;
};

//implementation to store MultiGen V2 informations (50x50cm^2; 1037 strips)
class MGv2_Event: public Event{
	public:
		MGv2_Event();
		MGv2_Event(const MGv2_Event& other);
		MGv2_Event& operator=(const MGv2_Event& other);
		MGv2_Event(Tanalyse_R * treeObject,const MGv2_Detector * const det, long entry);
		MGv2_Event(const Tanalyse_R * const treeObject,const MGv2_Detector * const det);
		MGv2_Event(const MGv2_Detector * const detector_, vector<vector<double> > strip_ampl_, int evn_, double evttime_);
		void set_strip_ampl(vector<vector<double> > strip_ampl_);
		~MGv2_Event();
		void MultiCluster();
		void ConvCluster();
		void HoughCluster(int hole_nb);
		TH1D * get_ampl_hist() const;
		TH1D * get_TOT_hist() const;
		Event * Clone() const;
};

//Group events objects of a same event for the whole cosmic bench
class CosmicBenchEvent{
	friend class Analyse;
	friend class Carac;
	public:
		CosmicBenchEvent();
		//copy constructor
		CosmicBenchEvent(const CosmicBenchEvent& other);
		//copy assignment
		CosmicBenchEvent& operator=(const CosmicBenchEvent& other);
		//build event for the given cosmic bench by using the treeObject data at the given entry
		CosmicBenchEvent(const CosmicBench * const detectors_, Tanalyse_R * treeObject, long entry);
		//build event for the given cosmic bench by using the treeObject data at the already loaded entry
		CosmicBenchEvent(const CosmicBench * const detectors_, const Tanalyse_R * const treeObject);
		//build event for the given cosmic bench by explicitely giving the Event object for all the detectors
		CosmicBenchEvent(const CosmicBench * const detectors_, const vector<Event*> events_);
		~CosmicBenchEvent();
		//do the tracking in order to create tracks for scattered muons
		void createPairs();
		//display the current event in a canvas
		void EventDisplay(TCanvas * c1 = 0);
		//retrieve the scattered muon of the given index
		RayPair get_rayPair(unsigned int i) const;
		//retrieve all the scattered muons of this event
		vector<RayPair> get_rayPairs() const;
		//retrieve the number of reconstructer scattered muons
		unsigned int get_rayPairs_N() const;
		//retrieve the number of events which is equal to the number of detectors inside the bench
		unsigned int get_event_N() const;
		//retrieve the total number of cluster across all the detectors
		unsigned int get_clus_N() const;
		//retrieve the number of cluster inside the given detector
		unsigned int get_clus_N_by_det(const Detector * const det) const;
		//do the CosMulti demultiplexing
		void Demux_CM();
		//delete clusters which do not pass the detectors specific cuts
		void do_cuts();
		//compute the straight muon tracks
		vector<Ray> get_absorption_rays(double chiSquare_threshold = -1);
		//compute the straight muon tracks using the hough method (currently not implemented)
		vector<Ray> get_hough_rays(double chiSquare_threshold = -1);
		//helper method to compute combinaison of one element of each key type in the set of sizes[key]
		template<typename T>
		static vector<map<T,int> > combinaisons(map<T,int> sizes, bool allow_drop = false);
		//void MultiCluster();
		//retrieve event id
		int get_evn() const;
		//retrieve evttime
		double get_evttime() const;
	protected:
		int evn;
		double evttime;
		const CosmicBench * detectors;
		vector<Event*> events;
		//After processing
		vector<RayPair> rayPairs;
};

#endif
