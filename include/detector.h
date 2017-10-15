#ifndef detector_h
#define detector_h
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <set>

#include "tomography.h"

//Boost
#include <boost/property_tree/ptree.hpp>

using boost::property_tree::ptree;

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;

class Cluster;
class Event;
class Tanalyse_R;
class TLine;

//abstract class to store the detectors caracteristics and 3D position
//has to be extended for a given type of detector
class Detector{
	public:
		//getters
		//retrieve the altitude of the detector
		double get_z() const;
		//check if the detector is measuring the X or Y coordinate
		bool get_is_X() const;
		//return the detector layer id in the bench
		int get_layer() const;
		//check if the detector data has to be used for the tracking
		bool get_is_ref() const;
		//retrieve the detector type
		virtual Tomography::det_type get_type() const = 0;
		//retrieve the detector measuring the perpendicular coordinate in case of 2D strip readout detectors such as the MG2D
		virtual Tomography::det_type get_perp_type() const;
		//retrieve the detector size accross the strip direction
		virtual double get_size() const = 0;
		//get the offset of the detector position in the direction accross the strips
		double get_offset() const;
		//check if the detector is upside down or not (inverse the direction axis)
		bool get_direction() const;
		//get the angle of the detector plane with respect to the x axis
		double get_angle_x() const;
		//get the angle of the detector plane with respect to the y axis
		double get_angle_y() const;
		//get the angle of the detector plane with respect to the z axis
		double get_angle_z() const;
		//retrieve the get_n_in_tree of the detector measuring the perpendicular coordinate in case of 2D strip readout detectors such as the MG2D
		int get_perp_n() const;
		//retrieve how many holes (un-hit strips) allowed inside a cluster
		int get_clustering_holes() const;
		//get the RMSPed of the ith channel of the detector
		double get_RMS(int i) const;
		//get the asic id to which the detector is connected and if the connector is plugged in the right way or not
		vector<pair<int,bool> > get_asic_n() const;
		//set the RMSPed of all the channel of the detectors
		virtual void set_RMS(vector<double> RMS_) = 0;
		//retrieve the channel to which the given strip id is connected (usefull for multiplexed detectors)
		virtual unsigned int StripToChannel(unsigned int i) const = 0;
		//retrieve the detector id in the storage ROOT tree
		int get_n_in_tree() const;
		//retrieve the total number of channels needed to read the detector
		virtual int get_Nchannel() const = 0;
		//retrieve the total number of strips composing the detectors
		virtual int get_Nstrip() const = 0;
		//retrieve the strip pitch
		virtual double get_StripPitch() const = 0;
		//retrieve the number of part the detector has to be divided in to compute the common noise
		virtual int get_CMN_div() const = 0;
		//check if the given cluster pass the detector dependent cuts
		virtual bool is_suitable(const Cluster * const clus) const = 0;
		virtual ~Detector();
		//copy the current detector
		virtual Detector * Clone() const = 0;
		//build the event corresponding to this detector using the treeObject data at the given entry
		virtual Event * build_event(Tanalyse_R * treeObject, int entry) const = 0;
		//build the event corresponding to this detector using the treeObject data at already loaded entry
		virtual Event * build_event(const Tanalyse_R * const treeObject) const = 0;
		//build the event corresponding to this detetector by giving the pedestal and common noise substracted amplitude (by channel and by sample), event id and event timestamp
		virtual Event * build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const = 0;
		//build a detector of this type using the data in the config tree
		virtual Detector * build_det(const ptree::value_type& child) const = 0;
		//channel mapping function while using feminos electronics
		virtual int feminos_mapping(int channel, bool connector_direction) const = 0;
		//channel mapping function while using dream electronics
		virtual int dream_mapping(int channel, bool connector_direction) const = 0;
		//retrieve the name of the detector (usually contain its type and id)
		virtual string Name() const = 0;
		//retrieve the maximum number of cluster by event
		virtual int get_MaxNClus() const = 0;
		//return a line to display the detector (usually a line at altitude z from offset-size/2 to offset+size/2)
		virtual TLine * get_line_display() const = 0;
	protected:
		Detector();
		//copy constructor
		Detector(const Detector& other);
		//copy assignment
		Detector& operator=(const Detector& other);
		//constructor using the explicit parameters
		Detector(double z_, bool is_X_, int layer_,int det_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, vector<pair<int,bool> > asic_n_);
		double z; //altitude inside cosmic bench
		bool is_X;//coordinate measured by the detector
		int layer;//layer which the detector is part of
		bool is_ref; //used to test detectors
		double offset; //used for alignement
		bool direction; // direction of the axis
		//used for alignement
		double angle_x;
		double angle_y;
		double angle_z;
		vector<double> RMS;
		int perp_n;
		int clustering_holes;
		int n_in_tree;
		vector<pair<int,bool> > asic_n;

};

//detector implementation for a single full connector dummy detector
class dummy_Detector: public Detector{
	public:
		dummy_Detector();
		dummy_Detector(int det_n_, int asic_n_, bool connector_direction_ = true);
		dummy_Detector(const dummy_Detector& other);
		dummy_Detector& operator=(const dummy_Detector& other);
		~dummy_Detector();
		Tomography::det_type get_type() const;
		double get_size() const;
		void set_RMS(vector<double> RMS_);
		unsigned int StripToChannel(unsigned int i) const;
		int get_Nchannel() const;
		int get_Nstrip() const;
		double get_StripPitch() const;
		int get_CMN_div() const;
		bool is_suitable(const Cluster * const clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry) const;
		Event * build_event(const Tanalyse_R * const treeObject) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const;
		Detector * build_det(const ptree::value_type& child) const;
		int feminos_mapping(int channel, bool connector_direction) const;
		int dream_mapping(int channel, bool connector_direction) const;
		string Name() const;
		int get_MaxNClus() const;
		TLine * get_line_display() const;
		static constexpr const int CMN_div = 2;
};

//detector implementation for the CosMulti detectors
class CM_Detector: public Detector{
	public:	
		CM_Detector();
		CM_Detector(const CM_Detector& other);
		CM_Detector& operator=(const CM_Detector& other);
		CM_Detector(double z_, bool is_X_, int layer_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int asic_n_, bool connector_direction = true);
		~CM_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		//CosMulti general charac
		static constexpr const double thinStripPitch = 500./1024.; // distance between the middle of 2 adjacent thin strips
		static constexpr const double wideStripPitch = 500./32.; // distance between the middle of 2 adjacent wide strips
		static constexpr const int Nchannel = 64;
		static constexpr const int Nstrip = 1024;
		static constexpr const double size = 500.;
		static constexpr const int CMN_div = 2;
		static constexpr const int MaxNClus = 600;
		unsigned int StripToChannel(unsigned int i) const;
		double get_size() const;
		int get_Nchannel() const;
		int get_Nstrip() const;
		double get_StripPitch() const;
		int get_CMN_div() const;
		//Cut
		void set_ClusMaxStripAmplCut_Min_Wide(double cut);
		void set_ClusSizeCut_Max_Wide(double cut);
		void set_ClusTOTCut_Min(double cut);
		void set_ClusMaxSampleCut_Min(double cut);
		void set_ClusMaxSampleCut_Max(double cut);
		//---
		bool get_use_thin_strip() const;
		Tomography::det_type get_type() const;
		void set_RMS(vector<double> RMS_);
		bool is_suitable(const Cluster * const clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry) const;
		Event * build_event(const Tanalyse_R * const treeObject) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const;
		int feminos_mapping(int channel, bool connector_direction) const;
		int dream_mapping(int channel, bool connector_direction) const;
		string Name() const;
		int get_MaxNClus() const;
		TLine * get_line_display() const;
	protected:
		bool use_thin_strip;
		//Detector dependent Cuts
		double ClusMaxStripAmplCut_Min_Wide;
		double ClusSizeCut_Max_Wide;
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;
};

//detector implementation for the MultiGen v1 detectors (50x50cm^2; 1024 strips)
class MG_Detector: public Detector{
	public:
		MG_Detector();
		MG_Detector(const MG_Detector& other);
		MG_Detector& operator=(const MG_Detector& other);
		MG_Detector(double z_, bool is_X_, int layer_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_, bool connector_direction = true);
		~MG_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		unsigned int StripToChannel(unsigned int i) const;
		static const vector<unsigned int> StripToChannel_a;
		static vector<unsigned int> ChannelToStrip(unsigned int channel_nb);
		//MultiGen general charac
		static constexpr const double StripPitch = 500./1024.; // distance between the middle of 2 adjacent strips
		static constexpr const int Nchannel = 61;
		static constexpr const int Nstrip = 1024;
		static constexpr const double size = 500.;
		static constexpr const int CMN_div = 2;
		static constexpr const int MaxNClus = 300;
		double get_size() const;
		int get_Nchannel() const;
		int get_Nstrip() const;
		double get_StripPitch() const;
		int get_CMN_div() const;
		//Cut
		void set_ClusSizeCut_Min(double cut);
		void set_ClusTOTCut_Min(double cut);
		void set_ClusMaxSampleCut_Min(double cut);
		void set_ClusMaxSampleCut_Max(double cut);
		double get_ClusSizeCut_Min() const;
		//---
		Tomography::det_type get_type() const;
		void set_RMS(vector<double> RMS_);
		void set_SRF(double offset_, double gauss, double lorentz, double ratio);
		double SRF_fit(double * x, double * p);
		bool is_suitable(const Cluster * const clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry) const;
		Event * build_event(const Tanalyse_R * const treeObject) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const;
		int feminos_mapping(int channel, bool connector_direction) const;
		int dream_mapping(int channel, bool connector_direction) const;
		string Name() const;
		int get_MaxNClus() const;
		TLine * get_line_display() const;
	protected:
		//Detector dependant cuts
		double ClusSizeCut_Min;
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;
		//Strip Response Function Parameters
		double srf_offset;
		double srf_gauss_width;
		double srf_lorentz_width;
		double srf_ratio;
};

//detector implementation for the MultiGen v2 detectors (50x50cm^2; 1037 strips)
class MGv2_Detector: public Detector{
	public:
		MGv2_Detector();
		MGv2_Detector(const MGv2_Detector& other);
		MGv2_Detector& operator=(const MGv2_Detector& other);
		MGv2_Detector(double z_, bool is_X_, int layer_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_, bool connector_direction = true);
		~MGv2_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		unsigned int StripToChannel(unsigned int i) const;
		static const vector<unsigned int> StripToChannel_a;
		static vector<unsigned int> ChannelToStrip(unsigned int channel_nb);
		//MultiGen general charac
		static constexpr const double StripPitch = 500./1037.; // distance between the middle of 2 adjacent strips
		static constexpr const int Nchannel = 61;
		static constexpr const int Nstrip = 1037;
		static constexpr const double size = 500.;
		static constexpr const int CMN_div = 2;
		static constexpr const int MaxNClus = 300;
		double get_size() const;
		int get_Nchannel() const;
		int get_Nstrip() const;
		double get_StripPitch() const;
		int get_CMN_div() const;
		//Cut
		void set_ClusSizeCut_Min(double cut);
		void set_ClusTOTCut_Min(double cut);
		void set_ClusMaxSampleCut_Min(double cut);
		void set_ClusMaxSampleCut_Max(double cut);
		double get_ClusSizeCut_Min() const;
		//---
		Tomography::det_type get_type() const;
		void set_RMS(vector<double> RMS_);
		void set_SRF(double offset_, double gauss, double lorentz, double ratio);
		double SRF_fit(double * x, double * p);
		bool is_suitable(const Cluster * const clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry) const;
		Event * build_event(const Tanalyse_R * const treeObject) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_, double evttime_) const;
		int feminos_mapping(int channel, bool connector_direction) const;
		int dream_mapping(int channel, bool connector_direction) const;
		string Name() const;
		int get_MaxNClus() const;
		TLine * get_line_display() const;
	protected:
		//Detector dependant cuts
		double ClusSizeCut_Min;
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;
		//Strip Response Function Parameters
		double srf_offset;
		double srf_gauss_width;
		double srf_lorentz_width;
		double srf_ratio;
};

//class to store all the bench detectors information
class CosmicBench{
	public:
		CosmicBench();
		//build a cosmic bench using the data in the config tree
		CosmicBench(ptree config_tree);
		//copy constructor
		CosmicBench(const CosmicBench& other);
		//copy assignment
		CosmicBench& operator=(const CosmicBench& other);
		~CosmicBench();
		//void add_MM(Detector * det);
		//retrieve the number of detector of given type inside the bench
		int get_det_N(Tomography::det_type det_t) const;
		//retrieve the number of detectors by type inside the bench
		map<Tomography::det_type,unsigned short> get_det_N() const;
		//retrieve the total number of detector inside the bench
		int get_det_N_tot() const;
		//retrive the total number of detector exclluded from the tracking
		int get_non_ref_N() const;
		//retrieve the number of different layers of detector
		int get_layers_n() const;
		//retrieve the detector at the given index (does not correspond to dector id (which is type dependant) but is the index in the vector number)
		Detector * get_detector(unsigned int i) const;
		//retrieve the detector index (does not correspond to dector id (which is type dependant) but is the index in the vector number) of the given cluster
		unsigned int find_det(const Cluster * const clus) const;
		//retrieve the detector index (does not correspond to dector id (which is type dependant) but is the index in the vector number) of the given detector type and tree id
		unsigned int find_det(Tomography::det_type det_t, unsigned int tree_n) const;
		//read the given Ped.dat file containing information about the detector list described by det_n_
		static map<Tomography::det_type,vector<vector<double> > > read_pedfile(string filename, map<Tomography::det_type,unsigned short> det_n_);
		//retrieve the altitude of the higher detector of the bench
		double get_z_Up() const;
		//retrieve the altitude of the lower detector of the bench
		double get_z_Down() const;
	protected:
		//some part of the construction is left in this init function
		void Init(ptree config_tree);
		vector<Detector*> detectors;
		map<Tomography::det_type,unsigned short> det_n;
		unsigned short non_ref_n;
		set<unsigned short> layers;
};

bool operator==(Detector const &det1, Detector const &det2);
bool operator!=(Detector const &det1, Detector const &det2);

#endif
