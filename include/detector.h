#ifndef detector_h
#define detector_h
#include <string>
#include <vector>
#include <map>

#include "tomography.h"
#include "cluster.h"
#include "event.h"

//Boost
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

using std::string;
using std::vector;
using std::map;

class Cluster;
class Event;
class Tanalyse_R;

class Detector{
	public:
		//getters
		double get_z() const;
		bool get_is_X() const;
		bool get_is_up() const;
		bool get_is_ref() const;
		virtual Tomography::det_type get_type() const = 0;
		virtual Tomography::det_type get_perp_type() const;
		virtual double get_size() const = 0;
		double get_offset() const;
		bool get_direction() const;
		double get_angle_x() const;
		double get_angle_y() const;
		double get_angle_z() const;
		int get_perp_n() const;
		int get_clustering_holes() const;
		double get_RMS(int i) const;
		int get_asic_n() const;
		virtual void set_RMS(vector<double> RMS_) = 0;
		virtual unsigned int StripToChannel(unsigned int i) const = 0;
		int get_n_in_tree() const;
		virtual int get_Nchannel() const = 0;
		virtual int get_Nstrip() const = 0;
		virtual double get_StripPitch() const = 0;
		virtual int get_CMN_div() const = 0;
		virtual bool is_suitable(Cluster * clus) const = 0;
		virtual ~Detector();
		virtual Detector * Clone() const = 0;
		virtual Event * build_event(Tanalyse_R * treeObject, int entry = -1) const = 0;
		virtual Event * build_event(vector<vector<double> > strip_ampl_, int evn_) const = 0;
		virtual Detector * build_det(const ptree::value_type& child) const = 0;
		virtual int feminos_mapping(int channel) const = 0;
		virtual int dream_mapping(int channel) const = 0;
	protected:
		Detector();	
		Detector(const Detector& other);
		Detector& operator=(const Detector& other);
		Detector(double z_, bool is_X_, bool is_up_,int det_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_);	
		double z; //altitude inside cosmic bench
		bool is_X;//coordinate measured by the detector
		bool is_up;//bloc (up|down) which the detector is part of
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
		int asic_n;

};

class CM_Detector: public Detector{
	public:	
		CM_Detector();
		CM_Detector(const CM_Detector& other);
		CM_Detector& operator=(const CM_Detector& other);
		CM_Detector(double z_, bool is_X_, bool is_up_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int asic_n_);
		~CM_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		//CosMulti general charac
		static const double thinStripPitch = 500./1024.; // distance between the middle of 2 adjacent thin strips
		static const double wideStripPitch = 500./32.; // distance between the middle of 2 adjacent wide strips
		static const int Nchannel = 64;
		static const int Nstrip = 1024;
		static const double size = 500.;
		static const int CMN_div = 2;
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
		bool is_suitable(Cluster * clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry = -1) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_) const;
		int feminos_mapping(int channel) const;
		int dream_mapping(int channel) const;
	protected:
		bool use_thin_strip;
		//Detector dependent Cuts
		double ClusMaxStripAmplCut_Min_Wide;
		double ClusSizeCut_Max_Wide;
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;
};

class MG_Detector: public Detector{
	public:
		MG_Detector();
		MG_Detector(const MG_Detector& other);
		MG_Detector& operator=(const MG_Detector& other);
		MG_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_);
		~MG_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		unsigned int StripToChannel(unsigned int i) const;
		static const vector<unsigned int> StripToChannel_a;
		static vector<unsigned int> ChannelToStrip(unsigned int channel_nb);
		//MultiGen general charac
		static const double StripPitch = 500./1024.; // distance between the middle of 2 adjacent strips
		static const int Nchannel = 61;
		static const int Nstrip = 1024;
		static const double size = 500.;
		static const int CMN_div = 2;
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
		bool is_suitable(Cluster * clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry = -1) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_) const;
		int feminos_mapping(int channel) const;
		int dream_mapping(int channel) const;
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

class MGv2_Detector: public Detector{
	public:
		MGv2_Detector();
		MGv2_Detector(const MGv2_Detector& other);
		MGv2_Detector& operator=(const MGv2_Detector& other);
		MGv2_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_, int asic_n_);
		~MGv2_Detector();
		Detector * build_det(const ptree::value_type& child) const;
		unsigned int StripToChannel(unsigned int i) const;
		static const vector<unsigned int> StripToChannel_a;
		static vector<unsigned int> ChannelToStrip(unsigned int channel_nb);
		//MultiGen general charac
		static const double StripPitch = 500./1037.; // distance between the middle of 2 adjacent strips
		static const int Nchannel = 61;
		static const int Nstrip = 1037;
		static const double size = 500.;
		static const int CMN_div = 2;
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
		bool is_suitable(Cluster * clus) const;
		Detector * Clone() const;
		Event * build_event(Tanalyse_R * treeObject, int entry = -1) const;
		Event * build_event(vector<vector<double> > strip_ampl_, int evn_) const;
		int feminos_mapping(int channel) const;
		int dream_mapping(int channel) const;
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

class CosmicBench{
	public:
		CosmicBench();
		CosmicBench(ptree config_tree);
		CosmicBench(const CosmicBench& other);
		CosmicBench& operator=(const CosmicBench& other);
		~CosmicBench();
		//void add_MM(Detector * det);
		int get_det_N(Tomography::det_type det_t) const;
		map<Tomography::det_type,unsigned short> get_det_N() const;
		int get_det_N_tot() const;
		Detector * get_detector(unsigned int i) const;
		static map<Tomography::det_type,vector<vector<double> > > read_pedfile(string filename, map<Tomography::det_type,unsigned short> det_n_);
	protected:
		void Init(ptree config_tree);
		vector<Detector*> detectors;
		map<Tomography::det_type,unsigned short> det_n;
};

bool operator==(Detector const &det1, Detector const &det2);
bool operator!=(Detector const &det1, Detector const &det2);

#endif