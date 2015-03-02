#ifndef detector_h
#define detector_h
#include <string>
#include <vector>

#include "tomography.h"

//Boost
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

using std::string;
using std::vector;

class Detector{
	public:
		//getters
		double get_z() const;
		bool get_is_X() const;
		bool get_is_up() const;
		bool get_is_ref() const;
		virtual Tomography::det_type get_type() const = 0;
		virtual double get_size() const = 0;
		double get_offset() const;
		bool get_direction() const;
		double get_angle_x() const;
		double get_angle_y() const;
		double get_angle_z() const;
		int get_perp_n() const;
		int get_clustering_holes() const;
		//seters
		void set_ClusTOTCut_Min(double cut);
		void set_ClusMaxSampleCut_Min(double cut);
		void set_ClusMaxSampleCut_Max(double cut);
		bool test_ClusTOT(double TOT) const;
		bool test_ClusMaxSample(double maxSample) const;
		double get_RMS(int i) const;
		virtual void set_RMS(vector<double> RMS_) = 0;
		virtual ~Detector();
	protected:
		Detector();	
		Detector(const Detector& other);
		Detector& operator=(const Detector& other);
		Detector(double z_, bool is_X_, bool is_up_, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_);	
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
		//Detector dependent Cuts
		double ClusTOTCut_Min;
		double ClusMaxSampleCut_Min;
		double ClusMaxSampleCut_Max;
		vector<double> RMS;
		int perp_n;
		int clustering_holes;

};

class CM_Detector: public Detector{
	public:	
		CM_Detector();
		CM_Detector(const CM_Detector& other);
		CM_Detector& operator=(const CM_Detector& other);
		CM_Detector(double z_, bool is_X_, bool is_up_, int cm_n, bool use_thin_strip_, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_);
		~CM_Detector();
		//CosMulti general charac
		static const double thinStripPitch = 500./1024.; // distance between the middle of 2 adjacent thin strips
		static const double wideStripPitch = 500./32.; // distance between the middle of 2 adjacent wide strips
		static const int Nstrip = 64;
		static const double size = 500.;
		double get_size() const;
		//Cut setters
		void set_ClusMaxStripAmplCut_Min_Wide(double cut);
		void set_ClusSizeCut_Max_Wide(double cut);
		bool test_ClusMaxStripAmpl_Wide(double maxStripAmpl) const;
		bool test_ClusSize_Wide(double size) const;
		int get_cm_n_in_tree() const;
		bool get_use_thin_strip() const;
		Tomography::det_type get_type() const;
		void set_RMS(vector<double> RMS_);
	protected:
		int cm_n_in_tree;
		bool use_thin_strip;
		//Detector dependent Cuts
		double ClusMaxStripAmplCut_Min_Wide;
		double ClusSizeCut_Max_Wide;
};

class MG_Detector: public Detector{
	public:
		MG_Detector();
		MG_Detector(const MG_Detector& other);
		MG_Detector& operator=(const MG_Detector& other);
		MG_Detector(double z_, bool is_X_, bool is_up_, int mg_n, bool is_ref_, double offset_, bool direction_, double angle_x_, double angle_y_, double angle_z_, int perp_n_, int clustering_holes_);
		~MG_Detector();
		static unsigned int StripToChannel_f(unsigned int strip_nb);
		static const vector<unsigned int> StripToChannel;
		static vector<unsigned int> ChannelToStrip(unsigned int channel_nb);
		//MultiGen general charac
		static const double StripPitch = 500./1024.; // distance between the middle of 2 adjacent strips
		static const int Nstrip = 61;
		static const double size = 500.;
		double get_size() const;
		//Cut setters
		void set_ClusSizeCut_Min(double cut);
		bool test_ClusSize(double size);
		int get_mg_n_in_tree() const;
		Tomography::det_type get_type() const;
		void set_RMS(vector<double> RMS_);
		void set_SRF(double offset_, double gauss, double lorentz, double ratio);
		double SRF_fit(double * x, double * p);
	protected:
		int mg_n_in_tree;
		//Detector dependant cuts
		double ClusSizeCut_Min;
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
		~CosmicBench();
		//void add_MM(Detector * det);
		int get_CM_N() const;
		int get_MG_N() const;
		Detector * get_detector(unsigned int i) const;
	protected:
		void Init(ptree config_tree);
		vector<Detector*> detectors;
		int CM_N;
		int MG_N;
};

bool operator==(Detector const &det1, Detector const &det2);
bool operator!=(Detector const &det1, Detector const &det2);

#endif