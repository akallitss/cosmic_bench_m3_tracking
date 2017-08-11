#ifndef cluster_h
#define cluster_h
#include <string>
#include "tomography.h"
#include "ray.h"

using std::string;

class Detector;
class Tanalyse_R;

//abstract class to store cluster informations
//has to be extended for a given type of detector
class Cluster{
	public:
		virtual ~Cluster();
		//return type of detector of said cluster
		Tomography::det_type get_type() const;
		//compute the cluster position taking into account the alignment parameters (angle and offset)
		virtual double get_pos_mm() const = 0;
		//compute a virtual strip id corresponding to the cluster position and taking into account the alignment parameters (angle and offset)
		virtual double correct_strip_nb(int strip_nb) const = 0;
		//return true if the cluster bears information about the X coordinate and false for the Y coordinate
		bool get_is_X() const;
		//return the detector layer id in the bench
		int get_layer() const;
		//return the altitude of the cluster taking into account the alignment parameters (angle and offset)
		virtual double get_z() const = 0;
		//return the cluster amplitude
		double get_ampl() const;
		//return the cluster size
		double get_size() const;
		//return the cluster position in strip
		double get_pos() const;
		//return the cluster Time Over Threshold
		double get_TOT() const;
		//return the cluster time
		double get_t() const;
		//return the maximum amplitude sample
		double get_maxSample() const;
		//return the amplitude of the maximum amplitude strip
		double get_maxStripAmpl() const;
		//return the strip id which carries the maximum amplitude
		int get_maxStrip() const;
		//for strip detector, give information about the position of the cluster along the strip direction (from another detector)
		void set_perp_pos_mm(double coord);
		//for strip detector, give information about the position of the cluster along the strip direction (from the muon track inducing this cluster)
		void set_perp_pos_mm(Ray ray);
		//get the external information set previously
		double get_perp_pos_mm() const;
		//return the index in the given array pointing to the detector in which the detector is sitting
		int find_det(const vector<Detector*> det_array) const;
		//check if the cluster sit in the given detector
		bool is_in_det(const Detector * const det) const;
		//get the detector index inside the storing tree
		int get_n_in_tree() const;
		//output cluster information inside a string
		string print() const;
		//make a copy of the current cluster
		virtual Cluster * Clone() const = 0;
	protected:
		Cluster();
		//copy constructor
		Cluster(const Cluster& other);
		//copy assignment
		Cluster& operator=(const Cluster& other);
		//build cluster of id number_ of the given detector from treeObject data at given entry
		Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry);
		//build cluster of id number_ of the given detector from treeObject data at already loaded entry
		Cluster(const Tanalyse_R * const treeObject,int number_,const Detector * const det);
		//build the cluster by giving explicitely its caracteristics
		Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		int evn;
		double evttime;
		int number;
		double ampl;
		double size;
		double pos;
		double maxStripAmpl;
		double maxSample;
		double TOT;
		double t;
		int maxStrip;
		Tomography::det_type type;
		double z;
		bool is_X;
		int layer;
		double offset;
		bool direction;
		double angle_x;
		double angle_y;
		double angle_z;
		double perp_pos_mm;
		int n_in_tree;
};

//implementation for a dummy detector of 64 strips (1 connector)
class dummy_Cluster: public Cluster{
	public:
		dummy_Cluster();
		dummy_Cluster(const dummy_Cluster& other);
		dummy_Cluster& operator=(const dummy_Cluster& other);
		~dummy_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

//implementation for CosMulti detector to store still multiplexed informations
class CM_Cluster: public Cluster{
	friend class CM_Demux_Cluster;
	public:
		CM_Cluster();
		CM_Cluster(const CM_Cluster& other);
		CM_Cluster& operator=(const CM_Cluster& other);
		CM_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry);
		CM_Cluster(const Tanalyse_R * const treeObject,int number_,const Detector * const det);
		CM_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~CM_Cluster();
		Tomography::strip_type get_strip_type() const;
		virtual double get_pos_mm() const;
		virtual double correct_strip_nb(int strip_nb) const;
		virtual double get_z() const;
		virtual Cluster * Clone() const;
	protected:
		Tomography::strip_type strip_type;
};

//implementation for CosMulti detector to store demultiplexed informations
class CM_Demux_Cluster: public CM_Cluster{
	public:
		CM_Demux_Cluster();
		CM_Demux_Cluster(const CM_Demux_Cluster& other);
		CM_Demux_Cluster& operator=(const CM_Demux_Cluster& other);
		CM_Demux_Cluster(const CM_Cluster& thinStrip_clus, const CM_Cluster& wideStrip_clus);
		CM_Demux_Cluster(const CM_Cluster& wideStrip_clus);
		~CM_Demux_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

//implementation to store MultiGen V1 informations (50x50cm^2; 1024 strips)
class MG_Cluster: public Cluster{
	public:
		MG_Cluster();
		MG_Cluster(const MG_Cluster& other);
		MG_Cluster& operator=(const MG_Cluster& other);
		MG_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry);
		MG_Cluster(const Tanalyse_R * const treeObject,int number_,const Detector * const det);
		MG_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~MG_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

//implementation to store MultiGen V2 informations (50x50cm^2; 1037 strips)
class MGv2_Cluster: public Cluster{
	public:
		MGv2_Cluster();
		MGv2_Cluster(const MGv2_Cluster& other);
		MGv2_Cluster& operator=(const MGv2_Cluster& other);
		MGv2_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry);
		MGv2_Cluster(const Tanalyse_R * const treeObject,int number_,const Detector * const det);
		MGv2_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~MGv2_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

#endif
