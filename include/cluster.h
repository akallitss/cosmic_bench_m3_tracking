#ifndef cluster_h
#define cluster_h
#include <string>
#include "cluster.h"
#include "detector.h"
#include "tomography.h"
#include "ray.h"

using std::string;

class Cluster;
class CM_Cluster;
class CM_Demux_Cluster;
class MG_Cluster;
class Detector;
class CM_Detector;
class MG_Detector;
class CosmicBench;
class Tanalyse_R;
class Ray;

class Cluster{
	public:
		virtual ~Cluster();
		Tomography::det_type get_type() const;
		virtual double get_pos_mm() const = 0;
		virtual double correct_strip_nb(int strip_nb) const = 0;
		bool get_is_X() const;
		bool get_is_up() const;
		virtual double get_z() const = 0;
		double get_ampl() const;
		double get_size() const;
		double get_pos() const;
		double get_TOT() const;
		double get_t() const;
		double get_maxSample() const;
		double get_maxStripAmpl() const;
		int get_maxStrip() const;
		void set_perp_pos_mm(double coord);
		void set_perp_pos_mm(Ray ray);
		double get_perp_pos_mm() const;
		int find_det(const vector<Detector*> det_array) const;
		bool is_in_det(const Detector * det) const;
		int get_n_in_tree() const;
		virtual Cluster * Clone() const = 0;
	protected:
		Cluster();
		Cluster(const Cluster& other);
		Cluster& operator=(const Cluster& other);
		Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry = -1);
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
		bool is_up;
		double offset;
		bool direction;
		double angle_x;
		double angle_y;
		double angle_z;
		double perp_pos_mm;
		int n_in_tree;
};

class CM_Cluster: public Cluster{
	friend class CM_Demux_Cluster;
	public:
		CM_Cluster();
		CM_Cluster(const CM_Cluster& other);
		CM_Cluster& operator=(const CM_Cluster& other);
		CM_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry = -1);
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

class MG_Cluster: public Cluster{
	public:
		MG_Cluster();
		MG_Cluster(const MG_Cluster& other);
		MG_Cluster& operator=(const MG_Cluster& other);
		MG_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry = -1);
		MG_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~MG_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

class MGv2_Cluster: public Cluster{
	public:
		MGv2_Cluster();
		MGv2_Cluster(const MGv2_Cluster& other);
		MGv2_Cluster& operator=(const MGv2_Cluster& other);
		MGv2_Cluster(Tanalyse_R * treeObject,int number_,const Detector * const det, long entry = -1);
		MGv2_Cluster(const Detector * const det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~MGv2_Cluster();
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		double get_z() const;
		Cluster * Clone() const;
};

#endif