#ifndef cluster_h
#define cluster_h
#include <string>
#include "T.h"
#include "cluster.h"
#include "detector.h"
#include "tomography.h"

using std::string;

class Cluster;
class CM_Cluster;
class CM_Demux_Cluster;
class MG_Cluster;
class Detector;
class CM_Detector;
class MG_Detector;
class CosmicBench;
class T;

class Cluster{
	public:
		virtual ~Cluster();
		Tomography::det_type get_type() const;
		virtual double get_pos_mm() const = 0;
		virtual double correct_strip_nb(int strip_nb) const = 0;
		bool get_is_X() const;
		bool get_is_up() const;
		double get_z() const;
		double get_ampl() const;
		double get_size() const;
		double get_pos() const;
		double get_TOT() const;
		double get_t() const;
		double get_maxSample() const;
		double get_maxStripAmpl() const;
		int get_maxStrip() const;
		void set_perp_pos_mm(double coord);
		double get_perp_pos_mm() const;
		int find_det(const vector<Detector*> det_array) const;
		virtual bool is_in_det(Detector * det) const = 0;
		virtual int get_n_in_tree() const = 0;
	protected:
		Cluster();
		Cluster(const Cluster& other);
		Cluster& operator=(const Cluster& other);
		Cluster(T * treeObject, long entry = -1);
		Cluster(double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
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
		double angle;
		double perp_pos_mm;
};

class CM_Cluster: public Cluster{
	friend class CM_Demux_Cluster;
	public:
		CM_Cluster();
		CM_Cluster(const CM_Cluster& other);
		CM_Cluster& operator=(const CM_Cluster& other);
		CM_Cluster(T * treeObject,int number_,CM_Detector * det, long entry = -1);
		CM_Cluster(CM_Detector * det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~CM_Cluster();
		static bool is_suitable(T * treeObject,int number_,CM_Detector * detector, long entry = -1);
		bool is_suitable(CM_Detector * detector);
		bool is_in_det(Detector * det) const;
		Tomography::strip_type get_strip_type() const;
		virtual double get_pos_mm() const;
		virtual double correct_strip_nb(int strip_nb) const;
		int get_n_in_tree() const;
	protected:
		Tomography::strip_type strip_type;
		int cm_n_in_tree;
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
};

class MG_Cluster: public Cluster{
	public:
		MG_Cluster();
		MG_Cluster(const MG_Cluster& other);
		MG_Cluster& operator=(const MG_Cluster& other);
		MG_Cluster(T * treeObject,int number_,MG_Detector * det, long entry = -1);
		MG_Cluster(MG_Detector * det, int number_, double pos_, double size_, double ampl_, double maxSample_, double maxStripAmpl_, double TOT_, double t_, int maxStrip_);
		~MG_Cluster();
		static bool is_suitable(T * treeObject,int number_,MG_Detector * detector, long entry = -1);
		bool is_suitable(MG_Detector * detector);
		bool is_in_det(Detector * det) const;
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
		int get_n_in_tree() const;
	protected:
		int mg_n_in_tree;
};

#endif