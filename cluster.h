#ifndef cluster_h
#define cluster_h
#include <string>
#include "T.h"
#include "ray.h"
#include "cluster.h"
#include "detector.h"

using std::string;

class Event;
class CM_Event;
class CM_Demux_Event;
class MG_Event;
class CosmicBenchEvent;
class Ray_2D;
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

class Cluster{
	friend class Ray;
	friend class Ray_2D;
	friend class RayPair;
	friend class CosmicBenchEvent;
	public:
		virtual ~Cluster();
		string get_type() const;
		virtual double get_pos_mm() const = 0;
		virtual double correct_strip_nb(int strip_nb) const = 0;
		bool get_is_X() const;
		bool get_is_up() const;
		double get_z() const;
		double get_ampl() const;
		double get_size() const;
		double get_pos() const;
		void set_perp_pos_mm(double coord);
		double get_perp_pos_mm() const;
		int find_det(const vector<Detector*> det_array) const;
	protected:
		Cluster();
		Cluster(const Cluster& other);
		Cluster& operator=(const Cluster& other);
		Cluster(T * treeObject, int entry = -1);
		virtual bool is_in_det(Detector * det) const = 0;
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
		string type;
		double z;
		bool is_X;
		bool is_up;
		double offset;
		bool direction;
		double angle;
		double perp_pos_mm;
};

class CM_Cluster: public Cluster{
	friend class Ray;
	friend class Ray_2D;
	friend class RayPair;
	friend class CM_Demux_Cluster;
	friend class CosmicBenchEvent;
	public:
		CM_Cluster();
		CM_Cluster(const CM_Cluster& other);
		CM_Cluster& operator=(const CM_Cluster& other);
		CM_Cluster(T * treeObject,int number_,CM_Detector * det, int entry = -1);
		~CM_Cluster();
		static bool is_suitable(T * treeObject,int number_,CM_Detector * detector, int entry = -1);
		bool is_in_det(Detector * det) const;
		string get_strip_type() const;
		virtual double get_pos_mm() const;
		virtual double correct_strip_nb(int strip_nb) const;
	protected:
		int maxStrip;
		string strip_type;
		int cm_n_in_tree;
};

class CM_Demux_Cluster: public CM_Cluster{
	friend class Ray;
	friend class Ray_2D;
	friend class RayPair;
	friend class CosmicBenchEvent;
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
	friend class Ray;
	friend class Ray_2D;
	friend class RayPair;
	friend class CosmicBenchEvent;
	public:
		MG_Cluster();
		MG_Cluster(const MG_Cluster& other);
		MG_Cluster& operator=(const MG_Cluster& other);
		MG_Cluster(T * treeObject,int number_,MG_Detector * det, int entry = -1);
		~MG_Cluster();
		static bool is_suitable(T * treeObject,int number_,MG_Detector * detector, int entry = -1);
		bool is_in_det(Detector * det) const;
		double get_pos_mm() const;
		double correct_strip_nb(int strip_nb) const;
	protected:
		int mg_n_in_tree;
};

#endif