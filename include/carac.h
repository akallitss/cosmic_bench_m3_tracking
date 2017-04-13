#ifndef carac_h
#define carac_h
#include "Tanalyse_R.h"
#include "tomography.h"
#include <boost/property_tree/ptree.hpp>
#include <string>

using std::string;

using boost::property_tree::ptree;

class TFile;

class Carac: public Tanalyse_R{
	public:
		Carac(string configFilePath);
		Carac(ptree config_tree_);
		~Carac();
		//void Efficacity();
		//void Residus();
		void Residus_ref(const double theta_max = 10000000);
		//double Residus_ref_cost();
		//void Residus_ref_2D();
	protected:
		void Init();
		long max_event;
		TFile * f;
		ptree config_tree;
};

#endif
