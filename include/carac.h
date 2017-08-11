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
		//build object given the file name containing the json config tree
		Carac(string configFilePath);
		//build object given the config tree
		Carac(ptree config_tree_);
		~Carac();
		//void Efficacity();
		//void Residus();
		//Do the same as analyse::Residus_ref (cf. analyse.h) but for every detectors in parrallel
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
