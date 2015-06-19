#ifndef carac_h
#define carac_h
#include "Tanalyse_R.h"
#include <TFile.h>
#include <boost/property_tree/ptree.hpp>
#include <string>

using std::string;

using boost::property_tree::ptree;

class Carac: public Tanalyse_R{
	public:
		Carac(string configFilePath);
		Carac(ptree config_tree_);
		~Carac();
		//void Efficacity();
		//void Residus();
		void Residus_ref();
		//double Residus_ref_cost();
		//void Residus_ref_2D();
	protected:
		long max_event;
		TFile * f;
		ptree config_tree;
		int MG_N;
		int CM_N;
};

#endif