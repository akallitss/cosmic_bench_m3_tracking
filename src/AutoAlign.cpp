#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <TMath.h>

#include <string>
#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <sstream>

#include "analyse.h"
#include "detector.h"

using boost::property_tree::ptree;

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::map;
using std::ostringstream;

using TMath::Sqrt;

class Minimizer{
public:
	Minimizer(ptree config, string path_);
	double operator()(const double * x);
	ptree get_tree_start() const;
	ptree get_tree_current() const;
private:
	ptree config_tree_start;
	ptree config_tree_current;
	int iteration_nb;
	string path;
	vector<pair<int,int> > MG2D_pair;
};

Minimizer::Minimizer(ptree config, string path_){
	config_tree_start = config;
	config_tree_current = config;
	iteration_nb = 0;
	path = path_;
	map<int,int> MG2D_pair_;
	BOOST_FOREACH(const ptree::value_type& child, config_tree_start.get_child("CosmicBench.MultiGens")){
		int current_n = child.second.get<int>("mg_n");
		int current_perp_n = child.second.get<int>("2D_perp_n");
		if(current_perp_n<0){
			cout << "error 1D detector !" << endl;
			return;
		}
		map<int,int>::iterator map_it = MG2D_pair_.find(current_n);
		if(map_it != MG2D_pair_.end()){
			cout << "detector appears more than 1 time" << endl;
			return;
		}
		map_it = MG2D_pair_.find(current_perp_n);
		if(map_it != MG2D_pair_.end()){
			if(map_it->second != current_n){
				cout << "problem in 2D association" << endl;
				return;
			}
			else continue;
		}
		MG2D_pair_[current_n] = current_perp_n;
	}
	for(map<int,int>::iterator map_it = MG2D_pair_.begin();map_it!=MG2D_pair_.end();++map_it){
		MG2D_pair.push_back(pair<int,int>(map_it->first,map_it->second));
	}
}

double Minimizer::operator()(const double * x){
	config_tree_current = config_tree_start;
	for(ptree::iterator tree_it = config_tree_current.get_child("CosmicBench.MultiGens").begin(); tree_it != config_tree_current.get_child("CosmicBench.MultiGens").end(); tree_it++){
		tree_it->second.put<bool>("is_ref",true);
	}
	for(unsigned int i=0;i<MG2D_pair.size();i++){
		for(ptree::iterator tree_it = config_tree_current.get_child("CosmicBench.MultiGens").begin(); tree_it != config_tree_current.get_child("CosmicBench.MultiGens").end(); tree_it++){
			if(tree_it->second.get<int>("mg_n") == MG2D_pair[i].first || tree_it->second.get<int>("mg_n") == MG2D_pair[i].second){
				tree_it->second.put<double>("z",tree_it->second.get<double>("z") + x[(6*i)+2]);
				tree_it->second.put<double>("angle_x",x[(6*i)+3]);
				tree_it->second.put<double>("angle_y",x[(6*i)+4]);
				tree_it->second.put<double>("angle_z",x[(6*i)+5]);
				if(tree_it->second.get<bool>("is_X")){
					tree_it->second.put<double>("offset",x[(6*i)]);
				}
				else{
					tree_it->second.put<double>("offset",x[(6*i)+1]);
				}
			}
		}
		/*
		BOOST_FOREACH(const ptree::value_type& child, config_tree_current.get_child("CosmicBench.MultiGens")){
			if(child.second.get<int>("mg_n") == MG2D_pair[i].first || child.second.get<int>("mg_n") == MG2D_pair[i].second){
				child.second.put<double>("z",child.second.get<double>("z") + x[(6*i)+2]);
				child.second.put<double>("angle_x",x[(6*i)+3]);
				child.second.put<double>("angle_y",x[(6*i)+4]);
				child.second.put<double>("angle_z",x[(6*i)+4]);
				if(child.second.get<bool>("is_X")){
					child.second.put<double>("offset",x[(6*i)]);
				}
				else{
					child.second.put<double>("offset",x[(6*i)+1]);
				}
			}
		}
		*/
	}
	ostringstream filename;
	filename << path << "_" << iteration_nb << ".cfg";
	write_json(filename.str().c_str(), config_tree_current);
	iteration_nb++;
	double cost = 0;
	for(unsigned int i=0;i<MG2D_pair.size();i++){
		for(ptree::iterator tree_it = config_tree_current.get_child("CosmicBench.MultiGens").begin(); tree_it != config_tree_current.get_child("CosmicBench.MultiGens").end(); tree_it++){
			if(tree_it->second.get<int>("mg_n") == MG2D_pair[i].first || tree_it->second.get<int>("mg_n") == MG2D_pair[i].second){
				tree_it->second.put<bool>("is_ref",false);
			}
			else{
				tree_it->second.put<bool>("is_ref",true);
			}
		}
		Analyse * current_analyse = new Analyse(config_tree_current);
		cost += current_analyse->Residus_ref_cost();
		delete current_analyse;
	}
	return Sqrt(cost);
}

ptree Minimizer::get_tree_start() const{
	return config_tree_start;
}
ptree Minimizer::get_tree_current() const{
	return config_tree_current;
}

int main(int argc, char ** argv){
	if(argc<3){
		cout << "You must indicate a config file which contains the Cosmic Bench caracs and a path were you want to store configs" << endl;
		return 1;
	}
	ptree configTree;
	read_json(argv[1], configTree);
	Minimizer current_minimizer(configTree, argv[2]);
	ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

	min.SetMaxFunctionCalls(1000000);
	min.SetMaxIterations(100000);
	min.SetTolerance(0.001);

	ROOT::Math::Functor f(current_minimizer,18); 
	CosmicBench bench_carac(configTree);
	if(bench_carac.get_CM_N() > 0){
		cout << "not implemented for CM" << endl;
		return 0;
	}
	if(bench_carac.get_MG_N() != 8){
		cout << "only implemented with 4 MG2D setup" << endl;
	}
	double step[18];
	double variable[18];
	string name[18];
	for(int i=0;i<3;i++){
		ostringstream common_name;
		common_name << "_MG_" << i;
		name[(6*i)] = ("dx" + common_name.str()).c_str();
		step[(6*i)] = 0.1;
		variable[(6*i)] = 0;
		name[(6*i)+1] = ("dy" + common_name.str()).c_str();
		step[(6*i)+1] = 0.1;
		variable[(6*i)+1] = 0;
		name[(6*i)+2] = ("dz" + common_name.str()).c_str();
		step[(6*i)+2] = 0.1;
		variable[(6*i)+2] = 0;
		name[(6*i)+3] = ("d_theta_x" + common_name.str()).c_str();
		step[(6*i)+3] = 0.001;
		variable[(6*i)+3] = 0;
		name[(6*i)+4] = ("d_theta_y" + common_name.str()).c_str();
		step[(6*i)+4] = 0.001;
		variable[(6*i)+4] = 0;
		name[(6*i)+5] = ("d_theta_z" + common_name.str()).c_str();
		step[(6*i)+5] = 0.001;
		variable[(6*i)+5] = 0;
	}

	min.SetFunction(f);

	// Set the free variables to be minimized!
	for(int i=0;i<18;i++){
		min.SetVariable(i,name[i].c_str(),variable[i], step[i]);
	}

	min.Minimize(); 
	write_json((argv[2] + string("_final.cfg")).c_str(), current_minimizer.get_tree_current());
	return 1;
}