#include "datareader.h"
#include <string>
#include <map>
#include <sstream>
#include <iomanip>

using std::string;
using std::map;
using std::ostringstream;
using std::setw;
using std::setfill;

int main(int argc, char ** argv){
	string outName = "../new_root_files/test_live";
	/*
	string baseName = "RawTest_465_400_141029_09H01";
	string dataDir = "/media/data/Clas12/CosmicBench/2014/W44";
	int FEU_N = 1;
	int dataFile_n = 1;
	int dataFile_n_offset = 0;
	*/
	map<int,string> det_type;
	map<int,int> det_asic;
	for(int i=0;i<8;i++){
		det_type[i] = "MG";
	}
	det_asic[0] = 3;
	det_asic[1] = 4;
	det_asic[2] = 5;
	det_asic[3] = 2;
	det_asic[4] = 0;
	det_asic[5] = 1;
	det_asic[6] = 6;
	det_asic[7] = 7;
	FeminosDataReader blah(outName,det_type,det_asic);//,false,false,false,100);
	blah.add_file_to_process("/media/pc_daq_new/Feminos/2014-11-21/R2014_11_27-10_31_35-018.aqs");
	/*
	for(int i=0;i<dataFile_n;i++){
		ostringstream dataFileName;
		dataFileName << dataDir << "/" << baseName << "_" << setw(3) << setfill('0') << i+dataFile_n_offset << "_" << setw(2) << setfill('0') << FEU_N << ".fdf";
		blah.add_file_to_process(dataFileName.str());
	}
	*/
	blah.process();
	blah.compute_ped();
	blah.do_ped_sub();
	blah.do_common_noise_sub();
	blah.compute_RMSPed();
	return 0;
}