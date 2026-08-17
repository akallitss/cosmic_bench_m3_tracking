#!/bin/sh

pathnameped="/mnt/cosmic_data/MX17/Run/mx17_det1_det2_overnight_6-17-26/resist_450V_drift_1000V/raw_daq_data/"
pathname="/mnt/cosmic_data/MX17/Run/mx17_det1_det2_overnight_6-17-26/resist_450V_drift_1000V/raw_daq_data/"

file_num=0

declare aped="MX17_resist_450V_drift_1000V_pedthr_260618_05H10"
declare acos="MX17_resist_450V_drift_1000V_datrun_260618_05H10"

#make pedestals/RMS
echo "Pedestal files: $pathnameped${aped}"_000_01.fdf  
#make ending of config file
echo    "\"data_file_first\" : 0," > end.txt
echo    "\"data_file_last\" : 0," >> end.txt
echo	"\"data_file_basename\" : \"$pathnameped${aped}_\"}"  >> end.txt
cat config_ref_2022.json end.txt > config.json
  
#clean previous ped
rm root_files/test_signal.root root_files/test_RMSPed.dat root_files/test_Ped.dat
./DataReader config.json read
./DataReader config.json ped

#make cosmic config file
#clean previous cos
rm config_cos.json end.txt
num=$(printf "%03d" "$file_num")
echo "Cosmics data files: $pathname${acos}"_"$num"_01.fdf
echo    "\"data_file_first\" : "$file_num"," > end.txt  
echo    "\"data_file_last\" : "$file_num"," >> end.txt  
echo	"\"data_file_basename\" : \"$pathname${acos}_\"}"  >> end.txt
cat config_ref_2022.json end.txt > config_cos.json
./DataReader config_cos.json analyse
./tracking config_cos.json rays output_"$num".root
