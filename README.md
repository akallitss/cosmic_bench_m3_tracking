### Remarks

If you are using ROOT6 or C++11, apply cpp11.patch before build : `patch -p1 < cpp11.patch`

### Use this soft

* Datareader :

`Datareader <config_file> {ped,data,analyse,live,read}`

config_file : path to config file which contain the cosmic bench caracteristics and different file path  
ped : use this option to build signal file ped and corr branches, calculate the Ped.dat and RMS.dat using an existing signal file with filled raw branches  
data : use this option to build signal file ped and corr branches using existing Ped.dat and signal file with filled raw branches
analyse : use this option to build analyse file using existing Ped.dat and RMS.dat  
live : use this option to read online data from FeuUdpControl using option -q to build the raw branches of the signal file  
read : use this option to read electronic binary files and build the raw branches of the signal file

* Multicluster

`Multicluster <config_file> {multicluster,SoB,dispersion,display,hough} [additional options]`

config_file : path to config file which contain the cosmic bench caracteristics and different file path  
multicluster : use this option to build analyse file using existing Ped.dat and RMS.dat and signal file  
rawcluster : use this option to build analyse file using existing Ped.dat and RMS.dat and signal file with only raw branches filled  
SoB : use this option to display signal and noise amplitude  
dispersion : use this option to display some correlation between the detectors signal  
display : use this option to display the signal shape for each detector for a signal type (raw,ped or corr) from event 0 to n (n being the additional option)  
hough : use this option to study hough tracking style for event n (n being the additional option)
conv : use this option to study convolution clustering

* tracking

`tracking <config_file> {efficacity,eff2D,residus,fluxmap,tomoAbs,raypairs,srf,correlation,eventdisplay,eventdisplaymult,SoN,acceptance,watto} [additional options]`

config_file : path to config file which contain the cosmic bench caracteristics and different file path  
eff2D : compute 2D efficiency of 2D detectors (you need to switch is_ref to
false for both coordinates  
residus : compute efficiency, residuals and various alignement distribution for
the is_ref=false detectors  
fluxmap : compute muon flux for a given height in additional option  
raypairs : export deviation tomography rays in a file which name is given in
additional option  
correlation : show correlation between detectors/coordinates for time and
amplitude  
eventdisplay : show the detectors, cluster and reconstructed tracks for the
event number given in additional option  
eventdisplaymult : do the same as eventdisplay but scan from event 0 to the
given event

* wrapper

`wrapper <config_wrapper> [config_file]`

config_wrapper : path to config file which contain path to the pedestal run and the data run in order to output the cosmic rays root file  
config_file : use a custom config_file for the cosmic bench caracteristics (use config_default.json by default)  

### Extend this soft :

Add a new type of detector :

* Add your type to Tomogrphy::det_type enum
* update Tomography::Static_Detector map
* Implement your `<type>_Detector` class (must inherit from `Detector` class)
* Implement your `<type>_Cluster` (must inherit from `Cluster` class)
* Implement your `<type>_Event` (must inherit from `Event` class)
* Update the config file corresponding to your bench caracteristcs
