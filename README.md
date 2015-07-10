### Utiliser le soft

* Datareader :

`Datareader <config_file> {ped,data,analyse}`

config_file : path to config file which contain the cosmic bench caracteristics and different file path

ped : use this option to build signal file with raw, ped and corr branches, and calculate the Ped.dat and RMS.dat

data : use this option to build signal file with raw, ped and corr branches using existing Ped.dat

analyse : use this option to build analyse file using existing Ped.dat and RMS.dat

* Multicluster

`Multicluster <config_file> {multicluster,SoB,dispersion,display,hough}`

config_file : path to config file which contain the cosmic bench caracteristics and different file path

multicluster : use this option to build analyse file using existing Ped.dat and RMS.dat and signal file

SoB : use this option to display signal and noise amplitude

dispersion : use this option to display some correlation between the detectors signal

display : use this option to display the signal shpape for each detector

hough : use this option to study hough tracking style


### Etendre le soft :

Ajouter un détecteur :

* Ajouter votre type de detecteur à l'enum Tomogrphy::det_type
* Mettre à jour Tomography::operator<<(ostream,det_type)
* Mettre à jour Tomography::Static_Detector
* Implementer sa classe Detector (hérite de Detector)
* Implementer sa classe Cluster (hérite de Cluster)
* Implementer sa classe Event (hérite de Event)
* Mettre à jour les classes de stockage T{signal,analyse}_{R,W} pour utiliser le stockage disque par fichier ROOT (fonctions constructeur, init, fill et get)
* Mettre à jour le fichier de config pour commencer à utiliser le soft