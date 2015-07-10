### Utiliser le soft

* Datareader :

`Datareader <config_file> {ped,data,analyse}`

* Multicluster

`Multicluster <config_file> {multicluster,SoB,dispersion,display,hough}`

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