IDIR          = include
SDIR          = src
ODIR          = obj
OUTPUTDIR     = ./
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux with egcs
CXX           = g++
CXXFLAGS      = -O2 -Wall -Wno-deprecated -fexceptions -fPIC  $(ROOTCFLAGS) -I$(IDIR) -DUNIX -DLINUX
#CXXFLAGS      = -g -O -Wall -Wno-deprecated -fexceptions -fPIC  $(ROOTCFLAGS) -I$(IDIR) -DUNIX -DLINUX
LD            = g++
LIBS          = $(ROOTLIBS) -lNetx -lm -ldl -rdynamic 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 -lm -ldl -rdynamic -lpthread -lMinuit2 -lcaenhvwrapper
LDFLAGS       =  $(GLIBS)

DataReader_obj_tmp = NewDataReader.o datareader.o header.o dataline.o tomography.o
DataReader_obj = $(patsubst %, $(ODIR)/%, $(DataReader_obj_tmp))

absorptionMap_obj_tmp = absorptionMap.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o acceptanceFunction.o
absorptionMap_obj = $(patsubst %, $(ODIR)/%, $(absorptionMap_obj_tmp))

tracking_obj_tmp = tracking.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o acceptanceFunction.o
tracking_obj = $(patsubst %, $(ODIR)/%, $(tracking_obj_tmp))

MultiCluster_obj_tmp = MultiCluster.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o datareader.o dataline.o tomography.o
MultiCluster_obj = $(patsubst %, $(ODIR)/%, $(MultiCluster_obj_tmp))

testCapa_obj_tmp = testCapa.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o datareader.o dataline.o tomography.o
testCapa_obj = $(patsubst %, $(ODIR)/%, $(testCapa_obj_tmp))

live_obj_tmp = live.o liveDisplay.o datareader.o header.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o
live_obj = $(patsubst %, $(ODIR)/%, $(live_obj_tmp))

AutoAlign_obj_tmp = AutoAlign.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o acceptanceFunction.o
AutoAlign_obj = $(patsubst %, $(ODIR)/%, $(AutoAlign_obj_tmp))

HV_Monitor_obj_tmp = HV_Monitor.o CAEN_comm.o
HV_Monitor_obj = $(patsubst %, $(ODIR)/%, $(HV_Monitor_obj_tmp))

#------------------------------------------------------------------------------

all: dir exec

dir: $(ODIR)

$(ODIR):
	mkdir -p $(ODIR)

exec: tracking absorptionMap MultiCluster testCapa DataReader live AutoAlign HV_Monitor

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $<

DataReader: $(DataReader_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

absorptionMap: $(absorptionMap_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

tracking: $(tracking_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

MultiCluster: $(MultiCluster_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

testCapa: $(testCapa_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

live: $(live_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

AutoAlign: $(AutoAlign_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

HV_Monitor: $(HV_Monitor_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
