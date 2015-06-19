IDIR          = include
SDIR          = src
ODIR          = obj
OUTPUTDIR     = ./
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux with egcs
WARNINGS      = -Wall -Wextra -Wchar-subscripts -Wundef -Wshadow -Wwrite-strings -Wsign-compare -Wunused -Wunused-parameter -Wuninitialized -Winit-self -Wpointer-arith -Wredundant-decls -Wformat-nonliteral -Wformat-zero-length -Wmissing-format-attribute -Wsequence-point -Wparentheses -Wmissing-declarations
CXX           = g++
#CXXFLAGS      = -O2 -Wall -Wextra -fexceptions -fPIC  $(ROOTCFLAGS) -fopenmp -I$(IDIR) -DUNIX -DLINUX
CXXFLAGS      = -g -O $(WARNINGS) -fexceptions -fPIC  $(ROOTCFLAGS) -I$(IDIR) -DUNIX -DLINUX
LD            = g++
LIBS          = $(ROOTLIBS) -lNetx -lm -ldl -rdynamic 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 -lm -ldl -rdynamic -lpthread -lMinuit2
LDFLAGS       =  $(GLIBS)

DataReader_obj_tmp = NewDataReader.o datareader.o ElecReader.o Tsignal_W.o dataline.o tomography.o
DataReader_obj = $(patsubst %, $(ODIR)/%, $(DataReader_obj_tmp))

absorptionMap_obj_tmp = absorptionMap.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o
absorptionMap_obj = $(patsubst %, $(ODIR)/%, $(absorptionMap_obj_tmp))

tracking_obj_tmp = tracking.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o
tracking_obj = $(patsubst %, $(ODIR)/%, $(tracking_obj_tmp))

MultiCluster_obj_tmp = MultiCluster.o Signal.o detector.o event.o cluster.o Tanalyse_W.o ray.o point.o Tsignal_R.o tomography.o Tray.o
MultiCluster_obj = $(patsubst %, $(ODIR)/%, $(MultiCluster_obj_tmp))

testCapa_obj_tmp = testCapa.o Signal.o detector.o event.o cluster.o Tanalyse_W.o ray.o point.o Tsignal_R.o tomography.o Tray.o
testCapa_obj = $(patsubst %, $(ODIR)/%, $(testCapa_obj_tmp))

live_obj_tmp = live.o liveDisplay.o datareader.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o
live_obj = $(patsubst %, $(ODIR)/%, $(live_obj_tmp))

AutoAlign_obj_tmp = AutoAlign.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o
AutoAlign_obj = $(patsubst %, $(ODIR)/%, $(AutoAlign_obj_tmp))

HV_Monitor_obj_tmp = HV_Monitor.o CAEN_comm.o
HV_Monitor_obj = $(patsubst %, $(ODIR)/%, $(HV_Monitor_obj_tmp))

wrapper_obj_tmp = wrapper.o detector.o event.o cluster.o ray.o point.o datareader.o Tsignal_W.o dataline.o ElecReader.o tomography.o Tray.o
wrapper_obj = $(patsubst %, $(ODIR)/%, $(wrapper_obj_tmp))

carac_all_obj_tmp = carac_all.o carac.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o tomography.o
carac_all_obj = $(patsubst %, $(ODIR)/%, $(carac_all_obj_tmp))

$(shell mkdir -p $(ODIR))

#------------------------------------------------------------------------------

default: msg

all: exec

exec: tracking absorptionMap MultiCluster testCapa DataReader AutoAlign HV_Monitor wrapper carac_all

todo: live

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
	$(LD) $^ -o $@ $(LDFLAGS) -lcaenhvwrapper

wrapper: $(wrapper_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

carac_all: $(carac_all_obj)
	$(LD) $^ -o $@ $(LDFLAGS)

msg:
	@echo "Warning : you have to select the module you want to compile : tracking, carac_all, absorptionMap, MultiCluster, testCapa, DataReader, AutoAlign, HV_Monitor, wrapper or all"

.PHONY: clean cleanall default all exec todo

clean:
	rm -f $(ODIR)/*.o

cleanall: clean
	rm -f tracking absorptionMap MultiCluster testCapa DataReader AutoAlign HV_Monitor wrapper carac_all live
