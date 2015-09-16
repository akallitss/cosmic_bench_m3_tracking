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

#Colors
GREEN=\033[0;32m
# No Color
NC=\033[0m

ECHO=/bin/echo -e

DataReader_obj_tmp = NewDataReader.o datareader.o ElecReader.o Tsignal_W.o dataline.o tomography.o Tanalyse_W.o event.o detector.o ray.o cluster.o point.o
DataReader_obj = $(patsubst %, $(ODIR)/%, $(DataReader_obj_tmp))

absorptionMap_obj_tmp = absorptionMap.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o
absorptionMap_obj = $(patsubst %, $(ODIR)/%, $(absorptionMap_obj_tmp))

tracking_obj_tmp = tracking.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o
tracking_obj = $(patsubst %, $(ODIR)/%, $(tracking_obj_tmp))

MultiCluster_obj_tmp = MultiCluster.o Signal.o detector.o event.o cluster.o Tanalyse_W.o ray.o point.o Tsignal_R.o tomography.o Tray.o
MultiCluster_obj = $(patsubst %, $(ODIR)/%, $(MultiCluster_obj_tmp))

live_obj_tmp = live.o liveDisplay.o datareader.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o
live_obj = $(patsubst %, $(ODIR)/%, $(live_obj_tmp))

AutoAlign_obj_tmp = AutoAlign.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o
AutoAlign_obj = $(patsubst %, $(ODIR)/%, $(AutoAlign_obj_tmp))

HV_Monitor_obj_tmp = HV_Monitor.o CAEN_comm.o tomography.o
HV_Monitor_obj = $(patsubst %, $(ODIR)/%, $(HV_Monitor_obj_tmp))

wrapper_obj_tmp = wrapper.o detector.o event.o cluster.o ray.o point.o datareader.o Tsignal_W.o dataline.o ElecReader.o tomography.o Tray.o
wrapper_obj = $(patsubst %, $(ODIR)/%, $(wrapper_obj_tmp))

carac_all_obj_tmp = carac_all.o carac.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o tomography.o
carac_all_obj = $(patsubst %, $(ODIR)/%, $(carac_all_obj_tmp))

$(shell mkdir -p $(ODIR))

#------------------------------------------------------------------------------

default: msg

all: exec

exec: echo_things tracking_b absorptionMap_b MultiCluster_b DataReader_b AutoAlign_b HV_Monitor_b wrapper_b carac_all_b

todo: live

echo_things:
	@$(ECHO) "$(GREEN)Building with flags:$(NC) $(CXXFLAGS)"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(CXX) -o $@ $(CXXFLAGS) -c $<

DataReader: echo_things DataReader_b

absorptionMap: echo_things absorptionMap_b

tracking: echo_things tracking_b

MultiCluster: echo_things MultiCluster_b

live: echo_things live_b

AutoAlign: echo_things AutoAlign_b

HV_Monitor: echo_things HV_Monitor_b

wrapper: echo_things wrapper_b

carac_all: echo_things carac_all_b

DataReader_b: $(DataReader_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o DataReader $(LDFLAGS)

absorptionMap_b: $(absorptionMap_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o absorptionMap $(LDFLAGS)

tracking_b: $(tracking_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o tracking $(LDFLAGS)

MultiCluster_b: $(MultiCluster_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o MultiCluster $(LDFLAGS)

live_b: $(live_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o live $(LDFLAGS)

AutoAlign_b: $(AutoAlign_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o AutoAlign $(LDFLAGS)

HV_Monitor_b: $(HV_Monitor_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o HV_Monitor $(LDFLAGS) -lcaenhvwrapper

wrapper_b: $(wrapper_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o wrapper $(LDFLAGS)

carac_all_b: $(carac_all_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o carac_all $(LDFLAGS)

msg:
	@$(ECHO) "Warning : you have to select the module you want to compile : tracking, carac_all, absorptionMap, MultiCluster, testCapa, DataReader, AutoAlign, HV_Monitor, wrapper or all"

.PHONY: clean cleanall default all exec todo

clean:
	rm -f $(ODIR)/*.o

cleanall: clean
	rm -f tracking absorptionMap MultiCluster testCapa DataReader AutoAlign HV_Monitor wrapper carac_all live
