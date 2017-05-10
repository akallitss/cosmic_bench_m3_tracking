IDIR          = include
SDIR          = src
ODIR          = obj
OUTPUTDIR     = ./
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux with egcs
WARNINGS      = -Wall -Wextra -Wchar-subscripts -Wundef -Wshadow -Wwrite-strings -Wsign-compare -Wunused -Wunused-parameter -Wuninitialized -Winit-self -Wpointer-arith -Wredundant-decls -Wformat-nonliteral -Wformat-zero-length -Wmissing-format-attribute -Wsequence-point -Wparentheses -Wmissing-declarations
WARNINGS     += -Wcast-align -Wformat=2 -Winline -Wmissing-include-dirs -Wunknown-pragmas -Wno-long-long -Wno-unused-function -std=gnu++11 -pedantic -pipe -fPIE -fstack-protector -fvisibility=hidden
CXX           = g++
CGNU          = gcc
#CXXFLAGS      = -O2 -Wall -Wextra -fexceptions -fPIC  $(ROOTCFLAGS) -fopenmp -I$(IDIR) -DUNIX -DLINUX
CXXFLAGS      = -g -O $(WARNINGS) -fexceptions -fPIC  $(ROOTCFLAGS) -I$(IDIR) -DUNIX -DLINUX
CFLAGS        = -g -O2 $(WARNINGS) -I$(IDIR) -pthread
LD            = g++
LIBS          = $(ROOTLIBS) -lNetx -lm -ldl -rdynamic 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 -lm -ldl -rdynamic -lpthread -lMinuit2 -lboost_system -lboost_filesystem
LDFLAGS       =  $(GLIBS)

#Colors
GREEN=\033[1;32m
RED=\033[1;31m
# No Color
NC=\033[0m

ECHO=/bin/echo -e

DataReader_obj_tmp = NewDataReader.o datareader.o ElecReader.o Tsignal_W.o dataline.o tomography.o Tanalyse_W.o event.o detector.o ray.o cluster.o point.o MT_tomography.o task/read_elec_task.o task/ped_task.o task/multicluster_task.o task/write_analyse_task.o task/read_live_task.o task/write_signal_task.o task/write_rays_task.o task/tracking_task.o Tray.o Pipes.o Platform.o
DataReader_obj = $(patsubst %, $(ODIR)/%, $(DataReader_obj_tmp))

absorptionMap_obj_tmp = absorptionMap.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o MT_tomography.o task/read_analyse_task.o task/tracking_task.o
absorptionMap_obj = $(patsubst %, $(ODIR)/%, $(absorptionMap_obj_tmp))

tracking_obj_tmp = tracking.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o dataline.o MT_tomography.o task/read_analyse_task.o task/tracking_task.o
tracking_obj = $(patsubst %, $(ODIR)/%, $(tracking_obj_tmp))

MultiCluster_obj_tmp = MultiCluster.o Signal.o detector.o event.o cluster.o Tanalyse_W.o ray.o point.o Tsignal_R.o tomography.o Tray.o MT_tomography.o task/read_signal_task.o task/ped_task.o task/write_analyse_task.o task/multicluster_task.o
MultiCluster_obj = $(patsubst %, $(ODIR)/%, $(MultiCluster_obj_tmp))

live_obj_tmp = live.o liveDisplay.o datareader.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o MT_tomography.o
live_obj = $(patsubst %, $(ODIR)/%, $(live_obj_tmp))

AutoAlign_obj_tmp = AutoAlign.o analyse.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o Tsignal_R.o tomography.o acceptanceFunction.o Tray.o MT_tomography.o task/read_analyse_task.o task/tracking_task.o
AutoAlign_obj = $(patsubst %, $(ODIR)/%, $(AutoAlign_obj_tmp))

HV_Monitor_obj_tmp = HV_Monitor.o CAEN_comm.o
HV_Monitor_obj = $(patsubst %, $(ODIR)/%, $(HV_Monitor_obj_tmp))

wrapper_obj_tmp = wrapper.o detector.o event.o cluster.o ray.o point.o datareader.o Tsignal_W.o dataline.o ElecReader.o tomography.o Tray.o MT_tomography.o task/read_live_task.o Pipes.o Platform.o
wrapper_obj = $(patsubst %, $(ODIR)/%, $(wrapper_obj_tmp))

carac_all_obj_tmp = carac_all.o carac.o Tanalyse_R.o event.o ray.o cluster.o detector.o point.o tomography.o MT_tomography.o
carac_all_obj = $(patsubst %, $(ODIR)/%, $(carac_all_obj_tmp))

$(shell mkdir -p $(ODIR)/task)
#------------------------------------------------------------------------------

default: msg

all: exec

exec: tracking absorptionMap MultiCluster DataReader AutoAlign HV_Monitor wrapper carac_all

todo: live

echo_things:
	@$(ECHO) "$(GREEN)Building with flags:$(NC) $(CXXFLAGS)"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(CXX) -o $@ $(CXXFLAGS) -c $<

$(ODIR)/%.o: $(SDIR)/%.c
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(CGNU) -o $@ $(CFLAGS) -c $<


DataReader: $(DataReader_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o DataReader $(LDFLAGS)

absorptionMap: $(absorptionMap_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o absorptionMap $(LDFLAGS)

tracking: $(tracking_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o tracking $(LDFLAGS)

MultiCluster: $(MultiCluster_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o MultiCluster $(LDFLAGS)

live: $(live_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o live $(LDFLAGS)

AutoAlign: $(AutoAlign_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o AutoAlign $(LDFLAGS)

HV_Monitor: $(HV_Monitor_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o HV_Monitor $(LDFLAGS) -lcaenhvwrapper

wrapper: $(wrapper_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o wrapper $(LDFLAGS)

carac_all: $(carac_all_obj)
	@$(ECHO) "$(GREEN)Building $@$(NC)"
	@$(LD) $^ -o carac_all $(LDFLAGS)

msg:
	@$(ECHO) "$(RED)Warning !$(NC) You have to select the module you want to compile : tracking, carac_all, absorptionMap, MultiCluster, testCapa, DataReader, AutoAlign, HV_Monitor, wrapper or all"

.PHONY: clean cleanall default all exec todo

clean:
	rm -f $(ODIR)/*.o $(ODIR)/task/*.o

cleanall: clean
	rm -f tracking absorptionMap MultiCluster testCapa DataReader AutoAlign HV_Monitor wrapper carac_all live
