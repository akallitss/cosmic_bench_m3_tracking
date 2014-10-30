ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

OutPutOpt     = -o

MYINCLUDE     = ~/.
MYLIB         = $(PPATH)/lib
OUTPUTDIR     = ./
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)



# Linux with egcs
CXX           = g++
#CXXFLAGS      = -g -O -Wall -Wno-deprecated -fno-exceptions -fPIC  $(ROOTCFLAGS) -I$(MYINCLUDE) 
CXXFLAGS      = -g -O -Wall -Wno-deprecated -fexceptions -fPIC  $(ROOTCFLAGS) -I$(MYINCLUDE)
LD            = g++
LIBS          = $(ROOTLIBS) -lNetx -lm -ldl -rdynamic 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -L$(MYLIB) -lXpm -lX11 -lm -ldl -rdynamic -lpthread -lMinuit
LDFLAGS       =  $(GLIBS)
SOFLAGS       = -shared -fPIC

#------------------------------------------------------------------------------

all: exec

exec: tracking absorptionMap MultiCluster

execDict: trackingDict absorptionMapDict MultiClusterDict

lib: libAnalyse.so

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

absorptionMap: absorptionMap.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o
	$(LD) $^ -o $@ $(LDFLAGS)

tracking: tracking.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o
	$(LD) $^ -o $@ $(LDFLAGS)

MultiCluster: MultiCluster.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o
	$(LD) $^ -o $@ $(LDFLAGS)

absorptionMapDict: absorptionMap.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

trackingDict: tracking.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

MultiClusterDict: MultiCluster.o signal.o detector.o event.o cluster.o Tanalyse.o Tsignal.o ray.o point.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

libAnalyse.so: analyse.o T.o event.o ray.o cluster.o detector.o point.o Tanalyse.o Tsignal.o signal.o MyDict.o
	$(CXX) $(SOFLAGS) $(LDFLAGS) -o $@ $^ 

MyDict.cpp: analyse.h T.h event.h ray.h cluster.h detector.h point.h Tanalyse.h Tsignal.h signal.h Linkdef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

clean:
	rm -f *.o *.so *Dict* *dict* Linkdef absorptionMap tracking MultiCluster

