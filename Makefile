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
#CXXFLAGS      = -O2 -Wall -Wno-deprecated -fexceptions -fPIC  $(ROOTCFLAGS) -I$(MYINCLUDE) 
CXXFLAGS      = -g -O -Wall -Wno-deprecated -fexceptions -fPIC  $(ROOTCFLAGS) -I$(MYINCLUDE)
LD            = g++
LIBS          = $(ROOTLIBS) -lNetx -lm -ldl -rdynamic 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -L$(MYLIB) -lXpm -lX11 -lm -ldl -rdynamic -lpthread -lMinuit
LDFLAGS       =  $(GLIBS)
SOFLAGS       = -shared -fPIC

#------------------------------------------------------------------------------

all: exec

exec: tracking absorptionMap MultiCluster testCapa DataReader live

execDict: trackingDict absorptionMapDict MultiClusterDict testCapaDict liveDict

lib: libAnalyse.so

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

DataReader: NewDataReader.o datareader.o header.o dataline.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

absorptionMap: absorptionMap.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

tracking: tracking.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

MultiCluster: MultiCluster.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o datareader.o dataline.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

testCapa: testCapa.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o datareader.o dataline.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

live: live.o liveDisplay.o datareader.o header.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o
	$(LD) $^ -o $@ $(LDFLAGS)

absorptionMapDict: absorptionMap.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

trackingDict: tracking.o analyse.o T.o event.o ray.o cluster.o detector.o point.o Tsignal.o tomography.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

MultiClusterDict: MultiCluster.o signal.o detector.o event.o cluster.o Tanalyse.o Tsignal.o ray.o point.o datareader.o dataline.o tomography.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

testCapaDict: testCapa.o signal.o detector.o event.o cluster.o Tanalyse.o ray.o point.o Tsignal.o datareader.o dataline.o tomography.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

liveDict: live.o liveDisplay.o datareader.o header.o dataline.o detector.o event.o cluster.o ray.o point.o tomography.o MyDict.o
	$(LD) $^ -o $@ $(LDFLAGS)

libAnalyse.so: analyse.o T.o event.o ray.o cluster.o detector.o point.o Tanalyse.o Tsignal.o signal.o tomography.o MyDict.o
	$(CXX) $(SOFLAGS) $(LDFLAGS) -o $@ $^ 

MyDict.cpp: analyse.h T.h event.h ray.h cluster.h detector.h point.h Tanalyse.h Tsignal.h signal.h acceptanceFunction.h tomography.h Linkdef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

clean:
	rm -f *.o *.so *Dict* *dict* Linkdef absorptionMap tracking MultiCluster acceptanceFunction DataReader live
