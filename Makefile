OBJ = bigrips reco toff_slope
#GO4SYS = /media/bae/data1/research/RIKEN/201605~06_EURICA/go4-5.2.0

all: $(OBJ)

#include $(GO4SYS)/Makefile.config

CXX	= g++

#LIBS            = -L/usr/local/lib -lXMLParser
ROOTCFLAGS	= $(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
TARTLIBS    = -L$(TARTSYS)/lib -lanacore -lanaroot -lanabrips -lanaloop -L/usr/local/lib -lXMLParser
#TARTLIBS    = -L$(TARTSYS)/lib $(TARTSYS)/lib/libanacore.so.0.0.0 -lanaloop -lanabrips -lanacore -lXMLParser
#OtherLIBS 	= -L/usr/local/lib -lXMLParser
INCLUDES	= -I$(TARTSYS)/include -I$(ROOTSYS)/include


#ifdef GO4_WIN32
#   GO4SYS = ../go4
#endif

bigrips: MakeFullBigRIPSTree.C BigRIPS.h
	$(CXX) -fPIC -DPIC $< -o $@ $(INCLUDES) $(ROOTLIBS) $(ROOTFLAGS) $(ROOTGLIBS) $(TARTLIBS) -std=c++11

reco: reco_main.cpp reco.cpp reco.h BigRIPS.h
	$(CXX) $< -o $@ $(INCLUDES) $(ROOTLIBS) $(ROOTFLAGS) $(ROOTGLIBS) $(TARTLIBS) -std=c++11 

toff_slope: toff_slope.cpp reco.cpp reco.h BigRIPS.h
	$(CXX) $< -o $@ $(INCLUDES) $(ROOTLIBS) $(ROOTFLAGS) $(ROOTGLIBS) $(TARTLIBS) -std=c++11 

clean: 
	rm -f $(OBJ) *.o

cleanbig:
	rm -f bigrips *.o
