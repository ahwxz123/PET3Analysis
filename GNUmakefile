#OBJ =   PET3_Coins.o
#EXE =   PET3_Coins
OBJ =   PET3Analysis.o
EXE =   PET3Analysis
#OBJ =   PET3Calibration_V1.o
#EXE =   PET3Calibration_V1
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include   
LDFLAGS = -L$(ROOTSYS)/lib        

#CXX = gcc -m32 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE
#####CXX = g++ -g -m32 -D_FILE_OFFSET_BITS=64  -D_LARGEFILE_SOURCE  -D_LARGEFILE64_SOURCE
#CXX = g++ -m32 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE
CXX = g++ 
FLAGS =  -Wall -std=c++11 -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE) 

$(EXE): $(OBJ)
	$(CXX) -o $(EXE) $(OBJ) $(ROOTFCLAGS) $(ROOTLIBS) -lSpectrum

%.o: %.cxx
	$(COMPILE)  $< 

clean:
	rm -f $(EXE) 
	rm -f $(OBJ) 
