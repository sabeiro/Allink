############USER#DEFINITION###################
#select a compiler
CC = g++ #mpic++
#debugging options
FLAG_DEBUG = -O0 -g -DDEBUG #-O2 -DVAR_DEBUG -DDRAW_DEBUG 
FLAG_GPROF = -pg
#standard libraries
LIB_STD = -lm 
#boost 
FLAG_BOOST = #-DUSE_BOOST
LIB_BOOST = #-lboost_iostreams
#Gnu scientific libraries
FLAG_GSL = #-D__GSL__
LIB_GSL = #-lgsl #-lgslcblas
#the fastest Fourier transform in the west
LIB_FFTW = -lfftw3
FLAG_FFTW = -DUSE_FFTW
#Local paths
LOCAL_LIB = -L../../lib/ -L$(HOME)/lib/ #-L/usr/local/lib
LOCAL_INC = -I$(HOME)/include/ -I../../include/ #-I/usr/local/include
EXEC = $(HOME)/bin/
#extra
FLAG_EXTRA = -fexceptions
LDFLAGS = -i_dynamic
#############SETTING#THE#FLAGS#####################
#collecting the flags
MY_FLAGS = $(FLAG_DEBUG) $(FLAG_BOOST)  $(FLAG_GSL)  $(FLAG_GPROF) $(FLAG_FFTW)
#to compile Matematica
MY_LIB = $(LOCAL_LIB)  $(LIB_GSL) $(LIB_FFTW) $(LIB_STD) $(LIB_BOOST)


OGG =  soliniziale.o
HEADERS = stopwatch.hpp matrix.h

all: PathOpt


PathOpt: $(OGG) $(HEADERS) 
	$(CC) $(OGG) $(MY_FLAGS) $(LOCAL_INC) $(LIB_ELPOLY) -o $(EXEC)$@

.cpp.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.c.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.cc.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<


clean:
	rm *.o
