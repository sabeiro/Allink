LOCAL=$(HOME)/local

#
# Select compiler
#
CC=icc -openmp

#
# for Intel Compiler
#
WARNINGS=-Wall -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
OPT=-O2 #-funroll-loops
#OPT=-fast -funroll-loops
#REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
LDLIBS=-L$(HOME)/local/lib -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgslcblas -lgsl -liomp5 -lguide -lpthread -lm
INCLUDES=-I$(HOME)/local/include -I../../include

#
# general
#
DEBUG=-g -DDEBUG

#
# now put it all together
#
CFLAGS=-DMKL_ILP64 -std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic


.PHONY: clean

all: forcematch comforces friction acfriction gslfit

friction : friction.o error.o
	$(CC) $(LDFLAGS) -o $@ friction.o error.o $(LDLIBS)
	
acfriction : acfriction.o error.o
	$(CC) $(LDFLAGS) -o $@ acfriction.o error.o $(LDLIBS)
	
forcematch : forcematch.o error.o
	$(CC) -openmp $(LDFLAGS) -o $@ forcematch.o error.o $(LDLIBS)
	
comforces : comforces.o error.o
	$(CC) $(LDFLAGS) -o $@ comforces.o error.o $(LDLIBS)
	
gslfit : gslfit.o error.o
	$(CC) $(LDFLAGS) -o $@ gslfit.o error.o $(LDLIBS)
	
clean :
	rm -f core* forcematch forcematch.o error.o friction friction.o comforces comforces.o gslfit gslfit.o

