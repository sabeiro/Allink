LOCAL=$(HOME)/local

#
# Select compiler
#
#CC = kinst-pomp mpicc
CC=mpicc

#
# for Intel Compiler
#
WARNINGS=-Wall -wd981 -wd810 -wd869 -wd1572
OPT=-O2
#OPT=-O1 #-funroll-loops
#OPT=-fast -funroll-loops
REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
LDLIBS=-L$(LOCAL)/lib -lmpfr -lgmp -lgsl -lgslcblas

#
# for GNU Compiler
#
#WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
#OPT=-O2
#LDLIBS=-lm

#
# general
#
DEBUG=#-g -debug full -DDEBUG -wn5 #-mp # -p
INCLUDES=-I$(LOCAL)/include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic


.PHONY: clean

all: dwham

dwham : dwham.o
	$(CC) $(LDFLAGS) -o $@ dwham.o $(LDLIBS)

clean :
	rm -f core* dwham dwham.o

