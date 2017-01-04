LOCAL=$(HOME)/local

#
# Select compiler
#
#CC = kinst-pomp mpicc
CC=mpicc
#CC=gcc -lmpi
#CC=icc -lmpi

#
# for Intel Compiler
#
#WARNINGS=-Wall -Wshadow -Wcheck -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
#OPT=-O2 #-funroll-loops
#OPT=-fast -funroll-loops
#REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
#LDLIBS=#-p
#DEBUG=#-g -debug full -DDEBUG -wn5 -mp # -p

#
# for GNU Compiler
#
WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
OPT=-O2
OPT=-O0 -DTYPE_MAX=2
LDLIBS=-lm

#
# general
#
DEBUG=-g -debug #-DDEBUG -DNDEBUG#full # -wn5 -mp # -p
#INCLUDES=-I/usr/local/fftw/3.2.1/include -I$(LOCAL)/include -I../../include
INCLUDES=-I$(LOCAL)/include -I../../include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic
#LDLIBS=-L/usr/local/fftw/3.2.1/lib -lfftw3 -L$(LOCAL)/lib -lconfuse
LDLIBS=-L$(LOCAL)/lib -lconfuse -lfftw3 

OBJECTS= exforce.o dynamics.o error.o main.o mdpd.o nblist.o rand.o rnemd.o stress.o surfpot.o usv.o usp.o vmd.o vv.o Peptide.o Rigid.o Tens.o Extern.o Widom.o Diff.o

.PHONY: clean

#all: CoarseLipid
all: fdmdpd2

fdmdpd2 : $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)


CoarseLipid: $(OBJECTS) 
	$(CC) $(LDFLAGS) -o $(HOME)/bin/CoarseLipid $(OBJECTS) $(LDLIBS)	
#	$(CC) $(CFLAGS) $(INCLUDES) -c -DPOLYMER_SYS mdpd.c
#	$(CC) $(LDFLAGS) -o $(HOME)/bin/CoarseMinimal Rigid.minimal.o $(OBJECTS) $(LDLIBS)
#	$(CC) $(CFLAGS) $(INCLUDES) -c mdpd.c

clean :
	rm -f core* fdmdpd2 $(OBJECTS)

