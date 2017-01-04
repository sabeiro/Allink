LOCAL=$(HOME)/local

#
# Select compiler
#
CC=gcc

#
# for Intel Compiler
#
#WARNINGS=-Wall -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
#OPT=-O2 #-funroll-loops
#OPT=-fast -funroll-loops
#REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
#LDLIBS=-L$(HOME)/local/lib -lfftw3 -lgslcblas -lgsl #-p

#
# for GNU Compiler
#
WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
OPT=-O0
LDLIBS=-L$(HOME)/local/lib -lfftw3 -lgslcblas -lgsl -lm

#
# general
#
DEBUG=-g -DDEBUG # -p
INCLUDES=-I$(HOME)/local/include -I../../include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic


.PHONY: clean

all: normpot

normpot : normpot.o error.o
	$(CC) $(LDFLAGS) -o $@ normpot.o error.o $(LDLIBS)
	
clean :
	rm -f core* normpot normpot.o error.o

