LOCAL=$(HOME)/local

#
# Select compiler
#
CC=icc

#
# for Intel Compiler
#
WARNINGS=-Wall -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
OPT=-O2 #-funroll-loops
#OPT=-fast -funroll-loops
REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
LDLIBS=-L$(HOME)/local/lib -lfftw3 #-p

#
# for GNU Compiler
#
#WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
#OPT=-O2
#LDLIBS=-lm

#
# general
#
DEBUG=-g -debug full -DDEBUG -wn5 -mp # -p
INCLUDES=-I$(HOME)/local/include -I../../include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic


.PHONY: clean

all: corr

corr : corr.o
	$(CC) $(LDFLAGS) -o $@ corr.o $(LDLIBS)
	
clean :
	rm -f core* corr corr.o

