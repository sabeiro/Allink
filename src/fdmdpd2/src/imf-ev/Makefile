LOCAL=$(HOME)/local

#
# Select compiler
#
#CC = kinst-pomp mpicc
CC=gcc

#
# for Intel Compiler
#
#WARNINGS=-Wall -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
#OPT=-Os #-funroll-loops
#OPT=-fast -funroll-loops -openmp
#REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
#LDLIBS=-L$(LOCAL)/lib -lfftw3 #-lguide -lpthread

#
# for GNU Compiler
#
WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
OPT=-O0
LDLIBS=-L$(MKLROOT)/lib/em64t -lmkl -lmkl_solver -lguide -lm -L$(LOCAL)/lib -lgslcblas -lgsl -lpthread

#
# general
#
DEBUG=-g
INCLUDES=-I$(LOCAL)/include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)


.PHONY: clean

all: imf-ev sl-ev solver solver2 imf2

imf-ev : imf-ev.o
	$(CC) $(LDFLAGS) -o $@ imf-ev.o $(LDLIBS)
	
imf2 : imf2.o
	$(CC) $(LDFLAGS) -o $@ imf2.o $(LDLIBS)

sl-ev : sl-ev.o
	$(CC) $(LDFLAGS) -o $@ sl-ev.o $(LDLIBS)

solver : solver.o
	$(CC) $(LDFLAGS) -o $@ solver.o $(LDLIBS)

clean :
	rm -f core* imf-ev imf-ev.o sl-ev sl-ev.o solver.o

