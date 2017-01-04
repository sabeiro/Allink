LOCAL=$(HOME)/local

#
# Select compiler
#
CC=icc
#CC=gcc

#
# for Intel Compiler
#
WARNINGS=-Wall -Wshadow -Wcheck -wd981 -wd810 -wd869 -wd1572
#OPT=-O0
OPT=-O2 #-funroll-loops
#OPT=-fast -funroll-loops
REPORT=-vec_report1 #-opt_report -opt_report_levelmax -opt_report_phaseall
#DEBUG=-g -debug full# -DDEBUG -wn5 -mp # -p
#LDLIBS=#-p

#
# for GNU Compiler
#
#WARNINGS=-Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wshadow -Wpointer-arith
#OPT=-O2
#LDLIBS=-lm
#DEBUG=-g

#
# general
#
INCLUDES=-I$(LOCAL)/include -I../../include

#
# now put it all together
#
CFLAGS=-std=c99 $(DEBUG) $(OPT) $(WARNINGS) $(REPORT) $(INCLUDES)
LDFLAGS=-i_dynamic

OBJECTS=ideal-chain.o

.PHONY: clean

all: ic

ic : $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

clean :
	rm -f core* ic $(OBJECTS)

