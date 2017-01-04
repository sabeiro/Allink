INC = -I../include/
LIB=  -L../lib/ -lMatematica -lm -lfftw3 
LIBR= ../lib/libMatematica.a
CC = g++
FLAG = -O2
PROG = ../bin/Segnali
OGG = Segnali.o VarSegnali.o 
SORG = Segnali.cpp VarSegnali.ccp
HEADERS = Segnali.h VarSegnali.h 

all: $(PROG)

$(PROG): $(OGG) $(HEADERS) ../src/Matematica.cpp
	$(CC) $(FLAG) $(OGG) $(INC) $(LIB) -o $(PROG)

.cpp.o:
	$(CC) -c $<

.c.o:
	$(CC) -c $<

../src/Matematica.cpp:
	cd ../src/
	make
	cd ../Segnali/

clean:
	rm *.o
