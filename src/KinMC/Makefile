Scrive = ../bin/Rinato
Legge = ../bin/Elabora
Disegna = ../bin/Disegna
LIB = -lm
LIBGL = -lglut
CC = g++
OGG = Rinato.o ran1.o Cubo.o Particella.o
HEADER = Variabili.h Cubo.h Particella.h Rinato.h
SORG = Rinato.cpp ran1.c Cubo.cpp Particella.cpp
FLAG = -O2

all: $(Scrive) $(Legge) $(Disegna)

$(Scrive): $(OGG) $(HEADER)
	$(CC) $(FLAG) $(LIB) $(OGG) -o $(Scrive)

$(Legge): Elabora.o Elabora.h
	$(CC) $(FLAG) $(LIB) Elabora.o -o $(Legge)

$(Disegna): Disegna.o Disegna.h
	$(CC) $(FLAG) $(LIB) $(LIBGL) Disegna.o -o $(Disegna)

.c.o:
	$(CC) -c $<

.cpp.o:
	$(CC) -c $<

clean:
	rm *.o *~

