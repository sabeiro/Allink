all:
	cd Matematica && make
	cd VarData && make 
	cd Addraw && make 
	cd ElPoly && make
	cd Addyn && make
	cd Avvis && make

.cpp.o:
	$(CC) $(INC) $(CFLAGS) -c $<

.c.o:
	$(CC) $(INC) $(CFLAGS) -c $<

clean:
	cd Matematica && make clean
	cd VarData && make clean
	cd ElPoly && make clean
	cd Addyn && make clean
	cd Addraw && make clean
	cd Avvis && make clean
	
	rm *.o
