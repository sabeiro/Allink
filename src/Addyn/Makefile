include ../../include/MakeInclude.mk
HEADERS = ../../include/VarData.h ../../include/Matematica.h
OGG = Dinamica.o Forces.o ForcesDraw.o ForcesCreate.o ForcesForceField.o ForcesIntegration.o ForcesBoundary.o ForcesLoop.o ForcesTens.o ForcesCalcEnergies.o

all: Addyn

Addyn: $(OGG) $(HEADERS) 
	$(CC) $(OGG) $(MY_FLAGS) $(LOCAL_INC) $(LIB_DINAMICA) -o $(EXEC)$@

.cpp.o:
	$(CC) $(LOCAL_INC)  $(MY_FLAGS) -c $<

.c.o:
	$(CC) $(MY_FLAGS) -c $<

clean:
	rm *.o
