include ../../include/MakeInclude.mk
OGG =  rnorrexp.o Matematica.o MatematicaFunc.o MatematicaSign.o MatematicaInterp.o MatematicaFilter.o MatematicaAlgebra.o MatematicaVect.o MatematicaPlane.o MatematicaQuaternion.o randomlib.o MatematicaChar.o mt19937ar.o #nist_spblas.o

all: ../../lib/libMatematica.a UsaMatematica #VarData Addraw ElPoly Addyn Avvis

../../lib/libMatematica.a: $(OGG)
	ar rvs ../../lib/libMatematica.a $(OGG)

UsaMatematica: $(OGG)
	cd ../VarData/ && make UsaMatematica

VarData: $(OGG)
	cd ../VarData && make

Addraw: $(OGG)
	cd ../Addraw && make

ElPoly: $(OGG)
	cd ../ElPoly && make

Addyn: $(OGG)
	cd ../Addyn && make

.cpp.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.c.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.cc.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<


clean:
	rm *.o
