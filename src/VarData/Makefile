include ../../include/MakeInclude.mk
HEADERS = ../../include/VarData.h ../../include/Matematica.h
OGGVARDATA = VarData.o VarDataEl.o VarDataCreate.o VarDataInterp.o VarDataString.o VarDataRead.o VarDataBackFold.o VarDataExp.o VarDataProfile.o VarDataWrite.o VarDataCGAL.o VarDataComm.o Cubo.o VarDataMarchCubes.o
OGGVARDATFILE = VarDatFile.o

all: libVarData UsaMatematica  #ElPoly Addyn

libVarData: ../../include/VarData.h $(OGGVARDATA) 
	ar rvs $(ALLINK_LIB)libVarData.a $(OGGVARDATA) 

UsaMatematica: UsaMatematica.o VarDatFile.o
	$(CC) $(CFLAGS) UsaMatematica.o VarDatFile.o $(LOCAL_INC) $(LIB_VARDATA) -o $(EXEC)$@
#	ar rvs $(LIBR) $(OGG)

ElPoly: $(OGGVARDATA)
	cd ../ElPoly/ && make

Addyn: $(OGGVARDATA)
	cd ../Addyn/ && make

.cpp.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.c.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

clean:
	rm *.o
