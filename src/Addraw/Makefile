include ../../include/MakeInclude.mk
LIB = $(LOCAL_LIB) $(LIB_MPI) $(LIB_CGAL) -lMatematica $(LIB_FFTW) $(LIB_GSL) -lm
LIBGL = $(LIB_GL)
CFLAGS = $(FLAG_DEBUG) $(FLAG_GL) $(FLAG_CGAL) $(FLAG_MPI) $(FLAG_GSL)
HEADERS = ../../include/VarData.h ../../include/Matematica.h
LDFLAGS=-i_dynamic
OGGDRAW = Draw.o DrawFile.o DrawDefinition.o DrawControl.o DrawScene.o

all:  Addraw ProvaAddraw DrImage

ProvaAddraw: ../../include/Draw.h $(OGGDRAW) ProvaDraw.o
	$(CC) ProvaDraw.o $(CFLAGS) -lDraw $(LIBGL) $(LIB) -o $(EXEC)$@

Addraw: ../../include/Draw.h $(OGGDRAW) 
	ar rvs $(ALLINK_LIB)/libDraw.a $(OGGDRAW)

ReadWritePng: ReadWritePng.cpp
	$(CC) -lpng ReadWritePng.cpp -o ../bin/$@

DrImage: DrImage.o
	g++ DrImage.o -o ../../bin/$@ `freetype-config --cflags`  $(LIB) $(LIB_PNG) 

ElPoly:
	cd ../ElPoly && make

.cpp.o:
	$(CC) $(LOCAL_INC) $(CFLAGS) -c $<

.c.o:
	$(CC) $(LOCAL_INC) $(CFLAGS) -c $<

clean:
	rm *.o
