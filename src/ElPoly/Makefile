include ../../include/MakeInclude.mk
OGG= ElPoly.o ElPolyEl.o ElPolyDraw.o ElPolyDrawSurf.o ElPolyDrawControl.o ElPolyRepr.o ElPolyMeasure.o ElPolyProfile.o ElPolyProfDens.o ElPolyProfPre.o ElPolyDrawSkinCGAL.o ElPolyDrawTriangCGAL.o  ElPolyOutput.o  #ElPolyDrawSurfCGAL.o
CFLAGS = $(FLAG_DEBUG) $(FLAG_GL) $(FLAG_CGAL) $(FLAG_MPI)
HEADERS = ../../include/VarData.h ../../include/Matematica.h

all: ElPoly #libElPoly

ElPoly: $(OGG) $(HEADERS) ElPoly.h 
	$(CC) $(OGG) $(MY_FLAGS) $(LOCAL_INC) $(LIB_ELPOLY) -o $(EXEC)$@

libElPoly: ElPoly.h $(OGG) 
	ar rvs ../../lib/libElPoly.a $(OGG) 

.cpp.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.c.o:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -c $<

.h.gch:
	$(CC) $(LOCAL_INC) $(MY_FLAGS) -x $<

clean:
	rm *.o
