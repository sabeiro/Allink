############USER#DEFINITION###################
#select a compiler
CC = g++ #mpic++
#debugging options
FLAG_DEBUG = -O0 -g -DDEBUG #-O2 -DVAR_DEBUG -DDRAW_DEBUG 
FLAG_GPROF = -pg
#standard libraries
LIB_STD = -lm 
#boost 
FLAG_BOOST = #-DUSE_BOOST
LIB_BOOST = #-lboost_iostreams
#graphic libraries
FLAG_PNG = -DUSE_PNG `freetype-config --cflags`
LIB_PNG = -lpng -lpngwriter -lz -lfreetype
FLAG_TIFF = -DUSE_TIFF
LIB_TIFF = -ltiff
FLAG_GL = -DUSE_GL $(FLAG_PNG) $(FLAG_TIFF) -D__glut_h__
LIB_GL = -lDraw -lglut -lX11 -lGL -lGLU $(LIB_TIFF) $(LIB_PNG) 
#CGAL libraries
FLAG_CGAL = -DUSE_CGAL
LIB_CGAL = -lCGAL
#Gnu scientific libraries
FLAG_GSL = -D__GSL__
LIB_GSL = -lgsl #-lgslcblas
#the fastest Fourier transform in the west
LIB_FFTW = -lfftw3
FLAG_FFTW = -DUSE_FFTW
#Message passing interface
FLAG_MPI = #-DUSE_MPI
LIB_MPI = 
#Local paths
LOCAL_LIB = -L../../lib/ -L$(HOME)/lib/ #-L/usr/local/lib
LOCAL_INC = -I$(HOME)/include/ -I../../include/ #-I/usr/local/include
EXEC = $(HOME)/bin/
ALLINK_LIB = ../../lib/
ALLINK_INCLUDE = ../../include/
#extra
FLAG_EXTRA = -fexceptions
LDFLAGS = -i_dynamic
#############SETTING#THE#FLAGS#####################
#collecting the flags
MY_FLAGS = $(FLAG_DEBUG) $(FLAG_BOOST) $(FLAG_PNG)  $(FLAG_GL) $(FLAG_CGAL) $(FLAG_GSL) $(FLAG_MPI) $(FLAG_GPROF) $(FLAG_FFTW)
#to compile Matematica
LIB_MATEMATICA = $(LOCAL_LIB) $(LIB_MPI) $(LIB_CGAL) $(LIB_GSL) $(LIB_FFTW) $(LIB_STD) 
#to compile Matematica and VarData
LIB_VARDATA = $(LOCAL_LIB) $(LIB_MPI) $(LIB_CGAL) -lMatematica $(LIB_GSL) $(LIB_FFTW) $(LIB_STD) $(LIB_BOOST)
#to compile ElPoly
LIB_ELPOLY = $(LOCAL_LIB) $(LIB_GL) $(LIB_MPI) $(LIB_CGAL) -lVarData -lMatematica $(LIB_GSL) $(LIB_FFTW) $(LIB_STD) $(LIB_BOOST)
#to compile Dinamica
LIB_DINAMICA = $(LOCAL_LIB)  $(LIB_GL) $(LIB_MPI) $(LIB_CGAL) -lVarData -lMatematica $(LIB_GSL) $(LIB_FFTW) $(LIB_STD) $(LIB_BOOST)

