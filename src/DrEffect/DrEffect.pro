#TEMPLETE = app
SOURCES += DrEffect.cpp DrFinestra.cpp DrOpenGl.cpp Animation.cpp DrScript.cpp
SOURCES += DrDefinition.cpp
SOURCES +=  ../Visualizza/ElementiGrafici.cpp  ../Visualizza/ElementiGraficiDis.cpp  ../Visualizza/ElementiGraficiComm.cpp  ../Visualizza/ElementiGraficiEl.cpp ../VarData/VarDatFile.cpp
HEADERS += ../../include/Draw.h DrEffect.h DrScript.h
HEADERS += ../../include/Matematica.h 
HEADERS += ../Visualizza/ElementiGrafici.h ../../include/VarDatFile.h
TARGET = DrEffect
DEFINES += USE_GL
QT += qt3support
QT += opengl
SRCMOC = moc_DrEffect.cpp
OBJMOC = moc_DrEffect.o
#CONFIG += qt debug warn_on release
INCLUDEPATH += ../include/ /usr/include/qt4/QtGui/ /usr/include/qt4/Qt/ ../../include
DESTDIR = ../../bin/
unix:LIBS+= -L../../lib/ -L$(HOME)/lib/ 
unix:LIBS+= -lDraw -lglut -lX11 -lXi -lXmu -lGL -lGLU -ltiff -lpng -lpngwriter -lz -lfreetype 
unix:LIBS+= -lVarData -lMatematica -lfftw3 -lCGAL -lgsl -lgslcblas -lm 
