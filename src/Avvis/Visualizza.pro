#TEMPLETE = app
SOURCES = Visualizza.cpp ../VarData/VarDatFile.cpp
SOURCES += ElementiGrafici.cpp ElementiGraficiComm.cpp ElementiGraficiDis.cpp ElementiGraficiEl.cpp Finestra.cpp 
HEADERS += ../../include/VarDatFile.h ElementiGrafici.h
HEADERS += ../../include/Matematica.h ../../include/VarData.h
TARGET = Avvis 
QT += qt3support
#SRCMOC = moc_ElementiGrafici.cpp
#OBJMOC = moc_ElementiGrafici.o
#CONFIG += qt debug warn_on release
#INCLUDEPATH += ../include/ /usr/include/qt4/QtGui/ /usr/include/qt4/Qt/ ../../include
INCLUDEPATH += ../include/ ../../include
DESTDIR = $(HOME)/bin/
#unix:LIBS+= -L../lib/ -L../../lib -lMatematica -lVarData -lm -lfftw3 -lgsl 
unix:LIBS+= -L../lib/ -L../../lib -lMatematica -lVarData -lm -lfftw3
