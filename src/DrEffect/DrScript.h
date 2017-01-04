#ifndef DRSCRIPT_H
#define DRSCRIPT_H

#include <qapplication.h>
#include <qwidget.h>
#include <qtimer.h>
#include <qpainter.h>
#include <qpixmap.h>
#include <qdatetime.h>
#include <qslider.h>
#include <qlabel.h>
				 //#include <qvbox.h>
#include <Q3VBox>
#include <stdlib.h>
#include <math.h>
#include <qerrormessage.h>
#include <qmessagebox.h>
#include <qstatusbar.h>
				 //#include <qpointarray.h>
#include <Q3PointArray>
#include <qimage.h>
#include <qpixmap.h>
#include <qevent.h>
#include <qprinter.h>
				 //#include <qpaintdevicemetrics.h>
#include <Q3PaintDeviceMetrics>
#include <Q3ButtonGroup>
#include <QGLWidget>
#include <QLineEdit>
				 //#include <QMenuBar>
#include "../include/VarData.h"

class DrScript : public QWidget
{
  Q_OBJECT
    public slots:
  void IncrSlide();
  signals:
  public:
  DrScript();
  ~DrScript();
  void ExecScript(QPainter *p);
  void ReadSlide(QPainter *p);
 private:
  void DrMessage(const char * s, ...);
  int Slide;
  FILE *File2Read;
 protected:
  virtual void paintEvent(QPaintEvent *);
  //virtual void mouseMoveEvent(QMouseEvent *);
};

#endif//DRSCRIPT
