#ifndef DREFFECT_H
#define DREFFECT_H

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
#include "../include/Matematica.h"
#include "../include/Draw.h"
#include "../include/VarDatFile.h"

class DrEffect:public Draw{
 private:
 public:
  DrEffect();
  /// Applies a filter to pixel
  void EffectFilter();
  /// Shifting of pixels
  void EffectMotion();
  /// Applies a random disposition of blocks
  void EffectMC();
  /// Defocus the images in blocks
  void EffectCoarseGrain(int Grana);
  /// Increase the resolution
  void EffectIncrease();
  /// Create an animation
  void Run();
  void Initialize();
  void Histo();
  void NablaPhi();
  void BlackWhite();
  void DReshape(int weight,int height);
  void DrEkeyboard(unsigned char key);
  void PrintIntensity();
  int Contrast();
  int Binary(int Mode);
  int Noise(unsigned char *Picture,unsigned char *OutPicture,int width,int height);
  int SwapBlocks(unsigned char *Picture,unsigned char *OutPicture,int width,int height,PERMUTE *Perm,int NGridw,int NGridh,int *Sequence1,int NPartition);
  int ShiftBlocks(unsigned char *Picture,unsigned char *OutPicture,int width,int height,int *Sequence,int NGridw,int NGridh,int *Sequence1,int NPartition);
  int Discretize(unsigned char *Picture,unsigned char *OutPicture,int ImWidth,int ImHeight,int WWidth,int WHeight,int Blockw,int Blokh);
  int Edges(unsigned char *Picture,unsigned char *OutPicture,int width,int height);
  int IncreaseResolution(unsigned char *Picture,unsigned char *OutPicture,int ImWidth,int ImHeight,int Times);
  int BuffSize(){return ImWidth*ImHeight;};
  int NChar;
  Matematica *Mat;
  unsigned char Median[4];
  unsigned char Quart1[4];
  unsigned char Quart3[4];
  double **Hist;
  double **Phi;
  int **Patch;
};
class Animation : public QGLWidget
{
  Q_OBJECT
    public slots:
  void setXRotation(int angle);
  void setYRotation(int angle);
  void setZRotation(int angle);
  void NomeFile(const QString &);
  void NomeFile(char *ExtName);
  void Open();
  void Filter();
  void Motion();
  void Run();
  void FilterMC();
  void ImpGrana(int Grana);
  void FilterCoarseGrain();
  void FilterIncrease();
  void FilterContrast();
  void Histo();
  void NablaPhi();
  void BlackWhite();
  void Binary();
  void Punta();
  void IncrSlide();
 signals:
  void xRotationChanged(int angle);
  void yRotationChanged(int angle);
  void zRotationChanged(int angle);
  void PuntaCoord(VarDatFile *v1);
  void PuntaCoord(double **,int,int,int);
  public:
  Animation(QWidget *parent);
  void InitConstant();
  void Menu();
  void Particle();
  void keyboard(unsigned char key, int x, int y);
  VarDatFile *v1;
  double **st;
 private:
  void DrMessage(const char * s, ...);
  char *FileName;
  double ScaleUn;
  double Edge[3];
  int xRot;
  int yRot;
  int zRot;
  int Grana;
  int Slide;
  int NMass;
  int NVar;
  int Valori;
  QPoint lastPos;
  int MainWindow;
 protected:
  virtual void initializeGL();
  virtual void paintEvent(QPaintEvent (event));
  virtual void paintGL();
  virtual void resizeGL(int width,int height);
  void mousePressEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
};




#endif//DREFFECT
