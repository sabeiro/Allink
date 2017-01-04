/***********************************************************************
ElementiGrafici: Header file for the ElementiGrafici class
 Copyright (C) 2008 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
***********************************************************************/
#ifndef ELEMENTIGRAFICI_H
#define ELEMENTIGRAFICI_H

#include <qprinter.h>
#include <qmessagebox.h>
#include <qimage.h>
#include <qerrormessage.h>
#include <Q3VBox>
#include <qslider.h>
#include <qtextdocument.h>
//#include <qvbox.h>

/* #include <qwidget.h> */
/* #include <Q3PointArray> */

/* #include <qtimer.h> */
/* #include <qpixmap.h> */
/* #include <qdatetime.h> */
/* /\* /\\* #include <stdlib.h> *\\/ *\/ */
/* /\* /\\* #include <math.h> *\\/ *\/ */
/* /\* 				 //#include <qpointarray.h> *\/ */
/* #include <qpixmap.h> */
/* 				 //#include <qpaintdevicemetrics.h> */
/* #include <QLineEdit> */
//#include <QMenuBar>
#include <VarDatFile.h>

#define LINEA_NIENTE 0x00000
#define LINEA_PUNTO  0x00001
#define LINEA_TRATTO 0x00002
#define LINEA_PUNTOETRATTO 0x00003
#define LINEA_ERRORE 0x00008

;
/// Drawing flags
enum VisDis{
  /// Single signal
  DIS_SEGNALE =  0x00002,
  /// All the signals
  DIS_TUTTI   =  0x00004,
  /// Fit label
  DIS_MOMENTI =  0x00008,
  /// Fit
  DIS_BLU     =  0x00010,
  /// Linear fit
  DIS_RETTA   =  0x00020,
  /// Parabola fit
  DIS_PARABOLA=  0x00040,
  /// Esponential fit
  DIS_EXP     =  0x00080,
  /// Gaussian
  DIS_GAUSS   =  0x00100,
  /// Points
  DIS_PUNTI   =  0x00200,
  /// Error bar
  DIS_ERR     =  0x00400,
};
#define RIS_UNO       0x00002
#define RIS_TUTTI     0x00004

#define INIT_NUM  15
;
/// Contains the function to interacts with the data stored in VarDatFile
class ElementiGrafici:public QWidget{
  Q_OBJECT
    public:
  /// General constructor
  ElementiGrafici( QWidget *parent=0,const char *name=0);
  /// Destructor
  ~ElementiGrafici();
  QSizePolicy Dimensionamento() const;
  //vector <double> st;
  /// Number of points per array
  int NMax;
  /// Current number of visualisation (debugging)
  int NDis;
  /// Class for the file handle
  VarDatFile *v1;
  /// Name of the current file opened
  char *nomeFile;
  /// Name of the config file
  char *nomeConf;
  /// Name of the output file
  char *nomeSalva;
  /// Title name
  QString nomeTit;
  /// X axis label
  QString nomeEtX;
  /// Y axis label
  QString nomeEtY;
 public:
  /// Choose the input file
  void ChooseDataFile(char* FileName);
  /// Choose the configuration file
  void ChooseConfFile(char* FileName);
  public slots:
  /// Set NBin
  void ImpNBin(int n);
  /// Set the bins for the running average
  void ImpNMobile(int n);
  /// Set the bins for the point correlation
  void ImpNCorrela(int n);
  /// Set the minimum visualisation point
  void ImpNVisMin(int n);
  /// Set the maximum visualisation point
  void ImpNVisMax(int n);
  /// Set the minimum elaboration point
  void ImpNElMin(int n);
  /// Set the maximum elaboration point
  void ImpNElMax(int n);
  void ImpNElMinY(int n);
  void ImpNElMaxY(int n);
  void ImpInter(int n);
  void ImpInterEl(int n);
  void ImpCoordX(int n);
  void ImpCoordY(int n);
  /// Set the y error bar
  void ImpCoordDY(int n);
  /// Type for coordinate change
  /* typedef void(ElementiGrafici::*IMP_COORD)(int n); */
  /// Pointer to a coordinate change function
  /* IMP_COORD Imp_Coord; */
  /* /// Pointer to a generic function */
  /* void ImpCoord(int n){return (*this.*Imp_Coord)(n);} */
  void ImpNVarX(int n,int m);
  void ImpNVarY(int n,int m);
  void ImpInter(int n,int m);
  void ImpInterEl(int n,int m);
  void ImpInterY(int n,int m);
  void Apri();
  void Apri(char **argv,int *FileList,int NFile);
  void ApriExt(double **ExtSt,int ExtNMax,int ExtNVar,int ExtNBin);
  void ApriExt(VarDatFile *Extv1);
  void Aggiungi();
  void PuntiSegnale();
  void PuntiDistribuzione();
  void PuntiSommaSegnali();
  void PuntiInterRett();
  void PuntiInterExp();
  void PuntiInterGauss();
  void PuntiParabola();
  void PuntiMediaMob();
  void PuntiCorrelaADue();
  void PuntiAutosimilarita();
  void PuntiIntegrale();
  void PuntiDerivata();
  void PuntiVarie();
  void PuntiSpettro();
  void PuntiNormalizza();
  void PuntiModulo();
  void PuntiRadice();
  void PuntiAutocor();
  void PuntiSum();
  void DisegnaLinee();
  void DisegnaPunti();
  void DisegnaGriglia();
  void DisegnaLogx();
  void DisegnaLogy();
  void SeRiscala();
  void SeRiscalaTutto();
  void NSet();
  void StampaFile();
  void PuntaExt(double **ExtSt,int nVar);
  void SulSegnale(bool);
  void SulGrafico(bool);
  void SuX(bool);
  void SuY(bool);
  void SuDX(bool);
  void SuDY(bool);
  void Ridisegna();
  void InfoStato(const char *);
  void InfoSequenza(const char *);
  void Salva();
  void NomeFile(const QString &);
  void NomeConf(const QString &);
  void NomeEntrata(const QString &);  
  void NomeSalva(const QString &);
  void NomeUscita(const QString &);
  void NomeTitolo(const QString &);
  void NomeEtichettaX(const QString &);
  void NomeEtichettaY(const QString &);
  void NomeTit(const QString &);
  void NomeEtX(const QString &);
  void NomeEtY(const QString &);
  
  private slots:

  signals:
  void NBinCambiati(int);
  void NMobileCambiato(int);
  void NCorrelaCambiato(int);
  void NVisMinCambiato(int);
  void NVisMaxCambiato(int);
  void NElMinCambiato(int);
  void NElMaxCambiato(int);
  void NElMinYCambiato(int);
  void NElMaxYCambiato(int);
  void LogxCambiato(bool);
  void LogyCambiato(bool);
  void MediaCambiata(int);
  void ScartoCambiato(int);
  void InterVisCambiato(int);
  void InterElCambiato(int);
  void InterVisCambiato(int,int);
  void InterElCambiato(int,int);
  void InterYCambiato(int,int);
  void CoordXCambiata(int);
  void CoordYCambiata(int);
  void InterCoordXCambiato(int,int);
  void InterCoordYCambiato(int,int);
  void FineSimulazione();
  void SegnaleGrafico(bool);
  void Stato(const QString &);
  void StatoSequenza(const QString &);
  void TestoCambiato(const QString &);
  void ConfCambiato(const QString &);
  void SalvaCambiato(const QString &);
  void TitoloCambiato(const QString &);
  void EtichettaXCambiato(const QString &);
  void EtichettaYCambiato(const QString &);
  
 protected:
  virtual void paintEvent(QPaintEvent *);
  virtual void mouseMoveEvent(QMouseEvent *);
  
 private:
  void DisegnaPunti(QPainter *);
  void GrStampante(QPainter *);
  void GrMomenti(QPainter *);
  void GrBlu(QPainter *);
  void GrSet(QPainter *p,int s);
  void GrPunti(QPainter *);
  void GrGriglia(QPainter *);
  void GrBarre(QPainter *);
  void GrRiscala();
  void GrLegenda(QPainter *);
  void GrMessage(const char *s, ... );
  void DisegnaTesto(QPainter *);
  void GrScript(char *File2Read,QPainter *p);
  void GrConf(char *File2Read);
  void StringToUnicode(char *cLine,QString *Label);
  FILE *SEGNALE;
  QRect PuntiRett(QPoint) const;
  QColor *GrLinee;
  QString *GrLabel;
  //  vector <QPoint> Punti;
  //double *Punti;
  //double *SetPunti;
  //double *PuntiBlu;
  //double *Ascissa;
  QPoint Punto;
  //  QPointArray *PuntiLinea;
  QPoint Punto1;
  QPoint Punto2;
  QPoint Punto3;
  int PointSize;
  int ContaTempo;
  QTimer *TempoGrafica;
  char *stringa;
  char *cPos;
  char **LineaQuale;
  char FormatPrecX;
  char FormatPrecY;
  char XFormula[60];
  char YFormula[60];
  int NBin;
  int NMobile;
  int NCorrela;
  int NVisMin;
  int NVisMax;
  int NElMin;
  int NElMax;
  int NElMinY;
  int NElMaxY;
  int PuntiMax;
  int RisY;
  int RisX;
  int GrigliaX;
  int GrigliaY;
  int StampaX;
  int StampaY;
  int DigPrecX;
  int DigPrecY;
  double yMax;
  double yMin;
  double yMax1;
  double yMin1;
  double xMax;
  double xMin;
  double yBluMin;
  double yBluMax;
  double Delta;
  double ViewportX;
  double ViewportY;
  double ScalaTopoX;
  double ScalaTopoY;
  double PosLegenda[4];
  double PosInterp[2];
  double RatioWidthHeight;
  int NyMin;
  int NyMax;
  int NxMin;
  int NxMax;
  int NVar;
  int CoordX;
  int CoordY;
  int CoordDY;
  int Tempx;
  int Tempy;
  int NSequenza;
  int Coord;
  int Confini[2];
  int *LineaCome;
  int PrimaVolta;
  int VoltaBlu;
  int VoltaMomenti;
  int VoltaRetta;
  int VoltaParabola;
  int Linee;
  int Quadrati;
  int LogLog;
  int IfLogx;
  int IfLogy;
  int Griglia;
  int IfSegnale;
  int IfStampa;
  int IfRiscala;
  int IfRiscalaTutto;
  int IfDisegna;
  int IfNorm;
  int FontSize;
  QErrorMessage *ErrPrima;
  QMessageBox *Messaggio;
  QImage Immagine;
  QPrinter *Stampante;
  QFont Font;
  MOMENTI m1;
  RETTA r1;
  PARABOLA p1;
  QColor *Colors;
};
/// Defines a window 
class Finestra: public QWidget{
 private:
  //         Elementi Grafici
  ElementiGrafici *e1;
 public:
  Finestra(QWidget *parent=0,const char *name=0);
  void DataFile(char **argv,int *FileList,int NFIle);
  void ConfFile(char *FileName);
  //  char *FileDaScrivere;
  //  void ScriviNomeFile(char *);
};
/// Homemade slider
class ImpostaNBin : public Q3VBox{
  Q_OBJECT
    public:
  ImpostaNBin(QWidget *parent,const char *name =0);
  void Testo(const char *);
  int Numero() const;

  public slots:
  void ImpNumero( int );
  void ImpInter(int Max);
  void ImpInter(int Min,int Max);
 signals:
  void ValoreCambiato( int );

 private:
  QSlider *slider;
  QLabel *label;
  QLabel *title;
};

/* class Sets : public QMenuBar{ */
/*   Q_OBJECT */
/*     public: */
/*   Sets(QWidget *parent,const char *name =0); */
/*   QAction *ASet; */
/*   int NSet; */
/*   public slots: */
/*   void ImpCoord(int Coord); */
/*  signals: */
/*   void ChoosedSet(void); */
/*  }; */


#endif //ELEMENTIGRAFICI_H
