/***********************************************************************
ElementiGrafici: Class that interface the VarData, VarElementiGrafici 
and Matematica classes function for comunicating with main window of 
Visualizza. The main function is the plotting function DisegnaPunti() 
that provides the visualisation of two set of data. The bars at the 
bottom permits to change the visualisation and the limit for the analisys
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
#include "ElementiGrafici.h"
;
#include <qevent.h>
#include <qpainter.h>
#include <qlabel.h>

ElementiGrafici::ElementiGrafici(QWidget *parent,const char *name)
:QWidget(parent,name){
  PuntiMax = 1000;
  NVisMin=0;
  NVisMax=1000;
  NBin = 100;
  NMobile = 1;
  NCorrela = 2;
  Coord = 2;
  yMax = 1.;
  CoordX = 0;
  CoordY = 0;
  CoordDY = 0;
  NVar = 1;
  FontSize = 42;
  RatioWidthHeight = 1.3;
  PointSize = 10;
  IfNorm = 0;
  GrigliaX = 6;
  GrigliaY = 6;
  //Punti = NULL;
  //PuntiBlu = NULL;
  //Ascissa = NULL;
  stringa = new char[256];
  //  PuntiLinea = new QPointArray[3];
  cPos = new char[30];
  nomeFile = new char[160];
  sprintf(nomeFile,"Density.dat");
  nomeConf = new char[160];
  sprintf(nomeConf,"Vis.conf");
  nomeSalva = new char[160];
  sprintf(nomeSalva,"Grafico.ps");
  // nomeTit = new QString();
  // nomeEtX = new QString();
  // nomeEtY = new QString();
  NomeEntrata(QString(nomeFile));
  setPalette( QPalette( QColor( 250, 250, 200) ) );
  //  TempoGrafica = new QTimer(this, "tempo della grafica");
  //  connect(TempoGrafica,SIGNAL(timeout() ),this,SLOT());
  ErrPrima = new QErrorMessage(this);
  Messaggio = new QMessageBox(this,"Messaggio");
  Stampante = new QPrinter(QPrinter::PrinterResolution);
  Stampante->setFullPage(TRUE);
  Stampante->setOutputToFile(TRUE);
  Stampante->setOutputFileName(QString(nomeSalva));
  Stampante->setOrientation(QPrinter::Landscape);
  Stampante->setPageSize(QPrinter::A4);
  Font.setFamily("Helvetica");
  Font.setWeight(QFont::Bold);
  IfDisegna = 0;
  IfRiscala = RIS_TUTTI;
  PrimaVolta = 1;
  Linee = 1;
  Quadrati = 0;
  Griglia = 1;
  NSequenza = 0;
  IfStampa = 0;
  NDis = 0;
  m1.Min = 0.;
  DigPrecX = 3;
  DigPrecY = 3;
  FormatPrecX = 'g';
  FormatPrecY = 'g';
  NElMinY=0;
  NElMaxY=0;
  ViewportX = 0.;//0.05;
  ViewportY = 0.;//0.05;
  GrLinee = new QColor[INIT_NUM];
  GrLabel = new QString[INIT_NUM];
  Colors =  new QColor[INIT_NUM];
  LineaCome = new int[INIT_NUM];
  LineaQuale = new char*[INIT_NUM];
  for(int c=0;c<10;c++) LineaQuale[c] = new char[INIT_NUM];
  VoltaMomenti = 0;
  QColor Colors1[INIT_NUM] = {QColor(157,128,43,255),
			      QColor(180,40,43,255),
			      QColor(20,180,70,255),
			      QColor(60,100,180,255),
			      QColor(180,50,180,255),
			      QColor(200,200,50,255),
			      QColor(50,180,180,255),
			      QColor(120,70,120,255),
			      QColor(180,10,240,255),
			      QColor(125,180,113,255),
			      QColor(160,110,110,255),
			      QColor(230,60,220,255),
			      QColor(132,20,220,255),
			      QColor(180,110,170,255),
			      QColor(225,0,143,255)};
  //QColor(0,140,210,255)};
  for(int c=0;c<10;c++){
    Colors[c] = Colors1[c];
  }
  PosLegenda[0] = 2.0;
  PosLegenda[1] = 2.0;
  PosLegenda[2] = 3.0;
  PosLegenda[3] = 3.0;
  PosInterp[0] = .5;
  PosInterp[1] = .2;
}
ElementiGrafici::~ElementiGrafici(){
  delete [] stringa;
  delete [] cPos;
  delete [] nomeFile;
  delete [] nomeSalva;
  delete [] nomeTit;
  delete [] nomeConf;
  delete [] nomeEtX;
  delete [] nomeEtY;
  delete [] GrLinee;
  delete [] GrLabel;
  delete [] Colors;
  delete [] LineaCome;
  delete [] LineaQuale;
  delete ErrPrima;
  delete Messaggio;
  delete Stampante;
}
void ElementiGrafici::ChooseDataFile(char *FileName){
  NomeEntrata(QString(FileName));
}
void ElementiGrafici::ChooseConfFile(char *FileName){
  NomeConf(QString(FileName));
}
void ElementiGrafici::Apri(){
  int NFile = 1;
  char **argv = (char **)calloc(NFile,sizeof(char));
  int FileList = 0;
  for(int f=0;f<NFile;f++){
    argv[f] = (char *)calloc(120,sizeof(char));
  }
  sprintf(argv[0],"%s",nomeFile);
  Apri(argv,&FileList,NFile);
  for(int f=0;f<NFile;f++){
    free(argv[f]);
  }
  free(argv);
}
void ElementiGrafici::Apri(char **argv,int *FileList,int NFile){
  if(!PrimaVolta){
    delete v1;
    //    delete Punti;
  }
  NomeFile(QString(argv[FileList[0]]));
  v1 = new VarDatFile(argv,FileList,NFile,NBin);
  //v1 = new Variabili(argv[FileList[0]],NBin);
  //for(int f=1;f<NFile;f++)
  //v1->Aggiungi(argv[FileList[f]]);
  NMax = v1->pNMax();
  PuntiMax = NMax;
  LogLog = DIS_NOLOG;
  NVar = v1->pNVar();
  GrLinee = (QColor *)realloc(GrLinee,NVar*sizeof(QColor));
  delete [] GrLabel;
  GrLabel = new QString[NVar];
  //GrLabel = (QString *)realloc(GrLabel,NVar*sizeof(QString));
  LineaCome = (int *)realloc(LineaCome,NVar*sizeof(int));
  //LineaQuale = (char **)realloc(LineaQuale,NVar*sizeof(char));
  //for(int c=0;c<NVar;c++)LineaQuale[c] = (char *)realloc(LineaQuale[c],60*sizeof(char));
  for(int c=0,cc=0,init=0;c<NVar;c++,cc++){
    if(!(c%10)) cc = init++;
    if(init > 10) init = 0;
    GrLinee[c] = Colors[cc];
    LineaCome[c] = LINEA_TRATTO;
    //sprintf(LineaQuale[c],"");
  }
  ImpNVarX(REF_SEQ,NVar);
  ImpNVarY(0,NVar);
  CoordX = 0;
  CoordY = 0;
  ImpInter(NVisMax=NMax);
  ImpInterEl(NElMax=NMax);
  ImpNVisMax(NVisMax=NMax);
  ImpNVisMin(NVisMin=0);
  ImpNElMax(NElMax=NMax);
  ImpNElMin(NElMin=0);
  if(NVar > 1){
    CoordY = 1;
  }
  else {
    CoordY = 0;
    CoordX = REF_SEQ;
  }
  sprintf(stringa,"Caricati %d vettori da %d Punti da %s",NVar,NMax,nomeFile);
  sprintf(stringa,"%s",v1->PrintHeader());
  ImpCoordX(CoordX);
  ImpCoordY(CoordY);
  PrimaVolta = 0;
  IfDisegna = 0;
  DIS_ADD_TYPE(IfDisegna,DIS_TUTTI);
  IfRiscala = 0;
  DIS_ADD_TYPE(IfRiscala,RIS_TUTTI);
  repaint();
}
void ElementiGrafici::ApriExt(double **ExtSt,int ExtNMax,int ExtNVar,int ExtNBin){
  if(!PrimaVolta){
    delete v1;
    //    delete Punti;
  }
  NBin = ExtNBin;
  NMax = ExtNMax;
  NVar = ExtNVar;
  NomeFile(QString(nomeFile));
  v1 = new VarDatFile(ExtSt,NMax,NVar,NBin);
  LogLog = DIS_NOLOG;
  NVar = v1->pNVar();
  GrLinee = (QColor *)realloc(GrLinee,NVar*sizeof(QColor));
  delete [] GrLabel;
  GrLabel = new QString[NVar];
  //  GrLabel = (QString *)realloc(GrLinee,NVar*sizeof(QString));
  LineaCome = (int *)realloc(LineaCome,NVar*sizeof(int));
  //LineaQuale = (char **)realloc(LineaQuale,NVar*sizeof(char));
  //  for(int c=0;c<NVar;c++)LineaQuale[c] = (char *)realloc(LineaQuale[c],60*sizeof(char));
  for(int c=0;c<NVar;c++){
    if(c<10)
      GrLinee[c] = Colors[c];
    else 
      GrLinee[c] = Colors[0];
    LineaCome[c] = LINEA_TRATTO;
    //sprintf(LineaQuale[c],"");
  }
  ImpNVarX(-1,NVar);
  ImpNVarY(0,NVar);
  ImpInter(NVisMax=NMax);
  ImpInterEl(NElMax=NMax);
  ImpNVisMax(NVisMax=NMax);
  ImpNVisMin(NVisMin=0);
  ImpNElMax(NElMax=NMax);
  ImpNElMin(NElMin=0);
  sprintf(stringa,"Caricati %d vettori da %d Punti da %s",NVar,NMax,nomeFile);
  sprintf(stringa,"%s",v1->PrintHeader());
  ImpCoordX(CoordX);
  ImpCoordY(CoordY);
  PrimaVolta = 0;
  //for(int i=0;i<NMax;i++) printf("%lf \n",ExtSt[0][i]);
  printf("%d %d %d\n",NMax,NVar,NBin);
  repaint();
}
void ElementiGrafici::ApriExt(VarDatFile *v1Ext){
  if(!PrimaVolta){
    delete v1;
    //    delete Punti;
  }
  v1 = new VarDatFile(0,0,0);
  v1 = v1Ext;
  NBin = v1->NBin;
  NMax = v1->pNMax();
  NVar = v1->pNVar();
  NomeFile(QString(nomeFile));
  LogLog = DIS_NOLOG;
  GrLinee = (QColor *)realloc(GrLinee,NVar*sizeof(QColor));
  delete [] GrLabel;
  GrLabel = new QString[NVar];
  //  GrLabel = (QString *)realloc(GrLabel,NVar*sizeof(QString));
  LineaCome = (int *)realloc(LineaCome,NVar*sizeof(int));
  //LineaQuale = (char **)realloc(LineaQuale,NVar*sizeof(char));
  //for(int c=0;c<NVar;c++)LineaQuale[c] = (char *)realloc(LineaQuale[c],60*sizeof(char));
  for(int c=0;c<NVar;c++){
    if(c<10)
      GrLinee[c] = Colors[c];
    else 
      GrLinee[c] = Colors[0];
    LineaCome[c] = LINEA_TRATTO;
    //    sprintf(LineaQuale[c],"");
  }
  ImpNVarX(-1,NVar);
  ImpNVarY(0,NVar);
  CoordX = 0;
  CoordY = 0;
  ImpInter(NVisMax=NMax);
  ImpInterEl(NElMax=NMax);
  ImpNVisMax(NVisMax=NMax);
  ImpNVisMin(NVisMin=0);
  ImpNElMax(NElMax=NMax);
  ImpNElMin(NElMin=0);
  if(NMax > 10000) Griglia = 0;
  sprintf(stringa,"Caricati %d vettori da %d Punti da %s",NVar,NMax,nomeFile);
  sprintf(stringa,"%s",v1->PrintHeader());
  IfSegnale = 1;
  ImpCoordX(CoordX);
  ImpCoordY(CoordY);
  IfRiscala = RIS_TUTTI;
  PrimaVolta = 0;
  //for(int i=0;i<NMax;i++) printf("%lf \n",ExtSt[0][i]);
  printf("%d %d %d\n",NMax,NVar,NBin);
  repaint();
}
void ElementiGrafici::Aggiungi(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file .dat prima"));
    return;
  }
  NomeFile(QString(nomeFile));
  if( v1->Aggiungi(nomeFile) == 0) return;
  NMax = v1->pNMax(); 
  PuntiMax = NMax;
  // printf("Caricato un vettore da %d Punti da %s",NMax,nomeFile);
  sprintf(stringa,"Caricato un vettore da %d Punti da %s",NMax,nomeFile);
  sprintf(stringa,"%s",v1->PrintHeader());
  NVar = v1->pNVar();
  GrLinee = (QColor *)realloc(GrLinee,NVar*sizeof(QColor));
  delete [] GrLabel;
  GrLabel = new QString[NVar];
  //  GrLabel = (QString *)realloc(GrLabel,NVar*sizeof(QString));
  LineaCome = (int *)realloc(LineaCome,NVar*sizeof(int));
  //LineaQuale = (char **)realloc(LineaQuale,NVar*sizeof(char));
  //  for(int c=0;c<NVar;c++)LineaQuale[c] = (char *)realloc(LineaQuale[c],60*sizeof(char));
  for(int c=0;c<NVar;c++){
    if(c<10)
      GrLinee[c] = Colors[c];
    else 
      GrLinee[c] = Colors[0];
    LineaCome[c] = LINEA_TRATTO;
  }
  ImpNVarX(-1,NVar);
  ImpNVarY(0,NVar);
  ImpCoordX(CoordX);
  ImpCoordY(CoordY);
  PuntiSegnale();
}
void ElementiGrafici::mouseMoveEvent(QMouseEvent *Topo){
  sprintf(cPos,"%.3g : %.3g",
	  (double)(Topo->x())/(double)width()*ScalaTopoX+xMin,
	  (height()-Topo->y())/(double)height()*ScalaTopoY+yMin);
  QString Info = QString(cPos);
  emit StatoSequenza(QString(cPos));  
}
void ElementiGrafici::DisegnaTesto(QPainter *p){
  p->setBrush( Qt::blue );
  p->setPen( Qt::blue );
  p->drawText(width()/2,height()/2,QString("L'autocorrelazione di un rumore non bianco\n `e lorentziano e decade come \nr^-2"));
  p->drawImage(QPoint(width()/2,height()/2),QImage(QString("Prova.png"),0),QRect(0,0,300,200),0);
  p->drawRect(QRect(width()/2-100,height()/2-50,200,100));
  p->drawText(0,height()/2,QString("L'autocorrelazione di un rumore non bianco\n `e lorentziano e decade come \nr^-2"));
}
QRect ElementiGrafici::PuntiRett(QPoint P) const{
  QRect r(0,0,2*PointSize,2*PointSize);
  r.moveCenter(P);
  return r;  
}
ImpostaNBin::ImpostaNBin(QWidget *parent,const char *name)
  :Q3VBox(parent,name){
  title = new QLabel(name,this,"title");
  title->setAlignment( Qt::AlignCenter );
  Qt::Orientation Orientazione = Qt::Horizontal;
  //  if(Modo == 1) Orientazione = Qt::Vertical;
  slider = new QSlider(Orientazione, this ,"slider");
  slider->setRange(0,100);
  slider->setValue(0);

    slider->setTickPosition(QSlider::TicksBelow);
  //  else slider->setTickPosition(QSlider::TicksRight);
    

  label = new QLabel(" ",this,"label");
  label->setAlignment( Qt::AlignCenter );
  connect(slider,SIGNAL(valueChanged(int)),label,SLOT(setNum(int)) );
  connect(slider,SIGNAL(valueChanged(int)),SIGNAL(ValoreCambiato(int)) );
}
void ImpostaNBin::ImpNumero(int n){
  if(n >= 0){
    slider->setValue(n);
  }
}
void ImpostaNBin::ImpInter(int Max){
  if(Max > 0){
    slider->setRange(0,Max);
  }
}
void ImpostaNBin::ImpInter(int Min,int Max){
  if(Min < Max){
    slider->setRange(Min,Max);
  }
}
int ImpostaNBin::Numero() const{
  return slider->value();
}
// Sets::Sets(QWidget *parent,const char *name )
//   :QMenuBar(parent,name){
//   ASet = NULL;
//   ASet = (QAction *)realloc(ASet,NSet*sizeof(QAction)); 
//   char *stringa;
//   stringa = (char *)malloc(100*sizeof(char));
//   for(int v=0;v<NSet;v++){
//     sprintf(stringa,"Set%d",v);
//     //ASet[v] = addAction(QString(stringa));
//   }
//   //connect(slider,SIGNAL(valueChanged(int)),label,SLOT(setNum(int)) );
//   //connect(slider,SIGNAL(valueChanged(int)),SIGNAL(ValoreCambiato(int)) );
// }
// void Sets::ImpCoord(int Coord){
  
// }
