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
;// α-a β-b γ-c δ-d ε-e ζ-f η-g θ-h ι-i κ-l λ-k μ-l ν-m ξ-n ο-o π-p ρ-q σ-r ς-r τ-s υ-t φ-u ψ-v ω-x χ-y  
#include "ElementiGrafici.h"
#include <qpainter.h>
#include <Q3PaintDeviceMetrics>

void ElementiGrafici::GrMessage(const char * s, ...){
#ifdef DEBUG
  va_list args;
  va_start(args, s);
  fprintf(stderr,"%d) ",NDis);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else  
  return;
#endif
}
void ElementiGrafici::paintEvent(QPaintEvent *e){
  GrMessage("Visualizza.paintEvent(begin)");
  //QRect DaNuovo = e->rect();
  //  if(DaNuovo.intersects( PuntiRett(0)) )
  QPainter p(this);
  DisegnaPunti(&p);
  if(PrimaVolta){
    Immagine.load("Logo.png");
    p.drawPixmap(200,300,QPixmap(Immagine),0,0,-1,-1);
  }
  GrMessage("Visualizza.paintEvent(end)");
}
void ElementiGrafici::DisegnaPunti(QPainter *p){
  if(PrimaVolta){
    return;
  }
  if(NVisMin < 0 || NVisMax > PuntiMax || NVisMin > NVisMax){
    sprintf(stringa,"Non `e corretto l'ordine 0<=%d<%d<=%d",NVisMin,NVisMax,PuntiMax);
    printf("Non `e corretto l'ordine 0<=%d<%d<=%d\n",NVisMin,NVisMax,PuntiMax);
    ErrPrima->message(stringa);
    return;
  }
  if( IfRiscala != 0 ) GrRiscala();
  GrStampante(p);
  GrConf(nomeConf);
  GrScript(nomeConf,p);
  if( DIS_IF_TYPE(IfDisegna,DIS_TUTTI) ){
    for(int s=0,sColor=0;s<NVar;s++){
      if(v1->IsAbscissa(s)) continue;
      sColor = s;
      p->setBrush( GrLinee[sColor] );
      p->setPen( QPen( GrLinee[sColor],2 ) );
      Quadrati = DIS_IF_TYPE(LineaCome[sColor],LINEA_PUNTO);
      Linee = DIS_IF_TYPE(LineaCome[sColor],LINEA_TRATTO);
      GrSet(p,s);
      sColor++;
    }
  }
  if( DIS_IF_TYPE(IfDisegna,DIS_SEGNALE) ){
    p->setBrush( GrLinee[CoordY] );
    p->setPen( QPen( GrLinee[CoordY],2 ) );
    GrSet(p,CoordY);
  }
  if( DIS_IF_TYPE(IfDisegna,DIS_MOMENTI) )
    GrMomenti(p);
  if( DIS_IF_TYPE(IfDisegna,DIS_PUNTI) ){
    p->setBrush( GrLinee[CoordY] );
    p->setPen( QPen( GrLinee[CoordY],2 ) );
    GrPunti(p);
  }
  if( DIS_IF_TYPE(IfDisegna,DIS_BLU) ) GrBlu(p);
  if( Griglia ) GrGriglia(p);
  if(!IfStampa) GrBarre(p);
  GrLegenda(p);
  GrMessage("Disegna %d = Tutti %d Set %d Momenti %d Punti %d Blu %d Riscala %d = RisUno %d RisTutti %d",IfDisegna,
	    DIS_IF_TYPE(IfDisegna,DIS_TUTTI),
	    DIS_IF_TYPE(IfDisegna,DIS_SEGNALE),
	    DIS_IF_TYPE(IfDisegna,DIS_MOMENTI),
	    DIS_IF_TYPE(IfDisegna,DIS_PUNTI),
	    DIS_IF_TYPE(IfDisegna,DIS_BLU),
	    IfRiscala,
	    DIS_IF_TYPE(IfRiscala,RIS_UNO),
	    DIS_IF_TYPE(IfRiscala,RIS_TUTTI));
  NDis++;
}
void ElementiGrafici::GrSet(QPainter *p,int s){
  GrMessage("Visualizza.Dis.GrSet");
  if(v1->IsAbscissa(s)) return;
  //printf("%d %d\n",s,v1->pRefAbsc(s));
  int IfDrawLine = 0;
  for(int i=MIN(NVisMin,v1->pNRow(s));i<MIN(NVisMax,v1->pNRow(s));i++){
    //printf("%lf %lf\n",v1->Ascissa(i),v1->Val(i));
    if(!DIS_IF_TYPE(LogLog,DIS_LOGX))
      Tempx = (int) ( (double)((v1->Abscissa(s,i)-xMin)*RisX)/(double)(xMax-xMin) ); 
    else {
      if(v1->Abscissa(s,i)<=0.) continue;
      Tempx = (int) ( (log10(v1->Abscissa(s,i))-xMin)/(xMax-xMin)*RisX);
    }
    if(!DIS_IF_TYPE(LogLog,DIS_LOGY))
      Tempy = RisY-(int)( (v1->Val(s,i)-yMin)*RisY/ (yMax-yMin) );
    else {
      if(v1->Val(s,i) <= 0.) continue;
      Tempy = RisY-(int)( (log10(v1->Val(s,i))-yMin)*RisY/ (yMax-yMin) );
      //printf("%d %d %lf\n",Tempy,RisY,log10(v1->Val(s,i)));
    }
    // if(Tempx > RisX){IfDrawLine=0;continue;}
    // if(Tempy > RisY){IfDrawLine=0;continue;}
    // if(Tempx < 0 ){IfDrawLine=0;continue;}
    // if(Tempy < 0 ){IfDrawLine=0;continue;}
    if(Tempx > RisX){Tempx = RisX;}
    if(Tempy > RisY){Tempy = RisY;}
    if(Tempx < 0 ){Tempx = 0;}
    if(Tempy < 0 ){Tempy = 0;}
    Punto = QPoint(Tempx, Tempy);
    //if(Quadrati) p->drawRect( PuntiRett(Punto) );
    if(Quadrati){
      Punto2 = QPoint(Punto.x()-PointSize,Punto.y()-PointSize);
      Punto3 = QPoint(Punto.x()+PointSize,Punto.y()+PointSize);
      p->drawLine(Punto2,Punto3);
      Punto2 = QPoint(Punto.x()+PointSize,Punto.y()-PointSize);
      Punto3 = QPoint(Punto.x()-PointSize,Punto.y()+PointSize);
      p->drawLine(Punto2,Punto3);
    }
    // if(Quadrati){
    //   Punto2 = QPoint(Punto.x()+PointSize,Punto.y()+PointSize);
    //   Punto3 = QPoint(Punto.x()+PointSize,Punto.y()-PointSize);
    //   p->drawLine(Punto2,Punto3);
    //   Punto2 = QPoint(Punto.x()+PointSize,Punto.y()-PointSize);
    //   Punto3 = QPoint(Punto.x()-PointSize,Punto.y()-PointSize);
    //   p->drawLine(Punto2,Punto3);
    //   Punto2 = QPoint(Punto.x()-PointSize,Punto.y()-PointSize);
    //   Punto3 = QPoint(Punto.x()-PointSize,Punto.y()+PointSize);
    //   p->drawLine(Punto2,Punto3);
    //   Punto2 = QPoint(Punto.x()-PointSize,Punto.y()+PointSize);
    //   Punto3 = QPoint(Punto.x()+PointSize,Punto.y()+PointSize);
    //   p->drawLine(Punto2,Punto3);
    // }
    if(Linee && IfDrawLine) p->drawLine(Punto,Punto1);
    IfDrawLine = 1;
    Punto1 = Punto;
  }
}
void ElementiGrafici::GrPunti(QPainter *p){
  GrMessage("Visualizza.Dis.GrPunti");
  for(int i=NVisMin;i<NVisMax*NMobile;i+=NMobile){// FIXME
    if(xMin > v1->Abscissa(CoordY,i)) xMin = v1->Abscissa(CoordY,i);
    if(xMax < v1->Abscissa(CoordY,i)) xMax = v1->Abscissa(CoordY,i);
  }
  ScalaTopoX = xMax-xMin;
  //FILE *Apri = fopen("Scritto.dat","w");
  for(int i=NVisMin;i<NVisMax;i++){
    //printf("%lf %lf\n",v1->Ascissa(i*NMobile),v1->pPunti(i));
    if(!DIS_IF_TYPE(LogLog,DIS_LOGX))
      Tempx = (int) ( (double)((v1->Abscissa(CoordY,i*NMobile)-xMin)*RisX)/(double)(xMax-xMin) ); 
    else{
      if(v1->Abscissa(CoordY,i)<=0.) continue;
      Tempx = (int) ( (log10(v1->Abscissa(CoordY,i))-xMin)/(xMax-xMin)*RisX);
    }
    if(!DIS_IF_TYPE(LogLog,DIS_LOGY))
      Tempy = RisY-(int)( (v1->pPunti(i)-yMin)*RisY/ (yMax-yMin) );
    else {
      if(v1->pPunti(i) <= 0.) continue;
      Tempy = RisY-(int)( (log10(v1->pPunti(i))-yMin)*RisY/ (yMax-yMin) );
      //printf("%d %d %lf\n",Tempy,RisY,log10(v1->Val(s,i)));
    }
    if(Tempx > RisX) Tempx=RisX;
    if(Tempy > RisY) Tempy=RisY;
    if(Tempx < 0 ) Tempx=0;
    if(Tempy < 0 ) Tempy=0;
    Punto = QPoint(Tempx, Tempy);
    if(Quadrati)    p->drawRect( PuntiRett(Punto) );
    if( Linee && i > NVisMin) p->drawLine(Punto,Punto1);
    Punto1 = Punto;
    int ErrY = (int)( .5*(v1->pPuntiErr(i))*RisY/ (yMax-yMin) );
    QPoint PuntoErr1(Tempx,Tempy+ErrY);
    QPoint PuntoErr2(Tempx,Tempy-ErrY);
    //p->setPen( QPen(Qt::SolidPattern,1,Qt::SolidLine,Qt::SquareCap,Qt::BevelJoin) );
    //fprintf(Apri,"%lf %lf %lf\n",v1->Ascissa(i*NMobile)*.2,v1->pPunti(i),v1->pPuntiErr(i));
    p->drawLine(PuntoErr1,PuntoErr2);
  }
  //fclose(Apri);
}
void ElementiGrafici::GrMomenti(QPainter *p){
  GrMessage("Visualizza.Dis.GrMomenti");
  p->setBrush( Qt::red );
  p->setPen( Qt::red );
  for(int i=NVisMin;i<NVisMax;i++){
    Tempx = (int) ( (double)((i-NVisMin)*RisX)/(double)(NVisMax-NVisMin) );
    int Tempx2 = (int) ( (double)((i+1-NVisMin)*RisX)/(double)(NVisMax-NVisMin) );
    Tempy = RisY-(int)( (v1->pInter(i)-yMin)*RisY/ (yMax-yMin) );
    Punto = QPoint(Tempx,Tempy);
    QPoint Punto2 = QPoint(Tempx2,Tempy);
    if( Quadrati ) p->drawRect( PuntiRett(Punto) );
    if( Linee && i>NVisMin){
      p->drawLine(Punto2,Punto);
      p->drawLine(Punto,Punto1);
    }
    Punto1 = Punto2;
  }
  GrMessage("Visualizza.Dis.GrMomenti.Blue");
  p->setBrush( Qt::blue );
  p->setPen( Qt::blue );
  yBluMin = v1->pInter1(NVisMin);
  yBluMax = v1->pInter1(NVisMin);
  for(int i=NVisMin;i<NVisMax;i++){
    Tempx = (int) ( (double)((i-NVisMin)*RisX)/(double)(NVisMax-NVisMin) );
    Tempy = RisY-(int)( (v1->pInter1(i)-yMin)*RisY/ (yMax-yMin) );
    Punto = QPoint(Tempx,Tempy);
    if( Quadrati ) p->drawRect( PuntiRett(Punto) );
    if( Linee && i>NVisMin) p->drawLine(Punto,Punto1);
    Punto1 = Punto;
    if( yBluMin > v1->pInter1(i) )
      yBluMin = v1->pInter1(i);
    if( yBluMax < v1->pInter1(i) )
      yBluMax = v1->pInter1(i);
  }
  GrMessage("Visualizza.Dis.GrMomenti.String");
  sprintf(stringa,"%.2g",yBluMax);
  p->drawText(RisX-70,15,QString(stringa));
  sprintf(stringa,"%.2g",yBluMax/2.);
  p->drawText(RisX-70,(int)RisY/2,QString(stringa));
  sprintf(stringa,"%.2g",yBluMin);
  p->drawText(RisX-70,RisY-15,QString(stringa));
  int PosMomX=(int)(PosInterp[0]*RisX);
  int PosMomY=(int)((1.-PosInterp[1])*RisY);
  sprintf(stringa,"m %.4g s %g N %d",m1.Uno,m1.Due,m1.Num);
  p->drawText(PosMomX,PosMomY,QString(stringa));
  if(!IfStampa){
    sprintf(stringa,"t %g Chi %.2g",m1.Tre,m1.Chi);
    p->drawText(RisX/2-70,RisY/2+120,QString(stringa));
  }
}
void ElementiGrafici::GrBlu(QPainter *p){
  GrMessage("Visualizza.Dis.GrBlu");
  int IfDrawLine = 0;
  p->setBrush( Qt::blue );
  p->setPen( Qt::blue );
  if(!IfStampa)
    p->setPen(QPen(Qt::blue,1,Qt::DashLine,Qt::RoundCap,Qt::RoundJoin));
  else 
    p->setPen(QPen(Qt::blue,3,Qt::SolidLine,Qt::RoundCap,Qt::RoundJoin));
  for(int i=NElMin;i<NElMax;i++){
    if(!DIS_IF_TYPE(LogLog,DIS_LOGX))
      Tempx = (int) ( (double)((v1->Abscissa(CoordY,i)-xMin)*RisX)/(double)(xMax-xMin) ) ; 
    else
      //Tempx = (int) ( log10((v1->Abscissa(CoordY,i))/(xMin))/log10((xMax)/(xMin))*RisX );
      Tempx = (int) ( (log10(v1->Abscissa(CoordY,i))-xMin)/(xMax-xMin)*RisX );
    if( DIS_IF_TYPE(IfDisegna,DIS_RETTA) ){
      if(!DIS_IF_TYPE(LogLog,DIS_LOGY)){
	Tempy = RisY-(int)( (v1->Abscissa(CoordY,i)*r1.m+r1.q-yMin)*RisY/ (yMax-yMin) );
      }
      else
	//Tempy = RisY-(int)( (log10(v1->Abscissa(CoordY,i))*r1.m+r1.q-log10(yMin))*RisY/ (log10(yMax/yMin)) );	
	Tempy = RisY-(int)( (log10(v1->Abscissa(CoordY,i))*r1.m+r1.q-yMin)*RisY/(yMax-yMin) );	
    }
    else if( DIS_IF_TYPE(IfDisegna,DIS_PARABOLA) ){
      if(!DIS_IF_TYPE(LogLog,DIS_LOGY)){
	Tempy = RisY-(int)( (v1->Abscissa(CoordY,i)*p1.a1+v1->Abscissa(CoordY,i)*v1->Abscissa(CoordY,i)*p1.a2+p1.a0-yMin)*RisY/ (yMax-yMin) );
      }
      else{
	double y = log10(v1->Abscissa(CoordY,i));
	// Tempy = RisY-(int)( (log10(v1->Abscissa(CoordY,i))*log10(v1->Abscissa(CoordY,i))*p1.a2+log10(v1->Abscissa(CoordY,i))*p1.a1+p1.a0-log10(yMin))*RisY/ (log10(yMax/yMin)) );	
	Tempy = RisY-(int)((y*y*p1.a2+y*p1.a1+p1.a0-yMin)*RisY/(yMax-yMin));
      }
    }
    else if( DIS_IF_TYPE(IfDisegna,DIS_EXP) ){
      if(!DIS_IF_TYPE(LogLog,DIS_LOGY)){
	Tempy = RisY-(int)((exp(r1.q+r1.m*v1->Abscissa(CoordY,i))-yMin)*RisY/(yMax-yMin));
      }
      else 
	Tempy = RisY-(int)((exp(r1.q+r1.m*log10(v1->Abscissa(CoordY,i)))-yMin)*RisY/(yMax-yMin));
    }
    else if( DIS_IF_TYPE(IfDisegna,DIS_GAUSS) ){
      Tempy = RisY-(int)( (v1->Gauss(m1.Uno,m1.Due,v1->Abscissa(CoordY,i))-yMin)*RisY/(yMax-yMin));
    }
    else
      Tempy = RisY-(int)( (v1->pInter1(i)-yMin)*RisY/ (yMax-yMin) );
    //printf("%d %d %d %d %lf %lf %lf %lf\n",Tempx,Tempy,RisX,RisY,yMin,yMax,v1->Abscissa(CoordY,i),v1->Val(CoordY,i));
    Punto = QPoint(Tempx,Tempy);
    if(Tempx > RisX){IfDrawLine=0;continue;}
    if(Tempy > RisY){IfDrawLine=0;continue;}
    if(Tempx < 0 ){IfDrawLine=0;continue;}
    if(Tempy < 0 ){IfDrawLine=0;continue;}
    //if( Quadrati ) p->drawRect( PuntiRett(Punto) );
    if( /*Linee && */ IfDrawLine) p->drawLine(Punto,Punto1);
    Punto1 = Punto;
    IfDrawLine = 1;
  }
  //----------------------LABEL-------------------------
  int PosLegX=(int)(PosInterp[0]*RisX);
  int PosLegY=(int)((1.-PosInterp[1])*RisY);
  if( DIS_IF_TYPE(IfDisegna,DIS_RETTA) ){
    //    NumLabel.setNum(r1.m,'g',DigPrec);
    if(!IfStampa){
      sprintf(stringa,"m %.3g+-%.3g q %.3g+-%.3g",r1.m,r1.ErrM,r1.q,r1.ErrQ);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
      sprintf(stringa,"r %.5g Corr %.3g",r1.Corr,r1.Cov);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY+15,QString(stringa));
    }
    else {
      sprintf(stringa,"y=%.3f x + %.3f",r1.m,r1.ErrM,r1.q,r1.ErrQ);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
      //sprintf(stringa,"r %.5e Corr %.3e",r1.Corr,r1.Cov);
      //p->drawText(RisX/2-RisX/8,RisY-15,QString(stringa));
    }      
  }
  else if( DIS_IF_TYPE(IfDisegna,DIS_PARABOLA) ){
    int RefY = (int)(.1*RisY);
    if(!IfStampa){
      sprintf(stringa,"%.2g + %.2g x + %.2g x^2",p1.a0,p1.a1,p1.a2);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
      sprintf(stringa,"Minimum (%.3g %.3g) ",p1.Minimo,p1.MinimoY);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY+25,QString(stringa));
    }
    else {
      p->setFont(QFont("Palatino",FontSize-6));
      sprintf(stringa,"%.2g + %.2g x + %.2g x^2",p1.a0,p1.a1,p1.a2);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
      sprintf(stringa,"min %.2g ",p1.Minimo);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY+25,QString(stringa));
      p->setFont(QFont("Palatino",FontSize));
    }
  }
  else if( DIS_IF_TYPE(IfDisegna,DIS_GAUSS) ){
    int RefY = (int)(.1*RisY);
    sprintf(stringa,"x_0 = %.3g s = %.3g ",m1.Uno,m1.Due);
    p->drawText(PosLegX,PosLegY,QString(stringa));
  }
  else if( DIS_IF_TYPE(IfDisegna,DIS_EXP) ){
    //    NumLabel.setNum(r1.m,'g',DigPrec);
    if(!IfStampa){
      sprintf(stringa,"y=(%.3g+-%.3g)exp(x/(%.3g+-%.3g))",log10(r1.q),log10(r1.ErrQ),1./r1.m,r1.ErrM/SQR(r1.m));
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
    }
    else {
      sprintf(stringa,"y=%.3g exp(%.3g x)",log10(r1.q),r1.m);
      p->drawText(PosLegX-strlen(stringa)*3.4,PosLegY,QString(stringa));
      //sprintf(stringa,"r %.5e Corr %.3e",r1.Corr,r1.Cov);
      //p->drawText(RisX/2-RisX/8,RisY-15,QString(stringa));
    }      
  }
  //printf("%d %d %d %d\n",NElMin,NElMax,NxMin,NxMax);
}
//     if( (i%3)==0 ){
//       for(int j=0;j<3;j++){
// 	PuntiLinea[j].setPoint(j,Mas(i,(int) ( i*RisX/(double)(NVisMax-NVisMin)) ),RisY-(int)( (double) Punti[i]*(double)(RisY)/yMax ) + (int) ( RisY/yMin) );
//       }
//     }
void ElementiGrafici::GrBarre(QPainter *p){
  GrMessage("Visualizza.Dis.GrBarre");
  p->setBrush( Qt::black );//              Barre
  p->setPen( QPen(Qt::black,1) );
  Tempx = 0;
  GrMessage("Visualizza.Dis.GrBarre.Ascissa");
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)){
    if(v1->Abscissa(CoordY,NElMin) > 0. ) 
      Tempx = (int)( (log10(v1->Abscissa(CoordY,NElMin))-xMin)/(xMax-xMin)*RisX);
  }
  else
    Tempx = (int) ( (double)((v1->Abscissa(CoordY,NElMin)-xMin)*RisX)/(double)(xMax-xMin) );
  p->drawLine( QPoint(Tempx ,0), QPoint(Tempx ,RisY));
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)){
    if(v1->Abscissa(CoordY,NElMax) > 0. ) 
      Tempx = (int)( (log10(v1->Abscissa(CoordY,NElMax))-xMin)/(xMax-xMin)*RisX);
  }
  else
    Tempx = (int) ( (double)((v1->Abscissa(CoordY,NElMax)-xMin)*RisX)/(double)(xMax-xMin) );
  p->drawLine( QPoint(Tempx ,0), QPoint(Tempx ,RisY));
  p->drawLine( QPoint(0,NElMinY), QPoint(RisX,NElMinY));
  p->drawLine( QPoint(0 ,NElMaxY), QPoint(RisX,NElMaxY));
  //emit SegnaleGrafico(1);
}
void ElementiGrafici::GrStampante(QPainter *p){
  GrMessage("Visualizza.Dis.GrStampante");
  Font.setFamily("Helvetica");
  if(!IfStampa){
    RisX = width();
    RisY = height();
    QRect view(0,0,width(),height());
    p->setViewport(view);
    ImpInterY(0,RisY);
    StampaX = 0;
    StampaY = 0;
    PointSize = 2;
    Font.setPointSize(12);
    p->setFont(Font);
    p->setFont(QFont("Palatino",12 ));
    int PosLabel = (int)(RisX/2. - 2.*nomeEtX.length());
    //int PosLabel = 0;
    p->save();
    p->translate(PosLabel,height()-40);
    QTextDocument doc;
    doc.setDefaultFont(QFont("Palatino",18 ));
    doc.setHtml(nomeEtX); // or populate using QTextCursor
    doc.drawContents(p);
    p->restore();
    p->save();
    PosLabel = RisY - (int)(RisY/2. - 2.*nomeEtY.length());
    p->translate(width()-40,PosLabel);
    p->rotate(270);
    doc.setHtml(nomeEtY); // or populate using QTextCursor
    doc.drawContents(p);
    p->restore();
  }
  else{
    Q3PaintDeviceMetrics metrics( p->device() );
    int dpiy = metrics.logicalDpiY();
    int margin = (int)  ((2/2.54)*dpiy);
    //QRect view(3*margin,margin,metrics.width() - 3*margin, metrics.height() - 2*margin );
    QRect window(0.,0.,RatioWidthHeight*metrics.height(),metrics.height());
    p->setWindow(window);
    QRect view(3*margin,margin,metrics.height() - 2*margin,metrics.height() - 2*margin );
    p->setViewport(view);
    Font.setPointSize(32);
    p->setFont(Font);
    p->setFont(QFont("Helvetica",FontSize));
    PointSize = 5;
    RisX = (int)(.9*metrics.width());
    RisY = (int)(.9*metrics.height());
    StampaX = 43;
    StampaY = -60;
    NElMinY = (int) (NElMinY*metrics.height()/(double)RisY);
    NElMaxY = (int) (NElMaxY*metrics.height()/(double)RisY);
    p->drawText(RisX/2-nomeTit.length(),-15,nomeTit);
    int PosLabel = (int)(RisX/2. - 2.*nomeEtX.length());
    //int PosLabel = 0;
    QString Label1(nomeTit);
    // p->drawText((int)(RisX/2. - nomeEtX.length()*7.),metrics.height()+40,nomeEtX,-1,Qt::AlignTop);
    p->save();
    p->translate(PosLabel,metrics.height()-10);
    QTextDocument doc;
    doc.setDefaultFont(QFont("Helvetica",FontSize));
    doc.setHtml(nomeEtX); // or populate using QTextCursor
    doc.drawContents(p);
    p->restore();
    p->save();
    PosLabel = RisY - (int)(RisY/2. - 2.*nomeEtY.length());
    p->translate(-(int)(3.5*margin),PosLabel);
    p->rotate(270);
    doc.setHtml(nomeEtY); // or populate using QTextCursor
    doc.drawContents(p);
    // PosLabel = RisY - (int)(RisY/2. - nomeEtY.length()*7.);
    //p->drawText(0,0,nomeEtY);//,-1,Qt::AlignRight);
    p->restore();
  }
}
void ElementiGrafici::GrGriglia(QPainter *p){
  GrMessage("Visualizza.Dis.GrGriglia");
  QString NumLabel = QString();
  if(!IfStampa)
    p->setPen(QPen(Qt::black,1,Qt::DashLine,Qt::RoundCap,Qt::RoundJoin));
  else 
    p->setPen(QPen(Qt::black,3,Qt::SolidLine,Qt::RoundCap,Qt::RoundJoin));
  //p->setPen( QPen(Qt::black,1) );
  //p->setBrush( Qt::black );
  //p->setPen( Qt::DashLine );
  p->drawLine(QPoint(0,RisY),QPoint(RisX,RisY) );	
  p->drawLine(QPoint(RisX,0),QPoint(RisX,RisY) );	
  p->drawLine(QPoint(0,0),QPoint(RisX,0) );
  p->drawLine(QPoint(0,0),QPoint(0,RisY) );
  if(!DIS_IF_TYPE(LogLog,DIS_LOGY)){
    //---------------------Griglia-y---------------------
    for(int i=0;i<=GrigliaY;i++){
      Tempy = (int) (RisY/(double)GrigliaY*i);
      if(!IfStampa)
	p->drawLine(QPoint(0,Tempy),QPoint(RisX,Tempy) );
      else{
	p->drawLine(QPoint(0,Tempy),QPoint(10,Tempy) );	
	p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-10,Tempy) );
	if(i<GrigliaY-1) {
	  int RiTempY = Tempy+(int)(.5*RisY/(double)GrigliaY);
	  p->drawLine(QPoint(0,RiTempY),QPoint(5,RiTempY) );
	  p->drawLine(QPoint(RisX,RiTempY),QPoint(RisX-5,RiTempY) );
	}
      }
      NumLabel.setNum((double)(yMax-yMin)/(double)GrigliaY*i+yMin,FormatPrecY,DigPrecY);
      if(i==GrigliaY && !IfStampa) Tempy -= 10;
      if(i==0 && !IfStampa) Tempy += 10;
      if(IfStampa) Tempy += 10;
      if(!IfStampa) p->drawText(StampaY,RisY-Tempy,NumLabel);	
      else p->drawText(StampaY,RisY-Tempy-15,50,20,Qt::AlignRight|Qt::TextDontClip,NumLabel);	
    }
  }
  if(!DIS_IF_TYPE(LogLog,DIS_LOGX)){
    //---------------------Griglia-x---------------------
    for(int i = 0;i<=GrigliaX;i++){
      Tempx = (int) (RisX/(double)GrigliaX*i);
      if(!IfStampa)
	p->drawLine(QPoint(Tempx,0),QPoint(Tempx,RisY) );
      else{
	p->drawLine(QPoint(Tempx,RisY-10),QPoint(Tempx,RisY) );
	p->drawLine(QPoint(Tempx,0),QPoint(Tempx,10) );
	if(i<GrigliaX-1) {
	  int RiTempX = Tempx+(int)(.5*RisX/(double)GrigliaX);
	  p->drawLine(QPoint(RiTempX,RisY-5),QPoint(RiTempX,RisY) );
	  p->drawLine(QPoint(RiTempX,0),QPoint(RiTempX,5) );
	}
      }
      NumLabel.setNum((double)(xMax-xMin)/(double)GrigliaX*i+xMin,FormatPrecX,DigPrecX);
      //p->drawText(Tempx,RisY+StampaX,QString(stringa));
      if(IfStampa) Tempx -= (int)(9.*NumLabel.length());
      p->drawText(Tempx,RisY+StampaX,NumLabel);    
      //      p->drawText(Tempx-(int)(.5*RisX/(double)GrigliaX),RisY+StampaX,NumLabel);    
    }
    if(!IfStampa){
      sprintf(stringa,"%% %.4g",(double)(yMax-yMin)/(double)((yMax+yMin)*.5)*100.);
      p->drawText(RisX/2,10,QString(stringa));
    }
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)){
    //--------LogLog-y-----------------------
    NyMin = (int)yMin;
    NyMax = (int)yMax;
    Delta = 1;
    if(NyMin < 0){
      for(int i=NyMin;i<0;i++){
	Delta /= 10;
      }
    }
    else if(NyMin > 0){
      for(int i=0;i<NyMin;i++){
	Delta *= 10;
      }
    }
    if( NyMin < NyMax){
      for(int i = NyMin;i<=NyMax;i++){
	Tempy = (int)( (double)(i-yMin)/((yMax-yMin))*RisY);
	if(Tempy < 0 || Tempy >= RisY) continue;
	sprintf(stringa,"%.2g",Delta);
	if(!IfStampa) p->drawText(StampaY,RisY-Tempy,stringa);	
	else p->drawText(StampaY,RisY-Tempy+5,50,20,Qt::AlignRight|Qt::TextDontClip,stringa);	
	//if(i==NyMax) Tempy -= 10;
	//if(i==NyMin) Tempy += 10;
	if(!IfStampa)
	  p->drawLine( QPoint(0,Tempy),QPoint(RisX,Tempy));
	else{
	  p->drawLine(QPoint(0,Tempy),QPoint(15,Tempy) );
	  p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-15,Tempy) );
	}
	Tempy = (int)( (double)(i-(yMin)+.30)/((yMax-yMin))*RisY);
	if(Tempy < 0 || Tempy >= RisY) continue;
	p->drawLine(QPoint(0,Tempy),QPoint(10,Tempy) );
	p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-10,Tempy) );
	Tempy = (int)( (double)(i-(yMin)+.60)/((yMax-yMin))*RisY);
	if(Tempy < 0 || Tempy >= RisY) continue;
	p->drawLine(QPoint(0,Tempy),QPoint(10,Tempy) );
	p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-10,Tempy) );
	Tempy = (int)( (double)(i-(yMin)+.77)/((yMax-yMin))*RisY);
	if(Tempy < 0 || Tempy >= RisY) continue;
	p->drawLine(QPoint(0,Tempy),QPoint(10,Tempy) );
	p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-10,Tempy) );
	Tempy = (int)( (double)(i-(yMin)+.90)/((yMax-yMin))*RisY);
	if(Tempy < 0 || Tempy >= RisY) continue;
	p->drawLine(QPoint(0,Tempy),QPoint(10,Tempy) );
	p->drawLine(QPoint(RisX,Tempy),QPoint(RisX-10,Tempy) );
	Delta *= 10;
      }
    }
    if(NyMin == NyMax){
      Tempy = (int)(RisY/2.);
      p->drawLine(QPoint(0,RisY-Tempy),QPoint(RisX,RisY-Tempy));
      sprintf(stringa,"%.2g",yMin/2.);
      p->drawText(0,RisY- Tempy -7,QString(stringa));
    }
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)){
    //--------------LogLog-x--------------------------
    NxMin = (int)xMin;
    NxMax = (int)xMax;
    Delta = 1;
    if(NxMin < 0){
      for(int i=NxMin;i<0;i++){
	Delta /= 10;
      }
    }
    else if(NxMin > 0){
      for(int i=0;i<NxMin;i++){
	Delta *= 10;
      }
    }
    double Ratio = RisX /(xMax-xMin);
    for(int i = NxMin;i<=NxMax;i++){
      sprintf(stringa,"%.2g",Delta);
      double x = log10(Delta);
      Tempx = (int)( (x - xMin)*Ratio);
      if(Tempx < 0 || Tempx >= RisX) continue;
      if(IfStampa) p->drawText(Tempx-strlen(stringa)*9., RisY+35,QString(stringa));
      else p->drawText( Tempx, RisY ,QString(stringa));
      if(!IfStampa)
	p->drawLine( QPoint(Tempx,RisY),QPoint(Tempx,0));
      else{
	p->drawLine(QPoint(Tempx,RisY-15),QPoint(Tempx,RisY) );
	p->drawLine(QPoint(Tempx,0),QPoint(Tempx,15) );
      }
      Tempx = (int)( (double)(x+.30-xMin)*Ratio);
      if(Tempx < 0 || Tempx >= RisX) continue;
      p->drawLine(QPoint(Tempx,RisY-10),QPoint(Tempx,RisY) );
      p->drawLine(QPoint(Tempx,0),QPoint(Tempx,10) );
      Tempx = (int)( (double)(x+.60-xMin)*Ratio);
      if(Tempx < 0 || Tempx >= RisX) continue;
      p->drawLine(QPoint(Tempx,RisY-10),QPoint(Tempx,RisY) );
      p->drawLine(QPoint(Tempx,0),QPoint(Tempx,10) );
      Tempx = (int)( (double)(x+.77-xMin)*Ratio);
      if(Tempx < 0 || Tempx >= RisX) continue;
      p->drawLine(QPoint(Tempx,RisY-10),QPoint(Tempx,RisY) );
      p->drawLine(QPoint(Tempx,0),QPoint(Tempx,10) );
      Tempx = (int)( (double)(x+.90-xMin)*Ratio);
      if(Tempx < 0 || Tempx >= RisX) continue;
      p->drawLine(QPoint(Tempx,RisY-10),QPoint(Tempx,RisY) );
      p->drawLine(QPoint(Tempx,0),QPoint(Tempx,10) );
      Delta *= 10;
    }
  }
}
void ElementiGrafici::GrRiscala(){
  GrMessage("Visualizza.Dis.GrRiscala.x");
  if(DIS_IF_TYPE(IfRiscala,RIS_UNO) ){
    if(v1->IsAbscissa(CoordX)){
      xMin = v1->pMin(CoordX,NVisMin,NVisMax);
      xMax = v1->pMax(CoordX,NVisMin,NVisMax);
    }
    else{
      xMin = NVisMin;
      xMax = NVisMax;
    }
    GrMessage("Visualizza.Dis.GrRiscala.y");
    yMin = v1->pMin(CoordY,NVisMin,NVisMax);
    yMax = v1->pMax(CoordY,NVisMin,NVisMax);
  }
  else if(DIS_IF_TYPE(IfRiscala,RIS_TUTTI) ){
    GrMessage("Visualizza.Dis.GrRiscala.Tutto %d",NVar);
    v1->pMinMaxGlob(NVisMin,NVisMax);
    v1->pGlobBorder(&xMin,&xMax,&yMin,&yMax);
    if(!v1->IsAbscissa(CoordX)){
      xMin = NVisMin;
      xMax = NVisMax;
    }
  }
  if( DIS_IF_TYPE(IfDisegna,DIS_PUNTI) ){
    GrMessage("Visualizza.Dis.GrRiscala.Punti");
    yMin = v1->PuntiMin();
    yMax = v1->PuntiMax();
  }
  if( DIS_IF_TYPE(IfDisegna,DIS_MOMENTI) ){
    GrMessage("Visualizza.Dis.GrRiscala.Momenti");
    yMin = m1.yMin;
    yMax = m1.yMax;
    xMin = m1.Min;
    xMax = m1.Max;
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGX)){// LogLog
    GrMessage("Visualizza.Dis.GrRiscala.Logx");
    if(v1->IsAbscissa(CoordX)){
      xMin = v1->pMinLog(CoordX,NVisMin,NVisMax);
      xMax = v1->pMaxLog(CoordX,NVisMin,NVisMax);
    }
    else{
      xMin = log10(NVisMin);
      xMax = log10(NVisMax);
    }
  }
  if(DIS_IF_TYPE(LogLog,DIS_LOGY)){// LogLog
    GrMessage("Visualizza.Dis.GrRiscala.Logy");
    if( DIS_IF_TYPE(IfRiscala,RIS_UNO) ){
      yMin = v1->pMinLog(CoordY,NVisMin,NVisMax);
      yMax = v1->pMaxLog(CoordY,NVisMin,NVisMax);
    }
    else if( DIS_IF_TYPE(IfRiscala,RIS_TUTTI) ){
      yMin = v1->pMinGlobLog(NVisMin,NVisMax);
      yMax = v1->pMaxGlobLog(NVisMin,NVisMax);
    }
    // if(yMin <= 0.){
    //   char SErr[120];
    //   sprintf(SErr,"Ordinates negative or null %f, no log log possible",yMin);
    //   ErrPrima->message(QString(SErr));
    //   //QMessageBox::information(this,"Visualizza","Numeri negativi o nulli\n Non posso visualizzare LogLog","Va Be!");
    //   DIS_REM_TYPE(LogLog,DIS_LOGY);
    //   emit LogyCambiato(LogLog);
    // }
  }
  GrMessage("xBound %lf-%lf yBound %lf-%lf\n",xMin,xMax,yMin,yMax);
  ScalaTopoX = xMax-xMin;
  ScalaTopoY = yMax-yMin;
  // yMax = (yMax*1.1);
  // yMin = (yMin*0.9);
  // xMax = (xMax*1.1);
  // xMin = (xMin*0.9);
  yMax += (yMax-yMin)*ViewportY;
  yMin -= (yMax-yMin)*ViewportY;
  xMax += (xMax-xMin)*ViewportX;
  xMin -= (xMax-xMin)*ViewportX;
}
void ElementiGrafici::GrLegenda(QPainter *p){
  GrMessage("Visualizza.Dis.GrLegenda");
  //QRect Leg((int)(PosLegenda[0]*width()),(int)((1.-PosLegenda[1])*height()),
  //	    (int)(PosLegenda[2]*width()),(int)((1.-PosLegenda[3])*height()));
  // QRect Leg((int)(.8*width()),(int)((1.-.8)*height()),
  //  	    (int)(.9*width()),(int)((1.-.9)*height()));
  //p->drawText(Leg,Qt::AlignCenter,tr("Ciccia"));
  //p->drawText((int)(PosLegenda[0]*width()),(int)((1.-PosLegenda[1])*height()),tr("Chem Pot"));
  //p->drawText((int)(.8*width()),(int)((1.-.8)*height()),tr("Chem Pot"));
  int PosLegX=(int)(PosLegenda[0]*RisX);
  int PosLegY=(int)((1.-PosLegenda[1])*RisY);
  int PosFinL=PosLegX, PosFinY=PosLegY;
  int IncrY = 15;
  if(IfStampa) IncrY = 30;
  for(int s=0;s<NVar;s++){
    if(v1->IsAbscissa(s)) continue;
    if(LineaCome[s]==0) continue;
    p->setBrush( GrLinee[s] );
    p->setPen( QPen( GrLinee[s],2 ) );
    Quadrati = DIS_IF_TYPE(LineaCome[s],LINEA_PUNTO);
    Linee = DIS_IF_TYPE(LineaCome[s],LINEA_TRATTO);
    Punto = QPoint(PosLegX,PosFinY-5);
    Punto1 = QPoint(PosLegX+15,PosFinY-5);
    p->drawLine(Punto,Punto1);
    p->drawText(PosLegX+25,PosFinY,GrLabel[s]);
    PosFinY += IncrY;
  }
}
void ElementiGrafici::GrScript(char *File2Read,QPainter *p){
  GrMessage("Visualizza.Dis.GrScript");
  FILE *File2Open;
  if( (File2Open = fopen(File2Read,"r"))==0){
    //    printf("The file %s does not exist\n");
    return ;
  }
  char cLine[256];
  double Color[4];
  double Pos[4];
  double Hue[4];
  int iHue[4];
  for(int k=0;!(fgets(cLine,256,File2Open)==NULL);k++){
    if(strncmp(cLine,"PutStr",6)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Color,Color+1,Color+2,Color+3);
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Pos,Pos+1,Pos+2,Pos+3);
      char String[256];
      fgets(cLine,256,File2Open);
      //for(int c=0;c<strlen(cLine);c++) String[c] = cLine[c];
      sprintf(String,"%s",cLine);
      QColor ColorLabel = QColor((int)(255.*Color[0]),(int)(255.*Color[1]),(int)(255.*Color[2]),(int)(255.*Color[3]));
      if(IfStampa) Pos[3] *= 2.;
      p->setFont(QFont("Palatino",Pos[3]));
      p->setBrush(ColorLabel);
      p->save();
      p->rotate(Pos[2]);
      p->drawText(Pos[0]*RisX,(1.-Pos[1])*RisY,QString(String));
      p->restore();
      if(IfStampa)
	p->setFont(QFont("Palatino",FontSize));
      else 
	p->setFont(QFont("Palatino",12));
    }
    else if(strncmp(cLine,"PutLine",7)==0){
      //Colore
      fgets(cLine,256,File2Open);
      p->setPen(QPen(Qt::blue,3,Qt::SolidLine,Qt::RoundCap,Qt::RoundJoin));
      sscanf(cLine,"%lf %lf %lf %lf\n",Color,Color+1,Color+2,Color+3);
      QColor ColorLabel = QColor((int)(255.*Color[0]),(int)(255.*Color[1]),(int)(255.*Color[2]),(int)(255.*Color[3]));
      p->setBrush(ColorLabel);
      //Posizione
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Pos,Pos+1,Pos+2,Pos+3);
      int Fromx = (int)(RisX*Pos[0]);
      int Fromy = (int)(RisY*(1.-Pos[1]));
      int Tox = (int)(RisX*Pos[2]);
      int Toy = (int)(RisY*(1.-Pos[3]));
      Punto = QPoint(Fromx,Fromy);
      Punto1 = QPoint(Tox,Toy);
      p->drawLine(Punto,Punto1);
    }
    else if(strncmp(cLine,"Text",4)==0){
      //Colore
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Hue,Hue+1,Hue+2,Hue+3);
      for(int i=0;i<4;i++) iHue[i] = (int)(255.*Hue[i]);
      p->setBrush( QColor(iHue[0],iHue[1],iHue[2],iHue[3]) );
      //Posizione
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf\n",Pos,Pos+1);
      int Fromx = (int)(width()*Pos[0]);
      int Fromy = (int)(height()*(1.-Pos[1]));
      //Testo
      fgets(cLine,256,File2Open);
      QChar ChType = QChar::NoCategory;
      QString String;
      int Posx = Fromx;
      int Posy = Fromy;
      int FontSize = 12;
      int IfSymbol = 0;
      int IfUpper = 0;
      for(int i=0,j=0;i<strlen(cLine);i++){
	if(cLine[i] == '\\'){
	  IfSymbol = 0x0350;
	  continue;
	}
	else if(cLine[i] == '^'){
	  IfUpper = 1;
	  continue;
	}
	else if(cLine[i] == '_'){
	  IfUpper = -1;
	  continue;
	}
	else if(cLine[i] == '%'){
	  IfSymbol = 0;
	  IfUpper = 1;
	  continue;
	}
	//String.insert(i,QChar(sTable(cLine+j)));
	String.insert(j,QChar((uint)(cLine[i]+IfSymbol)));
	j++;
      }
      p->drawText(Fromx,Fromy,String,-1,Qt::AlignTop);
    }
    else if(strncmp(cLine,"Picture",4)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Pos,Pos+1,Pos+2,Pos+3);
      int Fromx = (int)(width()*Pos[0]);
      int Fromy = (int)(height()*(1.-Pos[1]));
      int Tox = (int)(width()*Pos[2]);
      int Toy = (int)(height()*(1.-Pos[3]));
      char String[512];
      fgets(cLine,512,File2Open);
      for(int c=0;c<strlen(cLine)-1;c++) String[c] = cLine[c];
      p->drawImage(QRect(Fromx,Fromy,Tox,Toy),QImage(QString(String),0));
    }
    else if(strncmp(cLine,"Background",10)==0){
    }
  }
  fclose(File2Open);
}
void ElementiGrafici::GrConf(char *File2Read){
  GrMessage("Visualizza.Dis.GrConf");
  FILE *File2Open;
  if( (File2Open = fopen(File2Read,"r"))==0){
    //    printf("The file %s does not exist\n");
    return ;
  }
  char cLine[256];
  for(int k=0;!(fgets(cLine,256,File2Open)==NULL);k++){
    //printf("%s",cLine);
    //if(cLine[0] == '#') continue;
    if(strncmp(cLine,"NGrid",5)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d %d\n",&GrigliaX,&GrigliaY);
    }
    else if(strncmp(cLine,"XFormula",8)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%s\n",XFormula);
      v1->setXFormula(XFormula);
    }
    else if(strncmp(cLine,"YFormula",8)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%s\n",YFormula);
      v1->setYFormula(YFormula);
    }
    else if(strncmp(cLine,"XBound",6)==0){
      fgets(cLine,256,File2Open);
      if( DIS_IF_TYPE(IfRiscala,RIS_TUTTI) )
      	sscanf(cLine,"%lf %lf\n",&xMin,&xMax);
      if(DIS_IF_TYPE(LogLog,DIS_LOGX)){
	if(xMin < 0 ) printf("xMin can't be negative in log log\n");
	else xMin = log10(xMin);
	if(xMax < 0 ) printf("xMax can't be negative in log log\n");
	else xMax = log10(xMax);
      }
    }
    else if(strncmp(cLine,"YBound",6)==0){
      fgets(cLine,256,File2Open);
      if( DIS_IF_TYPE(IfRiscala,RIS_TUTTI) )
      	sscanf(cLine,"%lf %lf\n",&yMin,&yMax);
      if(DIS_IF_TYPE(LogLog,DIS_LOGY)){
	if(yMin < 0 ) printf("yMin can't be negative in log log\n");
	else yMin = log10(yMin);
	if(yMax < 0 ) printf("yMax can't be negative in log log\n");
	else yMax = log10(yMax);
      }
    }
    else if(strncmp(cLine,"PosLegend",6)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",(PosLegenda+0),(PosLegenda+1),(PosLegenda+2),(PosLegenda+3));
    }
    else if(strncmp(cLine,"PosInterp",6)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf\n",(PosInterp+0),(PosInterp+1));
    }
    else if(strncmp(cLine,"DigPrecX",8)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d %c\n",&DigPrecX,&FormatPrecX);
    }
    else if(strncmp(cLine,"DigPrecY",8)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d %c\n",&DigPrecY,&FormatPrecY);
    }
    else if(strncmp(cLine,"ViewportX",9)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf\n",&ViewportX);
    }
    else if(strncmp(cLine,"ViewportY",9)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf\n",&ViewportY);
    }
    else if(strncmp(cLine,"FontSize",8)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d\n",&FontSize);
    }
    else if(strncmp(cLine,"RatioWidthHeight",16)==0){
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf\n",&RatioWidthHeight);
    }
    else if(strncmp(cLine,"LabelX",6)==0){
      char *String = (char *)calloc(265,sizeof(char));
      fgets(cLine,256,File2Open);
      for(int c=0;c<strlen(cLine);c++) String[c] = cLine[c];
      nomeEtX.clear();
      StringToUnicode(String,&nomeEtX);
      free(String);
    }
    else if(strncmp(cLine,"LabelY",6)==0){
      char *String = (char *)calloc(265,sizeof(char));
      fgets(cLine,256,File2Open);
      for(int c=0;c<strlen(cLine);c++) String[c] = cLine[c];
      nomeEtY.clear();
      StringToUnicode(String,&nomeEtY);
      free(String);
    }
    else if(strncmp(cLine,"Title",5)==0){
      char *String = (char *)calloc(265,sizeof(char));
      fgets(cLine,256,File2Open);
      for(int c=0;c<strlen(cLine);c++) String[c] = cLine[c];
      nomeTit.clear();
      StringToUnicode(String,&nomeTit);
      free(String);
    }
    else if(strncmp(cLine,"Line",4)==0){
      int NLine;
      int Type;
      double Color[4];
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d\n",&NLine);
      if(NLine >= NVar) {
	//printf("The %d th column is not present\n",NLine);
	NLine = 0;
	continue;
      }
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%lf %lf %lf %lf\n",Color,Color+1,Color+2,Color+3);
      fgets(cLine,256,File2Open);
      sscanf(cLine,"%d\n",&Type);
      //printf("%d %d %lf %lf %lf\n",NLine,Type,Color[0],Color[1],Color[2],Color[3]);
      QColor ColorSet = QColor((int)(255.*Color[0]),(int)(255.*Color[1]),(int)(255.*Color[2]),(int)(255.*Color[3]));
      GrLinee[NLine] = ColorSet;
      LineaCome[NLine] = Type;
      fgets(cLine,256,File2Open);
      GrLabel[NLine].clear();
      StringToUnicode(cLine,GrLabel + NLine);
      GrMessage("Line n: %d type: %d name: %s color: %lf %lf %lf\n",NLine,Type,cLine,Color[0],Color[1],Color[2]);
    }
  }
  GrMessage("NGrid %d %d XBond %lf %lf YBound %lf %lf\n",GrigliaX,GrigliaY,xMin,xMax,yMin,yMax);
  fclose(File2Open);
}
void ElementiGrafici::StringToUnicode(char *cLine,QString *Label){
  QChar ChType = QChar::NoCategory;
  int UnicodeRow = 0;
  int CharOffset = 0;
  for(int i=0,j=0;i<strlen(cLine);i++){
    int Char = cLine[i];
    if(cLine[i] == '\\' && i < strlen(cLine)-1){
      if(cLine[i+1] == 's'){
	//ChType = QChar::Symbol_Math;
	CharOffset = 0x0350;
	UnicodeRow = 3;
      }
      else if(cLine[i+1] == 'l')
	ChType = QChar::Letter_Lowercase;
      else if(cLine[i+1] == 'u')
	ChType = QChar::Letter_Uppercase;
      else if(cLine[i+1] == 'n'){
	UnicodeRow = 0;
	CharOffset = 0;
	ChType = QChar::NoCategory;
      }
      i+=1;
      continue;
    }
    else 
      ChType = QChar::NoCategory;
    //QChar Character = QChar(QChar(Char).cell(),QChar(QChar::Symbol_Math).row());
    QChar Char1 = QChar(Char+CharOffset);
    //printf("%c %d %d %x\n",Char1.toAscii(),Char1.cell(),Char1.row(),Char1.unicode());
    //GrLabel[NLine].insert(i,Char1);
    Label->insert(j,QChar(Char1.cell(),UnicodeRow));//QChar(Char));
    j++;
  }
  //GrLabel[NLine].append(cLine);
  //sprintf(LineaQuale[NLine],"%s",cLine);
}
