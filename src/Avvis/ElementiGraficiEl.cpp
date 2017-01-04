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
#include <qpainter.h>
#include <Q3PaintDeviceMetrics>
;
void ElementiGrafici::PuntiSegnale(){
  PuntiMax = NMax;
  v1->ImpCoordY(CoordY);
  //v1->ImpCoordX(CoordX,CoordY);
  sprintf(stringa,"Mostro i %d Punti Coord %d %d/%d",PuntiMax,CoordX,CoordY,NVar);
  IfDisegna = DIS_SEGNALE;
  InfoStato(stringa);
  ImpInter(PuntiMax);
  ImpInter(0,PuntiMax);
  ImpNVisMin(0); ImpNVisMax(PuntiMax);
  ImpNElMin(0); ImpNElMax(PuntiMax);
  Ridisegna();
}
void ElementiGrafici::PuntiSommaSegnali(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  v1->SommaSegnali();
  sprintf(stringa,"Somma di tutti i segnali\n");
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiDistribuzione(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  m1 = v1->DistrGaussSegnale(NElMin,NElMax,IfNorm);
  sprintf(stringa,"Distribuzione per %d NBin di %d Punti",NBin,NElMax-NElMin);
  DIS_ADD_TYPE(IfDisegna,DIS_MOMENTI);
  ImpNVisMin(0);
  ImpNVisMax(NBin);
  ImpInter(NBin);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiSpettro(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  //  QMessageBox::information(this,"Visualizza","Sto calcolando\n","Attendere");
  PuntiMax = v1->SpettroSegnale(NElMin,NElMax);
  sprintf(stringa,"Spettro da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiNormalizza(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  int Campioni=0;
  //  QMessageBox::information(this,"Visualizza","Sto calcolando\n","Attendere");
  Campioni = v1->NormalizzaSegnale(NElMin,NElMax);
  sprintf(stringa,"Normalizzazione del segnale da %d campioni",Campioni);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiModulo(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  int Campioni=0;
  //  QMessageBox::information(this,"Visualizza","Sto calcolando\n","Attendere");
  v1->ModuloSegnale();
  Campioni = v1->ElabSegnale(NElMin,NElMax);
  sprintf(stringa,"Modulo del segnale da %d campioni",Campioni);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiRadice(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  //  QMessageBox::information(this,"Visualizza","Sto calcolando\n","Attendere");
  v1->RadiceSegnale();
  PuntiMax = NMax;
  sprintf(stringa,"Radice di %d NBin",PuntiMax);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiAutocor(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  PuntiMax = v1->AutocorSegnale(NElMin,NElMax);
  sprintf(stringa,"Autocorrelazione da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiSum(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  double Resp = v1->SumSegnale(CoordY,NElMin,NElMax);
  sprintf(stringa,"Somma %lf da %d a %d di %d Punti",Resp,NElMin,NElMax,NElMax-NElMin);
  //  IfDisegna = DIS_SEGNALE;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiInterRett(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  r1 = v1->InterRettSegnale(CoordY,NElMin,NElMax,LogLog);
  PuntiMax = NMax;
  sprintf(stringa,"Interpolazione lineare da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_ADD_TYPE(IfDisegna,DIS_RETTA);
  DIS_REM_TYPE(IfDisegna,DIS_GAUSS);
  DIS_REM_TYPE(IfDisegna,DIS_PARABOLA);
  DIS_REM_TYPE(IfDisegna,DIS_EXP);
  DIS_ADD_TYPE(IfDisegna,DIS_BLU);
  InfoStato(stringa);
  repaint();
}  
void ElementiGrafici::PuntiInterExp(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  r1 =  v1->InterExpSegnale(CoordY,NElMin,NElMax,LogLog);
  PuntiMax = NMax;
  sprintf(stringa,"Interpolazione lineare da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_REM_TYPE(IfDisegna,DIS_RETTA);
  DIS_REM_TYPE(IfDisegna,DIS_GAUSS);
  DIS_REM_TYPE(IfDisegna,DIS_PARABOLA);
  DIS_ADD_TYPE(IfDisegna,DIS_EXP);
  DIS_ADD_TYPE(IfDisegna,DIS_BLU);
  InfoStato(stringa);
  repaint();
}  
void ElementiGrafici::PuntiParabola(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  p1 = v1->ParabolaSegnale(CoordY,NElMin,NElMax,LogLog);
  printf("%.2g + %.2g x + %.2g x^2",p1.a0,p1.a1,p1.a2);
  PuntiMax = NMax;
  sprintf(stringa,"Interpolazione parabolica da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_REM_TYPE(IfDisegna,DIS_RETTA);
  DIS_REM_TYPE(IfDisegna,DIS_GAUSS);
  DIS_ADD_TYPE(IfDisegna,DIS_PARABOLA);
  DIS_REM_TYPE(IfDisegna,DIS_EXP);
  DIS_ADD_TYPE(IfDisegna,DIS_BLU);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiInterGauss(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  m1 = v1->InterGaussSegnale(CoordY,NElMin,NElMax,LogLog);
  PuntiMax = NMax;
  sprintf(stringa,"Interpolazione gaussiana da %d a %d di %d Punti",NElMin,NElMax,NElMax-NElMin);
  DIS_REM_TYPE(IfDisegna,DIS_RETTA);
  DIS_ADD_TYPE(IfDisegna,DIS_GAUSS);
  DIS_REM_TYPE(IfDisegna,DIS_PARABOLA);
  DIS_REM_TYPE(IfDisegna,DIS_EXP);
  DIS_ADD_TYPE(IfDisegna,DIS_BLU);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiMediaMob(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  PuntiMax = NVisMax = v1->MediaMobSegnale(NMobile);
  sprintf(stringa,"MediaMobile di %d Punti divisi in %d Parti",PuntiMax,NMobile);
  InfoStato(stringa);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNVisMin(0);
  ImpInter(NVisMax);
  ImpNVisMax(NVisMax);
  repaint();
}
void ElementiGrafici::PuntiCorrelaADue(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  PuntiMax = v1->CorrelaADuePunti(NCorrela);
  ImpNCorrela(NCorrela);
  sprintf(stringa,"MediaMobile di %d Punti divisi in %d Parti",PuntiMax,NCorrela);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiAutosimilarita(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  PuntiMax = 10;
  v1->AutosimilaritaSegnale(PuntiMax);
  sprintf(stringa,"Autosimilarita di %d Punti su %d Potenze",NMax,PuntiMax);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiIntegrale(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  //void (Matematica::*Qui)(double *,double *,int) = &Matematica::Integrazione;
  //v1->Elab = &Matematica::Integrazione;
  //v1->ElabSegnale(NElMin,NElMax);
  double Int = v1->IntSegnale();
  PuntiMax = NMax;
  sprintf(stringa,"Integrale su %d Punti= %lf",PuntiMax,Int);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiDerivata(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  v1->DerivataSegnale();
  v1->ElabSegnale(NElMin,NElMax);
  PuntiMax = NElMax - NElMin - 2;
  sprintf(stringa,"Derivata su %d Punti",PuntiMax);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::PuntiVarie(){
  if(PrimaVolta){
    ErrPrima->message(QString("Aprire un file prima"));
    return;
  }
  //v1->VarieSegnale();
  //v1->ElabSegnale(NElMin,NElMax);
  double Integr = v1->VarieSegnale(NElMin,NElMax)*32./200.;
  PuntiMax = NElMax - NElMin -2;
  NVisMax = PuntiMax;
  NVisMin = 0;
  double Volume = NMax/(32.*32.*32.);
  for(int n=0;n<PuntiMax;n++){
    //Punti[n] *= Volume;
  }
  sprintf(stringa,"Surface tension %lf",Integr);
  DIS_ADD_TYPE(IfDisegna,DIS_PUNTI);
  IfDisegna = DIS_PUNTI;
  ImpNMobile(1);
  InfoStato(stringa);
  repaint();
}
void ElementiGrafici::StampaFile(){
  QPainter paint;
  InfoStato("Stampando...");
  Stampante->setOutputFileName(QString(nomeSalva));
  if(Stampante->setup(this))
    if(!paint.begin(Stampante) )
      return;
  IfStampa = 1;
  DisegnaPunti(&paint);
  InfoStato("Stampato");
  IfStampa = 0;
}
void ElementiGrafici::PuntaExt(double **ExtSt,int NVar){
  v1->Punta(ExtSt,NVar);
}
