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
***********************************************************************/#include "ElementiGrafici.h"

void ElementiGrafici::Ridisegna(){
//   InfoStato(stringa);
//   ImpInter(PuntiMax);
//   ImpInter(0,PuntiMax);
//   ImpNVisMin(0); ImpNVisMax(PuntiMax);
//   ImpNElMin(0); ImpNElMax(PuntiMax);
  repaint();
  //  if(!VoltaMomenti) emit InterVisCambiato(NVisMin,NVisMax);
}
void ElementiGrafici::DisegnaLinee(){
  for(int s=0;s<NVar;s++){
    if(DIS_IF_TYPE(LineaCome[s],LINEA_TRATTO))
      LineaCome[s] = LINEA_PUNTO;
    else if(DIS_IF_TYPE(LineaCome[s],LINEA_PUNTO))
      LineaCome[s] = LINEA_TRATTO;
  }
  repaint();
}
void ElementiGrafici::DisegnaPunti(){
  for(int s=0;s<NVar;s++){
    if(DIS_IF_TYPE(LineaCome[s],LINEA_TRATTO))
      LineaCome[s] = LINEA_PUNTOETRATTO;
    else if(DIS_IF_TYPE(LineaCome[s],LINEA_PUNTOETRATTO))
      LineaCome[s] = LINEA_TRATTO;
  }
  repaint();
}
void ElementiGrafici::DisegnaGriglia(){
  Griglia = !Griglia;
  repaint();
}
void ElementiGrafici::DisegnaLogx(){
  if(DIS_IF_TYPE(LogLog,DIS_LOGX))
    DIS_REM_TYPE(LogLog,DIS_LOGX);
  else
    DIS_ADD_TYPE(LogLog,DIS_LOGX);
  repaint();
}
void ElementiGrafici::DisegnaLogy(){
  if(DIS_IF_TYPE(LogLog,DIS_LOGY))
    DIS_REM_TYPE(LogLog,DIS_LOGY);
  else
    DIS_ADD_TYPE(LogLog,DIS_LOGY);
  repaint();
}
void ElementiGrafici::SulSegnale(bool Segnale){
  if(PrimaVolta)
    return;
  DIS_ADD_TYPE(IfDisegna,DIS_TUTTI);
  DIS_REM_TYPE(IfDisegna,DIS_SEGNALE);
  if(Segnale)
    v1->Punta(CoordY);
}
void ElementiGrafici::SeRiscala(){
  if( DIS_IF_TYPE(IfRiscala,RIS_UNO) ){
    DIS_REM_TYPE(IfRiscala,RIS_UNO);
    DIS_ADD_TYPE(IfRiscala,RIS_TUTTI);
  }
  else{
    DIS_ADD_TYPE(IfRiscala,RIS_UNO);
    DIS_REM_TYPE(IfRiscala,RIS_TUTTI);
  }
  repaint();
}
void ElementiGrafici::SeRiscalaTutto(){
  if( DIS_IF_TYPE(IfRiscala,RIS_TUTTI) ){
    DIS_REM_TYPE(IfRiscala,RIS_TUTTI);
    DIS_ADD_TYPE(IfRiscala,RIS_UNO);
  }
  else {
    DIS_ADD_TYPE(IfRiscala,RIS_TUTTI);
    DIS_REM_TYPE(IfRiscala,RIS_UNO);
  }
  repaint();
}
void ElementiGrafici::NSet(){
  if(IfDisegna == DIS_SEGNALE){
    DIS_ADD_TYPE(IfDisegna,DIS_TUTTI);
    DIS_REM_TYPE(IfDisegna,DIS_SEGNALE);
  }
  else{
    DIS_REM_TYPE(IfDisegna,DIS_TUTTI);
    DIS_ADD_TYPE(IfDisegna,DIS_SEGNALE);
  }
  repaint();
}
void ElementiGrafici::SulGrafico(bool Grafico){
  if(PrimaVolta)
    return;
  IfSegnale = 0;
  if(Grafico)
    v1->CambiaPunti();
}
void ElementiGrafici::SuX(bool niente){
  if(niente)
    Coord = 0;
}
void ElementiGrafici::SuY(bool niente){
  if(niente)
    Coord = 1;
}
void ElementiGrafici::SuDX(bool niente){
  if(niente)
    Coord = 2;
}
void ElementiGrafici::SuDY(bool niente){
  if(niente)
    DIS_ADD_TYPE(IfDisegna,DIS_ERR);
  else 
    DIS_REM_TYPE(IfDisegna,DIS_ERR);
}
void ElementiGrafici::InfoStato(const char *testo){
  //  printf("Barra: %s\n",testo);
  QString Info = QString(testo);
  emit Stato(QString(testo));
}
void ElementiGrafici::InfoSequenza(const char *testo){
  QString Info = QString(testo);
  emit StatoSequenza(QString(testo));
}  
QSizePolicy ElementiGrafici::Dimensionamento() const{
  return QSizePolicy( QSizePolicy::Expanding,QSizePolicy::Expanding);
}
void ElementiGrafici::ImpNBin(int n){
  if(n>0){
    NBin = n;
    v1->CambiaNBin(n);
    emit NBinCambiati(n);
  }
}
void ElementiGrafici::ImpNMobile(int n){
  if(n>0){
    NMobile = n;
    emit NMobileCambiato(n);
  }
}
void ElementiGrafici::ImpNCorrela(int n){
  if(n>0){
    NCorrela = n;
    emit NCorrelaCambiato(n);
  }
}
void ElementiGrafici::ImpNVisMin(int n){
  if(n<0)
    return;
  else if(n < NVisMax){
    NVisMin = n;
    if(!VoltaMomenti) ImpNElMin(NVisMin);
    emit NVisMinCambiato(n);
  }
  else if(n>NVisMax){
    if(n < PuntiMax + 10){
      ImpNVisMax(n+10);
      NVisMin = n;
      if(!VoltaMomenti) ImpNElMin(NVisMin);
      emit NVisMinCambiato(n);
    }
    else
      return ;
  }
  if(!VoltaMomenti) ImpNElMin(NVisMin);
  if(PuntiMax < 1000 && !PrimaVolta)
    repaint();
}
void ElementiGrafici::ImpNVisMax(int n){
  if( n > NVisMin +5 && n<=PuntiMax){
    NVisMax = n;
    if(!VoltaMomenti) ImpNElMax(NVisMax);
    emit NVisMaxCambiato(n);
  }
  else if(n < NVisMin + 5 && n <= PuntiMax){
    if(NVisMin > 5){
      NVisMax = n;
      ImpNVisMin( n - 5);
      if(!VoltaMomenti) ImpNElMax(NVisMax);
      emit NVisMaxCambiato(n);
    }
    else 
      return ;
  }
  if(PuntiMax < 1000 && !PrimaVolta)
    repaint();
}
void ElementiGrafici::ImpNElMin(int n){
  if(n<0 )
    return;
  else if(n < NElMax){
    NElMin = n;
    emit NElMinCambiato(n);
    sprintf(cPos,"%.3g : 0.",NElMin*ScalaTopoX/(double)PuntiMax+xMin);
    emit StatoSequenza(QString(cPos));  
  }
  else if(n>NElMax){
    if(n < PuntiMax + 5){
      ImpNElMax(n+5);
      NElMin = n;
      emit NElMinCambiato(n);
      sprintf(cPos,"%.3g : 0.",NElMin*ScalaTopoX/(double)PuntiMax+xMin);
      emit StatoSequenza(QString(cPos));  
    }
    else
      return ;
  }
  if(NMax < 10000 && !PrimaVolta) repaint();
}
void ElementiGrafici::ImpNElMax(int n){
  if( n > NElMin +5 && n<=PuntiMax){
    NElMax = n;
    emit NElMaxCambiato(n);
    sprintf(cPos,"%.3g : 0.",NElMax*ScalaTopoX/(double)PuntiMax+xMin);
    emit StatoSequenza(QString(cPos));  
  }
  else if(n < NElMin + 5 && n <= PuntiMax){
    if(NElMin > 10){
      NElMax = n;
      ImpNElMin( n - 10);
      emit NElMaxCambiato(n);
      sprintf(cPos,"%.3g : 0.",NElMax*ScalaTopoX/(double)PuntiMax+xMin);
      emit StatoSequenza(QString(cPos));  
    }
    else {
      return ;
    }
  }
  if(NMax < 10000 && !PrimaVolta) repaint();
}
void ElementiGrafici::ImpNElMinY(int n){
  if(n<0 || n > RisY)
    return;
  NElMinY = n;
  emit NElMinYCambiato(n);
  sprintf(cPos,"%.3g : %.3g",0.,(height()-NElMinY)*ScalaTopoY/(double)height()+yMin);
  emit StatoSequenza(QString(cPos));  
  if(!PrimaVolta) repaint();
}
void ElementiGrafici::ImpNElMaxY(int n){
  if(n<0 || n > RisY)
    return;
  NElMaxY = n;
  emit NElMaxYCambiato(n);
  sprintf(cPos,"%.3g : %.3g",0.,(height()-NElMaxY)*ScalaTopoY/(double)height()+yMin);
  emit StatoSequenza(QString(cPos));  
  if(!PrimaVolta)   repaint();
}
void ElementiGrafici::NomeFile(const QString &S){
  sprintf(nomeFile,"%s",S.ascii());
  sprintf(stringa,"Carico il file %s",nomeFile);
  InfoStato(stringa);
}
void ElementiGrafici::NomeEntrata(const QString &S){
  sprintf(nomeFile,"%s",S.ascii());
  emit TestoCambiato(S);
}
void ElementiGrafici::NomeConf(const QString &S){
  sprintf(nomeConf,"%s",S.ascii());
  emit ConfCambiato(S);
}
void ElementiGrafici::NomeSalva(const QString &S){
  sprintf(nomeSalva,"%s",S.ascii());
  sprintf(stringa,"Salva il file %s",nomeSalva);
  InfoStato(stringa);
}
void ElementiGrafici::NomeTit(const QString &S){
  nomeTit.replace(0,S.length(),S);
}
void ElementiGrafici::NomeEtX(const QString &S){
  nomeEtX.replace(0,S.length(),S);
}
void ElementiGrafici::NomeEtY(const QString &S){
  nomeEtY.replace(0,S.length(),S);
}
void ElementiGrafici::NomeUscita(const QString &S){
  emit SalvaCambiato(S);
}
void ElementiGrafici::NomeTitolo(const QString &S){
  emit TitoloCambiato(S);
}
void ElementiGrafici::NomeEtichettaX(const QString &S){
  emit EtichettaXCambiato(S);
}
void ElementiGrafici::NomeEtichettaY(const QString &S){
  emit EtichettaYCambiato(S);
}
void ElementiGrafici::Salva(){
  NomeSalva(QString(nomeSalva));
  sprintf(stringa,"Salvo su %s",nomeSalva);
  InfoStato(stringa);
  if(!IfSegnale) v1->ScriviFile(nomeSalva,CoordY,LogLog,NVisMin,NVisMax);
  if( DIS_IF_TYPE(IfDisegna,DIS_PUNTI) )
    v1->ScriviPunti(nomeSalva);
  if( DIS_IF_TYPE(IfDisegna,DIS_MOMENTI) )
    v1->ScriviPunti(nomeSalva);
  else  v1->ScriviTutto(nomeSalva,LogLog,1.,1.);
}
void ElementiGrafici::ImpInter(int n){
  emit InterVisCambiato(n);
}
void ElementiGrafici::ImpInterEl(int n){
  emit InterElCambiato(n);
}
void ElementiGrafici::ImpInter(int Min,int Max){
  emit InterVisCambiato(Min,Max);
}
void ElementiGrafici::ImpInterEl(int Min,int Max){
  emit InterElCambiato(Min,Max);
}
void ElementiGrafici::ImpInterY(int Min,int Max){
  emit InterYCambiato(Min,Max);
}
void ElementiGrafici::ImpCoordX(int n){
  if(n < REF_SEQ || n > NVar) return;
  if(CoordX == n) return;
  CoordX = n;
  if(DIS_IF_TYPE(IfRiscala,RIS_UNO))
    v1->ImpCoordX(CoordY,CoordX);
  else if (DIS_IF_TYPE(IfRiscala,RIS_TUTTI)){
    v1->ImpCoordX(CoordX);
  }
  if(!PrimaVolta) repaint();
  emit CoordXCambiata(MIN(n,CoordX));
}
void ElementiGrafici::ImpCoordY(int n){
  if(n < 0 || n >= NVar) return;
  CoordY = n;
  v1->ImpCoordY(CoordY);
  if(!PrimaVolta) repaint();
  emit CoordYCambiata(CoordY);
}
void ElementiGrafici::ImpCoordDY(int n){
  if(n >= 0 && n < NVar){
    CoordDY = n;
  }
  if(!PrimaVolta) repaint();
  emit CoordYCambiata(CoordY);
}
void ElementiGrafici::ImpNVarY(int n, int m){
  emit InterCoordYCambiato(n,m-1);
}
void ElementiGrafici::ImpNVarX(int n, int m){
  emit InterCoordXCambiato(n,m-1);
}
