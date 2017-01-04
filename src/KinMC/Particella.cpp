#include "Particella.h"
Particella::Particella(char *cDis,char *cVal,int NPart){
  float N=0;
  nPart = NPart;
  nEffet = nPart;
  Part = new PART1[nPart];
  Pr0=0;Pr=0;E0=0;nDis=0;Norm=0;
  FDis = fopen(cDis,"r");
  if(FDis==NULL)printf("Non s'apre %s\n",cDis);
  Norm=0;
  PrDis.insert(PrDis.begin(),3);
  for(int i=1;fscanf(FDis,"%f",&Dis)==1;i++,nDis++){
    PrDis.push_back(Dis + PrDis[i-1]);
  }
  Norm = PrDis[nDis];
  for(int i=0;i<=nDis;i++){
    PrDis[i]/=Norm;
    //    printf("PrDis[%d]=%f nDis=%d\n",i,PrDis[i],nDis);
  }
  fclose(FDis);
  FVal=fopen(cVal,"r");
  if(FVal==NULL)printf("Non si apre %s\n",cVal);
  PrMoto[0]=0;
  PrSup[0]=0;
  for(int i=0;fscanf(FVal,"%f",&Dis)==1;i++){
    if(i<6){
      PrMoto[i+1] = Dis + PrMoto[i];
    }
    else if(i>=6 && i<9){
      PrSup[i+1-6] = Dis + PrSup[i-6];
    }
    else if(i==9)Pr0=Dis;
    else if(i==10)E0=Dis;
    else if(i==11)Tipo=(int)Dis;
  }
  N=PrMoto[6];
  for(int i=1;i<7;i++)PrMoto[i]/=N;
  N=PrSup[3];
  for(int i=1;i<4;i++)PrSup[i]/=N;
  T=300.;
  Pr=Pr0*exp(-E0/.00861/T);
  PrCum=Pr*nPart;
  printf("nPart=%d Pr0=%.3f E0=%.3f Pr=%.3g T=%.0f\n",nPart,Pr0,E0,Pr,T);
}
void Particella::CambiaTemperatura(float DeltaT){
  T += DeltaT;
  Pr=Pr0*exp(-E0/.00861/T);
}
int Particella::ScegliP(){
  float casuale=0.;
  int risp=0;
  if(nEffet==0){
    printf("Non ci sono pi`u particelle del tipo %d\n",Tipo);
    return -1;
  }
  else {
    casuale = ran1()*(nEffet);
    risp=(int)casuale;
    if(comm==1)printf("Al tempo %.3f scelgo n{%d}[%d] ",t,Tipo,Part[risp].Tocc);
    Deltat = -log(ran1())/PrCum;
    t += Deltat;
    if(risp>=nEffet){
      printf("La scelta %d `e oltre le possibilit`a %d\n",risp,nEffet);
      return -1;
    }
    return risp;
  }
}
void Particella::AggiornaTempo(){
  Deltat = -log(ran1())/PrCum;
  t += Deltat;
}
void Particella::Elimina(int Tocc){
  if(nEffet<2){
    nEffet=0;
    delete Part;
  }
  else{
    PART1 *p1 = new PART[nEffet];
    for(int i=0,j=0;i<nEffet;i++){
      if(Part[i].Tocc==Tocc){
	j=1;
	//      printf("Hop!\n");
      }
      p1[i].Tocc=Part[i+j].Tocc;
      for(int l=0;l<3;l++){
	p1[i].Pos[l]=Part[i+j].Pos[l];
      }
    }
    delete Part;
    nEffet--;
    Part = new PART1[nEffet];
    for(int i=0;i<nEffet;i++){
      Part[i].Tocc=p1[i].Tocc;
    //    printf("%d) Part[%d]",i,Part[i].Tocc);
      for(int l=0;l<3;l++){
	Part[i].Pos[l]=p1[i].Pos[l];
	//      printf("%d,",Part[i].Pos[l]);
      }
      //    printf("\b)\n");
    }
    delete p1;
  }
  PrCum=Pr*nEffet;
}
// void Particella::Elimina(int nTocc){
//   nEffet--;
//   for(int i=0,j=0;i<nEffet;i++){
//     if(Part[i].Tocc==nTocc)j++;
//     Part[i].Tocc=Part[i+j].Tocc;
//     for(int l=0;l<3;l++){
//       Part[i].Pos[l]=Part[i+j].Pos[l];
//     }
//   }
// }
int Particella::ScegliModo(){
  int modo=0;
  float casuale=0;
  casuale = ran1();
  for(int i=0;i<6;i++){
    if(casuale>PrMoto[i]&&casuale<=PrMoto[i+1]){
      modo=i;
    }
  }
  if(comm==1)printf("in modo %d\n",modo);
  return modo;
}
int Particella::ScegliSup(){
  float casuale=0.;
  int risp=0;
  casuale = ran1();
  for(int i=0;i<3;i++){
    if(casuale>PrSup[i]&&casuale<PrSup[i+1])
      risp=i;
  }
  //  printf("Scelgo %d\n",risp);
  return risp;
}
void Particella::Stato(){
  for(int i=0;i<nEffet;i++){
    printf("n{%d}[%d](",Tipo,Part[i].Tocc);
    for(int l=0;l<3;l++){
      printf("%d,",Part[i].Pos[l]);
    }
    printf("\b)\n");
  }
}
void Particella::Stato(int Tocc){
  for(int i=(Tocc>nEffet?Tocc:nEffet);i>0;i--){
    if(Part[i].Tocc==Tocc){
      printf("n{%d}[%d](",Tipo,Part[i].Tocc);
      for(int l=0;l<3;l++){
	printf("%d,",Part[i].Pos[l]);
      }
      printf("\b)\n");
    }
  }
}    
void Particella::fStato(int tempo){
   char *file;
   file = new char[30];
   sprintf(file,"Posizioni/Pos_t%d.dat",tempo);
   if((Posizioni = fopen(file,"w"))!=NULL){
     fprintf(Posizioni,"%d\n",nEffet);
     for(int i=0;i<nEffet;i++){
       fprintf(Posizioni,"%d\t%d\t%d\t%d\t%d\n",Part[i].Tocc,Tipo,Part[i].Pos[0],Part[i].Pos[1],Part[i].Pos[2]);
     }
     fclose(Posizioni);
   }
   else printf("Non s'apre %s\n",file);
}
int Particella::MTocc(int i){
  if(i<= nEffet)
    return Part[i].Tocc;
}
int Particella::MPos(int i,int l){
  if(l<3 && i<= nEffet)
    return Part[i].Pos[l];
}
int Particella::MTipo(void){
  return Tipo;
}
int Particella::MNEffet(void){
  return nEffet;
}
