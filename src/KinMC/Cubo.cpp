#include "Cubo.h"
Cubo::Cubo(int LatoCubo[3],int Partizioni[3]){
  //  if((Uscite = fopen("Situazione/Uscite.dat","w"))==0)
      printf("Non s'apre Situazione/Uscite.dat\n");
  if((Ricombinate = fopen("Situazione/Ricombinate.dat","w"))==0)
    printf("Non s'apre Situazione/Ricombinate.dat\n");
  for(int i=0;i<3;i++)Lato[i]=LatoCubo[i];
  for(int i=0;i<3;i++){
    nCelle[i] =(int) Lato[i]/Partizioni[i]; 
    nPartizioni[i]=Partizioni[i];
  }
  yMod10=nPartizioni[0];
  zMod10=nPartizioni[1]*nPartizioni[0];
  nMass=zMod10*(nPartizioni[2]+1);
  //  printf("nMass %d\n",nMass);
  Cella = new CELLA[nMass];
}
int Cubo::nCella(const int pos[3]){
  int ris[3];
  int risp=0;
  ris[0]= ((int)pos[0]/nCelle[0]);
  ris[1]= ((int)pos[1]/nCelle[1])*yMod10;
  ris[2]= ((int)pos[2]/nCelle[2])*zMod10;
  risp= ris[0]+ris[1]+ris[2];
  if(comm==1)printf("nCella (%d,%d,%d)->%d\n",pos[0],pos[1],pos[2],risp);
  if(risp>nMass)
    printf("Cella (%d,%d,%d)->%d fuori dal limite %d\n",pos[0],pos[1],pos[2],risp,nMass);
  return risp;
}
int Cubo::PosRet(const int Pos[3]){
  int pos[3]={0};
  int ret[3]={0};
  int risp=0;
  for(int l=0;l<3;l++){
    pos[l]=Pos[l]%nCelle[l];
  }
  for(int l=0;l<3;l++){
    if(pos[l]<=DisRet)ret[l]=-1;
    else if(pos[l]>DisRet && pos[l]<nCelle[l]-DisRet)
      ret[l]=0;
    else if(pos[l]>=nCelle[l]-DisRet)ret[l]=1;
  }
  risp = 9*ret[2]+3*ret[1]+ret[0];
  //  printf("(%d,%d,%d)->%d\n",pos[0],pos[1],pos[2],risp);
  return risp;
}
bool Cubo::Posiziona(Particella &p){
  //  int Posi[Lato[0]][Lato[1]][Lato[2]];
  int pos[3]={0};
  int ret[3]={0};
  int Temp=0;
  float casuale=0;
  OCCUP Occup={0};
  for(int i=0;i<p.nPart;i++){
    if(comm==1)printf("%d)",i);
    for(int l=0;l<2;l++){
      p.Part[i].Pos[l] = (int) (ran1()*(Lato[l]));
    }
    casuale=ran1();
    for(int k=0;k<p.nDis;k++){
      if(p.PrDis[k]<casuale&&casuale<p.PrDis[k+1]){
	p.Part[i].Pos[2]=(int)((k*Lato[2])/p.nDis);
	//	printf("k %d n=(%d,%d,%d) i %d\n",k,p.Part[i].Pos[0],p.Part[i].Pos[1],p.Part[i].Pos[2],i);
	break;
      }
    }
    //    if(Posi[dove[0]][dove[1]][dove[2]]!=1){
    p.Part[i].Tocc=i;
      //    }
      //    else{
      //      i--;
      //printf("Am'o gi:u\n");
      //    }
      //    Posi[dove[0]][dove[1]][dove[2]]=1;
    Temp=nCella(p.Part[i].Pos);
    //    printf("Temp %d\n",Temp);
    Occup.Ret=PosRet(p.Part[i].Pos);
    Occup.Tocc=p.Part[i].Tocc;
    Occup.Tipo=p.Tipo;
    //    Cella.insert(Cella.begin()+dove,Temp);
    Cella[Temp].Part.push_back(Occup);
  }
  return 1;
}
int Cubo::ScegliTipo(const Particella p[],const int tipi){//Da sistemare, va una sola classe
  float casuale=0.;
  float PrCum[tipi+1];
  float N=0;
  int Scelgo=0;
  PrCum[0]=0.;
  for(int h=0;h<tipi;h++){
    PrCum[h+1]=p[h].PrCum+PrCum[h];
    printf("PrCum %f PrCum[%d]=%f\n",p[h].PrCum,h,PrCum[h+1]);
  }
  N=PrCum[tipi+1];
  for(int h=0;h<tipi;h++){
    PrCum[h+1]/=N;
  }
  casuale=ran1();
  for(int h=1;h<tipi+1;h++){
    if(casuale>PrCum[h-1]&&casuale<=PrCum[h])Scelgo=h;
  }
  return Scelgo;
}
int Cubo::ScegliT(const Particella &p1,const Particella &p2){
  float PrCum[3];
  float N=0.;
  int Scelgo=0;
  float casuale=0.;
  PrCum[0]=0;
  PrCum[1]=p1.PrCum;
  PrCum[2]=p1.PrCum+p2.PrCum;
  N=PrCum[2];
  PrCum[1]/=N;
  PrCum[2]/=N;
  casuale=ran1();
  if(PrCum[0]<casuale&&casuale<PrCum[1])
    Scelgo=0;
  else if(PrCum[1]<casuale&&casuale<PrCum[2])
    Scelgo=1;
  if(comm==3)printf("Scelgo n{%d}\n",Scelgo);
  return Scelgo;
}
COORD Cubo::Coord(const int Ret){
  COORD Tre={0};
  Tre.x=(Ret+13+2)%3;
  if(Tre.x==2)Tre.x=-1;
  Tre.y=((int)((Ret+13+6)/3))%3;
  if(Tre.y==2)Tre.y=-1;
  Tre.z=((int)((Ret+13+18)/9))%3;
  if(Tre.z==2)Tre.z=-1;
  //  printf("%d->(%d,%d,%d)\n",Ret,Tre.x,Tre.y,Tre.z);
  return Tre;
}
COORD Cubo::NCella(const int mCella){
  COORD Tre={0};
  Tre.z=((int)(mCella/zMod10));
  Tre.y=((int)((mCella-Tre.z*zMod10)/yMod10));
  Tre.x=(int)(mCella-Tre.y*yMod10-Tre.z*zMod10);
  //  printf("%d->(%d,%d,%d)\n",mCella,Tre.x,Tre.y,Tre.z);
  return Tre;
}
int Cubo::Questa(const int Tipo,const int Tocc,const int dove){
  int questa=0;
  for(int i=0;i<Cella[dove].Part.size();i++){
    //    printf("n{%d}[%d], Cella[%d].Tocc %d\n",Cella[dove].Part[i].Tipo,Tocc,dove,Cella[dove].Part[i].Tocc);
    if(Cella[dove].Part[i].Tocc==Tocc){
      if(Cella[dove].Part[i].Tipo==Tipo){
	questa=i;
	break;
      }
    }
  }
  if(Cella[dove].Part[questa].Tocc!=Tocc)
    printf("La particella %d non `e nella cella %d\n",Tocc,dove);
  return questa;
}
bool Cubo::Sposta(Particella &p1,const int modo,const int Quale){
  int cella1=0;
  int cella2=0;
  int questa=0;
  int Scegli=0;
  OCCUP Occup={0};
  Occup.Tocc=p1.Part[Quale].Tocc;
  Occup.Tipo=p1.Tipo;
  cella1=nCella(p1.Part[Quale].Pos);
  if(cella1<0||cella1>nMass){
    printf("`E stato richiesto n{%d}[%d](%d,%d,%d) in modo %d su %d Part\n",p1.Tipo,p1.Part[Quale].Tocc,p1.Part[Quale].Pos[0],p1.Part[Quale].Pos[1],p1.Part[Quale].Pos[2],modo,p1.nEffet);
    Stato();
    return 0;
  }
  questa= Questa(p1.Tipo,p1.Part[Quale].Tocc,cella1);
  //  printf("Pos (%d,%d,%d)->(%d)\n",p1.Part[Quale].Pos[0],p1.Part[Quale].Pos[1],p1.Part[Quale].Pos[2],cella1);
  //Sposta secondo le condizioni periodiche
  if(modo==0){
    p1.Part[Quale].Pos[0]++;
    if(p1.Part[Quale].Pos[0]>=Lato[0])
      p1.Part[Quale].Pos[0]=0;
  }
  else if(modo==1){
    p1.Part[Quale].Pos[0]--;
    if(p1.Part[Quale].Pos[0]<0)
      p1.Part[Quale].Pos[0]=Lato[0]-1;
  }
  else if(modo==2){
    p1.Part[Quale].Pos[1]++;
    if(p1.Part[Quale].Pos[1]>=Lato[1])
      p1.Part[Quale].Pos[1]=0;
  }
  else if(modo==3){
    p1.Part[Quale].Pos[1]--;
    if(p1.Part[Quale].Pos[1]<0)
      p1.Part[Quale].Pos[1]=Lato[1]-1;
  }
  else if(modo==4){
    p1.Part[Quale].Pos[2]++;
    if(p1.Part[Quale].Pos[2]>=Lato[2])
      p1.Part[Quale].Pos[2]=0;
  }
  else if(modo==5){
    p1.Part[Quale].Pos[2]--;
    if(p1.Part[Quale].Pos[2]<0){
      Scegli=p1.ScegliSup();
      if(Scegli==0){
	p1.Part[Quale].Pos[2]+=2;
	if(comm==1) printf("n{%d}[%d] rimbalza\n",p1.Tipo,p1.Part[Quale].Tocc);
      }
      else if(Scegli==1){
	if(comm==1)printf("n{%d}[%d] si ferma\n",p1.Tipo,p1.Part[Quale].Tocc);
	p1.Part[Quale].Pos[2]++;
      }
      else if(Scegli==2){
	//	printf("n{%d}[%d](%d,%d) contro n{%d}[%d]\n",Cella[cella1].Part[questa].Tipo,Cella[cella1].Part[questa].Tocc,cella1,questa,p1.Tipo,p1.Part[Quale].Tocc);
	Cella[cella1].Part.erase(Cella[cella1].Part.begin()+questa);
 	p1.Elimina(p1.Part[Quale].Tocc);
	//	fprintf(Uscite,"%.0f\t%d\n",t,p1.Tipo);
	if(comm==1||comm==2)printf("n{%d}[%d] esce\n",p1.Tipo,p1.Part[Quale].Tocc);
	return 0;
      }
    }
  }
  cella2=nCella(p1.Part[Quale].Pos);
  if(cella1<0||cella1>nMass){
    printf("`E stato richiesto, su %d Part, n{%d}[%d](%d,%d,%d) in modo %d\n",p1.nEffet,p1.Tipo,p1.Part[Quale].Tocc,p1.Part[Quale].Pos[0],p1.Part[Quale].Pos[1],p1.Part[Quale].Pos[2],modo);
      return 0;
  }
  Occup.Ret=PosRet(p1.Part[Quale].Pos);
  if(cella1!=cella2){
    if(comm==1)printf("n{%d}[%d] Cambia Cella %d->%d\n",Cella[cella1].Part[questa].Tipo,Cella[cella1].Part[questa].Tocc,cella1,cella2);
    Cella[cella1].Part.erase(Cella[cella1].Part.begin()+questa);
    Cella[cella2].Part.push_back(Occup);
  }
  else {
    Cella[cella1].Part[questa].Ret=Occup.Ret;
  }
  return 1;
}
bool Cubo::nmTocc(Particella &p1,Particella &p2,const int dove1,const int dove2,const int Quale){
  int quale2=0;//p2 in particella
  int quella=0;//p2 nella cella
  int questa=0;//p1 nella cella
  int dist=0;//distanza tra le due
  //Contolla che non tocchi con altre 
  questa= Questa(p1.Tipo,p1.Part[Quale].Tocc,dove1);
  if(Cella[dove1].Part[questa].Tocc!=p1.Part[Quale].Tocc)
    return 0;
  else {
    for(int i=0;i<Cella[dove2].Part.size();i++){
      if(Cella[dove2].Part[i].Tipo==p2.Tipo){
	for(int j=Min(Cella[dove2].Part[i].Tocc,p2.nEffet);j>=0;j--){
	  if(p2.Part[j].Tocc==Cella[dove2].Part[i].Tocc){
	  quale2=j;
	  quella=i;
	  break;
	  }
	}
	dist=0;
	for(int l=0;l<3;l++){
	  dist+=quad(p1.Part[Quale].Pos[l]-p2.Part[quale2].Pos[l]);
	}
	if(dist<=quad(DisRet)){
	  if(comm==2 || comm==1)
	    printf("n{%d}[%d](%d,%d) su %d e n{%d}[%d](%d,%d) su %d toccano dist %d\n",p1.Tipo,p1.Part[Quale].Tocc,dove1,questa,p1.nEffet,p2.Tipo,p2.Part[quale2].Tocc,dove2,quella,p2.nEffet,dist);
	  // printf("n{%d}[%d](%d)-n{%d}[%d](%d)\n",Cella[dove1].Part[questa].Tipo,Cella[dove1].Part[questa].Tocc,questa,Cella[dove2].Part[quella].Tipo,Cella[dove2].Part[quella].Tocc,quella);
	  Cella[dove1].Part.erase(Cella[dove1].Part.begin()+questa);
	  if(dove1==dove2){
	    if(questa<quella){
	      quella--;
	    }
	  }
	  //printf("%d)%d\n",quella,Cella[dove2].Part[quella].Tocc);
	  Cella[dove2].Part.erase(Cella[dove2].Part.begin()+quella);
	  p1.Elimina(p1.Part[Quale].Tocc);
	  p2.Elimina(p2.Part[quale2].Tocc);
	  //fprintf(Ricombinate,"%.0f\n",t);
	  return 1;
	}
      }
    }
  }
  return 0;
}
bool Cubo::Controlla(Particella &p1,const int Quale,Particella &p2){
  int dove1=0;
  int questa=0;
  int dove=0;
  COORD TreRet={0};
  COORD TreCel={0};
  dove=nCella(p1.Part[Quale].Pos);
  if(dove <0||dove>nMass){
    printf("`E stato richiesto n{%d}[%d](%d,%d,%d) (quale %d) su %d Part\n",p1.Tipo,p1.Part[Quale].Tocc,p1.Part[Quale].Pos[0],p1.Part[Quale].Pos[1],p1.Part[Quale].Pos[2],Quale,p1.nEffet);
    return 0;
  }
  //Controlla che non tocchi altre nella cella
  if(nmTocc(p1,p2,dove,dove,Quale)==1){
    return 1;
  }
  //Controlla che non sia troppo vicino ad altre celle
  TreCel = NCella(dove);
  TreRet = Coord(Cella[dove].Part[questa].Ret);
  if(TreRet.x==0);
  else if(TreRet.x==1){
    if(TreCel.x<nCelle[0]){
      dove1=dove+1;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;
    }
  }
  else if(TreRet.x==-1){
    if(TreCel.x>0){
      dove1=dove-1;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;
    }
  }    
  if(TreRet.y==0);
  else if(TreRet.y==1){
    if(TreCel.y<nCelle[1]){
      dove1=dove+yMod10;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;
    }
  }
  else if(TreRet.y==-1){
    if(TreCel.y>0){
      dove1=dove-yMod10;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;      
    }
  }    
  if(TreRet.z==0);
  else if(TreRet.z==1){
    if(TreCel.z<nCelle[2]){
      dove1=dove+zMod10;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;
    }
  }
  else if(TreRet.z==-1){
    if(TreCel.z>0){
      dove1=dove-zMod10;
      if(nmTocc(p1,p2,dove,dove1,Quale)==1)
	return 1;
    }
  }
  return 1;
}
void Cubo::Stato(){
  for(int i=0;i<nMass;i++){
    for(int j=0;j<Cella[i].Part.size();j++){
      if(Cella[i].Part.size()>0)
	printf("%d)Part{%d}[%d](%d|%d)\n",j,Cella[i].Part[j].Tipo,Cella[i].Part[j].Tocc,i,Cella[i].Part[j].Ret);
    }
  }
}
void Cubo::Stato(int nCella){
  for(int i=0;i<Cella[nCella].Part.size();i++){
    printf("Part{%d}[%d](%d)\n",Cella[nCella].Part[i].Tipo,Cella[nCella].Part[i].Tocc,nCella);
  }
}
