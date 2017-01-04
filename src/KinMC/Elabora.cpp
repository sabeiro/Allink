#include "Elabora.h"

int main(int argc,char **argv){
  const int arc=argc-1;
  time_t ti,tf;
  (void) time(&ti);
  printf("------------------------Legge------------------------------\n\
*        Questo programma richiede in entrata un file      *\n\
*                formattato nel seguente modo:             *\n\
* %%d\t%%d\t%%d                                         *\n\
* Num\tt[ps]\ttipo                                      *\n\
* %%d\t%%d\t%%d\t%%d                                 *\n\
* tocc\tx\ty\tz                                  *\n\
* Richiede come argomento almeno un file                   *\n\
------------------------------------------------------------\n");
  FILE *Pos;
  char *nome,*nm,*tipo[arc];
  int *dec;
  int *tocc[arc],Tipo;
  double *pos0,*pos1,*pos2;
  float *xp[arc],*yp[arc],*zp[arc],*t,*tp;
  int Num[arc];
  dec  = (int *)malloc(sizeof(int));
  pos0 = (double *)malloc(sizeof(double));
  pos1 = (double *)malloc(sizeof(double));
  pos2 = (double *)malloc(sizeof(double));  
  tp =   (float *)malloc(sizeof(float));
  nome = (char *)malloc(30*sizeof(char));
  nm = (char *)malloc(2*sizeof(char));
  t  = (float *)malloc(arc*sizeof(int));
  sprintf(nm,"nm");
  if(argc<2)printf("Mi serve come argomento un file!\n");  
  for(int i=0;i<arc;i++){
    sprintf(nome,"%s",argv[i+1]);
    //    printf("Apro %s su %d\n",nome,arc);
    if((Pos=fopen(nome,"r"))==0)printf("Non s'apre %s\n",nome);
    fscanf(Pos,"%d",dec);
    Num[i]=*dec;
    fscanf(Pos,"%f",tp);
    t[i]=*tp;
    //    printf("Num[%d]=%d t=%.1f\n",i,Num[i],t[i]);
    tipo[i] = (char *)malloc(sizeof(char));
    fscanf(Pos,"%d",&Tipo);
    if(Tipo==0)
      sprintf(tipo[i],"n");
    else if(Tipo==1)
      sprintf(tipo[i],"m");
    //    printf("Tipo %d ->%s\n",Tipo,tipo[i]);
    *(xp+i) = (float *)malloc(Num[i]*sizeof(float));
    *(yp+i) = (float *)malloc(Num[i]*sizeof(float));
    *(zp+i) = (float *)malloc(Num[i]*sizeof(float));
    *(tocc+i) = (int *)malloc(Num[i]*sizeof(int));
    for(int j=0;fscanf(Pos,"%d",dec)==1&&dec!=0;j++){
      fscanf(Pos,"%lf",pos0);
      xp[i][j]=*pos0;
      fscanf(Pos,"%lf",pos1);
      yp[i][j]=*pos1;
      fscanf(Pos,"%lf",pos2);
      zp[i][j]=*pos2;
      tocc[i][j]=*dec;
      //    printf("toccc=%d pos0=%.0f pos1=%.0f pos2=%.0f\n",*dec,*pos0,*pos1,*pos2);
      //      printf("NUM[%d]=%d %c[%d]=(%.0f,%.0f,%.0f)\n",i,Num[i],tipo[i][0],tocc[i][j],xp[i][j],yp[i][j],xp[i][j]);
    }
  }
  printf("Elabora...\n");
  printf("  1) Particelle rimaste al tempo t\n");
  FILE *nPos,*mPos;
  system("rm -f Derivati/*.dat");
  if((nPos= fopen("Derivati/nDiff.dat","w+"))==0)
    printf("Non si apre Derivati/nDiff.dat");
  if((mPos= fopen("Derivati/mDiff.dat","w+"))==0)
    printf("Non si apre Derivati/mDiff.dat");
  float quadr[arc];
  int ini=0,imi=0;
  for(int i=0;i<arc;i++){
    if(*tipo[i]==*nm){
      if(t[i]==0.0){
	ini=i;
      }
    }
    else if(*tipo[i]==*(nm+1)){
      if(t[i]==0.0){
	imi=i;
      }
    }    
  }
  for(int i=0;i<arc;i++){
    quadr[i]=0;
    if(*tipo[i]==*nm){
      if(i!=ini){
	for(int k=0;k<Num[ini];k++){
	  for(int j=0;j<Num[i];j++){
	    if(tocc[i][j]==tocc[ini][k]){
	      quadr[i]+=pow((xp[ini][k]-xp[i][j]),2);
	      quadr[i]+=pow((yp[ini][k]-yp[i][j]),2);
	      quadr[i]+=pow((zp[ini][k]-zp[i][j]),2);
	      //printf("tocc[%d][%d]=%d tocc[%d][%d] =%d\n",ini,k,tocc[ini][k],i,j,tocc[i][j]);
	    }
	  }
	}
	//	printf("ini=%d Num[%d]=%d \t\t",ini,i,Num[i]);
	//	printf("quadr[%d]=%.0f\n",t[i],quadr[i]);
	fprintf(nPos,"%.1f\t%g\t%d\n",t[i],quadr[i]/(6*(t[i]-t[ini])*10000)/Num[i],Num[i]);
      }
    }
    else if(*tipo[i]==*(nm+1)){
      if(i!=imi){
	for(int k=0;k<Num[imi];k++){
	  for(int j=0;j<Num[i];j++){
	    if(tocc[i][j]==tocc[imi][k]){
	      quadr[i]+=quad((xp[imi][k]-xp[i][j]));
	      quadr[i]+=quad((yp[imi][k]-yp[i][j]));
	      quadr[i]+=quad((zp[imi][k]-zp[i][j]));
	    }
	  }
	}
	//printf("quadr[%d]=%.0f\n",t[i],quadr[i]);
	fprintf(mPos,"%.1f\t%g\t%d\n",t[i],quadr[i]/(6*(t[i]-t[imi])*10000)/Num[i],Num[i]);
      }
    }
    //printf("NUM[%d]=%d %c[%d]=(%.0f,%.0f,%.0f)\n",i,Num[i],tipo[i][0],*(tocc[i]+j),*(xp[i]+j),*(yp[i]+j),*(zp[i]+j));
  }
  fclose(nPos);
  fclose(mPos);
  printf("  2) Profondit`a media al tempo t\n");
  system("rm -f Fetta/*.dat");
  float zMedio[arc];
  FILE *nMedia,*mMedia,*nRim,*mRim;
  if((nMedia=fopen("Derivati/nMedia.dat","w"))==0)
    printf("Non si apre Derivati/nMedia.dat");
  if((mMedia=fopen("Derivati/mMedia.dat","w"))==0)
    printf("Non si apre Derivati/mMedia.dat");
  if((nRim  =fopen("Derivati/nRim.dat","w"))==0)
    printf("Non si apre Derivati/nRim.dat");
  if((mRim  =fopen("Derivati/mRim.dat","w"))==0)
    printf("Non si apre Derivati/mRim.dat");
  for(int i=0;i<arc;i++){
    zMedio[i]=0;
    for(int m=0;m<116;m++)Bin[m]=0;
    for(int j=0;j<zLim;j++)Fetta[j]=0;
    if(*tipo[i]==*nm){
      fprintf(nRim,"%.1f\t%d\n",t[i],Num[i]);
      for(int j=0;j<Num[i];j++){
 	zMedio[i]+=zp[i][j];
	*dec=(int)zp[i][j];
	if(*dec<zLim&&*dec>0){
	  Fetta[*dec]++;
	}
	else if(*dec>zLim)
	  printf("n Particella t=%.1f,%d oltre %d\n",t[i],tocc[i][j],zLim);
	else if(*dec<0)
	  printf("n Particella t=%.1f,%d oltre %d\n",t[i],tocc[i][j],0);
      }
      fprintf(nMedia,"%.1f\t%.0f\n",t[i],zMedio[i]/Num[i]);
      sprintf(nome,"Fetta/nFetta_%.1f.dat",t[i]);      
      if((Pos=fopen(nome,"w"))==0)
	printf("Non si apre %s",nome);
      for(int j=0,m=0;j<zLim;j++){
	Bin[m]+=Fetta[j];
	if(j%35==0)m++;
      }
      for(int m=0;m<116;m++){
 	fprintf(Pos,"%d\t%d\n",m,Bin[m]);
      }
      fclose(Pos);
    }

    if(*tipo[i]==*(nm+1)){
      fprintf(mRim,"%.1f\t%d\n",t[i],Num[i]);
      for(int j=0;j<Num[i];j++){
	zMedio[i]+=zp[i][j];
	*dec=(int)zp[i][j];
	if(*dec<zLim&&*dec>0){
	  Fetta[*dec]++;
	}
	else if(*dec>zLim)
	  printf("m Particella t=%.1f,%d oltre %d\n",t[i],tocc[i][j],zLim);
	else if(*dec<0)
	  printf("m Particella t=%.1f,%d oltre %d\n",t[i],tocc[i][j],0);
      }
      fprintf(mMedia,"%.1f\t%.0f\n",t[i],zMedio[i]/Num[i]);
      sprintf(nome,"Fetta/mFetta_%.1f.dat",t[i]);      
      if((Pos=fopen(nome,"w"))==0)
	printf("Non si apre %s",nome);
      for(int j=0,m=0;j<zLim;j++){
      	Bin[m]+=Fetta[j];
      	if(j%35==0)m++;
      }
      for(int m=0;m<116;m++){
	fprintf(Pos,"%d\t%d\n",m,Bin[m]);
      }
      fclose(Pos);
    }
  }
  fclose(nMedia);
  fclose(mMedia);
  fclose(nRim);
  fclose(mRim);
  FILE *Ricombinate,*nRico,*mRico;
  int iBin=1000;
  int DBin=1000;
  if(1==0){
    printf("  3) Ricombinate per tempo fisso\n");
    //nel file dev'esserci solo il momento di ricomb
    if((Ricombinate=fopen("Situazione/Ricombinate.dat","r"))=0)
      printf("Non si apre Situazione/Ricombinate.dat");
    if((nRico=fopen("Derivati/Ricombinate.dat","w"))==0)
      printf("Non si apre Derivati/Ricombinate.dat");
    int Ric=0;
    for(int i=0;fscanf(Ricombinate,"%d\n",dec)==1;i++){
      Ric=i;
      if(*dec>iBin){
	fprintf(nRico,"%d\t%d\n",iBin,Ric);
	Ric=0;
	iBin+=DBin;
	i=0;
      }
      //    printf("%d\t%d\t%d\n",iBin,*Ric,*dec);
    }
    fclose(Ricombinate);
    fclose(nRico);
    //free(Ric);
  }
  int nRic=0,mRic=0;
  int quale;
  FILE *Uscite;
  if(1==0){
    printf("  4) Uscite\n");
    if((Uscite=fopen("Situazione/Uscite.dat","r"))==0)
      printf("Non si apre Derivati/Uscite.dat");
    if((nRico=fopen("Derivati/nUscite.dat","w"))==0)
      printf("Non si apre Derivati/nUscite.dat");
    if((mRico=fopen("Derivati/mUscite.dat","w"))==0)
      printf("Non si apre Derivati/nUscite.dat");
    iBin=1000;
    for(int i=0,j=0;fscanf(Uscite,"%d\n",dec)==1;){
      fscanf(Uscite,"%d",&quale);
      if(quale==0){
	nRic=i;
	i++;
      }
      else if(quale==1){
	mRic=j;
	j++;
      }
      if(*dec>iBin){
	fprintf(nRico,"%d\t%d\n",iBin,nRic);
	fprintf(mRico,"%d\t%d\n",iBin,mRic);
	nRic=0;
	mRic=0;
	iBin+=DBin;
	i=0;
	j=0;
      }
      //printf("Bin=%d nRic=%d mRic=%d dec=%d quale=%c\n",iBin,*nRic,*mRic,*dec,*quale);
    }
    fclose(Uscite);
    fclose(nRico);
    fclose(mRico);
  }
  (void) time(&tf);
  printf("Il tutto dopo %d s\n",(int)tf-ti);
  printf("Te ve be te e?\n");
  
  return 0;
}
