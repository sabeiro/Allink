#include "DrImage.h"

#ifndef USE_PNG
int main(int agrc,char **argv){
  printf("No pngwriter libs used.\nExit.\n");
}
#else
int main(int argc, char **argv){
  DrImage *DrI = new DrImage(argc,argv); 
  for(int i=1;i<argc;i++){
    if(!strncmp(argv[i],"-ga",3)){
      DrI->Gauss();
    }
    else if(!strncmp(argv[i],"-c",2)){
      DrI->Cut();
    }
    else if(!strncmp(argv[i],"-g",2)){
      DrI->Gravity();
    }
    else if(!strncmp(argv[i],"-f",2)){
      DrI->Fourier();
    }
    else if(!strncmp(argv[i],"-i",2)){
      DrI->Ising();
    }
    else if(!strncmp(argv[i],"-l",2)){
      DrI->LennardJones();
    }
    else if(!strncmp(argv[i],"-m",2)){
      DrI->ConvMatrix();
    }
    else if(!strncmp(argv[i],"-s",2)){
      DrI->Shear();
    }
    else if(!strncmp(argv[i],"-r",2)){
      DrI->Rotor();
    }
    else if(!strncmp(argv[i],"-d",2)){
      DrI->Difference();
    }
    else if(!strncmp(argv[i],"-sl",2)){
      DrI->SlitScan();
    }
    else if(!strncmp(argv[i],"-y",2)){
    }
  }
  return 0;
}
DrImage::DrImage(int argc, char *argv){
  NFile = 0;
  NWidth = 0;
  NHeight = 0;
  for(int c=0;c<argc;c++){
    if(!strncmp(argv[i],"-",1)){
      continue;
    }

    NFile++;
  }
  FileList = (char **)calloc(NFile,sizeof(char));
  for(int c=0,f=0;c<argc;c++){
    if(!strncmp(argv[i],"-",1)){
      continue;
    }
    FileList[f] = (char *)calloc(120,sizeof(char));
    sprintf(FileList[f],argv[c]);
    f++;
  }
  Load(FileList[0]);
  Load2(FileList[NFile-1]);
}
DrImage::~DrImage(){
  for(int l=0;l<3;l++){
    free(data[l]);
    free(data1[l]);
  free(data);
  free(data1);
  for(int f=0;f<NFile;f++){
    free(FileList[f]);
  }
  free(FileList);
}
void DrImage::Load(char *FName){
  printf("Load %s\n",FName);
  ImageIn = new pngwriter(0,0,1.0,FName);
  ImageIn->readfromfile(FName);
  NWidth = ImageIn->getwidth();
  NHeight = ImageIn->getheight();
  if(NWidth*NHeight == 0){
    printf("The image is empty\n");
    exit(0);
  }
  data = (double **)calloc(3,sizeof(double));
  for(int l=0;l<3;l++){
    data[l] = (double *)calloc(NWidth*NHeight,sizeof(double));
  }
  for(int l=0;l<3;l++){
    for(int w=0;w<NWidth;w++){
      for(int h=0;h<NHeight;h++){
	data[l][h*NWidth + w] = ImageIn->dread(w,h,l+1);
      }
    }
  }
  //ImageIn.close();
}
void DrImage::ReLoad(char *FName){
  ImageIn->readfromfile(FName);
  for(int l=0;l<3;l++){
    for(int w=0;w<NWidth;w++){
      for(int h=0;h<NHeight;h++){
	data[l][h*NWidth + w] = ImageIn->dread(w,h,l+1);
      }
    }
  }
}
void DrImage::Load2(char *FName){
  printf("Load %s\n",FName);
  ImageIn->readfromfile(FName);
  data1 = (double **)calloc(3,sizeof(double));
  for(int l=0;l<3;l++){
    data1[l] = (double *)calloc(NWidth*NHeight,sizeof(double));
  }
  for(int l=0;l<3;l++){
    for(int w=0;w<NWidth;w++){
      for(int h=0;h<NHeight;h++){
	data1[l][h*NWidth + w] = ImageIn->dread(w,h,l+1);
      }
    }
  }
  //ImageIn.close();
}
void DrImage::Gravity(){
  int NStep = (int)(25.*60.*.2);
  double InvNStep = 1./(double)NStep;
  double InvNHei = 1./(double)NHeight;
  double Acc = -1.;
  double Vel = fabs(Acc)*NStep*.75;
  int Pos = 0.;
  for(int s=0;s<NStep;s++){
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      int h1 = h + Pos;
      if(h1 >= NHeight) h1 -= NHeight;
      if(h1 < 0) h1 += NHeight;
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h1,data[0][h*NWidth + w],data[1][h*NWidth + w],data[2][h*NWidth + w]);
      }
    }
    Vel += Acc;
    if(Vel > NHeight) Vel = NHeight;
    Pos += (int)Vel;
    Pos -= (int)(floor(Pos*InvNHei)*NHeight);
    ImageOut.close();
  }
  printf("\n");
}
void DrImage::Ising(){
  int NStill = 50;
  int NStep = (int)(25.*60.*.2) - NStill;
  int NChange = 60;
  int NSquare = 4;
  double InvNStep = 1./(double)NStep;
  Matematica *Mate = new Matematica();
  for(int s=NStill+NStep;s>=NStep;s--){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  for(int s=NStep;s>=0;s--){
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int c=0;c<NChange;c++){
      int w = (int)(NWidth*Mate->Casuale());
      int h = (int)(NHeight*Mate->Casuale());
      w = w - w%NSquare;
      h = h - h%NSquare;
      if(Mate->Casuale() < .5){
	for(int ws=0;ws<NSquare;ws++){
	  if(ws+w > NWidth) continue;
	  for(int hs=0;hs<NSquare;hs++){
	    if(hs+h > NHeight) continue;
	    for(int l=0;l<3;l++){
	      data[l][(h+hs)*NWidth+(w+ws)] = 0.;
	    }
	  }
	}
      }
      else {
	for(int ws=0;ws<NSquare;ws++){
	  if(ws+w > NWidth) continue;
	  for(int hs=0;hs<NSquare;hs++){
	    if(hs+h > NHeight) continue;
	    for(int l=0;l<3;l++){
	      data[l][(h+hs)*NWidth+(w+ws)] = 1.;
	    }
	  }
	}
      }
    }
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  printf("\n");
}
void DrImage::ConvMatrix(){
  double Average = 0.;
  double Count = 0.;
  Matrice ImIn(NWidth,NHeight);
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      double Sum = 1.- .3333*(data[0][h*NWidth+w]+data[1][h*NWidth+w]+data[2][h*NWidth+w]);
      //double Sum = 1. - (data[0][h*NWidth+w]+data[1][h*NWidth+w]+data[2][h*NWidth+w]);
      ImIn.Set(w,h,Sum);
      Average += Sum;
      Count += 1.;
    }
  }
  Average /= Count;
  //---------Canny---------------------
  Matematica *Mat = new Matematica;
  Matrice Canny(5,5);
  Canny.FillCanny();
  Canny.Print();
  // Mat->ApplyFilter(&ImIn,&Canny);
  // Mat->ApplyFilter(&ImIn,&Canny);
  //-----------Edge-------------------
  SPLINE Weight;
  Weight.a0 = 0.;
  Weight.a1 = 1.; Weight.a2 = 0.;
  Weight.a3 = 0.; Weight.a4 = 0.;
  Matrice *Mask = new Matrice(Weight,3);
  Mask->Print();
  //Mat->ApplyFilter(&ImIn,Mask);
  // Mask->Transpose();
  // Mat->ApplyFilter(ImIn,Mask);
  //-------------Smooth----------------
  const int NMatSize = 5;
  Matrice GaussEdge(NMatSize,NMatSize);
  GaussEdge.FillGaussian(.5,1.);
  //double LatPos[5] = {.125,-.25,0.,.25,-.125};
  // double LatPos[3] = {-1.,0.,1.};
  // for(int w=0;w<NMatSize;w++){
  //   for(int h=0;h<NMatSize;h++){
  //     GaussEdge.Set(w,h,GaussEdge.Val(w,h)*LatPos[w]);
  //   }
  // }
  GaussEdge.Print();
  //Mat->ApplyFilter(ImIn,&GaussEdge);
  Matrice GaussSmooth(5,5);
  GaussSmooth.FillGaussian(.5,3.);
  // Mat->ApplyFilter(ImIn,&GaussSmooth);
  //------------PixDev------------------
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      ImIn.Set(w,h,Average-ImIn.Val(w,h));
    }
  }
  int PixelDev = 5;
  double ValStep[3] = {0.,0.,0.};
  int ValNStep[3] = {0,0,0};
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      ValNStep[0] = w - PixelDev;
      if(ValNStep[0] < 0 ) continue;
      if(ValNStep[0] >= NWidth) continue;
      ValNStep[1] = w;
      ValNStep[2] = w + PixelDev;
      if(ValNStep[2] < 0 ) continue;
      if(ValNStep[2] >= NWidth) continue;
      for(int d=0;d<3;d++){
	ValStep[d] = ImIn.Val(ValNStep[d],h);
 	if(d == 1){
	  ValStep[d] = ValStep[d] > 0. ? ValStep[d] : 0.;
	  continue;
	}
	ValStep[d] = ValStep[d] < 0. ? -ValStep[d] : 0.;
      }
      double Resp = ValStep[0]*ValStep[1]*ValStep[2];
      //double Resp = ValStep[1];
      //ImIn.Set(w,h,Resp);
    } 
  }
  char cImage[160];
  sprintf(cImage,"Canny.png");
  pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
  FILE *Ciccia = fopen("Pos3d.dat","w");
  double NormH = 1./(double)NHeight;
  double NormW = 1./(double)NWidth;
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      fprintf(Ciccia,"%lf %lf %lf\n",w*NormH,h*NormH,ImIn.Val(w,h));
      ImageOut.plot(w,h,ImIn.Val(w,h),ImIn.Val(w,h),ImIn.Val(w,h));
    }
  }
  fclose(Ciccia);
  ImageOut.close();
}
void DrImage::Cut(){
  char cImage[160];
  sprintf(cImage,"Small.png");
  double pixel[4];
  double xBound[2] = {.5,.75};
  double yBound[2] = {.2,.45};
  int xNBound[2];
  int yNBound[2];
  for(int d=0;d<2;d++){
    xNBound[d] = (int)(xBound[d]*NWidth);
    yNBound[d] = (int)(yBound[d]*NHeight);
    printf("%d %d %d %d\n",xNBound[d],NWidth,yNBound[d],NHeight);
  }
  int Dx = xNBound[1]-xNBound[0];
  int Dy = yNBound[1]-yNBound[0];
  double *Plot = new double[Dx*Dy];
  Matematica *Mat = new Matematica;
  //VarData *Var = new VarData;
  Matrice *ImIn = new Matrice(Dx,Dy);
  pngwriter ImageOut(xNBound[1]-xNBound[0],yNBound[1]-yNBound[0],1.0,cImage);
  double Average = 0.;
  double Count = 0.;
  for(int h=yNBound[0];h<yNBound[1];h++){
    for(int w=xNBound[0];w<xNBound[1];w++){
      int hh = h -yNBound[0];
      int ww = w -xNBound[0];
      double Sum = data[0][h*NWidth+w]+data[1][h*NWidth+w]+data[2][h*NWidth+w];
      //Sum *= .5;
      // ImageOut.plot(ww,hh,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      ImIn->Set(ww,hh,Sum);
      Plot[hh*Dx+ww] = Sum;
    }
  }


  for(int h=0;h<Dy;h++){
    for(int w=0;w<Dx;w++){
      //if(ImIn->Val(w,h) >= 0.)
      ImageOut.plot(w,h,10.*ImIn->Val(w,h),ImIn->Val(w,h),ImIn->Val(w,h));
    }
  }
  ImageOut.close();
  delete [] Plot;
}
void DrImage::Difference(){
  char cImage[160];
  sprintf(cImage,"Difference.png");
  pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
  double pixel[4];
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      for(int l=0;l<3;l++){
	pixel[l] = data[l][h*NWidth+w] - data1[l][h*NWidth+w];
      }
      ImageOut.plot(w,h,pixel[0],pixel[1],pixel[2]);
    }
  }
  ImageOut.close();
  printf("\n");
  return;
  printf("Difference\n");
  //char cImage[60];
  double **data2 = (double **)calloc(3,sizeof(double));
  for(int l=0;l<3;l++){
    data2[l] = (double *)calloc(NHeight*NWidth,sizeof(double));
  }
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      for(int l=0;l<4;l++){
  	data2[l][h*NWidth+w] = data[l][h*NWidth+w];// -  data1[l][h*NWidth+w];
      }
    }
  }
  sprintf(cImage,"Difference.png");
  //pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
  for(int h=0;h<NHeight;h++){
    for(int w=0;w<NWidth;w++){
      ImageOut.plot(w,h,data2[0][h*NWidth+w],data2[1][h*NWidth+w],data2[2][h*NWidth+w]);
    }
  }
  ImageOut.close();
}
void DrImage::Shear(){
  int NStill = 50;
  int NStep = (int)(25.*60.*.2) - NStill;
  double InvNStep = 1./(double)NStep;
  double InvNHei = 1./(double)NHeight;
  int NStrip = 50;
  double InvNStrip = 1./(double)NStrip;
  double *Acc = (double *)calloc(NStrip,sizeof(double));
  double *Vel = (double *)calloc(NStrip,sizeof(double));
  int *Pos = (int *)calloc(NStrip,sizeof(int));
  double NotRational = .1*exp(1)*sqrt(M_PI)/(double)(NStrip);
  for(int s=0;s<NStrip;s++){
    Acc[s] = (s*NotRational);
    Vel[s] = 0.;
    Pos[s] = 0;
  }
  for(int s=NStill+NStep;s>=NStep;s--){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  for(int s=NStep;s>=0;s--){
    //for(int s=0;s<NStep;s++){
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      int st = (int)(h*InvNHei*NStrip);
      for(int w=0;w<NWidth;w++){
	int w1 = w + Pos[st];
	if(w1 >= NWidth) w1 -= NWidth;
	if(w1 < 0) w1 += NWidth;
    	ImageOut.plot(w1,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    for(int st=0;st<NStrip;st++){
      Vel[st] += Acc[st];
      //if(Vel[st] > NHeight) Vel[st] = NHeight;
      Pos[st] += (int)Vel[st];
      Pos[st] -= (int)(floor(Pos[st]*InvNHei)*NHeight);
    }
    ImageOut.close();
  }
  printf("\n");
}
void DrImage::EffectOnDataR(double *data2,int l,int w,int h,int ws,int hs,int NSquare){
  data[l][(h+hs)*NWidth+(w+ws)] = data2[(NSquare-1-ws)*NSquare+hs];
}
void DrImage::EffectOnDataT(double *data2,int l,int w,int h,int ws,int hs,int NSquare){
  data[l][(h+hs)*NWidth+(w+ws)] = data2[ws*NSquare+hs];
}
void DrImage::EffectOnDataM(double *data2,int l,int w,int h,int ws,int hs,int NSquare){
  data[l][(h+hs)*NWidth+(w+ws)] = data2[(NSquare-1-hs)*NSquare+ws];
}
void DrImage::Rotor(){
  int NStill = 50;//50;
  int NStep = (int)(25.*60.*.2) - NStill;
  const int NSquare = 6;
  int NDecr = (int)((NStep)/(double)NSquare);
  int NChange = 4;
  int NSq[7] = {64,49,36,25,16,9,4};
  int NIncrement = 2;
  int NSkip = 15;
  double *data2 = (double *)calloc(NSq[0]*NSq[0],sizeof(double));
  double InvNStep = 1./(double)NStep;
  Matematica *Mate = new Matematica();
  for(int s=NStill+NStep;s>=NStep;s--){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  for(int s=NStep,sq=0;s>=0;s--){
    // if( (s%NDecr) == 0 ) sq++;
    // if(sq < 0) sq = 0;
    // if(sq > NSquare) sq = NSquare-1;
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    if((s%NSkip) == 0){
      for(int c=0;c<NChange;c++){
	int w = (int)(NWidth*Mate->Casuale());
	int h = (int)(NHeight*Mate->Casuale());
	for(int l=0;l<3;l++){
	  for(int ws=0;ws<NSq[sq];ws++){
	    if(ws+w > NWidth) continue;
	    for(int hs=0;hs<NSq[sq];hs++){
	      if(hs+h > NHeight) continue;
	      data2[hs*NSq[sq]+ws] = data[l][(h+hs)*NWidth+(w+ws)];
	    }
	  }
	  for(int ws=0;ws<NSq[sq];ws++){
	    if(ws+w > NWidth) continue;
	    for(int hs=0;hs<NSq[sq];hs++){
	      if(hs+h > NHeight) continue;
	      EffectOnDataR(data2,l,w,h,ws,hs,NSq[sq]);
	    //data[l][(h+hs)*NWidth+(w+ws)] = data2[(NMin-hs)*NSq[sq]+ws];
	    }
	  }
	}
      }
    }
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  free(data2);
  printf("\n");
}
void DrImage::LennardJones(){
  //---------------------alloc----------------------
  int NStill = 50;
  int NStep = (int)(25.*60.*.2) - NStill;
  int NRow = 20;
  int NCol = 20;
  double NPRow = 1.;//(double)NRow;
  double NPCol = 1.;//(double)NCol;
  double NStepPRow = NRow/(double)NStep;
  double NStepPCol = NCol/(double)NStep;
  double Eps = 0.1;
  double Sigma = 3.0;
  double Dt = 2.;//.1;
  double InvNStep = 1./(double)NStep;
  double InvWid = 1./(double)NWidth;
  double InvHei = 1./(double)NHeight;
  double InvRow = 1./(double)NRow;
  double InvCol = 1./(double)NCol;
  //int HeiSize = (int)(NHeight*InvRow);
  //int WidSize = (int)(NWidth*InvCol);
  double **data2 = (double **)calloc(3,sizeof(double));
  for(int l=0;l<3;l++){
    data2[l] = (double *)calloc(NHeight*NWidth,sizeof(double));
  }
  double *PPos = (double *)calloc(2*NRow*NCol,sizeof(double));
  int *PBound  = (int *)calloc(4*NRow*NCol,sizeof(int));
  double *PVel = (double *)calloc(2*NRow*NCol,sizeof(double));
  double *PFor = (double *)calloc(2*NRow*NCol,sizeof(double));
  double Edge[2] = {NWidth,NHeight};
  double InvEdge[2] = {1./(double)NWidth,1./(double)NHeight};
  Matematica *Mate = new Matematica();
  PERMUTE *Perm = (PERMUTE *)calloc(NRow*NCol,sizeof(PERMUTE));
  //-------------------------initializing-----------------------
  for(int l=0;l<3;l++){
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
	data2[l][h*NWidth+w] = data[l][h*NWidth+w];
      }
    }
  }
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      double x = (c)*InvCol*NWidth;
      double y = (r)*InvRow*NHeight;
      int p1 = (r*NCol+c)*2;
      PPos[p1  ] = x;
      PPos[p1+1] = y;
      PVel[p1  ] = (2.*Mate->Casuale()-1.);
      PVel[p1+1] = (2.*Mate->Casuale()-1.);
      PBound[(r*NCol+c)*4  ] = (int)(c*InvCol*NWidth);
      PBound[(r*NCol+c)*4+1] = (int)((c+1)*InvCol*NWidth);
      PBound[(r*NCol+c)*4+2] = (int)(r*InvRow*NHeight);
      PBound[(r*NCol+c)*4+3] = (int)((r+1)*InvRow*NHeight);
      Perm[r*NCol+c].n = r*NCol+c;
      Perm[r*NCol+c].m = r*NCol+c;      
    }
  }
  //squares of different sizes
  if(1==1){
    for(int c=0;c<NCol-1;c++){
      int NBorder = (int)((2.*Mate->Casuale()-1.)*10.);
      for(int r=0;r<NRow;r++){
	int p1 = (r*NCol+c)*4;
	int p2 = (r*NCol+(c+1))*4;
	PBound[p1+1] += NBorder;
	PBound[p2  ] += NBorder;
      }
    }
    for(int r=0;r<NRow-1;r++){
      int NBorder = (int)((2.*Mate->Casuale()-1.)*10.);
      for(int c=0;c<NCol;c++){
	int p1 = (r*NCol+c)*4;
	int p2 = ((r+1)*NCol+c)*4;
	PBound[p1+3] += NBorder;
	PBound[p2+2] += NBorder;
      }
    }
    for(int r=0;r<NRow;r++){
      for(int c=0;c<NCol;c++){
	int p1 = (r*NCol+c)*4;
	double x = PBound[p1  ];
	double y = PBound[p1+2];
	int p2 = (r*NCol+c)*2;
	PPos[p2  ] = x;
	PPos[p2+1] = y;
	PPos[p2  ] -= floor(PPos[p2  ]*InvEdge[0])*Edge[0];
	PPos[p2+1] -= floor(PPos[p2+1]*InvEdge[1])*Edge[1];
      }
    }
  }
  // for(int c=0;c<NCol;c++){
    //   int p1 = (0*NCol+c)*4;
  //   printf("x %d) %d %d %d\n",c,PBound[p1],PBound[p1+1],PBound[p1+1]-PBound[p1],NWidth);
  // }
  // for(int r=0;r<NRow;r++){
  //   int p1 = (r*NCol+0)*4;
  //   printf("y %d) %d %d %d\n",r,PBound[p1+2],PBound[p1+3],PBound[p1+3]-PBound[p1+2],NHeight);
  // }
  Mate->PermuteRandomAll(Perm,NRow*NCol);
  //------------------------loop----------------------------
  for(int s=NStill+NStep;s>=NStep;s--){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  for(int s=NStep;s>=0;s--){
    NPCol += NStepPCol;
    if(NPCol > NCol) NPCol = NCol;
    NPRow += NStepPRow;
    if(NPRow > NRow) NPRow = NRow;
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    for(int r=0;r<NRow;r++){
      for(int c=0;c<NCol;c++){
	int p1 = (r*NCol+c)*2;
	PFor[p1  ] = 0.;
	PFor[p1+1] = 0.;
      }
    }
    for(int r=0;r<NRow;r++){
      for(int c=0;c<NCol;c++){
	if(c > (int)NPCol) continue;
	if(r > (int)NPRow) continue;
	int p1 = Perm[r*NCol+c].m*2;
	if(p1 < 0 || p1 >= NRow*NCol*2){
	  continue;
	}
	for(int r1=r+1;r1<NRow;r1++){
	  for(int c1=c+1;c1<NCol;c1++){
	    //int p2 = (r1*NCol+c1)*2;
	    int p2 = Perm[r1*NCol+c1].m*2;
	    if(p2 < 0 || p2 >= NRow*NCol) continue;
	    double Dist[3];
	    for(int d=0;d<2;d++){
	      Dist[d] = PPos[p1+d] - PPos[p2+d];
	      Dist[d] -= floor(Dist[d]/Edge[d] + .5)*Edge[d];
	    }
	    double Dist2 = SQR(Dist[0]) + SQR(Dist[1]);
	    if(Dist2 > SQR(3.*Sigma))continue;
	    Dist[2] = sqrt(Dist2);
	    double a = (Sigma/sqrt(Dist2));
	    //double Force = Eps*12.*pow(a,13.) - Eps*6.*pow(a,7.);
	    // if(Force > 20.) Force = 20.;
	    double Force = -Eps/Dist2;
	    for(int d=0;d<2;d++){
	      PFor[p1+d] -= Force*Dist[d]/Dist[2];
	      PFor[p2+d] += Force*Dist[d]/Dist[2];
	    }
	  }
	}
      }
    }
    //-------------integration--------------
    for(int r=0;r<NRow;r++){
      for(int c=0;c<NCol;c++){
	if(c > (int)NPCol) continue;
	if(r > (int)NPRow) continue;
	int p1 = Perm[r*NCol+c].m*2;
	if(p1 < 0 || p1 >= NRow*NCol*2){
	  continue;
	}
	for(int d=0;d<2;d++){
	  PVel[p1+d] += PFor[p1+d]*Dt;
	  PPos[p1+d] += PVel[p1+d]*Dt;
	  PPos[p1+d] -= floor(PPos[p1+d]*InvEdge[d])*Edge[d];
	}
	//	printf("%d %lf %lf\n",r*NCol+c,PPos[p1],PPos[p1+1]);
      }
    }
    for(int l=0;l<3;l++){
      for(int h=0;h<NHeight;h++){
    	for(int w=0;w<NWidth;w++){
    	  data[l][h*NWidth+w] = 0.;
    	}
      }
    }
   //---------------updating position-------------
    for(int r=0;r<NRow;r++){
      for(int c=0;c<NCol;c++){
	int p1 = (r*NCol+c)*2;
	int xn = (int)(PPos[p1  ]);
	int yn = (int)(PPos[p1+1]);
	int xo = PBound[(r*NCol+c)*4  ];
	int yo = PBound[(r*NCol+c)*4+2];
	int WidSize = PBound[(r*NCol+c)*4+1]-PBound[(r*NCol+c)*4  ]; 
	int HeiSize = PBound[(r*NCol+c)*4+3]-PBound[(r*NCol+c)*4+2];
	//printf("%d %d) %d %d -> %d %d |%d %d| %f %f\n",s,r*NCol+c,xo,yo,xn,yn,WidSize,HeiSize,PPos[p1],PPos[p1+1]);
	for(int l=0;l<3;l++){
	  for(int ws=0;ws<WidSize;ws++){
	    int wo = xo+ws;
	    int wn = xn+ws;
	    if(wo >= NWidth){
	      printf("x out of bound %d > %d | %d\n",wo,NWidth,WidSize);
	      continue;
	    }
	    if(wn >= NWidth) wn -= NWidth;
	    for(int hs=0;hs<HeiSize;hs++){
	      int ho = yo+hs;
	      int hn = yn+hs;
	      if(ho >= NHeight){
	      	printf("y out of bound %d+%d = %d > %d |%d\n",yo,hs,ho,NHeight,HeiSize);
	      	continue;
	      }
	      if(hn > NHeight) hn -= NHeight;
	      //printf("%d= %d %d) %d->%d (%d %d) %d->%d (%d %d)\n",p1,r,c,xo,xn,WidSize,NWidth,yo,yn,HeiSize,NHeight);
	      if(hn*NWidth+wn < 0 || hn*NWidth+wn >= NWidth*NHeight){
		printf("input out of border hn %d wn %d > %d %d\n",hn,wn,yn,xn,hn*NWidth+wn,NWidth,NHeight);
		continue;
	      }
	      if(ho*NWidth+wo < 0 || ho*NWidth+wo >= NWidth*NHeight){
		printf("output out of border %d %d %d > %d %d\n",hn,wn,hn*NWidth+wn,NWidth,NHeight);
		continue;
	      }
	      data[l][hn*NWidth+wn] = data2[l][ho*NWidth+wo];
	    }
	  }
	}
      }
    }
    //-------------saving images
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  printf("\n");
  for(int l=0;l<3;l++){
    free(data2[l]);
  }
  free(data2);
  free(PPos);
  free(PFor);
  free(PVel);
  free(PBound);
}
#include <fftw3.h>
void DrImage::Fourier(){
  int NStill = 0;//50;
  int NStep = 10;//(int)(25.*60.*.2) - NStill;
  int NChange = 60;
  int NSquare = 2;
  double Norm = 1./(double)(NWidth*NHeight);
  double InvNStep = 1./(double)NStep;
  Matematica *Mate = new Matematica();
  fftw_complex *out = (fftw_complex *)fftw_malloc(NWidth*NHeight*sizeof(fftw_complex));
  fftw_complex *in  = (fftw_complex *)fftw_malloc(NWidth*NHeight*sizeof(fftw_complex));
  fftw_plan direct = fftw_plan_dft_2d(NWidth,NHeight,
  				     in,out,FFTW_FORWARD,FFTW_PATIENT);
  fftw_plan reverse = fftw_plan_dft_2d(NWidth,NHeight,
  				     out,in,FFTW_BACKWARD,FFTW_PATIENT);
  int NhInit = NHeight;
  int NwInit = NWidth;
  double **data2 = (double **)calloc(3,sizeof(double));
  for(int l=0;l<3;l++){
    data2[l] = (double *)calloc(NHeight*NWidth,sizeof(double));
  }
  for(int l=0;l<3;l++){
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
  	data2[l][h*NWidth+w] = data[l][h*NWidth+w];
      }
    }
  }
  // for(int s=NStill+NStep;s>=NStep;s--){
  //   char cImage[160];
  //   sprintf(cImage,"Image%05u.png",s);
  //   pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
  //   for(int h=0;h<NHeight;h++){
  //     for(int w=0;w<NWidth;w++){
  //   	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
  //     }
  //   }
  //   ImageOut.close();
  // }
  for(int s=NStep;s>0;s--){
    // NhInit--;
    // if(NhInit < 0) NhInit = 0;
    // NwInit--;
    // if(NwInit < 0) NwInit = 0;
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    for(int l=0;l<3;l++){
      for(int h=0;h<NHeight;h++){
	for(int w=0;w<NWidth;w++){
	  in[h*NWidth+w][0] = data2[l][h*NWidth+w];
	  out[h*NWidth+w][0] = 0.;
	  out[h*NWidth+w][1] = 0.;
	}
      }
      fftw_execute(direct);
      // for(int h=NhInit;h<NHeight;h++){
      // 	for(int w=NwInit;w<NWidth;w++){
      // 	  printf("%d %d\n",h,w);
      // 	  out[h*NWidth+w][0] = 0.;
      // 	  out[h*NWidth+w][1] = 0.;
      // 	}
      // }
      fftw_execute(reverse);
      for(int h=0;h<NHeight;h++){
      	for(int w=0;w<NWidth;w++){
      	  data[l][h*NWidth+w] = in[h*NWidth+w][0]*Norm;
      	}
      }
    }
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
    // FILE *FOut = fopen("Result1.dat","w");
    // for(int h=0;h<NHeight;h++){
    //   for(int w=0;w<NWidth;w++){
    // 	double x = w/(double)NWidth;
    // 	double y = h/(double)NHeight;
    // 	fprintf(FOut,"%lf %lf %lf\n",x,y,data[2][h*NWidth+w]);
    //   }
    // }
    // fclose(FOut);
    // FOut = fopen("Result2.dat","w");
    // for(int h=0;h<NHeight;h++){
    //   for(int w=0;w<NWidth;w++){
    // 	double x = w/(double)NWidth;
    // 	double y = h/(double)NHeight;
    // 	fprintf(FOut,"%lf %lf %lf\n",x,y,in[h*NWidth+w][0]);
    //   }
    // }
    // fclose(FOut);
  }
  fftw_destroy_plan(direct);
  fftw_destroy_plan(reverse);
  fftw_free(out);
  fftw_free(in);
  printf("\n");
}
void DrImage::Gauss(){
  int NStill = 0;//50;
  int NStep = (int)(25.*60.*.2) - NStill;
  int NSkip = 15;
  double *Plot2 = (double *)calloc(NWidth*NHeight,sizeof(double));
  double InvNStep = 1./(double)NStep;
  Matrice Mask(3,3);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NMat2 = (int)Mask.pNCol()/2;
  for(int s=NStill+NStep;s>=NStep;s--){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  for(int s=NStep,sq=0;s>=0;s--){
    // if( (s%NDecr) == 0 ) sq++;
    // if(sq < 0) sq = 0;
    // if(sq > NSquare) sq = NSquare-1;
    fprintf(stderr,"Elaborating %lf %%\r",s*InvNStep*100.);
    char cImage[160];
    sprintf(cImage,"Image%05u.png",s);
    pngwriter ImageOut(NWidth,NHeight,1.0,cImage);
    //if((s%NSkip) == 0)
      {
      //convolute filter
      for(int l=0;l<3;l++){
	for(int gx=0;gx<NWidth;gx++){
	  for(int gy=0;gy<NHeight;gy++){
	    Plot2[gy*NWidth+gx] = 0.;
	    for(int mx=0;mx<Mask.pNRow();mx++){
	      int g1x = gx + mx - NMat2;
	      if(g1x >= NWidth) continue;//g1x -= NGrid;
	      if(g1x < 0) continue;//g1x + NGrid;
	      for(int my=0;my<Mask.pNCol();my++){
	    	int g1y = gy + my - NMat2;
	    	if(g1y >= NHeight) continue;//g1y -= NGrid;
	    	if(g1y < 0) continue;//g1y + NGrid;
	    	Plot2[gy*NWidth+gx] += data[l][g1y*NWidth+g1x]*Mask.Val(mx,my);
	      }
	    }
	  }
	}
	for(int gx=0;gx<NWidth;gx++){
	  for(int gy=0;gy<NHeight;gy++){
	    data[l][gy*NWidth+gx] = Plot2[gy*NWidth+gx];
	  }
	}
      }
    }
    for(int h=0;h<NHeight;h++){
      for(int w=0;w<NWidth;w++){
    	ImageOut.plot(w,h,data[0][h*NWidth+w],data[1][h*NWidth+w],data[2][h*NWidth+w]);
      }
    }
    ImageOut.close();
  }
  free(Plot2);
  printf("\n");
}
void DrImage::SlitScan(){
  int NStep = 1;
  pngwriter ImageOut[NStep];
  for(int f=0;f<NStep;f++){
    char cImage[160];
    sprintf(cImage,"Image%05u.png",f);
    ImageOut[f] = new pngwriter(NWidth,NHeight,1.0,cImage);
  }
  for(int f=0;f<NFile;f++){
    fprintf(stderr,"Elaborating %lf %%\r",f/(double)NFile*100.);
    ReLoad(FileList[f]);
    int wInit = f;
    if(wInit >= NWidth) break;//wInit -= NWidth;
    for(int w=wInit;w<NWidth;w++){
      int f1 = 0;
      for(int h=0;h<NHeight;h++){
	//int set = f + h*NFile;
    	ImageOut[f1].plot(w,h,data[0][h*NWidth + w],data[1][h*NWidth + w],data[2][h*NWidth + w]);
      }
    }
  }
  for(int s=0;s<NStep;s++){
    ImageOut[s].close();
  }
  printf("\n");
}
#endif //USE_PNG
