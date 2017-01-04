#include "DrEffect.h"
#ifdef __glut_h__
DrEffect::DrEffect(){
  Mat = new Matematica();
  NChar = 256;
  Hist = (double **)calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++){
    *(Hist + l) = (double *)calloc(NChar,sizeof(double));
    Quart1[l] =  (unsigned char ) 0;
    Median[l] =  (unsigned char ) 0;
    Quart3[l] =  (unsigned char ) 0;
  }
  Initialize();
};
void DrEffect::Run(){
  int NStep = 1000;
  int Acc = -1;
  int Vel = 100;
  int Pos = 0.;
  Matrice *ImIn[3];
  for(int d=0;d<3;d++){
    ImIn[d] = new Matrice(ImHeight,ImWidth);
  }
  Matrice *ImOut = new Matrice(ImHeight,ImWidth);
  for(int l=0;l<NLevel-1;l++){
    for(int r=0;r<ImHeight;r++)
      for(int c=0;c<ImWidth;c++)
  	ImIn[l]->Set(r,c,(double)pixel[(r*ImWidth+c)*NLevel+l]);
  }
  for(Step=NStep;Step>0;Step--){
    for(int l=0;l<NLevel-1;l++){
  //     for(int r=0;r<ImHeight;r++){
  // 	int r1 = r + Pos;
  // 	if(r1 >= In->pNRow()) r1 -= ImHeight;
  // 	if(r1 < 0) r1 += ImHeight;
  // 	for(int c=0;c<Imwidth;c++){
	  
  // 	  Out->Set(r1,c,In->Val(r,c));
  // 	}
  // }
      Mat->BackFold(ImIn[l],ImOut,Pos);
      for(int r=0;r<ImHeight;r++)
	for(int c=0;c<ImWidth;c++)
	  pixel[(r*ImWidth+c)*NLevel+l] = (unsigned char)ImOut->Val(r,c);
    }
    Vel += Acc;
    Pos += Vel;
    //WritePixel();
    WritePngwriter();
  }
  delete ImOut;
  delete ImIn;
}
void DrEffect::EffectIncrease(){
  int Times = 2;
  int SubDiv = Times*Times;
  unsigned char *raster = (unsigned char *)malloc(SubDiv*NLevel*ImWidth*ImHeight*sizeof(char));
  //IncreaseResolution(pixel,raster,ImWidth,ImHeight,Times);
  free(pixel);
  ImWidth *= Times;
  ImHeight *= Times;
  pixel = (unsigned char *)malloc(NLevel*ImWidth*ImHeight*sizeof(char));
  memcpy(pixel,raster,NLevel*ImWidth*ImHeight);
  free(raster);
  //SubWindow1 = glutCreateSubWindow(MainWindow, 0, 0,width/2,height/2);
  return;
}
void DrEffect::EffectMC(){
  unsigned char *raster = (unsigned char *)malloc(NLevel*ImWidth*ImHeight*sizeof(char));
  int NGridw = 32;
  int NGridh = 32;//24;
  int NStep = NGridw*NGridh/2+30;
  //int NBlock = 4;
  int NPartition = 1;//(int)(NGridw*NGridh/(double)NStep);
  int *Sequence = (int *)calloc(NGridw*NGridh,sizeof(int));
  PERMUTE *Perm = (PERMUTE *)calloc(NGridw*NGridh,sizeof(PERMUTE));
  Mat->PermuteRandomAll(Perm,NGridw*NGridh);
  Mat->PermuteRandomAll(Sequence,NGridw*NGridh/2);
  //for(Step=NStep;Step>0;Step--)
    {
      SwapBlocks(pixel,raster,ImWidth,ImHeight,Perm,NGridw,NGridh,Sequence,NPartition);
      Sequence += NPartition;
      memcpy(pixel,raster,NLevel*ImWidth*ImHeight);
      glutPostRedisplay();
      printf("%d/%d %lf\n",Step,NStep,Step/(double)NStep);
      //WritePng();
      WritePixel();
    }
  free(Perm);
  free(raster);
  //free(Sequence);
  return;
}
void DrEffect::EffectCoarseGrain(int Grana){
  unsigned char *raster = (unsigned char *)malloc(4*ImWidth*ImHeight*sizeof(char));
  int NStep = 8;//NGridw*NGridh/2;
  int *Sequence = (int *)calloc(NStep*NStep,sizeof(int));
  int Discrx = 1;
  int Discry = 1;
  int TtWidth = (int)(ImWidth/(double)NStep);
  int TtHeight = (int)(ImHeight/(double)NStep);
  int TWidth = 0;
  int THeight= 0;
  for(int m=Grana;m>=0;m--){
    Discrx *= 2;
    Discry *= 2;
  }
  int m = 0;
    {
    Mat->PermuteRandomAll(Sequence,NStep*NStep);
    THeight = 0;
    TWidth = 0;
    int n=0;
    for(Step=NStep*NStep*(m+1);Step>NStep*NStep*m;Step--,n++)
      {
      THeight = (int)(((Sequence[n])/(double)NStep))*TtHeight;
      TWidth = ((Sequence[n]%NStep))*TtWidth;
      memcpy(raster,pixel,4*ImWidth*ImHeight);
      int nPos = ((THeight)*ImWidth + (TWidth))*NLevel;
      //int nPos = ((THeight-TtHeight)*ImWidth + (TWidth-TtWidth))*NLevel;
      //printf("(%d %d) (%d %d) %d %d %d %d %d\n",TWidth,THeight,ImWidth,ImHeight,Discrx,Discry,Sequence[n],nPos,(ImHeight*ImWidth + (ImWidth))*NLevel);
      Discretize(pixel+nPos,raster+nPos,ImWidth,ImHeight,TtWidth,TtHeight,Discrx,Discry);
      if((n+1)%NStep==0){
	THeight+= TtHeight;
	TWidth = 0;
      }
      TWidth += TtWidth;
      memcpy(pixel,raster,4*ImWidth*ImHeight);
      glutPostRedisplay();
      //printf("%d/%d %lf\n",Step,NStep,Step/(double)NStep);
      //WritePng();
      //WritePixel();
    }
  }
  free(raster);
  //free(Sequence);
  return;
}
void DrEffect::NablaPhi(){
  if(Median[0] == 0){printf("Create an histogram before\n");return;}
  SPLINE Weight;
  Weight.a0 = 0.;
  Weight.a1 = 1.; Weight.a2 = 0.;
  Weight.a3 = 0.; Weight.a4 = 0.;
  Matrice *Mask = new Matrice(Weight,3);
  Mask->Print();
  Matrice *ImIn = new Matrice(ImHeight,ImWidth);
  Phi = (double **)calloc(2,sizeof(double));
  for(int i=0;i<2;i++)
    *(Phi+i) = (double *)calloc(ImWidth*ImHeight,sizeof(double));
  for(int r=0;r<ImHeight;r++)
    for(int c=0;c<ImWidth;c++){
      int i = r*ImWidth+c;
      double Val = (pixel[i*NLevel]==0?1.:(double)pixel[i*NLevel])/(double)Median[0];
      Phi[0][i] = log10(Val);
    ImIn->Set(r,c,Phi[0][i]);
    }
  Mat->ApplyFilter(ImIn,Mask);
  Mask->Transpose();
  Mat->ApplyFilter(ImIn,Mask);
  for(int r=0;r<ImHeight;r++)
    for(int c=0;c<ImWidth;c++)
      Phi[1][r*ImWidth+c] = ImIn->Val(r,c);
  free(Mask);
  free(ImIn);
  printf("Fatto\n");
}
void DrEffect::EffectFilter(){
  //Mat->Binary(pixel,width,height,100);return;
  //Mat->Contrast(pixel,width,height);return;
  SPLINE Weight;
  Weight.a0 = 0.;
  Weight.a1 = 1.; Weight.a2 = 0.;
  Weight.a3 = 0.; Weight.a4 = 0.;
  Matrice *Mask = new Matrice(Weight,3);
  double Sigma = .5;
  double CutOff = 3.;
  //  Mask->FillGaussian(Sigma,CutOff);
  Mask->Print();
  Matrice *ImIn = new Matrice(ImHeight,ImWidth);
  //int l=1;
  for(int l=0;l<NLevel-1;l++){
    for(int r=0;r<ImHeight;r++)
      for(int c=0;c<ImWidth;c++)
  	ImIn->Set(r,c,(double)pixel[(r*ImWidth+c)*NLevel+l]);
    Mat->ApplyFilter(ImIn,Mask);
    // Mask->Transpose();
    // Mat->ApplyFilter(ImIn,Mask);
    for(int r=0;r<ImHeight;r++)
      for(int c=0;c<ImWidth;c++)
	pixel[(r*ImWidth+c)*NLevel+l] = (unsigned char)ImIn->Val(r,c);
  }
  delete Mask;
  delete ImIn;
  return;
}
void DrEffect::EffectMotion(){
  Matrice *ImIn = new Matrice(ImHeight,ImWidth);
  Matrice *ImOut = new Matrice(ImHeight,ImWidth);
  for(int l=0;l<NLevel-1;l++){
    for(int r=0;r<ImHeight;r++)
      for(int c=0;c<ImWidth;c++)
  	ImIn->Set(r,c,(double)pixel[(r*ImWidth+c)*NLevel+l]);
    Mat->BackFold(ImIn,ImOut,120);
    for(int r=0;r<ImHeight;r++)
      for(int c=0;c<ImWidth;c++)
	pixel[(r*ImWidth+c)*NLevel+l] = (unsigned char)ImOut->Val(r,c);
  }
  delete ImOut;
  delete ImIn;
  return;
}  
void DrEffect::Histo(){
  double Area[4] = {0.,0.,0.,0.};
  for(int l=0;l<NLevel;l++){
    for(int i=0;i<NChar;i++){
      Hist[l][i] = 0.000001;
    }
  }
  for(int l=0;l<3;l++){
    for(int h=0;h<ImHeight; h++) {
      for(int w=0;w<ImWidth; w++) {
	int Intens = (int)pixel[(h*ImWidth+w)*NLevel +l];
	//if( Intens <= 0 || Intens >= NChar -1){ continue;}
	Hist[l][Intens] += 1.;
	Area[l] += 1.;
      }
    }
  }
  double Max[4] = {0.,0.,0.,0.};
  for(int l=0;l<NLevel;l++)
    for(int i=0;i<NChar;i++){
      if(Max[l] < Hist[l][i])
	Max[l] =  Hist[l][i] ;
    }
  for(int l=0;l<NLevel;l++){
    double dMedian = 0.;
    int Fatto[3] = {0,0,0};
    for(int i=0;i<NChar;i++){
      //      Hist[l][i] /= Max[l] > 0. ? Max[l] : 1.;
      Hist[l][i] /= Area[l] > 0. ? Area[l] : 1.;
      dMedian += Hist[l][i];
      if(dMedian > .25 && !Fatto[0]){
	Quart1[l] = i;
	Fatto[0] = 1;
      }
      if(dMedian > .5 && !Fatto[1]){
	Median[l] = i;
	Fatto[1] = 1;
      }
      if(dMedian > .75 && !Fatto[2]){
	Quart3[l] = i;
	Fatto[2] = 1;
      }
      //printf("%d %d %lf Quartiles %d %d %d\n",l,i,Hist[l][i],Quart1[l],Median[l],Quart3[l]);
    }
  }
  //PrintIntensity();
}
void DrEffect::PrintIntensity(){
  FILE *Ciccia = fopen("Intensity.dat","w");
  for(int l=0;l<3;l++)
    for(int h=0;h<ImHeight; h++) 
      for(int w=0;w<ImWidth; w++) 
	fprintf(Ciccia,"%d %lf %lf %lf",(h*ImWidth+w),(double)pixel[(h*ImWidth+w)*NLevel +0],(double)pixel[(h*ImWidth+w)*NLevel +1],(double)pixel[(h*ImWidth+w)*NLevel +2]);
  fclose(Ciccia);
}
void DrEffect::BlackWhite(){
  for(int h=0;h<ImHeight; h++) {
    for(int w=0;w<ImWidth; w++) {
      double Temp = 0.;
      for(int l=0;l<3;l++)
	Temp += (double)pixel[(h*ImWidth+w)*NLevel +l];
      for(int l=0;l<3;l++)
	pixel[(h*ImWidth+w)*4 +l] = (GLubyte)(Temp*.333);
    }
  }
}
int DrEffect::Binary(int Mode){
  Histo();
  for(int l=0;l<3;l++){
    for(int i=0;i<ImWidth*ImHeight;i++){
      if(pixel[i*NLevel+l] < Median[l])
	pixel[i*NLevel+l] = (unsigned char)0;
      else 
	pixel[i*NLevel+l] = (unsigned char)255;
    }
  }
  return 1;
}
int DrEffect::Contrast(){
  Histo();
  double Factor[NLevel];
  double Min[NLevel];
  for(int l=0;l<NLevel;l++){
    Factor[l]  = .1*(double)NChar/(double)((Quart3[l] - Quart1[l]));
    Min[l] = (double) Quart1[l];
  }
  for(int l=0;l<NLevel;l++){
    for(int i=0;i<ImWidth*ImHeight;i++){
      double dNew = (double)pixel[i*NLevel+l]*Factor[l];// - Min[l];
      unsigned char New = dNew > 0. ? (unsigned char) dNew : 0;
      pixel[i*NLevel+l] = New;
    }
  }
  return 1;
}
int DrEffect::SwapBlocks(unsigned char *Picture,unsigned char *OutPicture,int width,int height,PERMUTE *Perm,int NGridw,int NGridh,int *Sequence,int NPartition){
  int Stepw = (int)(width/NGridw);
  int Steph = (int)(height/NGridh);
  memcpy(OutPicture,Picture,4*width*height);
  for(int p=0;p<NPartition;p++){
    int PSwap = Sequence[p];
    int vw1 = (int)(PSwap/(double)NGridw);
    int vh1 = PSwap%NGridh;
    int w1 = vw1*Stepw;
    int h1 = vh1*Steph;
    int vw2 = (int)(Perm[PSwap].m/(double)NGridw);//(int)(Sequence[PSwap]/(double)NGridw);
    int vh2 = Perm[PSwap].m%NGridh;//Sequence[PSwap]%NGridh;
    int w2 = vw2*Stepw;
    int h2 = vh2*Steph;
    for(int vw=0;vw<Stepw;vw++){
      for(int vh=0;vh<Steph;vh++){
	int Pos1 = (vw+w1) + (vh + h1)*width;
	int Pos2 = (vw+w2) + (vh + h2)*width;
	if(Pos1 < 0 || Pos1 >= width*height) {printf("Width out of range 0<=%d<%d\n",Pos1,width);continue;}
	if(Pos2 < 0 || Pos2 >= width*height) {printf("Height out of range 0<=%d<%d\n",Pos2,height);continue;}
	for(int l=0;l<NLevel;l++){
	  OutPicture[Pos1*NLevel+l] = Picture[Pos2*NLevel+l];
	  OutPicture[Pos2*NLevel+l] = Picture[Pos1*NLevel+l];
	}
      }
    }
  }
  return 1;
}
int DrEffect::ShiftBlocks(unsigned char *Picture,unsigned char *OutPicture,int width,int height,int *Sequence,int NGridw,int NGridh,int *Sequence1,int NPartition){
  int Stepw = (int)(width/NGridw);
  int Steph = (int)(height/NGridh);
  //memcpy(OutPicture,Picture,4*width*height);
  for(int wi=0;wi<width;wi++){
    int vw = (int) (wi/(double)width*NGridw);
    for(int hi=0;hi<height;hi++){
      int vh = (int) (hi/(double)height*NGridh);
      int PosIn = vw*NGridh + vh;
      int PosOut = Sequence[PosIn];
      int vvw = (int)(PosOut/(double)NGridw);
      int vvh = PosOut%NGridh;
      int wo = wi - (vw - vvw)*Stepw;
      int ho = hi - (vh - vvh)*Steph;
      //printf("%d %d -> %d %d / %d %d, %d %d -> %d %d / %d %d\n",wi,hi,wo,ho,width,height,vw,vh,vvw,vvh,Stepw,Steph);
      if(wo < 0 || wo >= width) {printf("Width out of range 0<=%d<%d\n",wo,width);continue;}
      if(ho < 0 || ho >= height) {printf("Height out of range 0<=%d<%d\n",ho,height);continue;}
      for(int l=0;l<NLevel;l++){
	OutPicture[(ho*width+wo)*NLevel+l] = Picture[(hi*width+wi)*NLevel+l];
      }
    }
  }
  free(Sequence);
  return 1;
}
int DrEffect::Discretize(unsigned char *Picture,unsigned char *OutPicture,int ImWidth,int ImHeight,int WWidth,int WHeight,int Blockw,int Blockh){
  for(int wo=0;wo<WWidth;wo+=Blockw){
    for(int ho=0;ho<WHeight;ho+=Blockh){
      for(int l=0;l<NLevel;l++){
	double ColAverage=0.;
	for(int bwi=0;bwi<Blockw;bwi++){
	  for(int bhi=0;bhi<Blockh;bhi++){
	    int wi = wo + bwi;
	    int hi = ho + bhi;
	    ColAverage += (double)Picture[(hi*ImWidth+wi)*NLevel+l];
	  }
	}
	ColAverage /= (double)(Blockw*Blockh);
	for(int bwi=0;bwi<Blockw;bwi++){
	  for(int bhi=0;bhi<Blockh;bhi++){
	    int wi = wo + bwi;
	    int hi = ho + bhi;
	    OutPicture[(hi*ImWidth+wi)*NLevel+l] = (unsigned char) ColAverage;
	  }
	}
      }
    }
  }
  return 1;
}
int DrEffect::IncreaseResolution(unsigned char *Picture,unsigned char *OutPicture,int ImWidth,int ImHeight,int Times){
  int OutWidth = ImWidth*Times;
  int OutHeight = ImHeight*Times;
  double DeltaIn=1./(double)(ImWidth-1);
  double DeltaOut=1./(double)(OutWidth-1);
  int NOrder = 3+1;
  double *dArray = (double *)calloc(ImWidth+NOrder+1,sizeof(double));
  for(int p=0;p<=ImWidth+NOrder;p++){
    if(p<NOrder){
      dArray[p] = 0;
    }
    else if( (NOrder<=p) && (p<=ImWidth) ){
      dArray[p] = (p-NOrder+1)*DeltaIn;//Pm[p-NOrder].Pos[CLat1];//
    }
    else if( p>ImWidth){
      dArray[p] = (p-NOrder+2)*DeltaIn;
    }
  }
//   for(int vo=0;vo<ImWidth;vo++){
//    for(int vvo=0;vvo<ImHeight;vvo++){
//      PlOut[vo][vvo] = 0.;
//      double x = DeltaOut*vo;
//      //for(int vi=0;vi<NIn;vi++){
//      for(int vi=vo-1;vi<vo+NOrder+1;vi++){
//        if(vi < 0 || vi >= ImWidth) continue;
//        double Blendx = Blend(dArray,x,vi,NOrder);
//        double y = DeltaOut*vvo;
//        //for(int vvi=0;vvi<NIn;vvi++){
//        for(int vvi=vvo-1;vvi<vvo+NOrder+1;vvi++){
// 	 if(vvi < 0 || vvi >= ImWidth) continue;
// 	 double Blendy = Mat->Blend(dArray,y,vvi,NOrder);
// 	 PlOut[vo][vvo] += Blendx*Blendy * PlIn[vi][vvi];
//        }
//      }
//    }
//   }
   for(int l=0;l<NLevel;l++){
    for(int w=0;w<ImWidth;w++){
      for(int h=0;h<ImHeight;h++){
	int ho = h*Times;
	int wo = w*Times;
	if(w>0 && h>0)
	  OutPicture[((ho+0)*OutWidth+(wo+0))*NLevel+l] = 
	    .5*(Picture[( (h+0)*ImWidth+(w-1) )*NLevel+l] + 
		Picture[( (h-1)*ImWidth+(w+0) )*NLevel+l]);
	if(h>0 && w<ImWidth-1)
	  OutPicture[( (ho+0)*OutWidth+(wo+1))*NLevel+l] = 
	    .5*(Picture[( (h+0)*ImWidth+(w+1) )*NLevel+l] + 
		Picture[( (h-1)*ImWidth+(w+0) )*NLevel+l]);
	if(h<ImHeight-1 && w<ImWidth-1)	
	  OutPicture[( (ho+1)*OutWidth+(wo+1))*NLevel+l] = 
	    .5*(Picture[( (h+0)*ImWidth+(w+1) )*NLevel+l] + 
		Picture[( (h+1)*ImWidth+(w+0) )*NLevel+l]);
	if(h<ImHeight-1 && w>0)
	  OutPicture[( (ho+1)*OutWidth+(wo+0))*NLevel+l] = 
	    .5*(Picture[( (h+0)*ImWidth+(w-1) )*NLevel+l] + 
		Picture[( (h+1)*ImWidth+(w+0) )*NLevel+l]);
     }
    }
  }
  return 1;
}
int DrEffect::Noise(unsigned char *Picture,unsigned char *OutPicture,int width,int height){
  int NMaskh = 3;
  int NMaskw = 3;
  double Mask[3][3];
  //unsigned char * raster = (unsigned char *) malloc(sizeof(unsigned char) * width * height * 4);
  double Det = 1./9.;
  Mask[0][0] = 1; Mask[0][1] = 1; Mask[0][2] = 1;
  Mask[1][0] = 1; Mask[1][1] = 1; Mask[1][2] = 1;
  Mask[2][0] = 1; Mask[2][1] = 1; Mask[2][2] = 1;
  Det = 1./10.;
  Mask[0][0] = 1; Mask[0][1] = 1; Mask[0][2] = 1;
  Mask[1][0] = 1; Mask[1][1] = 2; Mask[1][2] = 1;
  Mask[2][0] = 1; Mask[2][1] = 1; Mask[2][2] = 1;
  Det = 1./16.;
  Mask[0][0] = 1; Mask[0][1] = 2; Mask[0][2] = 1;
  Mask[1][0] = 2; Mask[1][1] = 4; Mask[1][2] = 2;
  Mask[2][0] = 1; Mask[2][1] = 2; Mask[2][2] = 1;
  //Gaussian   
  // Det = 1./249.;
  //   3 6  8   6 3
  //   6 14 19 14 6 
  //   8 19 25 19 8                
  //   6 14 19 14 6 
  //   3 6  8   6 3
  for(int l=0;l<NLevel;l++){
    for(int h=0;h<height-NMaskh; h++) {
      for(int w=0;w<width-NMaskw; w++) {
	OutPicture[(h*width+w)*NLevel +l] = 0;
	for(int lh=0;lh<NMaskh;lh++){
	  for(int lw=0;lw<NMaskw;lw++){
	    OutPicture[(h*width+w)*NLevel +l] += (unsigned char)(
								 Mask[lh][lw]*Picture[((h+lh)*width+(w+lw))*NLevel +l]*Det);
	  }
	}
      }
    }
  }
  return 1;
}
      
int DrEffect::Edges(unsigned char *Picture,unsigned char *OutPicture,int width,int height){
  int NMaskh = 3;
  int NMaskw = 3;
  //unsigned char *raster = (unsigned char *)malloc(width*height*NLevel*sizeof(unsigned char));
  double Mask[3][3];
  double Det = 1./4.;
  Mask[0][0] = 1; Mask[0][1] = 0; Mask[0][2] = -1;
  Mask[1][0] = 2; Mask[1][1] = 0; Mask[1][2] = -2;
  Mask[2][0] = 1; Mask[2][1] = 0; Mask[2][2] = -1;
  //Gradiente per riga
  for(int l=0;l<NLevel;l++)
    {
      //int l=2;
      for(int h=0;h<height-NMaskh; h++) {
	for(int w=0;w<width-NMaskw; w++) {
	  OutPicture[(h*width+w)*NLevel +l] = 0;
	  double Temp=0.;
	  for(int lh=0;lh<NMaskh;lh++){
	    for(int lw=0;lw<NMaskw;lw++){
	      if(h+lh-1 <0 ) continue;
	      if(w+lw-1<0 ) continue;
	      if(h+lh-1 > height-1) continue;
	      if(w+lw-1 > width-1) continue;
	      Temp +=  ASS((Mask[lh][lw]*Picture[((h+lh-1)*width+(w+lw-1))*NLevel +l]*Det));
	      //OutPicture[(h*width+w)*NLevel +l] += (unsigned char)(ASS((Mask[lh][lw]*Picture[((h+lh-1)*width+(w+lw-1))*NLevel +l]*Det)));
	    }
	  }
	  OutPicture[(h*width+w)*NLevel +l] += (unsigned char) Temp;
	}
      }
    }
  //Gradiente per colonna
  for(int l=0;l<NLevel;l++)
    {
      //int l=2;
      for(int h=0;h<height-NMaskh; h++) {
	for(int w=0;w<width-NMaskw; w++) {
	  double Temp=0.;
	  for(int lh=NMaskh;lh>0;lh--){
	    for(int lw=0;lw<NMaskw;lw++){
	      if(h+lh-1 <0 ) continue;
	      if(w+lw-1<0 ) continue;
	      if(h+lh-1 > height-1) continue;
	      if(w+lw-1 > width-1) continue;
	      //OutPicture[(h*width+w)*NLevel +l] += (unsigned char)(ASS((Mask[lw][lh]*Picture[((h+lh-1)*width+(w+lw-1))*NLevel +l]*Det)));
	      Temp +=  ASS((Mask[lw][lh]*Picture[((h+lh-1)*width+(w+lw-1))*NLevel +l]*Det));
	    }
	  }
	  OutPicture[(h*width+w)*NLevel +l] += (unsigned char) Temp;
	}
      }
    }
  return 1;
}
void DrEffect::DrEkeyboard(unsigned char key){
  switch (key){
  case 'a':
    glutPostRedisplay();
    break;
  case 'b':
    glutPostRedisplay();
    break;
  case 'c':
    glutPostRedisplay();
    break;
  case 'd':
    glutPostRedisplay();
    break;
  case 'e':
    glutPostRedisplay();
    break;
  case 27:
    exit(0);
    break;
  case 40:
    break;
  default:
    break;
  }
  keyboardDraw(key);
}

#endif// __glut_h__
