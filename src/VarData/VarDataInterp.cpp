#include "../include/VarData.h"

int VarData::InterParab(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Par1;
  SPLINE Par2;
  double Maxx=0.;
  double Minx=0.;
  int NShowOut=1;
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(Maxx-Minx)/(double)(NIn-1);
  double DeltaxOut=(Maxx-Minx)/(double)(NOut-1);
  double x = Minx;
  for(int p=0,pp=0;p<NIn;p++){
    if(p < NIn-2){
      Par1 = Mat->Parab(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
      Par2 = Mat->Parab(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    }
    for(;x < DeltaxIn*(p+1);x += DeltaxOut){
      if(pp == NOut-1) return NOut;
      PmOut[pp].Pos[0] = x;
      PmOut[pp].Pos[1] = Par1.a2*x*x+Par1.a1*x+Par1.a0;
      PmOut[pp].Pos[2] = Par2.a2*x*x+Par2.a1*x+Par2.a0;
      //      printf("%d) (%lf %lf %lf) \n",p,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
      pp++;
      NShowOut = pp;
      PmOut[pp].Idx = PmIn[p].Idx;
    }
  }
  return NShowOut;
}
int VarData::InterParab2(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Par1;
  SPLINE Par2;
  double Maxx=PmIn[0].Pos[0];
  double Minx=PmIn[0].Pos[0];
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(double)NIn/((Maxx-Minx)*(double)NOut);
  double DeltaxOut=0.;
  int NOutShow=NOut;
  double x = Minx;
  for(int p=0,pp=0;p<NIn;p++){
    if(p == 0){
      Par1 = Mat->Parab2(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
      Par2 = Mat->Parab2(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
      DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
      for(double x = PmIn[p].Pos[0], Stepx = 0.;
	  x <= PmIn[p+1].Pos[0];
	  x += DeltaxOut,pp++){
	Stepx = (x-PmIn[p].Pos[0]);
	PmOut[pp].Pos[0] = x;
	PmOut[pp].Pos[1] = Par1.a2*Stepx*Stepx+Par1.a1*Stepx+Par1.a0;
	PmOut[pp].Pos[2] = Par2.a2*Stepx*Stepx+Par2.a1*Stepx+Par2.a0;
      }
      continue;
    }
    else if(p < NIn-1){
      Par1 = Mat->Parab2(PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,0,1);
      Par2 = Mat->Parab2(PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,0,2);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Stepx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
	Stepx = (x-PmIn[p].Pos[0]);
	if(pp == NOut-1) return NOut;
	PmOut[pp].Pos[0] = x;
	PmOut[pp].Pos[1] = Par1.a2*Stepx*Stepx+Par1.a1*Stepx+Par1.a0;
	PmOut[pp].Pos[2] = Par2.a2*Stepx*Stepx+Par2.a1*Stepx+Par2.a0;
	NOutShow = pp;
	PmOut[pp].Idx = PmIn[p].Idx;
	//printf("%d/%d) %lf (%lf %lf %lf) \n",p,pp,Maxx,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
      }
    }
  }
  return NOutShow;
}
int VarData::InterCubica(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Cub1;
  SPLINE Cub2;
  double Maxx=PmIn[0].Pos[0];
  double Minx=PmIn[0].Pos[0];
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(double)NIn/((Maxx-Minx)*(double)NOut);
  double DeltaxOut=0.;
  int NOutShow=NOut;
  double x = Minx;
  for(int p=0,pp=0;p<NIn-2;p++){
    Cub1 = Mat->Cubica(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,1);
    Cub2 = Mat->Cubica(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,2);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Dx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      Dx = (x-PmIn[p+1].Pos[0]);
      if(pp == NOut-1) return NOutShow;
      PmOut[pp].Pos[0] = x;
      PmOut[pp].Pos[1] = Cub1.a0+Cub1.a1*Dx+Cub1.a2*Dx*Dx+Cub1.a3*Dx*Dx*Dx;
      PmOut[pp].Pos[2] = Cub2.a0+Cub2.a1*Dx+Cub2.a2*Dx*Dx+Cub2.a3*Dx*Dx*Dx;
      NOutShow = pp;
      PmOut[pp].Idx = PmIn[p].Idx;
      //printf("%d/%d) %lf (%lf %lf %lf) \n",p,pp,Maxx,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
    }
  }
  return NOutShow;
}
int VarData::InterForth(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Forth1;
  SPLINE Forth2;
  double Maxx=PmIn[0].Pos[0];
  double Minx=PmIn[0].Pos[0];
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(double)NIn/((Maxx-Minx)*(double)NOut);
  double DeltaxOut=0.;
  int NOutShow=NOut;
  double x = Minx;
  double Ref;
  for(int p=0,pp=0;p<NIn;p++){
    if(p < 1){
      Forth1 = Mat->Forth(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,PmIn[p+4].Pos,0,1);
      Forth2 = Mat->Forth(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,PmIn[p+4].Pos,0,2);
      Ref = PmIn[p+2].Pos[0]; 
    }
    else if(p < 2){
      Forth1 = Mat->Forth(PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,1);
      Forth2 = Mat->Forth(PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,2);
      Ref = PmIn[p+1].Pos[0]; 
    }
    else if(p < NIn-2){
      Forth1 = Mat->Forth(PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
      Forth2 = Mat->Forth(PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
      Ref = PmIn[p].Pos[0]; 
    }
    else if(p == NIn-2){
      Forth1 = Mat->Forth(PmIn[p-3].Pos,PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,0,1);
      Forth2 = Mat->Forth(PmIn[p-3].Pos,PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,PmIn[p+1].Pos,0,2);
      Ref = PmIn[p-1].Pos[0]; 
    }
    else if(p == NIn-1){
      Forth1 = Mat->Forth(PmIn[p-4].Pos,PmIn[p-3].Pos,PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,0,1);
      Forth2 = Mat->Forth(PmIn[p-4].Pos,PmIn[p-3].Pos,PmIn[p-2].Pos,PmIn[p-1].Pos,PmIn[p].Pos,0,2);
      Ref = PmIn[p-2].Pos[0]; 
    }
    //printf("%d/%d (%lf %lf %lf) Sp (%lf %lf %lf %lf)\n",p,pp,PmIn[p].Pos[0],PmIn[p].Pos[1],PmIn[p].Pos[2],Forth2.a0,Forth2.a1,Forth2.a2,Forth2.a3);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Dx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      Dx = (x-Ref);
      if(pp == NOut-1) return NOut;
      PmOut[pp].Pos[0] = x;
      PmOut[pp].Pos[1] = Forth1.a0+Forth1.a1*Dx+Forth1.a2*Dx*Dx+Forth1.a3*Dx*Dx*Dx+Forth1.a4*Dx*Dx*Dx*Dx;
      PmOut[pp].Pos[2] = Forth2.a0+Forth2.a1*Dx+Forth2.a2*Dx*Dx+Forth2.a3*Dx*Dx*Dx+Forth2.a4*Dx*Dx*Dx*Dx;;
      NOutShow = pp;
      PmOut[pp].Idx = PmIn[p].Idx;
      //printf("%d/%d) %lf (%lf %lf %lf) \n",p,pp,Maxx,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
    }
  }
  return NOutShow;
}int VarData::InterSpline3(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Sp1;
  SPLINE Sp2;
  double Maxx=PmIn[0].Pos[0];
  double Minx=PmIn[0].Pos[0];
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(double)(NIn-1)/((Maxx-Minx)*(double)NOut);
  double DeltaxOut=0.;
  int NOutShow=0;
  double RatioOutIn = DeltaxOut/DeltaxIn;
  for(int p=0,pp=0;p<NIn;p+=1){
    if( p == 0)
      Sp1 = Mat->Spline3Beg(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    else if( p < NIn-2 )
      Sp1 = Mat->Spline3(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    else if( p == NIn - 2 )
      Sp1 = Mat->Spline3End(PmIn[p].Pos,PmIn[p+1].Pos,0,2);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Stepx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      //printf("%d %d %lf %lf\n",p,pp,x,DeltaxIn*(p+1));
      if(pp == NOut-1) break;
      Stepx = (x-PmIn[p].Pos[0]);
      PmOut[pp].Pos[0] = x;
      PmOut[pp].Pos[2] = Sp1.a3*Stepx*Stepx*Stepx + 
	Sp1.a2*Stepx*Stepx + 
	Sp1.a1*Stepx + Sp1.a0;
    }
  }
  for(int p=0,pp=0;p<NIn;p+=1){
    //printf("%d/%d (%lf %lf %lf) Sp (%lf %lf %lf %lf)\n",p,pp,PmIn[p].Pos[0],PmIn[p].Pos[1],PmIn[p].Pos[2],Sp2.a0,Sp2.a1,Sp2.a2,Sp2.a3);
    if( p == 0)
      Sp2 = Mat->Spline3Beg(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
    else if( p < NIn-2 )
      Sp2 = Mat->Spline3(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
    else if( p == NIn - 2 )
      Sp2 = Mat->Spline3End(PmIn[p].Pos,PmIn[p+1].Pos,0,1);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Stepx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      if(pp == NOut-1)	return NOutShow;
      Stepx = (x-PmIn[p].Pos[0]);
      PmOut[pp].Pos[1] = Sp2.a3*Stepx*Stepx*Stepx + 
	Sp2.a2*Stepx*Stepx + 
	Sp2.a1*Stepx + Sp2.a0;
      //printf("%d/%d) %lf (%lf %lf %lf)\n",p,pp,Stepx,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
      NOutShow = pp;
      PmOut[pp].Idx = PmIn[p].Idx;
    }
  }
  return NOutShow;
}
int VarData::InterSpline4(PART *PmIn,PART *PmOut,int NIn,int NOut){
  SPLINE Sp1;
  SPLINE Sp2;
  double Maxx=PmIn[0].Pos[0];
  double Minx=PmIn[0].Pos[0];
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[0] < Minx)
      Minx = PmIn[p].Pos[0];
    if(PmIn[p].Pos[0] > Maxx)
      Maxx = PmIn[p].Pos[0];
  }
  double DeltaxIn=(double)(NIn-1)/((Maxx-Minx)*(double)NOut);
  double DeltaxOut=0.;
  int NOutShow=0;
  double RatioOutIn = DeltaxOut/DeltaxIn;
  for(int p=0,pp=0;p<NIn;p+=1){
    if( p == 0)
      Sp1 = Mat->Spline4Beg(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,1);
    //else if( p < NIn - 3 )
      //Sp1 = Mat->Spline4(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,1);
    //else if( p == NIn - 3 )
    //Sp1 = Mat->Spline3(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
    //Sp1 = Mat->Spline4PreEnd(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
    else if(p < NIn - 2)
      Sp1 = Mat->Spline4(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,1);
    else if( p == NIn - 2 )
      Sp1 = Mat->Spline3End(PmIn[p].Pos,PmIn[p+1].Pos,0,1);
      //Sp1 = Mat->Spline4End(PmIn[p].Pos,PmIn[p+1].Pos,0,1);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Dx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      //printf("%d %d %lf %lf\n",p,pp,x,DeltaxIn*(p+1));
      if(pp == NOut-1) break;
      Dx = (x-PmIn[p].Pos[0]);
      PmOut[pp].Pos[0] = x;
      PmOut[pp].Pos[1] = Sp1.a0+Sp1.a1*Dx+Sp1.a2*Dx*Dx+Sp1.a3*Dx*Dx*Dx+Sp1.a4*Dx*Dx*Dx*Dx;
    }
  }
  for(int p=0,pp=0;p<NIn;p+=1){
    //printf("%d/%d (%lf %lf %lf) Sp (%lf %lf %lf %lf)\n",p,pp,PmIn[p].Pos[0],PmIn[p].Pos[1],PmIn[p].Pos[2],Sp2.a0,Sp2.a1,Sp2.a2,Sp2.a3);
    if( p == 0)
    Sp2 = Mat->Spline4Beg(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,2);
    //else if( p < NIn - 3 )
    //Sp2 = Mat->Spline4(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,PmIn[p+3].Pos,0,2);
    //else if( p == NIn - 3 )
    //Sp2 = Mat->Spline3(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    //Sp2 = Mat->Spline4PreEnd(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    else if(p < NIn - 2){
      Sp2 = Mat->Spline4(PmIn[p].Pos,PmIn[p+1].Pos,PmIn[p+2].Pos,0,2);
    }
    else if( p == NIn - 2 )
      //Sp2 = Mat->Spline4End(PmIn[p].Pos,PmIn[p+1].Pos,0,2);
      Sp2 = Mat->Spline3End(PmIn[p].Pos,PmIn[p+1].Pos,0,2);
    DeltaxOut = (PmIn[p+1].Pos[0] - PmIn[p].Pos[0])*DeltaxIn;
    for(double x = PmIn[p].Pos[0], Dx = 0.;
	x <= PmIn[p+1].Pos[0];
	x += DeltaxOut,pp++){
      if(pp == NOut-1)	return NOutShow;
      Dx = (x-PmIn[p].Pos[0]);
      PmOut[pp].Pos[2] = Sp2.a0+Sp2.a1*Dx+Sp2.a2*Dx*Dx+Sp2.a3*Dx*Dx*Dx+Sp2.a4*Dx*Dx*Dx*Dx;
      //printf("%d/%d) %lf (%lf %lf %lf)\n",p,pp,Dx,PmOut[pp].Pos[0],PmOut[pp].Pos[1],PmOut[pp].Pos[2]);
      NOutShow = pp;
      PmOut[pp].Idx = PmIn[p].Idx;
    }
  }
  return NOutShow;
}
int VarData::InterBSpline(PART *PmIn,PART *PmOut,int NIn,int NOut){
  double Maxx=PmIn[0].Pos[CLat1];
  double Minx=PmIn[0].Pos[CLat1];
  int NOutShow=0;
  for(int p=0;p<NIn;p++){
    if(PmIn[p].Pos[CLat1] < Minx)
      Minx = PmIn[p].Pos[CLat1];
    if(PmIn[p].Pos[CLat1] > Maxx)
      Maxx = PmIn[p].Pos[CLat1];
  }
  double DeltaxIn=(Maxx-Minx)/(double)(NIn-1);
  double DeltaxOut=(Maxx-Minx)/(double)(NOut-1);
  int NOrder = 3+1;
  double *dArray = (double *)calloc(NIn+NOrder+1,sizeof(double));
  for(int p=0;p<=NIn+NOrder;p++){
    if(p<NOrder)
      dArray[p] = Minx;
    else if( (NOrder<=p) && (p<=NIn) )
      dArray[p] = (p-NOrder+1)*DeltaxIn+Minx;
    else if( p>NIn)
      dArray[p] = (p-NOrder+2)*DeltaxIn+Minx;
  }
  for(int pp=0;pp<NOut-1;pp++){
    for(int d=0;d<3;d++) PmOut[pp].Pos[d] = 0.;
    double x = DeltaxOut*pp+Minx;
    int vPos = (int)((x-Minx)/(Maxx-Minx)*NIn);
    for(int p=vPos-1;p<vPos+NOrder;p++){
      if(p >= NIn || p <0) continue;
      double Blend = Mat->Blend(dArray,x,p,NOrder);
      PmOut[pp].Pos[0] += Blend * PmIn[p].Pos[0];
      PmOut[pp].Pos[1] += Blend * PmIn[p].Pos[1];
      PmOut[pp].Pos[2] += Blend * PmIn[p].Pos[2];
    }
  }
  PmOut[NOut-1].Pos[0] = PmIn[NIn-1].Pos[0];
  PmOut[NOut-1].Pos[1] = PmIn[NIn-1].Pos[1];
  PmOut[NOut-1].Pos[2] = PmIn[NIn-1].Pos[2];
  NOutShow = NOut;
  return NOutShow;
}
int VarData::InterBSpline2D(double **PlIn,double **PlOut,int NIn,int NOut){
  double Max=Gen->Edge[CLat1];
  double Min=0.;
  int NOutShow=0;
  double DeltaIn=(Max-Min)/(double)(NIn-1);
  double DeltaOut=(Max-Min)/(double)(NOut-1);
  double RatioInOut = (double)(NIn/(double)NOut);
  int NOrder = 3+1;
  double *dArray = (double *)calloc(NIn+NOrder+1,sizeof(double));
  for(int p=0;p<=NIn+NOrder;p++){
    if(p<NOrder){
      dArray[p] = Min;
    }
    else if( (NOrder<=p) && (p<=NIn) ){
      dArray[p] = (p-NOrder+1)*DeltaIn+Min;
    }
    else if(p>NIn){
      dArray[p] = (p-NOrder+2)*DeltaIn+Min;
    }
  }
  for(int vo=0;vo<NOut;vo++){
   for(int vvo=0;vvo<NOut;vvo++){
     PlOut[vo][vvo] = 0.;
     double x = DeltaOut*vo+Min;
     int vxPos = (int)(vo*RatioInOut);
     for(int vi=vxPos-1,vn=0;vi<vxPos+NOrder+1;vi++){
       vn = vi;
       if(vi < 0) vn = NIn + vi;
       if(vi >= NIn) vn = vi - NIn;
       double Blendx = Mat->Blend(dArray,x,vn,NOrder);
       double y = DeltaOut*vvo+Min;
       int vyPos = (int)(vvo*RatioInOut);
       for(int vvi=vyPos-1,vvn=0;vvi<vyPos+NOrder+1;vvi++){
	 vvn = vvi;
	 if(vvi < 0) vvn = NIn + vvi;
	 if(vvi >= NIn) vvn = vvi - NIn;
	 double Blendy = Mat->Blend(dArray,y,vvn,NOrder);
	 PlOut[vo][vvo] += Blendx*Blendy*PlIn[vn][vvn];
       }
     }
   }
  }
  // for(int vo=0;vo<NOut;vo++){
  //   if(vo < NIn){//To be arranged
  //     PlOut[vo][NOut-1] = PlIn[vo][NIn-1];
  //     PlOut[NOut - 1][vo] = PlIn[NIn-1][vo];
  //   }
  // }
  //for(int vo=0;vo<NOut;vo++)for(int vvo=0;vvo<NOut;vvo++)printf("%d %d %lf %lf \n",vo,vvo,PlOut[vo][vvo],PlIn[vo][vvo]);
  NOutShow = NOut;
  free(dArray);
  return NOutShow;
}
/**
   Call InterBSpline2d to update the particle position;
 */
void VarData::SmoothGrid(int NSample){
  double *PlotIn = (double *)calloc(SQR(NSample),sizeof(double));
  double *PlotOut = (double *)calloc(SQR(NSample),sizeof(double));
  for(int nx=0;nx<NSample;nx++){
    for(int ny=0;ny<NSample;ny++){
      int n = nx*NSample+ny;
      PlotIn[n] = Pm[n].Pos[2];
    }
  }
  InterBSpline2D(PlotIn,PlotOut,NSample,NSample);
  for(int nx=0;nx<NSample;nx++){
    for(int ny=0;ny<NSample;ny++){
      int n = nx*NSample+ny;
      if(pType(n) != 0) continue;
      Pm[n].Pos[2] = PlotOut[n];
    }
  }
  free(PlotIn);
  free(PlotOut);
}
/**
   Perform a 2d BSpline interpolation on two square arrays.
 */
int VarData::InterBSpline2D(double *PlIn,double *PlOut,int NIn,int NOut){
  double Max=Gen->Edge[CLat1];
  double Min=0.;
  int NOutShow=0;
  double DeltaIn = (Max-Min)/(double)(NIn-1);
  double DeltaOut = (Max-Min)/(double)(NOut-1);
  double RatioInOut = (double)(NIn/(double)NOut);
  int NOrder = 3+1;
  double *dArray = (double *)calloc(NIn+NOrder+1,sizeof(double));
  for(int p=0;p<=NIn+NOrder;p++){
    if(p<NOrder){
      dArray[p] = Min;
    }
    else if( (NOrder<=p) && (p<=NIn) ){
      dArray[p] = (p-NOrder+1)*DeltaIn+Min;
    }
    else if(p>NIn){
      dArray[p] = (p-NOrder+2)*DeltaIn+Min;
    }
  }
  for(int vo=0;vo<NOut;vo++){
   for(int vvo=0;vvo<NOut;vvo++){
     PlOut[vo*NOut+vvo] = 0.;
     double x = DeltaOut*vo+Min;
     int vxPos = (int)(vo*RatioInOut);
     for(int vi=vxPos-1,vn=0;vi<vxPos+NOrder+1;vi++){
       vn = vi;
       if(vi < 0) vn = NIn + vi;
       if(vi >= NIn) vn = vi - NIn;
       double Blendx = Mat->Blend(dArray,x,vn,NOrder);
       double y = DeltaOut*vvo+Min;
       int vyPos = (int)(vvo*RatioInOut);
       for(int vvi=vyPos-1,vvn=0;vvi<vyPos+NOrder+1;vvi++){
	 vvn = vvi;
	 if(vvi < 0) vvn = NIn + vvi;
	 if(vvi >= NIn) vvn = vvi - NIn;
	 double Blendy = Mat->Blend(dArray,y,vvn,NOrder);
	 PlOut[vo*NOut+vvo] += Blendx*Blendy*PlIn[vn*NIn+vvn];
       }
     }
   }
  }
  // for(int vo=0;vo<NOut;vo++){
  //   if(vo < NIn){//To be arranged
  //     PlOut[vo][NOut-1] = PlIn[vo][NIn-1];
  //     PlOut[NOut - 1][vo] = PlIn[NIn-1][vo];
  //   }
  // }
  //for(int vo=0;vo<NOut;vo++)for(int vvo=0;vvo<NOut;vvo++)printf("%d %d %lf %lf \n",vo,vvo,PlOut[vo][vvo],PlIn[vo][vvo]);
  NOutShow = NOut;
  free(dArray);
  return NOutShow;
}
//minimum image convention does not seem to work
int VarData::InterBSpline1D(double *PlIn,double *PlOut,int NIn,int NOut){
  double Max=Gen->Edge[CLat1];
  double Min=0.;
  int NOutShow=0;
  double DeltaIn=(Max-Min)/(double)(NIn-1);
  double DeltaOut=(Max-Min)/(double)(NOut-1);
  int NOrder = 3+1;
  double *dArray = (double *)calloc(NIn+2*NOrder+2,sizeof(double));
  for(int vi=0;vi<=NIn+2*NOrder;vi++){
    if(vi<NOrder){
      dArray[vi] = Min;
    }
    else if(vi > NIn){
      dArray[vi] = (vi-NOrder+2)*DeltaIn+Min;
    }
    else {
      dArray[vi] = (vi-NOrder+1)*DeltaIn+Min;
    }
  }
  for(int vo=0;vo<NOut;vo++){
    //to avoid the minimum image convention
    if(vo > NOut-NOrder) continue;
    PlOut[vo] = 0.;
    double x = DeltaOut*vo+Min;
    for(int vi=vo-1,vn=0;vi<vo+NOrder+1;vi++){
      vn = vi;
      if(vi < 0){
	vn = NIn + vi;
      }
      if(vi >= NIn){
	vn = vi - NIn;
      }
      double Blendx = Mat->Blend(dArray,x,vn,NOrder);
      PlOut[vo] += Blendx*PlIn[vn];
    }
  }
  NOutShow = NOut;
  free(dArray);
  return NOutShow;
}
int VarData::InterPoly(PART *PmIn,PART *PmOut,int NIn,int NOut){
  double *Pz = (double *)calloc(NIn,sizeof(double));
  double *Px = (double *)calloc(NIn,sizeof(double));
  double *Py = (double *)calloc(NIn,sizeof(double));
  for(int p=0;p<NIn;p++){
    Px[p] = PmIn[p].Pos[CLat1];
    Py[p] = PmIn[p].Pos[CLat2];
    Pz[p] = PmIn[p].Pos[CNorm];
  }
  Spline *Sp1  = new Spline(NIn);
  Spline *Sp2  = new Spline(NIn);
  Mat->Polinomio(Px,Py,NIn,Sp1);
  Mat->Polinomio(Px,Pz,NIn,Sp2);
  double Length = PmIn[NIn-1].Pos[CLat1] - PmIn[0].Pos[CLat1];
  double Deltax = Length / (double)NOut;
  double Pos = PmIn[0].Pos[CLat1];
  for(int p=0;p<NOut;p++){
    Pos += Deltax;
    PmOut[p].Pos[CLat1] = Pos;
    PmOut[p].Pos[CLat2] = 0.;
    PmOut[p].Pos[CNorm] = 0.;
    for(int n=0;n<NIn;n++){
      PmOut[p].Pos[CLat2] += Sp1->GetCoe(n)*Mat->Elevato(Pos,n);
      PmOut[p].Pos[CNorm] += Sp2->GetCoe(n)*Mat->Elevato(Pos,n);
    }
    //printf("(%lf %lf %lf)\n",PmOut[p].Pos[CLat1],PmOut[p].Pos[CLat2],PmOut[p].Pos[CNorm]);
  }
  delete  Sp1;delete  Sp2;
  free(Pz);free(Py);free(Px);
  return NOut;
}
// int VarData::InterDerMatrix(PART *Pm,int NMass,SPLINE Wg,double Offset){
//   Matrice *Coeff = new Matrice(Wg,1);
//   Matrice *Resp = new Matrice(NMass,NMass);
//   for(int r=0;r<NMass;r++){
//     if(Pm[r].Typ != 0){ Resp->Set(r,r,1.);continue;}
//     if(r >= 2) Resp->Set(r,r-2,Coeff->Val(2,0));
//     if(r >= 1) Resp->Set(r,r-1,Coeff->Val(2,1));
//     if(r < NMass-1) Resp->Set(r,r+1,Coeff->Val(2,3));
//     if(r < NMass-2) Resp->Set(r,r+2,Coeff->Val(2,4));
//     Resp->Set(r,r,Coeff->Val(2,2));
//   }
//   double *Known = (double *) calloc(NMass,sizeof(double));
//   double *UnKnown = (double *) calloc(NMass,sizeof(double));
//   for(int p=0;p<NMass;p++)
//     Known[p] = Pm[p].Pos[CNorm] - Offset;
//   Resp->Solve(Known,UnKnown);
//   for(int p=0;p<NMass;p++){
//     if(Pm[p].Typ != 0){continue;}
//     Pm[p].Pos[CNorm] = UnKnown[p] + Offset;
//     //printf("%lf %lf\n",Known[p],UnKnown[p]);
//   }
//   free(Known);
//   free(UnKnown);
//   // Resp->invert();
//   // for(int r=0;r<NMass;r++){
//   //   double Agg=0.;
//   //   for(int c=0;c<NMass;c++){
//   //     if(Pm[r].Typ == 0)
//   // 	Agg += Resp->Val(r,c) * Pm[c].Pos[CNorm];	
//   //     else 
//   // 	Agg += Resp->Val(r,c) * Pm[c].Pos[CNorm];
//   //   }
//   //   Pm[r].Pos[CNorm] = Agg;
//   // }
//   delete Coeff;
//   delete Resp;
//   return NMass;  
// }
void VarData::SmoothGrid(int NSample,char *FWrite){
  double InvNSample = 1./(double)NSample;
  double **PlotIn = (double **)calloc(SQR(NSample),sizeof(double));
  double **PlotOut = (double **)calloc(SQR(NSample),sizeof(double));
  double **Count = (double **)calloc(SQR(NSample),sizeof(double));
  double Round = 0.00001;
  for(int v=0;v<NSample;v++){
    PlotIn[v] = (double *)calloc(NSample,sizeof(double));
    PlotOut[v] = (double *)calloc(NSample,sizeof(double));
    Count[v] = (double *)calloc(NSample,sizeof(double));
  }
  int Nx = 0;
  int Ny = 0;
  for(int p=0;p<pNPart();p++){
    if(Pm[p].Pos[CLat2] < Pm[p+1].Pos[CLat2]) Ny++;
    else break;
  }
  Ny++;
  Nx = (int)(pNPart()/(double)Ny);
  for(int p=0;p<pNPart();p++){
    int vx = (int)((Pm[p].Pos[CLat1]+Round)*pInvEdge(CLat1)*NSample);
    int vy = (int)((Pm[p].Pos[CLat2]+Round)*pInvEdge(CLat2)*NSample);
    if(vx < 0 || vx >= NSample)continue;
    if(vy < 0 || vy >= NSample)continue;
    PlotIn[vx][vy] += Pm[p].Pos[CNorm];
    Count[vx][vy] += 1.;
  }
  for(int vx=0;vx<NSample;vx++)
    for(int vy=0;vy<NSample;vy++)
      PlotIn[vx][vy] *= Count[vx][vy] > 0. ? 1./Count[vx][vy] : 1.;
  InterBSpline2D(PlotIn,PlotOut,NSample,NSample);
  if(1==0){
    FILE *F2Write = fopen(FWrite,"w");
    for(int vx=0;vx<NSample;vx++){
      for(int vy=0;vy<NSample;vy++){
	double x = vx*pEdge(CLat1)*InvNSample;
	double y = vy*pEdge(CLat2)*InvNSample;      
      fprintf(F2Write,"%lf %lf %lf\n",x,y,PlotOut[vx][vy]);
      }
    }
    fclose(F2Write);
  }
  for(int v=0;v<NSample;v++){
    free(PlotIn[v]);
    free(PlotOut[v]);
    free(Count[v]);
  }
  free(PlotIn);
  free(PlotOut);
  free(Count);
}
