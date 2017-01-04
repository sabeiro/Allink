#include "../include/Matematica.h"

double Matematica::QBezier(double *P1,double *P2,double *P3,double x,int y){
  double Resp;
//   int Order = 3;
//   double *Coeff;
//   Coeff = (double *)malloc(Order*sizeof(double));
//   for(int i=0;i<Order;i++){
//     Coeff[i] = Binomial(Order,i)*P1[y];
//   }//deCasteljau algorithm
//   //P(x) = sum_i Coeff[i]*t^i*(t-1)^(n-1)
  Resp = QUAD((1-x))*P1[y]+2*(1-x)*x*P2[y]+QUAD(x)*P3[y];
  return Resp;
}
SPLINE Matematica::Spline3Beg(double *P1,double *P2,double *P3,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  Sp.a0 = P1[y];
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a2 = 0.;
  //  printf("%lf %lf\n",Sp.a2,Par.A);
  Sp.a3 = (Par.a2 - Sp.a2)/(3.*Deltax);
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax;
  SpMem.a2 = Par.a2;
  SpMem.a1 = Sp.a1+2*Sp.a2*Deltax+3*Sp.a3*Deltax*Deltax;
  Sp.a4 = 0.;
  //printf("%lf %lf|%lf %lf %lf\n",Sp.a0,Sp.a1,SplineA1,Sp.a2,Sp.a3);
  return Sp;
}
SPLINE Matematica::Spline3(double *P1,double *P2,double *P3,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  Sp.a3 = (Par.a2 - Sp.a2)/(3.*Deltax);
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax;
  SpMem.a2 = Par.a2;
  //printf("%lf %lf|%lf %lf|%lf %lf\n",Sp.a0,SplineA1,Sp.a1,Sp.a2,SplineA2,Sp.a3);
  SpMem.a1 = Sp.a1;
  Sp.a4 = 0.;
  return Sp;
}
SPLINE Matematica::Spline3End(double *P1,double *P2,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  Sp.a3 = (0. - SpMem.a2)/(3.*Deltax);
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax;
  //  printf("%lf %lf %lf %lf\n",Sp.a0,Sp.a1,Sp.a2,Sp.a3);
  SpMem.a2 = 0.;
  Sp.a4 = 0.;
  return Sp;
}
SPLINE Matematica::Spline4Beg(double *P1,double *P2,double *P3,double *P4,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  SPLINE Cub = Cubica(P1,P2,P3,P4,x,y);
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a0 = P1[y];
  Sp.a2 = 0.;
  Sp.a4 = 3.*Cub.a3/Deltax + (Sp.a2 - Cub.a2)/QUAD(Deltax);
  Sp.a4 /= 12.;
  Sp.a3 = Cub.a3 - 6.*Sp.a4*Deltax;
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax - Sp.a4*Deltax*Deltax*Deltax;
  SpMem.a2 = Cub.a2;
  SpMem.a3 = Sp.a3;
  return Sp;
}
SPLINE Matematica::Spline4(double *P1,double *P2,double *P3,double *P4,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  SPLINE Cub = Cubica(P1,P2,P3,P4,x,y);
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  Sp.a4 = 3.*Cub.a3/Deltax + (Sp.a2 - Cub.a2)/QUAD(Deltax);
  Sp.a4 /= 12.;
  Sp.a3 = Cub.a3 - 6.*Sp.a4*Deltax;
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax - Sp.a4*Deltax*Deltax*Deltax;
  SpMem.a2 = Cub.a2;
  //  printf("%lf %lf|%lf %lf|%lf %lf %lf\n",Sp.a0,SpMem.a1,Sp.a1,Sp.a2,SpMem.a2,Sp.a3,Sp.a4);
  SpMem.a1 = Sp.a1;
  SpMem.a3 = Cub.a3;
  return Sp;
}
SPLINE Matematica::Spline4(double *P1,double *P2,double *P3,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  Sp.a1 = SpMem.a1;
//   Sp.a1 = Par.B - 13.*Sp.a2*Deltax/2. - 3.*(P2[y] - P1[y])/Deltax - Par.A*Deltax/2.;
//   Sp.a1 /= 13.;
  Sp.a4 = Par.a2 + 2.*Sp.a2 - 3.*(P2[y] - P1[y])/(Deltax*Deltax) + 3.*Sp.a1/Deltax;
  Sp.a4 /= 3.;
  Sp.a3 = (P2[y] - P1[y])/(Deltax*Deltax*Deltax) - Sp.a1/(Deltax*Deltax) - Sp.a2/Deltax - Sp.a4*Deltax;
  SpMem.a2 = Par.a2;
  //printf("%lf %lf|%lf %lf|%lf %lf %lf\n",Sp.a0,SpMem.a1,Sp.a1,Sp.a2,SpMem.a2,Sp.a3,Sp.a4);
  SpMem.a1 = Par.a1;
  return Sp;
}
SPLINE Matematica::Spline4PreEnd(double *P1,double *P2,double *P3,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  SPLINE Par = Parab2(P1,P2,P3,x,y);
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  Sp.a3 = (Sp.a2 - Par.a2)/Deltax + 4.*SpMem.a3;
  Sp.a4 = (SpMem.a3 - Sp.a3)/(6.*Deltax);
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax - Sp.a4*Deltax*Deltax*Deltax;
  SpMem.a2 = Par.a2;
  SpMem.a1 = Sp.a1;
  SpMem.a3 = Sp.a3;
  return Sp;
}
SPLINE Matematica::Spline4End(double *P1,double *P2,int x,int y){
  SPLINE Sp; Sp.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Sp;
  }
  double Deltax = P2[x] - P1[x];
  Sp.a0 = P1[y];
  Sp.a2 = SpMem.a2;
  //Sp.a3 = (Sp.a2 - 0.)/Deltax + 4.*SpMem.a3;
  Sp.a3 = 0.;
  Sp.a4 = (SpMem.a3 - Sp.a3)/(6.*Deltax);
  Sp.a1 = (P2[y] - P1[y])/Deltax - Sp.a2*Deltax - Sp.a3*Deltax*Deltax - Sp.a4*Deltax*Deltax*Deltax;
  SpMem.a2 = 0.;
  SpMem.a1 = Sp.a1;
  SpMem.a3 = Sp.a3;
  return Sp;
}
//BSplines
//   double w = x - floor(x) - 1.;
//   xSp.a3 = (1./6.)*w * w * w;
//   xSp.a0 = (1./6.) + (1./2.)*w*(w-1.)-xSp.a3;
//   xSp.a2 = w + xSp.a0 - 2*xSp.a3;
//   xSp.a1 = 1. - xSp.a0 - xSp.a2 - xSp.a3;
//   Resp = xSp.a0*P1[0]+xSp.a1*P2[0]+xSp.a2*P3[0];
  //  w = y - 3.;
//   ySp.a3 = (1./6.)*w * w * w;
//   ySp.a0 = (1./6.) + (1./2.)*w*(w-1.)-ySp.a3;
//   ySp.a2 = w + ySp.a0 - 2*ySp.a3;
//   ySp.a1 = 1. - ySp.a0 - ySp.a2 - ySp.a3;
  //  Resp *= ySp.a0*P1[1]+ySp.a1*P2[1]+ySp.a2*P3[1];
SPLINE Matematica::Parab(double *P1,double *P2,double *P3,int x,int y){
  SPLINE Par; Par.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Par;
  }
  double Deltax = P2[x] - P1[x];
  double Deltax2 = QUAD(P2[x]) - QUAD(P1[x]);
  double Dx = P3[x] - P1[x];
  double Dx2 = QUAD(P3[x]) - QUAD(P1[x]);
  double Deltay = P2[y] - P1[y];
  double Dy = P3[y] - P1[y];
  double b1= Dy - Deltay/Deltax2*Dx2;
  double b2= Dx - Deltax/Deltax2*Dx2;
  Par.a1 = b1/b2;
  Par.a2 = Deltay/Deltax2 - Par.a1*Deltax/Deltax2;
  Par.a0 = P1[y] - Par.a1*P1[x] - Par.a2*QUAD(P1[x]);
  Par.a3 = 0.;
  Par.a4 = 0.;
  return Par;
}
SPLINE Matematica::Parab2(double *PA,double *PB,double *PC,int x,int y){
  SPLINE Par; Par.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Par;
  }
  double DxAB = PA[x] - PB[x];
  double DxCB = PC[x] - PB[x];
  double DyCB = PC[y] - PB[y];
  double DyAB = PA[y] - PB[y];
  double a1 = DyCB * DxAB - DyAB * DxCB;
  double a2 = DxCB*DxCB*DxAB - DxCB*DxAB*DxAB;
  Par.a2 = a1/a2;
  Par.a1 = DyAB/DxAB - Par.a2*DxAB;
  Par.a0 = PB[y];
  Par.a3 = 0.;
  Par.a4 = 0.;
  return Par;
}
CIRCLE Matematica::Osculante(double *PA,double *PB,double *PC,int x,int y){
  CIRCLE Cir; Cir.yC = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Cir;
  }
  double DxCA = PC[x] - PA[x];
  double DxBA = PB[x] - PA[x];
  double DyAC = PA[y] - PC[y];
  double DyAB = PA[y] - PB[y];
  double DxAB2 = QUAD(PA[x]) - QUAD(PB[x]);
  double DxAC2 = QUAD(PA[x]) - QUAD(PC[x]);
  double DyBA2 = QUAD(PB[y]) - QUAD(PA[y]);
  double DyCA2 = QUAD(PC[y]) - QUAD(PA[y]);
  double a1 = DxCA*(DyBA2 - DxAB2) - DxBA*(DyCA2 - DxAB2);
  double a2 = 2.*(DyAC*DxBA - DyAB*DxCA);
  Cir.yC = a1/a2;
  Cir.xC = (DyBA2 - DxAB2 + 2.*Cir.yC*DyAB)/DxBA;
  Cir.Rad = sqrt( QUAD((PA[x] - Cir.xC)) + QUAD((PA[y] - Cir.yC)) );
  return Cir;
}
SPLINE Matematica::Cubica(double *PA,double *PB,double *PC,double *PD,int x,int y){
  SPLINE Cub; Cub.a0 = 0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Cub;
  }
  double DyDB = PD[y] - PB[y];
  double DyCB = PC[y] - PB[y];
  double DyAB = PA[y] - PB[y];
  double DxDB = PD[x] - PB[x];
  double DxCB = PC[x] - PB[x];
  double DxAB = PA[x] - PB[x];
  double DxCA = PC[x] - PA[x];
  double Num3 = DyDB/DxDB - DyAB/DxAB - (DyCB/DxCB - DyAB/DxAB)*(DxDB - DxAB)/DxCA;
  double Den3 = DxDB*DxDB - DxAB*DxAB - (DxCB*DxCB - DxAB*DxAB)*(DxDB - DxAB)/DxCA;
  Cub.a3 = Num3/Den3;
  Cub.a2 = DyCB/DxCB - DyAB/DxAB - Cub.a3*(DxCB*DxCB - DxAB*DxAB);
  Cub.a2 /= DxCA;
  Cub.a1 = DyAB/DxAB - Cub.a2*DxAB - Cub.a3*DxAB*DxAB;
  Cub.a0 = PB[y];
  Cub.a4 = 0.;
  return Cub;
}
SPLINE Matematica::Forth(double *PA,double *PB,double *PC,double *PD,double *PE,int x,int y){
  SPLINE Forth; Forth.a0=0.;
  if(x < 0 || x>3 || y < 0 || y>3){
    printf("x,y must be 0 <= %d , %d < 3\n",x,y);
    return Forth;
  }
  double Dx = ASS((PC[x] - PB[x]));
  Forth.a3 = PA[y]/12. - PB[y]/6. + PD[y]/6. - PE[y]/12.;
  Forth.a3 /= -Dx*Dx*Dx;
  Forth.a1 = - PA[y]/12. + 4.*PB[y]/6. - 4.*PD[y]/6. + PE[y]/12.;
  Forth.a1 /= -Dx;
  Forth.a4 = PA[y]/24. - PB[y]/6. + 1./4.*PC[y] - PD[y]/6. + PE[y]/24.;
  Forth.a4 /= Dx*Dx*Dx*Dx;
  Forth.a2 = -PA[y]/24. + 2./3.*PB[y] - 5./4.*PC[y] + 2./3.*PD[y] - PE[y]/24.;
  Forth.a2 /= Dx*Dx;
  Forth.a0 = PC[y];
  //printf("Dx %lf Coeff %lf %lf %lf %lf %lf \n",Dx,Forth.a0,Forth.a1,Forth.a2,Forth.a3,Forth.a4);
  return Forth;

  double DxCA = PC[x] - PA[x];
  double DxCB = PC[x] - PB[x];
  double DxCD = PC[x] - PD[x];
  double DxCE = PC[x] - PE[x];
  double DxAB = PA[x] - PB[x];
  double DxAD = PA[x] - PD[x];
  double DyAC = (PA[y] - PC[y])/DxCA;
  double DyBC = (PB[y] - PC[y])/DxCB;
  double DyDC = (PD[y] - PC[y])/DxCD;
  double DyEC = (PE[y] - PC[y])/DxCE;
  //if(DxCB == DxCD)
  if (1==1)
    {
      double DxDB = PD[x] - PB[x];
      double DxEA = PE[x] - PA[x];
      double DyBD = (PB[y] - PD[y])/DxDB;
      double DyAE = (PA[y] - PE[y])/DxEA;
      //      Forth.a4 = - Noto34/Num34;
      Forth.a4 = DyBC - DyDC - DyEC + DyAC;
      Forth.a4 /= (DxCE*DxCE*DxCE - DxCA*DxCA*DxCA) - (DxCB*DxCB*DxCB - DxCD*DxCD*DxCD);
      Forth.a2 = DyBC - DyDC - Forth.a4*(DxCB*DxCB*DxCB - DxCD*DxCD*DxCD);
      Forth.a3 = DyBD - DyAE;
      Forth.a3 /= (DxCB*DxCB*DxCB - DxCD*DxCD*DxCD)/DxDB - 
	(DxCA*DxCA*DxCA - DxCE*DxCE*DxCE)/DxEA;
      Forth.a1 = DyAE - Forth.a3*(DxCA*DxCA*DxCA - DxCE*DxCE*DxCE)/DxEA;
    }
  else
    {
      double DxAE = PA[x] - PE[x];
      double Den34 = DxAB*(DxCB*DxCB - DxCA*DxCA) - DxAB*(DxCD*DxCD - DxCA*DxCA);// 0.;
      double Num34 = DxAB*(DxCD*DxCD*DxCD - DxCA*DxCA*DxCA) - DxAD*(DxCB*DxCB*DxCB - DxCA*DxCA*DxCA);//21Dx
      double Noto34 = DyBC*DxAD - DyAC*DxAD - DyDC*DxAB + DyAC*DxAB;
      double Num4 = -DxAB*(Noto34*(DxCE*DxCE - DxCA*DxCA) - DyEC + DyAC) + 
	DxAE*(Noto34*(DxCB*DxCB - DxCA*DxCA) - DyBC + DyAC);
      double Den4 = DxAB*(Num34/Den34*(DxCE*DxCE - DxCA*DxCA) + (DxCE*DxCE*DxCE -  DxCA*DxCA*DxCA) ) - 
	DxAE*(Num34/Den34*(DxCB*DxCB - DxCA*DxCA) + (DxCB*DxCB*DxCB - DxCA*DxCA*DxCA));
      //printf("Dx %lf %lf %lf %lf %lf %lf %lf\n",DxCB,Den34,Num34,Noto34,Num4,Den4,Num34/Den34);
      Forth.a4 = Num4/Den4;
      Forth.a3 = Num34/Den34*Forth.a4 + Noto34;
      Forth.a2 = DyBC - DyAC - Forth.a3*(DxCB*DxCB - DxCA*DxCA) - Forth.a4* (DxCB*DxCB*DxCB - DxCA*DxCA*DxCA);
      Forth.a2 /= DxAB;
      Forth.a1 = DyAC - Forth.a2*DxCA - Forth.a3*DxCA*DxCA - Forth.a4*DxCA*DxCA*DxCA;
    }
  Forth.a0 = PC[y];
  //printf("%lf %lf %lf %lf %lf \n",Forth.a0,Forth.a1,Forth.a2,Forth.a3,Forth.a4);
  return Forth;
}
int Matematica::Polinomio(double *Px,double *Py,int NMass,Spline *Sp){
  Matrice *Coeff = new Matrice(NMass,NMass);
  for(int r=0;r<NMass;r++){
    Sp->SetCoe( 0. , r);
    for(int c=0;c<NMass;c++){
      Coeff->Set(r,c,Elevato(Px[r],c));
    }
  }
  Coeff->Invert();
  for(int c=0;c<NMass;c++){
    for(int r=0;r<NMass;r++){
      Sp->AddCoe( Coeff->Val(r,c)*Py[c] , r);
    }
  }
  delete  Coeff;
  return 1;
}
int Matematica::DerMatrix(double *Px,double *Py,int NMass,SPLINE Wg,Spline *Sp){
  Matrice *Coeff = new Matrice(NMass,NMass);
  double MenoDue = - .75*Wg.a3 + 1.5*Wg.a4;	
  //MenoDue += .125*Wg.a1 - .125*Wg.a2; //O(h^4)
  double MenoUno = -.5*Wg.a1 + Wg.a2 + 1.5*Wg.a3 - 6.*Wg.a4;
  //MenoUno += -.125*Wg.a1 - .5*Wg.a2;
  double Zero = Wg.a0 - 2.*Wg.a2 + 9.*Wg.a4;
  //Zero += 0.75*Wg.a2;
  double PiuUno = .5*Wg.a1 + Wg.a2 - 1.5*Wg.a3 - 6.*Wg.a4;
  //PiuUno += .25*Wg.a1 + .5*Wg.a2;
  double PiuDue = .75*Wg.a3 + 1.5*Wg.a4;
  //PiuDue += -.125*Wg.a1 -.125*Wg.a2;
  for(int r=0;r<NMass;r++){
    if(r> 1) Coeff->Set(r,r-2,MenoDue*Px[r-2]);
    if(r> 0) Coeff->Set(r,r-1,MenoUno*Px[r-1]);
    if(r< NMass-1) Coeff->Set(r,r+1,PiuUno*Px[r+1]);
    if(r< NMass-2) Coeff->Set(r,r+2,PiuDue*Px[r+2]);
    Coeff->Set(r,r,Zero*Px[r]);
  }
  //Coeff->Print();
  Coeff->Invert();
  for(int r=0;r<NMass;r++){
    for(int c=0;c<NMass;c++){
      Sp->AddCoe( Coeff->Val(r,c)*Py[c] , r);
    }
  }
  delete Coeff;
  return 0;
}
double Matematica::LinInterp(double Px1,double Px2,double Py1,double Py2,double x){
  double m = (Py2-Py1)/(Px2-Px1);
  double q = Py1 - m*Px1;
  return m*x+q;
}
RETTA Matematica::InterRett(double *Px,double *Py,int NMass){
  RETTA r1;
  double Uno=0.;double Due=0.;double Tre=0.;double Quattro=0.;
  double UnoUno=0.;double DueUno=0.;double ZeroUno=0.;double ZeroDue=0.;
  for(int i=0;i<NMass;i++){
    //printf("%lf %lf\n",Px[i],Py[i]);
    Uno += Px[i];
    Due += QUAD(Px[i]);
    ZeroUno += Py[i];
    UnoUno += Px[i]*Py[i];
    ZeroDue += QUAD(Py[i]);
  }
  double Mediax = Uno / (double) NMass;
  double Mediay = ZeroUno / (double) NMass;
  double Scartox = (Due - NMass*Uno*Uno)/(double)(NMass-0);
  double Scartoy = (ZeroDue - NMass*ZeroUno*ZeroUno)/(double)(NMass-0);
  r1.m = (NMass*UnoUno - Uno*ZeroUno) / (NMass*Due - Uno*Uno);
  r1.q = (ZeroUno - r1.m*Uno)/NMass;
  double Posteriori = 0.;
  r1.Cov = UnoUno/(double)NMass - Mediax*Mediay;
  for(int i=0;i<NMass;i++){
    Posteriori += QUAD(( Px[i]*r1.m + r1.q - Py[i] ));
  }
  r1.Corr = r1.Cov / ( sqrt(Scartox*Scartoy) );
  r1.ErrY = sqrt(Posteriori/(double)(NMass-2));
  r1.ErrM = r1.ErrY * sqrt( NMass / (NMass * Due - Uno*Uno)); 
  r1.ErrQ = r1.ErrY * sqrt( Due / (NMass*Due - Uno*Uno) );
  //  printf("m %lf q %lf r %lf sigma %lf Mediax %lf Mediay %lf\n",r1.m,r1.q,r1.r,r1.Corr,Mediax,Mediay);
  return r1;
}
// Vettore Matematica::Directive(double *Px,int NMass){
//   Vettore v1(3);
  

// }
RETTA Matematica::InterRett(double *Px,double *Py,double *Peso,int NMass){
  RETTA r1;
  double Mediax = 0.;
  double Mediay = 0.;
  double Scartox = 0.;
  double Scartoy = 0.;
  double Pesox = 0.;
  double Uno=0.;double Due=0.;double Tre=0.;double Quattro=0.;
  double UnoUno=0.;double DueUno=0.;double ZeroUno=0.;double ZeroDue=0.;
  for(int i=0;i<NMass;i++){
    Pesox += Peso[i];
    Uno += Px[i]*Peso[i];
    Due += QUAD(Px[i]*Peso[i]);
    ZeroUno += Py[i];
    UnoUno += Px[i]*Py[i]*Peso[i];
    ZeroDue += QUAD(Py[i]);
  }
  Mediax = Uno / (double) NMass/Pesox;
  Mediay = ZeroUno / (double) NMass;
  Scartox = (Due - Uno*Uno/QUAD(Pesox))/(double)NMass;
  Scartoy = (ZeroDue - ZeroUno*ZeroUno)/(double)NMass;
  r1.Corr = (UnoUno - Uno*ZeroUno) / (double) NMass/Pesox ;
  r1.r = r1.Corr / ( sqrt(Scartox*Scartoy) );
  r1.m = (NMass*UnoUno - Uno*ZeroUno) / (NMass*Due - Uno*Uno)/Pesox;
  r1.ErrM = Scartoy * sqrt( Due / (NMass * Due - Uno*Uno)/QUAD(Pesox)); 
  r1.q = (ZeroUno - r1.m*Uno)/NMass;
  r1.ErrQ = Scartoy * sqrt( NMass / (NMass * Due - Uno*Uno)/QUAD(Pesox) );
  return r1;
}
RETTA Matematica::InterExp(double *Px,double *Py,int NMass){
  RETTA r1;
  double Uno=0.;double Due=0.;double Tre=0.;double Quattro=0.;
  double UnoUno=0.;double DueUno=0.;double ZeroUno=0.;double ZeroDue=0.;
  for(int i=0;i<NMass;i++){
    if(Py[i] <= 0.){
      //printf("Negative number not allowed for an exponential interpolation %lf\n",Py[i]);
      continue;
    }
    Uno += Px[i];
    Due += QUAD(Px[i]);
    ZeroUno += log10(Py[i]);
    UnoUno += Px[i]*log10(Py[i]);
    ZeroDue += QUAD(log10(Py[i]));
  }
  double Mediax = Uno / (double) NMass;
  double Mediay = ZeroUno / (double) NMass;
  double Scartox = (Due - NMass*Uno*Uno)/(double)(NMass-1);
  double Scartoy = (ZeroDue - NMass*ZeroUno*ZeroUno)/(double)(NMass-1);
  r1.m = (NMass*UnoUno - Uno*ZeroUno) / (NMass*Due - Uno*Uno);
  r1.q = (ZeroUno - r1.m*Uno)/NMass;
  double Posteriori = 0.;
  r1.Cov = UnoUno/(double)NMass - Mediax*Mediay;
  for(int i=0;i<NMass;i++){
    Posteriori += QUAD(( Px[i]*r1.m + r1.q - log10(Py[i]) ));
  }
  r1.Corr = r1.Cov / ( sqrt(Scartox*Scartoy) );
  r1.ErrY = sqrt(Posteriori/(double)(NMass-2));
  r1.ErrM = r1.ErrY * sqrt( NMass / (NMass * Due - Uno*Uno)); 
  r1.ErrQ = r1.ErrY * sqrt( Due / (NMass*Due - Uno*Uno) );
  //  printf("m %lf q %lf r %lf sigma %lf Mediax %lf Mediay %lf\n",r1.m,r1.q,r1.r,r1.Corr,Mediax,Mediay);
  return r1;
}
MOMENTI Matematica::InterGauss(double *Px,double *Py,int NMax){
  MOMENTI m1; m1.Uno=0.; m1.Due=0.; m1.Tre=0.;m1.Delta=0.;m1.Num=0;
  m1.Min = Px[0];
  m1.Max = Px[0];
  double Count = 0.;
  double Sum2 = 0.;
  for(int i=0;i<NMax;i++){
#ifdef MAT_DEBUG
    if(Py[i] < 0.){printf("Invalid weight %lf < 0\n",Py[i]);}
#endif
    if(Px[i]<m1.Min) m1.Min = Px[i];
    if(Px[i]>m1.Max) m1.Max = Px[i];
    m1.Uno += Px[i]*Py[i];
    Sum2 += SQR(Px[i])*Py[i];
    Count += Py[i];
  }
  m1.Uno /= Count;
  m1.Due = sqrt(Sum2 - m1.Uno*NMax)/(double)(NMax-1);
  return m1;
}
PARABOLA Matematica::MinimoParabola(double a, double b,double *Px,double *Py,int NMass){
  PARABOLA Par;
  double Uno=0.;double Due=0.;double Tre=0.;double Quattro=0.;
  double UnoUno=0.;double DueUno=0.;double ZeroUno=0.;double ZeroDue=0.;
  for(int i=0;*(Px+i)<b || i<NMass;i++){
    //    printf("%.0f\t",*(Px+i));
    if(*(Px+i)>=a){
      Uno += *(Px+i);
      Due += QUAD(*(Px+i));
      Tre += *(Px+i)*QUAD(*(Px+i));
      Quattro += QUAD(QUAD(*Px+i));
      UnoUno += *(Px+i) * *(Py+i);
      DueUno += QUAD(*(Px+i)) * *(Py+i);
      ZeroUno += *(Py+i);
      ZeroDue += QUAD(*(Py+i));
      //      printf("%.2g\t%.2g\t%.2g\t%.2g\n",Uno,Due,Tre,Quattro);
    }
  }
  double X2 = NMass*Due - Uno*Uno;
  double X3 = NMass*Tre - Due*Uno;
  double X4 = NMass*Quattro - Due*Due;
  double XY = NMass *UnoUno - Uno*ZeroUno;
  double X2Y = NMass*DueUno - Due*ZeroUno;
  Par.a2 = (X4*X2 - X3*X3)/(X2Y*X2 - X3*XY);
  Par.a1 = (Par.a2*XY - X3)/X2;
  Par.a0 = (ZeroUno - Par.a1*Uno - Par.a2*Due)/(double)NMass;
  Par.Minimo = -Par.a1/(2*Par.a2);
  Par.MinimoY = Par.a2*SQR(Par.Minimo) + Par.a1*Par.Minimo + Par.a0;
  //printf("y = %.2f*x^2 + %.2f*x + %.2f\t a=0.075,b=-50.93,c=8623\n",Par.a2,Par.a1,Par.a0);
  return Par;
}
PARABOLA Matematica::MinimoParabola(double *Px,double *Py,int NMass){
  PARABOLA Par;
  double Uno=0.;double Due=0.;double Tre=0.;double Quattro=0.;
  double UnoUno=0.;double DueUno=0.;double ZeroUno=0.;double ZeroDue=0.;
  for(int i=0;i<NMass;i++){
    Uno += Px[i];
    Due += Px[i]*Px[i];
    Tre += Px[i]*Px[i]*Px[i];
    Quattro += Px[i]*Px[i]*Px[i]*Px[i];
    UnoUno += Px[i]*Py[i];
    DueUno += Px[i]*Px[i]*Py[i];
    ZeroUno += Py[i];
    ZeroDue += Py[i]*Py[i];
  }
  double X2 = NMass*Due - Uno*Uno;
  double X3 = NMass*Tre - Due*Uno;
  double X4 = NMass*Quattro - Due*Due;
  double XY = NMass *UnoUno - Uno*ZeroUno;
  double X2Y = NMass*DueUno - Due*ZeroUno;
  Par.a2 = (X2Y*X2 - X3*XY)/(X4*X2 - X3*X3);
  Par.a1 = (XY - Par.a2*X3)/X2;
  Par.a0 = (ZeroUno - Par.a1*Uno - Par.a2*Due)/(double)NMass;
  Par.Minimo = -Par.a1/(2.*Par.a2);
  Par.MinimoY = Par.a2*SQR(Par.Minimo) + Par.a1*Par.Minimo + Par.a0;
  //printf("y = %.2f*x^2 + %.2f*x + %.2f\t a=0.075,b=-50.93,c=8623\n",Par.a2,Par.a1,Par.a0);
  return Par;
}
double Matematica::Blend(const double *dPoint,double x,int nPoint,int nOrder){
  // if( ( x < dPoint[nPoint]) || ( x >= dPoint[nPoint+1]))
  //   return 0.;
  if(nOrder == 1){
    if( ( x >= dPoint[nPoint])&&( x < dPoint[nPoint+1]))
      return 1.;
    else 
      return 0.;
  }
  double Resp = 0.;
  double Primo = (x - dPoint[nPoint])/(dPoint[nPoint+nOrder] - dPoint[nPoint])*Blend(dPoint,x,nPoint,nOrder-1);
  double Secondo = (dPoint[nPoint+nOrder+1]-x)/(dPoint[nPoint+nOrder+1] - dPoint[nPoint+1])*Blend(dPoint,x,nPoint+1,nOrder-1);
  if( (dPoint[nPoint+nOrder-1] == dPoint[nPoint]) && (dPoint[nPoint+nOrder] == dPoint[nPoint+1]) )
    Resp = 0.;
  else if(dPoint[nPoint+nOrder] == dPoint[nPoint])
    Resp = Primo;
  else if(dPoint[nPoint+nOrder] == dPoint[nPoint+1])
    Resp = Secondo;
  else 
    Resp = Primo + Secondo;
  //printf("%d %d - (%lf-%lf) = %lf\n",nPoint,nOrder,x,dPoint[nPoint],Resp);
  return Resp;
}
double Matematica::Blend(double *dPoint,size_t Incr,double x,int nP,int nO){
/*********************************************************************

 Simple b-spline curve algorithm

 Copyright 1994 by Keith Vertanen (vertankd@cda.mrs.umn.edu)

 Released to the public domain (your mileage may vary)

**********************************************************************/
  double Resp;			
  if(nO == 1){
    if( (dPoint[nP*Incr] <= x)&&(x<dPoint[(nP+1)*Incr]))
      Resp = 1.;
    else 
      Resp = 0.;
  }
  else {
    if( (dPoint[(nP+nO-1)*Incr] == dPoint[nP*Incr]) && (dPoint[(nP+nO)*Incr] == dPoint[(nP+1)*Incr]) )
      Resp = 0.;
    else if(dPoint[(nP+nO)*Incr] = dPoint[(nP)*Incr])
      Resp = (dPoint[(nP+nO)*Incr] - x)/(dPoint[(nP+nO)*Incr] - dPoint[(nP+1)*Incr])*Blend(dPoint,Incr,x,nP,nO-1);
    else if(dPoint[(nP+nO)*Incr] == dPoint[(nP+1)*Incr])
      Resp = (x - dPoint[(nP)*Incr])/(dPoint[(nP+nO-1)*Incr] - dPoint[(nP)*Incr])*Blend(dPoint,Incr,x,nP+1,nO-1);
    else 
      Resp = (dPoint[(nP+nO)*Incr] - x)/(dPoint[(nP+nO)*Incr] - dPoint[(nP+1)*Incr])*Blend(dPoint,Incr,x,nP,nO-1) + (x - dPoint[nP*Incr])/(dPoint[(nP+nO-1)*Incr] - dPoint[nP*Incr])*Blend(dPoint,Incr,x,nP+1,nO-1);
  }
  //printf("%d %d - (%lf-%lf) %d = %lf\n",nP,nO,x,dPoint[nP],Incr,Resp);
  return Resp;
}
// int Matematica::InterBSpline2D(Matrice *MaIn,Matrice *MaOut){
//   double Max=1.;
//   double Min=0.;
//   int NOutShow=0;
//   double DeltaIn=(Max-Min)/(double)(NIn-1);
//   double DeltaOut=(Max-Min)/(double)(NOut-1);
//   int NOrder = 3+1;
//   double *dArray = (double *)calloc(NIn+NOrder+1,sizeof(double));
//   for(int p=0;p<=NIn+NOrder;p++){
//     if(p<NOrder){
//       dArray[p] = Min;
//     }
//     else if( (NOrder<=p) && (p<=NIn) ){
//       dArray[p] = (p-NOrder+1)*DeltaIn+Min;//Pm[p-NOrder].Pos[CLat1];//
//     }
//     else if( p>NIn){
//       dArray[p] = (p-NOrder+2)*DeltaIn+Min;
//     }
//   }
//   for(int vo=0;vo<NOut;vo++){
//    for(int vvo=0;vvo<NOut;vvo++){
//      MaOut->Set(vo,vvo,0.);
//      double x = DeltaOut*vo+Min;
//      //for(int vi=0;vi<NIn;vi++){
//      for(int vi=vo-1;vi<vo+NOrder+1;vi++){
//        if(vi < 0 || vi >= NIn) continue;
//        double Blendx = Mat->Blend(dArray,x,vi,NOrder);
//        double y = DeltaOut*vvo+Min;
//        //for(int vvi=0;vvi<NIn;vvi++){
//        for(int vvi=vvo-1;vvi<vvo+NOrder+1;vvi++){
// 	 if(vvi < 0 || vvi >= NIn) continue;
// 	 double Blendy = Mat->Blend(dArray,y,vvi,NOrder);
// 	 MaOut->Add(vo,vvo,Blendx*Blendy * PlIn->Val(vi,vvi));
//        }
//      }
//    }
//   }
//   for(int vo=0;vo<NOut;vo++){
//     if(vo < NIn){//To arrange
//       PlOut[vo][NOut-1] = PlIn[vo][NIn-1];
//       PlOut[NOut - 1][vo] = PlIn[NIn-1][vo];
//     }
//   }
//   NOutShow = NOut;
//   free(dArray);
//   return NOutShow;
// }
// 
