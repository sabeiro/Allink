#include "../include/Matematica.h"
#include <iostream>
#include <fstream>
void Matrice::Shout(const char * s, ...){
#ifdef MAT_DEBUG
  va_list args;
  va_start(args, s);
  fprintf(stderr, "Matrice] ");
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else
  return;
#endif
}
// Matrice::Matrice(int newNSize, int newNRow)  { // the only public ctor
//   if (newNSize <= 0) newNSize = 5;
//   NSize = newNSize; 
//   if ((newNRow <= newNSize)&&(newNRow>0))
//     NRow = newNRow;
//   else 
//     NRow = newNSize;
//   // since allocate() will first call delete[] on data:
//   //data = new double [NSize*NSize];
//   data = (double *)calloc(NSize*NSize,sizeof(double));
//   NCol = NRow = NSize;
//   //    allocate();
// }
//----------------Constructors------------------------
Matrice::Matrice(int ExtNRow, int ExtNCol){ 
  Shout("Matrice(row,col)");
  NRow = ExtNRow;
  NCol = ExtNCol;
  NZed = 0;
  NSize = NCol*NRow; 
  NDim = 2;
  data = (double *)calloc(NSize,sizeof(double));
}
Matrice::Matrice(int ExtNRow, int ExtNCol, double *pointer){ 
  Shout("Matrice(row,col)");
  NRow = ExtNRow;
  NCol = ExtNCol;
  NSize = NCol*NRow; 
  NDim = 2;
  data = (double *)malloc(sizeof(double));
  data = pointer;
}
Matrice::Matrice(int ExtNRow){ // the only public  data = (double *)calloc(NSize*NSize,sizeof(double));
  Shout("Matrice(row=col)");
  NRow = ExtNRow;
  NCol = 0;
  NZed = 0;
  NSize = NRow;
  NDim = 1;
  data = (double *)calloc(NSize,sizeof(double));
}
Matrice::Matrice(int ExtNRow,int ExtNCol,int ExtNZed){
  Shout("Matrice(row,col)");
  NRow = ExtNRow;
  NCol = ExtNCol;
  NZed = ExtNZed;
  NSize = NCol*NRow*NZed; 
  NDim = 3;
  data = (double *)calloc(NSize,sizeof(double));  
}
Matrice::Matrice(SPLINE Wg){
  Shout("Matrice(SPLINE)");
  NRow = NCol = 3;
  NZed = 0;
  NSize = NRow * NCol;
  data = (double *)calloc(NSize,sizeof(double));
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
  double Norm=0.;
  for(int r=0;r<NCol;r++){
    data[NRow*0+r] = MenoUno;
    data[NRow*1+r] = Zero;
    data[NRow*2+r] = PiuUno;
  }
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      //Norm += data[NSize*r+c] > 0. ? data[NSize*r+c] : 0;
      Norm += ASS(data[NRow*c+r]);
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      data[NRow*c+r] /= Norm*.5;
  NCol = NRow;
  NSize = NCol * NRow;
  NDim = 2;
}
Matrice::Matrice(double *Axis,double Angle,int ExtNRow){
  double c = cos(Angle);
  double s = sin(Angle);
  double t = 1. - c;
  double Norm = 0.;
  for(int d=0;d<3;d++){
    Norm += SQR(Axis[d]);
  }
  Norm = Norm > 0. ? sqrt(Norm) : 1.;
  double x = Axis[0]/Norm;
  double y = Axis[1]/Norm;
  double z = Axis[2]/Norm;
  NDim = 2;
  NZed = 0;
  if(ExtNRow == 4){
    NCol = NRow = 4;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0]  = t*x*x + c;
    data[NRow*0+1]  = t*x*y + z*s;
    data[NRow*0+2]  = t*x*z - y*s;
    data[NRow*0+3]  = 0.;
    
    data[NRow*1+0]  = t*x*y - z*s;
    data[NRow*1+1]  = t*y*y + c;
    data[NRow*1+2]  = t*y*z + x*s;
    data[NRow*1+3]  = 0.;
    
    data[NRow*2+0]  = t*x*z + y*s;
    data[NRow*2+1]  = t*y*z - x*s;
    data[NRow*2+2]  = t*z*z + c;
    data[NRow*2+3]  = 0.;
    
    data[NRow*3+0]  = 0.;
    data[NRow*3+1]  = 0.;
    data[NRow*3+2]  = 0.;
    data[NRow*3+3]  = 1.;
  }
  else{
    NCol = NRow = 3;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0]  = t*x*x + c;
    data[NRow*0+1]  = t*x*y + z*s;
    data[NRow*0+2]  = t*x*z - y*s;
    
    data[NRow*1+0]  = t*x*y - z*s;
    data[NRow*1+1]  = t*y*y + c;
    data[NRow*1+2]  = t*y*z + x*s;
    
    data[NRow*2+0]  = t*x*z + y*s;
    data[NRow*2+1]  = t*y*z - x*s;
    data[NRow*2+2]  = t*z*z + c;
  }  
}
Matrice::Matrice(Quadri q,int dim){
  //FIXME: the determinant is not zero!
  NDim = 2;
  NZed = 0;
  if(dim == 4){
    NCol = NRow = 4;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0]  = q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z;
    data[NRow*0+1]  = 2.*q.x*q.y + 2.*q.w*q.z;
    data[NRow*0+2]  = 2.*q.x*q.z - 2.*q.w*q.y;
    data[NRow*0+3]  = 0.;
    
    data[NRow*1+0]  = 2.*q.x*q.y - 2.*q.w*q.z;
    data[NRow*1+1]  = q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z;
    data[NRow*1+2]  = 2.*q.y*q.z + 2.*q.w*q.x;
    data[NRow*1+3]  = 0.;
    
    data[NRow*2+0]  = 2.*q.x*q.z + 2.*q.w*q.y;
    data[NRow*2+1]  = 2.*q.y*q.z - 2.*q.w*q.x;
    data[NRow*2+2]  = q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z;
    data[NRow*2+3]  = 0.;
    
    data[NRow*3+0]  = 0.;
    data[NRow*3+1]  = 0.;
    data[NRow*3+2]  = 0.;
    data[NRow*3+3]  = q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
  }
  else{
    NCol = NRow = 3;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0]  = 1. - 2.*SQR(q.y) - 2.*SQR(q.z);
    data[NRow*0+1]  = 2.*q.x*q.y + 2.*q.w*q.z;
    data[NRow*0+2]  = 2.*q.x*q.z - 2.*q.w*q.y;
    
    data[NRow*1+0]  = 2.*q.x*q.y - 2.*q.w*q.z;
    data[NRow*1+1]  = 1. - 2.*SQR(q.x) - 2.*SQR(q.z);
    data[NRow*1+2]  = 2.*q.y*q.z + 2.*q.w*q.x;
    
    data[NRow*2+0]  = 2.*q.x*q.z + 2.*q.w*q.y;
    data[NRow*2+1]  = 2.*q.y*q.z - 2.*q.w*q.x;
    data[NRow*2+2]  = 1. - 2.*SQR(q.x) - 2.*SQR(q.y);
  }  
}
Matrice::Matrice(double Roll,double Pitch,double Yaw,int ExtNRow){
  //Defining Rrpy(r,p,y) := RzrRypRxy
  double CR = cos(Roll);
  double SR = sin(Roll);
  double CP = cos(Pitch);
  double SP = sin(Pitch);
  double CY = cos(Yaw);
  double SY = sin(Yaw);
  NDim = 2;
  NZed = 0;
  if(ExtNRow == 4){
    NCol = NRow = 4;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0] = CY*CP;
    data[NRow*0+1] = SY*CP;
    data[NRow*0+2] = -SP;
    data[NRow*0+3] = 0.;

    data[NRow*1+0] = CY*SP*SR - SY*CR;
    data[NRow*1+1] = SY*SP*SR + CY*CR;
    data[NRow*1+2] = CP*SR;
    data[NRow*1+3] = 0.;
    
    data[NRow*2+0] = CY*SP*CR + SY*SR;
    data[NRow*2+1] = SY*SP*CR - CY*SR;
    data[NRow*2+2] = CP*CR;
    data[NRow*2+3] = 0.;

    data[NRow*3+0] = 0.;
    data[NRow*3+1] = 0.;
    data[NRow*3+2] = 0.;
    data[NRow*3+3] = 1.;
    // Other definition of the angles (to check)
    // data[NRow*0+0] = CR*CP*CY - SR*SY;
    // data[NRow*0+1] = CR*CP*SY - SR*CY;
    // data[NRow*0+2] = -CR*SP;
    // data[NRow*0+3] = 0.;

    // data[NRow*1+0] = -SR*CP*CY - CR*SY;
    // data[NRow*1+1] = -SR*CP*SY + CR*CY;
    // data[NRow*1+2] = SR*SP;
    // data[NRow*1+3] = 0.;
    
    // data[NRow*2+0] = SP*CY;
    // data[NRow*2+1] = SP*SY;
    // data[NRow*2+2] = CP;
    // data[NRow*2+3] = 0.;

    // data[NRow*3+0] = 0.;
    // data[NRow*3+1] = 0.;
    // data[NRow*3+2] = 0.;
    // data[NRow*3+3] = 1.;
  }
  if(ExtNRow == 3){
    NCol = NRow = 3;
    NSize = NCol * NRow;
    data = (double *)calloc(NSize,sizeof(double));
    data[NRow*0+0] = CR*CP;
    data[NRow*0+1] = SR*CP;
    data[NRow*0+2] = -SP;

    data[NRow*1+0] = CR*SP*SY - SR*CY;
    data[NRow*1+1] = SR*SP*SY + CR*SY;
    data[NRow*1+2] = CP*SY;
    
    data[NRow*2+0] = CR*SP*CY+SR*SY;
    data[NRow*2+1] = SR*SP*CY-CR*SY;
    data[NRow*2+2] = CR*CY;
  }
}
Matrice::Matrice(double *M,int Row,int Col){
  NRow = Row;
  NCol = Col;
  NZed = 0;
  NSize = NRow*NCol;
  //data = M;
  NDim = 2;
  data = (double *)calloc(NSize,sizeof(double));
  for(int c=0;c<NCol;c++)
    for(int r=0;r<NRow;r++)
      data[NRow*c+r] = M[NRow*c+r];
}
Matrice::Matrice(SPLINE Wg,int Dim){
  Shout("Matrice(SPLINE,dim)");
  NCol = NRow = 5;
  NZed = 0;
  NSize = NCol * NRow;
  data = (double *)calloc(NSize,sizeof(double));
  FillDiffOperator(Wg,Dim);
}
Matrice::~Matrice() { 
  //  delete[] data; 
  free(data);
}
void Matrice::RandomFill(double Max){
  Shout("RandomFill");
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      data[NRow*c+r] = Max*drand48();//Casuale();
}
void Matrice::FillDiffOperator(SPLINE Wg,int Dim){
  Shout("Matrice(SPLINE,dim)");
  double MenoDue = - .75*Wg.a3 + 2.*Wg.a4;//1.5	
  //MenoDue += .125*Wg.a1 - .125*Wg.a2; //O(h^4)
  MenoDue -= 1./12.*Wg.a2; //O(h^4)
  double MenoUno = -.5*Wg.a1 + Wg.a2 + 1.5*Wg.a3 - 8.*Wg.a4;//6
  //MenoUno += -.125*Wg.a1 - .5*Wg.a2;
  MenoUno += 1./3.*Wg.a2; //O(h^4)
  double Zero = Wg.a0 - 2.*Wg.a2 + 12.*Wg.a4;//9
  //Zero += 0.75*Wg.a2;
  Zero -= 0.5*Wg.a2; //O(h^4)
  double PiuUno = .5*Wg.a1 + Wg.a2 - 1.5*Wg.a3 - 8.*Wg.a4;
  //PiuUno += .25*Wg.a1 + .5*Wg.a2;
  PiuUno += 1./3.*Wg.a2;
  double PiuDue = .75*Wg.a3 + 2.*Wg.a4;
  //PiuDue += -.125*Wg.a1 -.125*Wg.a2;
  PiuDue -= 1./12.*Wg.a2; //O(h^4)
  double Norm=0.;
  if (Dim == 0){
    for(int r=0;r<NCol;r++){
      data[NRow*0+r] = MenoDue;
      data[NRow*1+r] = MenoUno;
      data[NRow*2+r] = Zero;
      data[NRow*3+r] = PiuUno;
      data[NRow*4+r] = PiuDue;
    }
  }      
  else if(Dim == 1){
    int r=2;
    data[NRow*0+r] = MenoDue;
    data[NRow*1+r] = MenoUno;
    data[NRow*2+r] = Zero;
    data[NRow*3+r] = PiuUno;
    data[NRow*4+r] = PiuDue;
  }
  else if (Dim == 2){
    int r=2;
    data[NRow*0+r] = MenoDue;
    data[NRow*1+r] = MenoUno;
    data[NRow*2+r] = Zero;
    data[NRow*3+r] = PiuUno;
    data[NRow*4+r] = PiuDue;
    int c=2;
    data[NRow*c+0] = MenoDue;
    data[NRow*c+1] = MenoUno;
    data[NRow*c+2] += Zero;
    data[NRow*c+3] = PiuUno;
    data[NRow*c+4] = PiuDue;
  }
  else if (Dim == 3){
    for(int r=0;r<NCol;r++){
      data[NRow*0+r] = MenoDue;
      data[NRow*1+r] = MenoUno;
      data[NRow*2+r] = Zero;
      data[NRow*3+r] = PiuUno;
      data[NRow*4+r] = PiuDue;
    }
  }
}
void Matrice::FillCanny(){
  if(NRow != 5){
    printf("I don't know the expresion for matrices of NRow != 5\n");
    return ;
  }
  if(NCol != 5){
    printf("I don't know the expresion for matrices of NRow != 5\n");
    return ;
  }
  Set(0,0,2.); Set(0,1,4.); Set(0,2,5.); Set(0,3,4.); Set(0,4,2.);
  Set(1,0,4.); Set(1,1,9.); Set(1,2,12.); Set(1,3,9.); Set(1,4,4.);
  Set(2,0,5.); Set(2,1,12.); Set(2,2,15.); Set(2,3,12.); Set(2,4,5.);
  Set(3,0,4.); Set(3,1,9.); Set(3,2,12.); Set(3,3,9.); Set(3,4,4.);
  Set(4,0,2.); Set(4,1,4.); Set(4,2,5.); Set(4,3,4.); Set(4,4,2.);
  Multiply(1./159.);
}
void Matrice::FillGaussian(double Sigma,double CutOff){
  int Half = NRow/2;
  double Norm = 0.;
  if(NDim == 3){
    for(int r=0;r<NRow;r++){
      double x = CutOff*(r-Half)/(double)NRow;
      for(int c=0;c<NCol;c++){
	double y = CutOff*(c-Half)/(double)NCol;
	for(int q=0;q<NCol;q++){
	  double z = CutOff*(q-Half)/(double)NCol;
	  double r2 = SQR(x) + SQR(y) + SQR(z);
	  double Gauss = exp(-r2*.5/SQR(Sigma));
	  data[(NRow*c+r)*NCol+q] = Gauss;
	  Norm += Gauss;
	}
      }
    }
  }
  else if(NDim == 2){
    for(int r=0;r<NRow;r++){
      double x = CutOff*(r-Half)/(double)NRow;
      for(int c=0;c<NCol;c++){
	double y = CutOff*(c-Half)/(double)NCol;
	double r2 = SQR(x) + SQR(y);
	double Gauss = exp(-r2*.5/SQR(Sigma));
	data[NRow*c+r] = Gauss;
	Norm += Gauss;
      }
    }
  }
  else if(NDim == 1){
    for(int r=0;r<NRow;r++){
      double x = CutOff*(r-Half)/(double)NRow;
      double r2 = SQR(x);
      double Gauss = exp(-r2*.5/SQR(Sigma));
      data[r] = Gauss;
      Norm += Gauss;
    }
  }
  Multiply(1./Norm);
}
void Matrice::FillGaussian5(){
  if(NRow != 5 && NCol != 5){
    printf("Invalid number of rows and cols %dx%d != 7x7\n",NRow,NCol);
    return;
  }
  data[NRow*0+0] = 1;
  data[NRow*0+1] = 4;
  data[NRow*0+2] = 7;
  data[NRow*0+3] = 4;
  data[NRow*0+4] = 1;
  
  data[NRow*1+0] = 4;
  data[NRow*1+1] = 16;
  data[NRow*1+2] = 26;
  data[NRow*1+3] = 16;
  data[NRow*1+4] = 4;
  
  data[NRow*2+0] = 7;
  data[NRow*2+1] = 26;
  data[NRow*2+2] = 41;
  data[NRow*2+3] = 26;
  data[NRow*2+4] = 7;

  data[NRow*3+0] = 4;
  data[NRow*3+1] = 16;
  data[NRow*3+2] = 26;
  data[NRow*3+3] = 16;
  data[NRow*3+4] = 4;

  data[NRow*4+0] = 1;
  data[NRow*4+1] = 4;
  data[NRow*4+2] = 7;
  data[NRow*4+3] = 4;
  data[NRow*4+4] = 1;

  Multiply(1./273.);
}


//---------------Manage-entries-------------------
double Matrice ::Val(int row){
#ifdef MATR_DEBUG
  if ( (row>=NRow) || (row<0) )
    {  printf("Values %d %d out of range %d %d\n",row,NRow);
      return 1.;    }
#endif
  return data[row];
} 
double Matrice ::Val(int row,int col){
#ifdef MATR_DEBUG
  if ( (row>=NRow) || (col>=NCol) 
       || (row<0) || (col<0) )
    {  printf("Values %d %d out of range %d %d\n",row,col,NRow,NCol);
      return 1.;    }
#endif
  return data[ NRow*col+row ];
} 
double Matrice ::Val(int row,int col,int zed){
  return data[ (NRow*col+row)*NCol+zed ];
} 
bool Matrice ::Set(int row, int column, double newvalue)  {
#ifdef MATR_DEBUG
  if ( (row >= NRow) || (column >= NCol) 
       || (row<0) || (column<0) ) return false;
#endif
  data[ NRow*column + row ] = newvalue;
  return true;
}
bool Matrice ::Add(int row, int column, double Value)  {
#ifdef MATR_DEBUG
  if ( (row >= NRow) || (column >= NCol) 
       || (row<0) || (column<0) ) return false;
#endif
  data[ NRow*column + row ] += Value;
  return true;
}
void Matrice ::Print(){
  if(NDim == 3){
    for(int q=0;q<NZed;q++){
      printf("%d)\n",q);
      for(int r=0;r<NRow;r++){
	printf("|");
	for(int c=0;c<NCol;c++){
	  printf("%.2g ",data[(NRow*c+r)*NCol+q]);
	}
	printf("|\n");
      }
      printf("\n");
    }
    return;
  }
  else if(NDim == 1){
    printf("|");
    for(int r=0;r<NRow;r++){
      printf("%.2g ",data[r]);
    }
    printf("|\n");
    printf("\n");
    return;
  }
  for(int r=0;r<NRow;r++){
    printf("|");
    for(int c=0;c<NCol;c++){
      printf("%.2g ",data[NRow*c+r]);
    }
    printf("|\n");
  }
  printf("\n");
}
//---------------inversion-(-stolen,-to-be-checked-)---------
void Matrice ::comparetoidentity()  {
  int worstdiagonal = 0;
  double maxunitydeviation = 0.0;
  double currentunitydeviation;
  for ( int i = 0; i < NRow; i++ )  {
    currentunitydeviation = data[i+i*NRow] - 1.;
    if ( currentunitydeviation < 0.0) currentunitydeviation *= -1.;
    if ( currentunitydeviation > maxunitydeviation )  {
      maxunitydeviation = currentunitydeviation;
      worstdiagonal = i;
    }
  }
  int worstoffdiagonalrow = 0;
  int worstoffdiagonalcolumn = 0;
  double maxzerodeviation = 0.0;
  double currentzerodeviation ;
  for ( int i = 0; i < NRow; i++ )  {
    for ( int j = 0; j < NRow; j++ )  {
      if ( i == j ) continue;  // we look only at non-diagonal terms
      currentzerodeviation = data[i+j*NRow];
      if ( currentzerodeviation < 0.0) currentzerodeviation *= -1.0;
      if ( currentzerodeviation > maxzerodeviation )  {
	maxzerodeviation = currentzerodeviation;
	worstoffdiagonalrow = i;
	worstoffdiagonalcolumn = j;
      }
      
    }
  }
  cout << "Worst diagonal value deviation from unity: " 
       << maxunitydeviation << " at row/column " << worstdiagonal << endl;
  cout << "Worst off-diagonal value deviation from zero: " 
       << maxzerodeviation << " at row = " << worstoffdiagonalrow 
       << ", column = " << worstoffdiagonalcolumn << endl;
}
void Matrice ::settoproduct(Matrice& left, Matrice& right)  {
//   NRow = left.getNRow();
//   if ( NSize < left.getNRow() )   {
//     NSize = left.getNRow();
//     allocate();
//   }
  for ( int i = 0; i < NRow; i++ )
    for ( int j = 0; j < NRow; j++ )  {
      double sum = 0.0;
      double leftvalue, rightvalue;
      bool success;
      for (int c = 0; c < NRow; c++)  {
	left.getvalue(i,c,leftvalue,success);
	right.getvalue(c,j,rightvalue,success);
	sum += leftvalue * rightvalue;
      }
      Set(i,j,sum);
    }
}
void Matrice ::copymatrix(Matrice&  source)  {
//   NRow = source.getNRow();
//   if ( NSize < source.getNRow() )  {
//     NSize = source.getNRow();
//     allocate();
//   }
  for ( int i = 0; i < NRow; i++ )
    for ( int j = 0; j < NRow; j++ )  {
      double value;
      bool success;
      source.getvalue(i,j,value,success);
      data[i+j*NRow] = value;
    }
}
void Matrice ::setNRow(int newNRow) {
//   if ( newNRow > NSize )
//     {
//       NSize = newNRow ; // * 2;  // wastes memory but saves
//       // time otherwise required for
//       // operation new[]
//       allocate();
//     }
//   if (newNRow >= 0) NRow = newNRow;
}
int Matrice ::getNRow() { return NRow; }
void Matrice ::getvalue(int row, int column, double& returnvalue, bool& success)   {
#ifdef MATR_DEBUG
  if ( (row>=NRow) || (column>=NCol) 
       || (row<0) || (column<0) )
    {  success = false;
      return;    }
#endif
  returnvalue = data[ row + column*NRow ];
  success = true;
}
void Matrice::Invert()  {
#ifdef MATR_DEBUG //CHECKIT
  if (NRow <= 0) return;  // sanity check
  if (NRow == 1) return;  // must be of dimension >= 2
#endif
  for (int i=1; i < NRow; i++) data[i] /= data[0]; // normalize row 0
  for (int i=1; i < NRow; i++)  { 
    for (int j=i; j < NRow; j++)  { // do a column of L
      double sum = 0.0;
      for (int k = 0; k < i; k++)  
	sum += data[j+k*NRow] * data[k+i*NRow];
      data[j+i*NRow] -= sum;
    }
    if (i == NRow-1) continue;
    for (int j=i+1; j < NRow; j++)  {  // do a row of U
      double sum = 0.0;
      for (int k = 0; k < i; k++)
	sum += data[i+k*NRow]*data[k+j*NRow];
      data[i+j*NRow] = 
	(data[i+j*NRow]-sum) / data[i+i*NRow];
    }
  }
  for ( int i = 0; i < NRow; i++ )  // invert L
    for ( int j = i; j < NRow; j++ )  {
      double x = 1.0;
      if ( i != j ) {
	x = 0.0;
	for ( int k = i; k < j; k++ ) 
	  x -= data[j+k*NRow]*data[k+i*NRow];
      }
      data[j+i*NRow] = x / data[j+j*NRow];
    }
  for ( int i = 0; i < NRow; i++ )   // invert U
    for ( int j = i; j < NRow; j++ )  {
      if ( i == j ) continue;
      double sum = 0.0;
      for ( int k = i; k < j; k++ )
	sum += data[k+j*NRow]*( (i==k) ? 1.0 : data[i+k*NRow] );
      data[i+j*NRow] = -sum;
    }
  for ( int i = 0; i < NRow; i++ )   // final inversion
    for ( int j = 0; j < NRow; j++ )  {
      double sum = 0.0;
      for ( int k = ((i>j)?i:j); k < NRow; k++ )  
	sum += ((j==k)?1.0:data[j+k*NRow])*data[k+i*NRow];
      data[j+i*NRow] = sum;
    }
};
//-------------Other-operations-----------------
int Matrice ::Solve(double *Known,double *UnKnown){
  Shout("Solve");
#ifdef MATR_DEBUG
#endif
  Matrice Lower(NCol,NCol);
  Matrice Upper(NCol,NCol);
  double *AlmostKnown = (double *)calloc(NCol,sizeof(double));
  if(!AlmostKnown) {printf("Not allocated\n");return 0;}
  /* assegna i valori della diagonale di L */
  for(int r = 0; r < NRow ; r ++){
    Lower.Set(r,r,1.);
  }
  /* calcola gli elementi fuori della diagonale */
  for(int j = 0; j < NCol ; j ++){
    for(int i = 0; i <= j ; i++){
      Upper.Set(i,j,Val(i,j));
      for (int k = 0; k <= i - 1; k++) {
	Upper.Add(i,j,-Lower.Val(i,k)*Upper.Val(k,j));
      }
    }
    for (int  i = j + 1; i < NCol ; i++) {
      Lower.Set(i,j,Val(i,j));
      for (int k = 0; k <= j - 1; k ++) {
	Lower.Add(i,j,-Lower.Val(i,k)*Upper.Val(k,j));
      }//FIXME
      double Diag = POS(Upper.Val(j,j)) > 0. ? Lower.Val(i,j)/Upper.Val(j,j) : 1.;
      Lower.Set(i,j,Diag);
    }
  }//Solve the Uz = y
  //Upper = NULL;
  for(int r=0;r<NRow;r++){
    AlmostKnown[r] = Known[r];
    for(int c=r-1;c>=0;c--){
      AlmostKnown[r] -= AlmostKnown[c]*Lower.Val(r,c);
      //printf("Coeff[%d][%d]=%lf %lf\n",r,c,UnKnown[r],Lower.Val(r,c));
    }
    if(POS(Lower.Val(r,r))> 0.)
      AlmostKnown[r] /=  POS(Lower.Val(r,r))>0. ? Lower.Val(r,r) : 1.;
    else
      AlmostKnown[r] = 0.;
    //printf("Coeff[%d][%d]=%lf %lf\n",r,r,AlmostKnown[r],Lower.Val(r,r));
  }
  //Lower.Print();
  for(int r=NRow-1;r>=0;r--){
    UnKnown[r]  = AlmostKnown[r];
   for(int c=r+1;c<NCol;c++){
     UnKnown[r] -= UnKnown[c]*Upper.Val(r,c);
     // printf("Coeff[%d][%d]=%lf %lf\n",r,c,AlmostKnown[r],Upper.Val(r,c));
    }
   if(POS(Upper.Val(r,r))> 0.)
      UnKnown[r] /=  Upper.Val(r,r);
    else
      UnKnown[r] = 0.;
   //printf("Coeff[%d][%d]=%lf %lf\n",r,r,UnKnown[r],Upper.Val(r,r));
  }//finish solving Lx = z
  //Upper.Print();
  free(AlmostKnown);
  return 1;
};
void Matrice::Apply(double *Known,double *UnKnown){
  Shout("Apply");
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      UnKnown[r] += Known[c]*data[r+c*NRow];
    }
  }
};
double Matrice::Det(){
  Shout("Det");
  Matrice Lower(NCol);
  Matrice Upper(NCol);
  /* assegna i valori della diagonale di L */
  for (int r = 0; r < NCol ; r ++) {
    Lower.Set(r,r,1.);
  }
  /* calcola gli elementi fuori della diagonale */
  for (int j = 0; j < NCol ; j ++) {
    for (int i = 0; i <= j ; i++) {
      Upper.Set(i,j,Val(i,j));
      for (int k = 0; k <= i - 1; k++) {
	Upper.Add(i,j,-Lower.Val(i,k)*Upper.Val(k,j));
      }
    }
    for (int  i = j + 1; i < NCol ; i++) {
      Lower.Set(i,j,Val(i,j));
      for (int k = 0; k <= j - 1; k ++) {
	Lower.Add(i,j,-Lower.Val(i,k)*Upper.Val(k,j));
      }
      Lower.Set(i,j,Lower.Val(i,j)/Upper.Val(j,j));
    }
  }
  double Det = 1.;
  for(int r=0;r<NCol;r++)
    Det *= Upper.Val(r,r);
  return Det;
};
void Matrice::Transpose(){
  Shout("Transpose");
  double *Temp = (double *)calloc(NCol*NRow,sizeof(double));
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      Temp[r+c*NRow] = data[r+c*NRow];
    }
  }
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      data[r+c*NRow] = Temp[c+r*NRow];
    }
  }
  free(Temp);
}
void Matrice::Clear(){
  Shout("Clear");
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      data[r+c*NRow] = 0.;
    }
  }
}
void Matrice::Normalize(){
  Shout("Normalize");
  double Max=0.;
  double Min=0.;
  for(int r=0;r<NRow;r++){
    for(int c=0;c<NCol;c++){
      if(Max < data[r+c*NRow])
	Max = data[r+c*NRow];
      if(Min > data[r+c*NRow])
	Min = data[r+c*NRow];
    }
  }
  if(fabs(Min) > 0.0e-5){
    for(int r=0;r<NRow;r++)
      for(int c=0;c<NCol;c++)
	data[r+c*NRow] = data[r+c*NRow]/(Max-Min)-Min;
 
  }
  else {
    for(int r=0;r<NRow;r++)
      for(int c=0;c<NCol;c++)
	data[r+c*NRow] = data[r+c*NRow]/Max;
  }
}  
void Matrice::Multiply(double Val){
  for(int s=0;s<NSize;s++)
    data[s] *= Val;
  // for(int r=0;r<NRow;r++)
  //   for(int c=0;c<NCol;c++)
  //     data[r+c*NRow] = data[r+c*NRow]*Val;
}
void Matrice::CopyOn(Matrice *B){
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      B->Set(r,c,data[r+c*NRow]);
}
Matrice Matrice::operator+ (Matrice& A){
  if(NCol != A.NCol || NRow != A.NCol){
    printf("Incompatible matrices. Dim: %d/%d\n",NSize,A.NSize); 
    return 0;
  }
  Matrice Resp(NRow,NCol);
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      Resp.Set(r,c,Val(r,c)+A.Val(r,c));
  return Resp;
}
Matrice Matrice::operator* (Matrice &A){
  if(NCol != A.NRow || NRow != A.NCol){
    printf("Incompatible matrices. Dim: %d/%d\n",NSize,A.NSize); 
    return 0;
  }
  Matrice Resp(NRow,A.NCol);
  for(int r=0;r<NRow;r++)
    for(int c=0;c<A.NCol;c++){
      double Temp=0.;
      for(int i=0;i<NRow;i++)
	//Temp += Val(r,i)*A.Val(i,c);
	Temp += Resp.data[r+i*NRow]*A.data[i+c*NRow];
      Resp.Set(r,c,Temp);
    }
  Resp.Print();
  return Resp;
}
Matrice Matrice::operator= (Matrice& A){
  if(NRow != A.NRow || NCol != A.NCol){
    printf("Incompatible matrices. Dim: %d/%d\n",NSize,A.NSize); 
    return 0;
  }
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      Set(r,c,A.Val(r,c));
  return *this;
}
Matrice Matrice::operator* (const double& Molt)const{
  for(int r=0;r<NRow;r++)
    for(int c=0;c<NCol;c++)
      data[r+c*NRow] *= Molt;
  return *this;
}
Matrice Matrice::operator^(Matrice& A)const{
  assert(  NCol == A.NRow || NRow == A.NCol);
  printf("Still to be written!!!\n");
  Matrice B(NCol,NRow);
  return B;
}
Vettore Matrice::operator* (const Vettore& u)const {
  if(NCol != u.NDim){
    printf("Incompatible matrices. Dim: %d/%d\n",NSize,u.NDim); 
    return 0;
  }
  Vettore Resp(NRow);
  for(int r=0;r<NRow;r++){
      double Temp=0.;
    for(int c=0;c<NCol;c++)
      Temp += data[r+c*NRow]*u.x[c];
    Resp.x[r] = Temp;
  }
  return Resp;
}
void Matrice::Mult(Matrice &A,Matrice &B){
  if(B.NCol != A.NRow || B.NRow != A.NCol){
    printf("Incompatible matrices. Dim: %d/%d\n",B.NSize,A.NSize); 
    return ;
  }
  for(int r=0;r<A.NRow;r++)
    for(int c=0;c<B.NCol;c++){
      double Temp=0.;
      for(int i=0;i<A.NRow;i++)
	Temp += A.Val(r,i)*B.Val(i,c);
      Set(r,c,Temp);
    }
}
void Matrice::Mult(Matrice &A){
  if(NCol != A.NRow || NRow != A.NCol){
    printf("Incompatible matrices. Dim: %d/%d\n",NSize,A.NSize); 
    return ;
  }
  for(int r=0;r<A.NRow;r++)
    for(int c=0;c<NCol;c++){
      double Temp=0.;
      for(int i=0;i<A.NRow;i++)
	Temp += A.Val(r,i)*Val(i,c);
      Set(r,c,Temp);
    }
}
Vettore Matrice::Mult(Matrice &A,Vettore &v){
  Vettore Resp(NRow);
  for(int r=0;r<NRow;r++){
    double Temp=0.;
    for(int c=0;c<NCol;c++)
      Temp += A.data[r+c*NRow]*v.x[c];
    Resp.x[r] = Temp;
  }
  return Resp;
}
Vettore Matrice::Mult(Vettore &v){
  Vettore Resp(NRow);
  for(int r=0;r<NRow;r++){
    double Temp=0.;
    for(int c=0;c<NCol;c++)
      Temp += data[r+c*NRow]*v.x[c];
    Resp.x[r] = Temp;
  }
  return Resp;
}
void Matrice::Mult(Vettore &v,Vettore &u){
  for(int r=0;r<NRow;r++){
    double Temp=0.;
    for(int c=0;c<NCol;c++)
      Temp += data[r+c*NRow]*v.x[c];
    u.x[r] = Temp;
  }
}
void Matrice::ConvoluteMatrix(double *Plot,int NGrid,int NDim,int IfMinImConv){
  if(NDim == 1){
    if(!IfMinImConv) ConvoluteMatrix1(Plot,NGrid);
    else ConvoluteMatrix1MinImConv(Plot,NGrid);
  }
  else if(NDim == 2){
    if(!IfMinImConv) ConvoluteMatrix2(Plot,NGrid);
    else ConvoluteMatrix2MinImConv(Plot,NGrid);
  }
  else if(NDim == 3){
    ConvoluteMatrix3(Plot,NGrid);
  }
}
//without minimum image convention
void Matrice::ConvoluteMatrix1(double *Plot,int NGrid){
  double *Plot2 = (double *)calloc(NGrid,sizeof(double));
  int NMat2 = (int)pNRow()/2;
  for(int gx=0;gx<NGrid;gx++){
    if(gx <= NMat2  || gx >= NGrid - NMat2){
      Plot2[gx] = Plot[gx];
      continue;
    }
    for(int mx=0;mx<pNRow();mx++){
      int g1x = gx + mx - NMat2;
      if(g1x >= NGrid){
	Plot2[gx] = Plot[gx];
	break;
      }	
      if(g1x < 0){
	Plot2[gx] = Plot[gx];
	break;
      }	
      Plot2[gx] += Plot[g1x]*Val(mx);
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    Plot[gx] = Plot2[gx];
  }
  free(Plot2);
}
//with minimum image convention
void Matrice::ConvoluteMatrix1MinImConv(double *Plot,int NGrid){
  double *Plot2 = (double *)calloc(NGrid,sizeof(double));
  int NMat2 = (int)pNRow()/2;
  for(int gx=0;gx<NGrid;gx++){
    if(gx <= NMat2  || gx >= NGrid - NMat2){
      Plot2[gx] = Plot[gx];
      continue;
    }
    for(int mx=0;mx<pNRow();mx++){
      int g1x = gx + mx - NMat2;
      if(g1x >= NGrid) g1x -= NGrid;
      if(g1x < 0) g1x + NGrid;
      Plot2[gx] += Plot[g1x]*Val(mx);
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    Plot[gx] = Plot2[gx];
  }
  free(Plot2);
}
///no minimum image convention
void Matrice::ConvoluteMatrix2(double *Plot,int NGrid){
  double *Plot2 = (double *)calloc(SQR(NGrid),sizeof(double));
  int NMat2 = (int)pNCol()/2;
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      for(int mx=0;mx<pNRow();mx++){
	int g1x = gx + mx - NMat2;
	if(g1x >= NGrid){
	  Plot2[gx*NGrid+gy] = Plot[gx*NGrid+gy];
	  break;
	}
	if(g1x < 0) {
	  Plot2[gx*NGrid+gy] = Plot[gx*NGrid+gy];
	  break;
	}
	for(int my=0;my<pNCol();my++){
	  int g1y = gy + my - NMat2;
	  if(g1y >= NGrid){
	    Plot2[gx*NGrid+gy] = Plot[gx*NGrid+gy];
	    break;
	  } 
	  if(g1y < 0){
	    Plot2[gx*NGrid+gy] = Plot[gx*NGrid+gy];
	    break;
	  }
	  Plot2[gx*NGrid+gy] += Plot[g1x*NGrid+g1y]*Val(mx,my);
	}
      }
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot[gx*NGrid+gy] = Plot2[gx*NGrid+gy];
    }
  }
  free(Plot2);
}
///with minimum image convention
void Matrice::ConvoluteMatrix2MinImConv(double *Plot,int NGrid){
  double *Plot2 = (double *)calloc(SQR(NGrid),sizeof(double));
  int NMat2 = (int)pNCol()/2;
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      for(int mx=0;mx<pNRow();mx++){
	int g1x = gx + mx - NMat2;
	if(g1x >= NGrid) continue;//g1x -= NGrid;
	if(g1x < 0) continue;//g1x + NGrid;
	for(int my=0;my<pNCol();my++){
	  int g1y = gy + my - NMat2;
	  if(g1y >= NGrid) continue;//g1y -= NGrid;
	  if(g1y < 0) continue;//g1y + NGrid;
	  Plot2[gx*NGrid+gy] += Plot[g1x*NGrid+g1y]*Val(mx,my);
	}
      }
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      Plot[gx*NGrid+gy] = Plot2[gx*NGrid+gy];
    }
  }
  free(Plot2);
}
//with minimum image convention
void Matrice::ConvoluteMatrix3(double *Plot,int NGrid){
  double *Plot2 = (double *)calloc(CUBE(NGrid),sizeof(double));
  int NMat2 = (int)pNCol()/2;
  printf("NMat2 %d\n",NMat2);
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      for(int gz=0;gz<NGrid;gz++){
	for(int mx=0;mx<pNRow();mx++){
	  int g1x = gx + mx - NMat2;
	  if(g1x >= NGrid) g1x -= NGrid;
	  if(g1x < 0) g1x += NGrid;
	  for(int my=0;my<pNCol();my++){
	    int g1y = gy + my - NMat2;
	    if(g1y >= NGrid) g1y -= NGrid;
	    if(g1y < 0) g1y += NGrid;
	    for(int mz=0;mz<pNZed();mz++){
	      int g1z = gz + mz - NMat2;
	      if(g1z >= NGrid) g1z -= NGrid;
	      if(g1z < 0) g1z += NGrid;
	      Plot2[(gx*NGrid+gy)*NGrid+gz] += Plot[(g1x*NGrid+g1y)*NGrid+g1z]*Val(mx,my,mz);
	    }
	  }
	}
      }
    }
  }
  for(int gx=0;gx<NGrid;gx++){
    for(int gy=0;gy<NGrid;gy++){
      for(int gz=0;gz<NGrid;gz++){
	Plot[(gx*NGrid+gy)*NGrid+gz] = Plot2[(gx*NGrid+gy)*NGrid+gz];
      }
    }
  }
  free(Plot2);
}
