#ifndef MATEMATICASTRUCT_H
#define MATEMATICASTRUCT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/// Parabolas coefficients
typedef struct {
  /// a0 + a1 x + a2 x^2
  double a0;
  /// Error on the first coefficient
  double Erra0;
  /// a0 + a1 x + a2 x^2
  double a1;
  /// Error on the second coefficient
  double Erra1;
  /// a0 + a1 x + a2 x^2
  double a2;
  /// Error on the third coefficient
  double Erra2;
  /// Minimum of the parabola
  double Minimo;
  /// Minimum of the parabola
  double MinimoY;
}PARABOLA;
/// Moments of a distribution
typedef struct {
  /// First moment
  double Uno;
  /// Second moment
  double Due;
  /// Third moment
  double Tre;
  /// Step interval
  double Delta;
  /// Chi square
  double Chi;
  /// Minimum value
  double Min;
  /// Maximum value
  double Max;
  /// Minimum of the ordinate
  double yMin;
  /// Maximum of the ordinate
  double yMax;
  /// Number of points considered
  int Num;
} MOMENTI;
/// Radius and center of a circle
typedef struct {
  /// Radius of the circle
  double Rad;
  /// x position of the center
  double xC;
  /// y position of the center
  double yC;}CIRCLE;
/// Linear interpolation
typedef struct {
  /// y = m*x + q
  double m;
  /// Error on the slope
  double ErrM;
  /// y = m*x + q
  double q;
  /// Error on the intercept
  double ErrQ;
  /// Correlation
  double Corr;
  /// Covariance
  double Cov;
  /// r factor
  double r;
  /// Error a posteriori
  double ErrY;
}RETTA;
/// Three dimentional vector
typedef struct {
  ///Position
  double x[3];
}VETT;
/// Four dimentional vector/quaternion
typedef struct {
  /// w,x,y,z
  double x[4];
}QUADRI;
/// Where a root was searched
typedef struct {
  /// Iferior limit
  double iLim;
  /// Superior limit
  double sLim;
  /// Point of the zero
  double Zero;
  /// If the zero was founds
  int IfRis;//Se `e andato a buon fine
} RADICE;
/// Indices for the permutation
typedef struct {
  /// First number of the couple 
  int n;
  /// Second number of the couple
  int m;
}PERMUTE;
/// Coefficient of a spline
typedef struct {
  /// a0 + a1*x + a2*x^2 + a3*x^3 + a4^4
  double a0;
  /// a0 + a1*x + a2*x^2 + a3*x^3 + a4^4
  double a1;
  /// a0 + a1*x + a2*x^2 + a3*x^3 + a4^4
  double a2;
  /// a0 + a1*x + a2*x^2 + a3*x^3 + a4^4
  double a3;
  /// a0 + a1*x + a2*x^2 + a3*x^3 + a4^4
  double a4;
}SPLINE;
//-----------------Spline-----------------
/**  \class Spline manages the coeficients of a N order spline */
class Spline{
 private:
  /// Vector of the coefficients
  double *Coe;
  /// Number of coefficients
  int NCoeff;
 public:
  /// Allocate the memory
  Spline(int N){
    NCoeff = N;
    Coe = NULL;
    Coe = (double*)calloc(N,sizeof(double));
    if(!Coe)printf("Spline: not allocated\n");
  };
  /// Free the memory
  ~Spline(){free(Coe);};
  /// Print the number of coefficients
  int GetN(){return NCoeff;};
  /// Print the n coefficient
  double GetCoe(int n){
    if(n<0||n>=NCoeff){printf("Spline: Wrong index\n");
      return 1.;} 
    return Coe[n];
  };
  /// Set the n coefficient
  void SetCoe(double Val,int n){
    if(n<0||n>=NCoeff){
      printf("Spline: Wrong index\n");
      return ;} 
    Coe[n] = Val;};
  /// Add to the n coefficient
  void AddCoe(double Val,int n){
    if(n<0||n>=NCoeff){
      printf("Spline: Wrong index\n");return ;} 
    Coe[n] += Val;
  };
};
/// Quaternion operations
class Quaternione{
  //Stolen from NeHe
private:
  double m_x;
  double m_y;
  double m_z;
  double m_w;
public:
  Quaternione();
  void AxisRotation(double x,double y,double z,double degrees);
  void CreateMatrix(double *pMatrix);
  Quaternione operator *(Quaternione q);

};



#endif //MATEMATICA_H
