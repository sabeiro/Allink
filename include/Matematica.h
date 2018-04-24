#ifndef MATEMATICA_H
#define MATEMATICA_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <fstream>
#include <iostream>
//#include <time.h>
//#include <complex>
//to look up!
//#include <boost/math/special_functions/fpclassify.hpp>
#include <assert.h>
#include "randomlib.h"
#include "MatematicaVect.h"
#include "MatematicaQuadri.h"
#include "MatematicaMatrix.h"
#include "MatematicaStruct.h"
#include "MatematicaPlane.h"
#include <mt19937ar.h>
#ifndef __GSL__
#else
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#endif
#ifndef FFTW3_H
#include  <fftw3.h>
#endif
#define STRSIZE 512
#define STRLEN  STRSIZE
#define QUAD(A) ((A)*(A))
#define SQUARE(A) QUAD(A)
#define SQR(A) QUAD(A)
#define CUB(A) (QUAD(A)*(A))
#define CUBE(A) (QUAD(A)*(A))
#define SIGN(A) (A>0?1.:-1.)
#define NTAB 32
#define DUE_PI (6.2831853071795864769252867665590057683943)
#define DEG_RAD (0.0174532925199433)
#define RAD_DEG (57.2957795130823)
#define INV_DUE_PI (0.159154943091895)
#define PI M_PI
#define SQRT_PI 1.77245385090552
#define MIN(a,b) ( ( (a)<=(b) ) ? (a) : (b) )
#define MIN3(a,b,c) ( MIN(MIN(a,b),c) )
#define MAS(a,b) ( ( (a)>=(b) ) ? (a):(b) )
#define MAX(a,b) MAS(a,b)
#define POS(a)   ( ( a>=0 ) ? (a):(-a) )
#define ASS(a)   POS(a)
#define ABS(a)   POS(a)
#define OP_IF(c,t) ( ((c)&(t))==(t) )
#define OP_INVERT     0x00001
#define OP_ROT_90     0x00002
#define OP_ROT_180    0x00004
#define OP_MIRROR     0x00008
#define OP_ROT_270    0x0000c
#define OP_TRANSPOSE  0x0000f

using namespace std;

//---------------------Miscellaneous------------------
int NumeroStringa(const char *NumeroC);
void abort_(const char * s, ...);
void SigErr(int IfWrong,const char * s, ...);
void Message(const char * s, ...);
double rnor();
//float nfix(void);
//float efix(void);
float zigset(unsigned long);
/* double Evalx(void (*pt2Function)(double x)){ */
/*   return pt2Function(.2); */
/* } */
//-------------------Matematica------------------------
/// Implementation of useful algorythms
class Matematica{
 private:
  /// Spline to store in the memory
  SPLINE SpMem;
  /// Minimum for the convergence
  double PrecMinimo;
  /// number of sub intervals
  int NPassi;
#ifdef __GSL__
  gsl_rng *rng;
#endif
 public:
  /// Constructor
  Matematica();
  /// 
  ~Matematica();
  //-----------------------Filter----------------------------
  /// Permutes a sequence without repeating
  int PermuteRandomAll(int *Sequence,int NMass);
  /// Permutes a sequence without repeating
  int PermuteRandomAll(PERMUTE *Sequence,int NMass);
  /// Applies the filter \param Mask on \param Point to \param Res
  int ApplyFilter(Matrice *Point,Matrice *Res,Matrice *Mask);
  /// Applies  the filter \param Mask on \param Res
  int ApplyFilter(Matrice *Res,Matrice *Mask);
  /// Boh
  int Transform(int *Out,int *In,int NEdge,int operation);
  /// Shift all the rows upwards/downwards
  void BackFold(Matrice *In,Matrice *Out,int NShift);
  /// Smooth the line with BSplines
  void Smooth(double *st,double *sw,int NIn,int NOut);
  //------------------------Func------------------------
  /// A tipical f(x) function
  typedef double(Matematica::*FUNC)(double x);
  /// Pointer to a function 
  FUNC Func;
  /// Pointer to a generic function
  double Evalx(double x){return (*this.*Func)(x);}
  /// Definition of the contact angle
  double ContactAngle(double x);
  /// Trial function
  double fProva(double x);
  /* /// Boh */
  /* double *fp; */
  /// Boh
  double F(double TD,double T);
  /// Boh
  double Df(double x,double Delta);
  //template <class Elabora> Elabora ElabFunc(double *st,double *sw,int NMass);
  //template <class Elabora> Elabora Elaboro(double *st,double *sw,int NMass){return (*this.*ElabFunc)(st,sw,NMass);};
  // NPunti = Elabora<int>(st,sw,NMass);
  typedef void(Matematica::*ELAB)(double *st,double *sw,int NMass);
  ELAB Elab;
  /// Pointer to a function which operates on \param sw with the data of \param st
  void ElabSt(double *st,double *sw,int NMass){return (*this.*Elab)(st,sw,NMass);};
  /// Derivate of \param st
  void Derivata(double *st,double *sw,int NMass);
  /// Derivate O(4) of \param st
  void DerO4(double *st,double *sw,int NMass);
  /// Integral of \param st
  double Integrazione(double *Punti,double *sw,int NMass);
  /// Perform a integration of a LJ6 Potential
  void IntegraA3();
  /// Square of the gradient
  void SquareGradient(double *st,double *sw,int NMass);
  /// Integrates the function within \param a and \param b
  double Integrazione(double a,double b);
  /// Itegrate a Gaussian
  double IntegrazioneGauss(double a,double b,double Scarto);
  /** Find the \param NRadici zeros of the pointed function between a\param a and \param b using  different algorithm */
  int Zeri(double a,double b,double *Radici,int NRadici);
  /// Use regula falsi algorithm to find the roots
  RADICE RegulaFalsi(double a,double b);
  /// Use Newton to find the roots
  RADICE Newton(double a);
  /// Other algorithm to find the roots
  double Estremo(double a,double b);
  /// Compute the factorial
  double Fattoriale(int n);
  /// Euler's gamma
  double Gamma(int n);
  /// Integer power
  double Elevato(double x,int Volte);
  /// Bessel function
  double Bessel(double Val,int Ord);
  /// Neumann function
  double Neumann(double Val,int Ord);
  /// Sign of -^n
  double Segno(int n);
  /// A faster Bessel
  double QuasiBessel(double Val,int Ord);
  /// A faster Neumann
  double QuasiNeumann(double Val,int Ord);
  /// Definition of a weighting function
  double WeightFunction(double x,double a);
  /// Definition of a weighting function
  double WeightFunction2(double x,double a);
  /// Integration of the LJ 6 term
  double LJHamaker(double r,double r_np,double theta);
  /// Integrate over r_np up to RadNp
  double LJHamakerCum(double Rad,double RadNpMin,double RadNpMax);
  /// Integrate over theta
  double LJHamaker(double Rad,double r_np);
  /// Integration of the LJ 6 term
  double LJ39(double r,double r_np);
  //-----------------------Sign---------------------
  /// Execute a command defined in string
  void ExecCommand(double *st,double *st1,int NMass,char *cmd);
  /// Execute a formula
  double ExecFormula(double x,double y,char *cmd);
  /// Execute a formula
  double ExecFormula(double **st,int n,char *cmd);
  /// Gaussian
  double Gauss(double Media,double Scarto,double x);
  /// Initialize the Gaussian number generator
  bool InizializzaGaussiano(double Scarto,int N);
  /// Gaussian random number
  double Gaussiano(double Media,double Scarto);
  /// Compute the spectrum
  void Spettro(double *st,double *sw,int NMass);
  /// Compute the 2d spectrum of \param st, return the 1d
  void Spettro2d(double *st,double *sw,int NMass);
  /// Compute the 2d spectrum of \param st, return the 2d 
  void Spettro2d(double *st,double **sw,int NMass);
  /// DFT implementation of the spectrum, slow
  void SpettroDFT(double *st,double *sw,int NMass);
  /// Compute the root of the signal
  void Radice(double *st,double *sw,int N);
  /// Compute the autocorrelation of a boolean signal
  void Autocor(bool *st,double *sAuto,int N);
  /// Compute the autocorrelation of the signal
  void Autocor(double *st,double *sAuto,int NMass);
  /// Norm of an array
  double Norm(double *st,int NMass);
  /// Normalize
  int NormalizeArea(double *st,int NMass);
  /// Normalize
  void NormalizeVect(double *st,int NMass);
  /// Normalize
  int Normalizza(double *st,int NMass);
  /// Normalize
  int  Normalizza(double *st,double *sw,int NMass);
  /// Compute the modulus
  void Modulo(double *st,double *sw,int NMass);
  /// Running average
  void MediaMobile(double *st,int NMass,double *sw,int Parti);
  /// Running average
  int MediaMobile(double *st,int NMass,double *sw,double *sErr,int Parti);
  /// Two points correlation
  int CorrelaDuePunti(double *st,int NMass,double *sw,int Punti);
  /// Self similarity
  void Autosimilarita(double *st,int NMass,double *sw,int Valori);
  /// Moments of a signal
  MOMENTI Distribuzione(const double *st,int NMass);
  /// Moments and histogram of a signal
  MOMENTI Distribuzione(const double *st,int NMass,double *Intervalli,int Valori,int IfNorm);
  /// Moments and histogram of a signal between two values
  MOMENTI Distribuzione(const double *st,int NMass,double *Intervalli,int Valori,double *Confine,int IfNorm);
  /// Moments and histogram of a signal between two values
  MOMENTI DistrErr(const double *st,int NMass,double *Intervalli,double *Err,int Valori,double *Confine,int IfNorm);
  /// Look for the Gaussian distribution
  MOMENTI DistribuzioneGauss(const double *st,int NMass,double *Intervalli,double *dInt,int Valori,int IfNorm);
  /// Look for the Maxwellian distribution
  MOMENTI DistribuzioneMaxwell(const double *st,int NMass,double *Intervalli,double *dInt,int Valori,int IfNorm);
  /// Compare the distribution of a sample of data
  void DistrSample(double *Px,double *Py,int NMax,double **Distr,int NBin,const int NSample,int IfNorm,double *xBound);
  /// Calculate the weighted average
  MOMENTI WeightAverage(const double *sx,const double *sy,int NMax);
  /// Weighted histogram analysis
  void WeightHisto(double **hist,double *Border,int NBin,int NHisto,double tolerance,double *OrPos,double *kSpring);
  /// Sort 
  void Sort(double *Sign,int NMass);
  /// Swap to indices
  void Swap(int i,int j,double *Sign);
  /// Swap to arrays
  void Swap(double *s,int si,double *t,int ti,const int NDim);
  /// Sort 
  void Sort(int *Sign,int NMass);
  /// Swap to indices
  void Swap(int i,int j,int *Sign);
  /// Create a file with a sign function
  void FileSin1d(char *FName);
  /// Create a file with a sign function in 2d
  void FileSin2d(char *FName);
  /// Convolute with a weight
  void ConvWeight(double *st,int NMax,double *sw,int *WIndex,int NWeight);
  /// Fill the weight array with a gaussian fuction
  void FillWeightGauss(double *st,int *WIndex,int NWeight,double CutOff,double Sigma);
  //--------------------Interp--------------------------------
  /// Linear interpolation between two points
  double LinInterp(double Px1,double Px2,double Py1,double Py2,double x);
  /// Linear interpolation
  RETTA InterRett(double *Px,double *Py,int NMass);
  /// Exponential interpolation
  RETTA InterExp(double *Px,double *Py,int NMass);
  /// Gaussian interpolation
  MOMENTI InterGauss(double *Px,double *Py,int NMass);
  /// Linear weighted interpolation
  RETTA InterRett(double *Px,double *Py,double *Peso,int NMass);
  /// Minimum of the Parabola between \param a and \param b
  PARABOLA MinimoParabola(double a, double b,double *Px,double *Py,int NMass);
  /// Global minimum interpolating via a Parabola
  PARABOLA MinimoParabola(double *Px,double *Py,int NMass);
  /// Three points parabolic interpolation
  SPLINE Parab(double *P1,double *P2,double *P3,int x,int y);
  /// Three points parabolic interpolation
  SPLINE Parab2(double *PA,double *PB,double *PC,int x,int y);
  /// Osculant circle
  CIRCLE Osculante(double *PA,double *PB,double *PC,int x,int y);
  /// Four point cubic interpolation
  SPLINE Cubica(double *PA,double *PB,double *PC,double *PD,int x,int y);
  /// Five points four order interpolation
  SPLINE Forth(double *PA,double *PB,double *PC,double *PD,double *PE,int x,int y);
  /// Three order spline
  SPLINE Spline3(double *P1,double *P2,double *P3,int x,int y);
  /// Three order spline first boundary
  SPLINE Spline3Beg(double *P1,double *P2,double *P3,int x,int y);  
  /// Three order spline last boundary
  SPLINE Spline3End(double *P1,double *P2,int x,int y);  
  /// Four order spline first boundary
  SPLINE Spline4Beg(double *P1,double *P2,double *P3,double *P4,int x,int y);
  /// Four order spline
  SPLINE Spline4(double *P1,double *P2,double *P3,double *P4,int x,int y);
  /// Four order spline
  SPLINE Spline4(double *P1,double *P2,double *P3,int x,int y);
  /// Four order spline just before the end
  SPLINE Spline4PreEnd(double *P1,double *P2,double *P3,int x,int y);
  /// Four order spline last boundary
  SPLINE Spline4End(double *P1,double *P2,int x,int y);
  /// Polinimial interpolation of \paran NMass order
  int Polinomio(double *P1,double *P2,int NMass,Spline *Sp);
  /// Boh
  int DerMatrix(double *Px,double *Py,int NMass,SPLINE Wg,Spline *Sp);
  /// Random uniform number
  double Casuale();
  /// Random number following a discrete probability
  double RandDiscrProb(double *Prob,int NBin);
  /// QBezier curve of three points
  double QBezier(double *P1,double *P2,double *P3,double x,int y);
  /// Computes \param times!
  int Factorial(int times);
  /// For the BSpline
  double Binomial(int times,int n);
  /// For the BSpline
  double Blend(const double *dPoint,double x,int nPoint,int nOrder);
  /// For the BSpline
  double Blend(double *dPoint,size_t Incr,double x,int nPoint,int nOrder);
  /// Computes the BSpline of a given \param MaIn
  int InterBSpline2D(Matrice *MaIn,Matrice *MaOut);
  /// Voronoi tassellation, in progress
  int Voronoi();
  //---------------------------Algebra------------------------------
  //  int Invert(Matrice A,Matrice B);
  /// External parameter to calculate the contact angle
  double Ypsilon;
  /// External parameter in the definition of the contact angle
  double PreFact;
};
//-----------------------Char---------------------
/// Trasform a string in the unicode number, aren't there any libaries?
uint sTable(char *String);
/// Trasform a character in the unicode number in the greek alphabet
uint sTable(char Char);
//--------------------Trials---------------------
class Funtore{
 public:
  virtual double operator()(double x)=0;
  virtual double Call(double x)=0;
};
template <class Classe> class SpecFuntore : public Funtore{
 private:
  double (Classe::*ftp)(const double x);
  Classe* pt2Object;
 public:
  SpecFuntore(Classe* _pt2Object,void (Classe::*_ftp)(double x))
    { pt2Object = _pt2Object; ftp = _ftp; };
  virtual double operator()(const double x)
  { return (*pt2Object.*ftp)(x);};
  virtual double Call(double x)
  { return (*pt2Object.*ftp)(x);};
};

#endif //MATEMATICA_H
