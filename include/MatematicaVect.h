#ifndef MATEMATICAVECT_H
#define MATEMATICAVECT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//-----------------Vettore-----------------
/** Geometrical operations on vectors  */
class Vettore{
 public:
  /// Allocates
  Vettore(int N);
  /// Allocates
  Vettore();
  /// Frees
  ~Vettore();
  /// Allocates and assigns
  Vettore(double *Pos,int N);
  /// Three dim vector
  Vettore(double x,double y,double z);
  /// Two dim vector
  Vettore(double x,double y);
  /// Dimension allocated
  int NDim;
  /// Where the data are stored
  double *x;
  //  Vettore *CrossProduct(Vettore *u,Vettore *v);
  /// Return the absolute value
  double Abs();
  /// Return the norm of a Vettore
  double Norm();
  /// Normlizes a Vettore
  double Normalize();
  //Vettore Normal(const Vettore u,const Vettore v);
  /// Computes the normal with respect to the Vettore u and v
  void Normal(const Vettore *u,const Vettore *v);
  /// Computes the normal to the plane described by u,v,w
  void NormalSurf(const Vettore *u,const Vettore *v,const Vettore *w);
  /// Project a point P on the point PS perpendicular to the surface described by the points S1, S2 and S3
  double ProjOnSurf(Vettore *S1,Vettore *S2,Vettore *S3,Vettore *P);
  /// multiply by a scalar
  void Mult(double Fact);
  /// substruct two Vettore
  void Subs(const Vettore *u,const Vettore *v);
  /// Multiplies the components of a Vettore for a scalar
  double ScalS(const Vettore *u,const Vettore *v);
  /// Computes the cosine with respect to \param u
  double CosAngle(Vettore *u);
  /// Computes the cosine between two Vettore
  double CosAngle(Vettore *u,Vettore *v);
  /// Computes the sine between two Vettore
  double SinAngle(Vettore *u,Vettore *v);
  /// Computes the sine with respect to \param u
  double SinAngle(Vettore *u);
  /// Computes the angle between two Vetttore
  double Angle(Vettore *u,Vettore *v);
  /// Computes the angle between two Vetttore
  double Angle(Vettore *u);
  /// Value of the N column
  double Col(int N);
  /// Value of the N column
  double Val(int N){return Col(N);};
  /// Set the N column
  void Set(double Val,int Col);
  /// Calculates the axis formed by two Vettore
  void Axis(Vettore *u,Vettore *v);
  /// Scalar product between two Vettore
  void ScalV(const Vettore *u,const Vettore *v);
  /// Vectorial product between two Vettore returns the area
  double VetV(const Vettore *u,const Vettore *v);
  /// Vectorial product between two Vettore in three dimension (faster) returns the area
  double VetV3(const Vettore *u,const Vettore *v);
  /// Vectorial product between two Vettore returns the area
  double VetV(const Vettore *u);
  /// Projects along the axis
  double ProjOnAxis(Vettore *a);
  /// The length of Pos on Axis
  double ProjOnAxis(Vettore *Pos,Vettore *Axis);
  /// Apply on a origin
  void ApplyOn(Vettore *o);
  /// Copy the vector 
  void Copy(Vettore *o);
  /// Export
  void Export(double *x);
  /// The vector perpendicolar 
  void PerpTo(Vettore *o);
  /// The vector perpendicolar 
  double PerpTo(Vettore *Pos,Vettore *Axis);
  /// The vector perpendicolar in three dimension (faster)
  double PerpTo3(Vettore *Pos,Vettore *Axis);
  /// Prints the components
  void Print();
  /// Rescale the total length of the vector
  void Rescale(double NewLength);
  /// Sum component by component
  Vettore operator+(const Vettore &Vet);
  /// Sum component by component
  Vettore operator+(const Vettore &Vet) const;
  /// Sum component by component
  Vettore& operator+=(const Vettore &Vet);
  /// Difference component by component
  Vettore operator-(const Vettore &Vet);
  /// Difference component by component
  Vettore operator-(const Vettore &Vet) const;
  /// Difference component by component
  Vettore& operator-=(const Vettore &Vet);
  /// Moltiplication component by component
  Vettore operator*(const Vettore &Vet); 
  /// Scalar product
  Vettore& operator*=( const Vettore& vec );
  /// Moltiplication component by a scalar
  Vettore operator*(const double Fact);
  /// Moltiplication component by a scalar
  Vettore& operator*=( const double f ); 
  /// Scalar product
  double operator%(const Vettore &Vet); 
  /// Scalar product
  double& operator%=(const Vettore &Vet); 
  /// Vectorial product
  Vettore operator^(const Vettore& vec);
  /// Vectorial product
  Vettore& operator^=(const Vettore& vec);
  /// Division component by component
  Vettore operator/(const Vettore& vec);
  /// Division component by component
  Vettore& operator/=(const Vettore& vec);
  /// Assigns the entries of the rhs Vettore to the lhs
  Vettore operator=(const Vettore &Vet);
  /// Returns a entry
  double * getPtr() {
    return (double*)&x;
  }
  /// Returns a entry
  const double * getPtr() const {
    return (const double *)&x;
  }
  /// Returns a entry
  double& operator[]( int col ){
    //return getPtr()[col];
    return x[col];
  }
  /// Returns a entry
  double operator[]( int col ) const {
    //return getPtr()[col];
    return x[col];
  }
};
  /// Moltiplication component by scalar
  Vettore operator*(double Fact, const Vettore& vec);
#endif //MATEMATICA_H
