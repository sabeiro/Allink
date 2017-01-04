#ifndef MATEMATICAQUADRI_H
#define MATEMATICAQUADRI_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/// Quaternion class
class Quadri{
 private:
 public:
  /// An empty quaternion
  Quadri();
  /// A quaternion generated from an axis and an angle
  Quadri(double *Axis,double Angle);
  /// A quaternion generated form Euler's angle
  Quadri(double Ang1,double Ang2,double Ang3);
  /// A quaternion generated specifying the basis
  Quadri(double ww,double xx,double yy,double zz);
  /// Print the component of the axis
  double *Axis();
  /// First basis component
  double x;
  /// Second basis component
  double y;
  /// Third basis component
  double z;
  /// Forth basis component
  double w;
  /// Norm of the quaternion
  double Norm();
  /// Inverse norm of a 3d-vector
  double NormInv(double *Vett);
  /// Square of a quaternion
  double Sqr();
  /// Normalize a quaternion
  double Normalize();
  /// Normalize a vector
  double Normalize(double *Vett);
  /// Return the rotation angle 
  double Angle();
  /// Product between two quaternions
  void Prod(Quadri q);
  /// Product between two quaternions
  Quadri Prod(Quadri q,Quadri p);
  /// Give the conjugate
  Quadri GetConj();
  /// Inverse
  void Inv();
  /// Create a 4x4 rotation matrix
  void Matrix4x4(double *M);
  /// Create a 3x3 rotation matrix
  void Matrix3x3(double *M);
  /// Alternative creation of a rotation matrix
  void RotMatrix(double *data,int dim);
  /// Boh
  void Basis(double a,double b,double c,double d,double *Matr);
  /// Conjugate
  void Conj();
  /// Print a rotation matrix
  void PrintMatrix(double *M);
  /// Scalar product with a quaternion
  Quadri operator* (const Quadri &rq) const;
  /// Scalar product with a vector
  double *operator* (const double *Vet) const;
  /// Equal operator
  Quadri operator= (const Quadri &rq) const;
};
#endif //MATEMATICA_H
