#ifndef MATEMATICAPLANE_H
#define MATEMATICAPLANE_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MatematicaVect.h"

//-----------------Piano-----------------
/** Define a plane  */
class Piano{
 public:
  /// Allocates
  Piano(Vettore *P1,Vettore *P2,Vettore *P3);
  /// Frees
  ~Piano();
  /// Distance
  double Distance(Vettore *P);
  /// Reflect velocity
  int Impact(Vettore *P,Vettore *V);
  /// Get vertex
  Vettore GetVertex(int i);
  /// Project on surface (point)
  Vettore ProjOnSurf(Vettore *Pos);
  /// Project on normal (vector)
  Vettore ProjOnNorm(Vettore *v);
  /// Is the orientation of the difference vectors on the same side?
  int SameSide(Vettore *P,Vettore *A,Vettore *B,Vettore *C);
  /// Calculate the inverse
  double Inv(double x);
  /// If the point is inside the triangle
  int IsOnSurf(Vettore *P);
  /// If the point is inside the triangle first method
  int IsOnSurf1(Vettore *P);
  /// If the point is inside the triangle second method
  int IsOnSurf2(Vettore *P);
  /// Reflect a vector by the normal
  Vettore Reflect(Vettore *V);
  /// Points defining the plane
  Vettore P1,P2,P3,P4;
  /// Direction vectors
  Vettore Dir21, Dir31, Dir23;
  /// Normal and inverse to the normal
  Vettore Norm,InvNorm;
  /// Boundaries
  double Bound[6];
  /// If the inverse to the normal is infinite
  int IsInf[3];
  /// d of ax+by+cz+d=0
  double dPar;
  /// Radius of the contact
  double Rad;
  /// Slope
  double mxy[3],mxz[3];
  /// Intercept
  double qxy[3],qxz[3];
};
#endif //MATEMATICA_H
