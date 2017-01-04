#ifndef MATEMATICAMATRIX_H
#define MATEMATICAMATRIX_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MatematicaStruct.h"
//---------------------Matrice--------------------
/** Matrice computes the algebric operations on matrices
 */
class Matrice{
 private:
  /// dangerous
  void allocate(){
    delete[] data;
    data = new double [NSize*NSize];
  };
  /// Allocated size of the matrix
  int NSize;
  /// Number of columns
  int NCol;
  /// Number of rows
  int NRow;
  /// Number of third indexx
  int NZed;
  /// Number of dimensions
  int NDim;
 public:
  /// creates a square Matrice
  Matrice(int newNSize);// {Matrice(newNSize,newNSize);};
  //Matrice(int newNSize, int newactualsize);
  /// Creates a NRow x NCol Matrice
  Matrice(int NRow, int NCol);
  /// Creates a NRow x NCol x NZed Matrice
  Matrice(int NRow, int NCol,int NZed);
  /** Creates a NRow x NCol Matrice which points 
      to an already allocated double */
  Matrice(int ExtNRow, int ExtNCol, double *pointer);
  /** A 3x3 Matrice defined as discreate two dimensional
     derivation operator */
  Matrice(SPLINE Wg);
  /** A 5x5 Matrice defined as discrete two dimensional
     derivation operator, Dim defines different options */
  Matrice(SPLINE Wg,int Dim);
  /** a 4x4 or 3x3 Matrice defined by a quaternion \param Quadri q */
  Matrice(Quadri q,int dim);
  /** The data array point to an alreay allocated array */
  Matrice(double *M,int Nr,int Nc);
  /** Define a rotation matrix from Euler angles */
  Matrice(double Roll,double Pitch,double Yaw,int NDim);
  /** Define a rotation matrix from a axis and a rotation angle */
  Matrice(double *Axis,double Angle,int NDim);
  /// Freeing
  ~Matrice();
  /// Set a coefficient
  bool Set(int row, int column,double Val);
  /// Add the value of the coefficient to the previous one
  bool Add(int row,int col,double Val);
  /// Size of Matrice
  int Size(){return NSize;};
  /// Number of columns
  int pNCol(){return NCol;};
  /// Number of rows
  int pNRow(){return NRow;};
  /// Number of zed
  int pNZed(){return NZed;};
  /// Apply Ax = y
  void Apply(double *Known,double *UnKnown);
  /// Solve a system A|b = y
  int Solve(double *Known,double *UnKnown);
  /// Return size
  int getNRow();
  /// Name of the last function called. For debugging
  void Shout(const char *s, ... );
  /// Compare to identity
  void comparetoidentity();
  /// Set to product
  void settoproduct(Matrice& left, Matrice& right);
  /// Copy matrix
  void copymatrix(Matrice&  source);
  /// Set new size
  void setNRow(int newNRow);
  /// Return value
  void getvalue(int row, int column, double& returnvalue, bool& success);
  /// Invert 
  void Invert();
  /// Multiplication between two matrices
  void Mult(Matrice &A,Matrice &B);
  /// Multiplication by matrices
  void Mult(Matrice &A);
  /// Multiplication by a vector
  Vettore Mult(Matrice &A,Vettore &v);
  /// Multiplication by a vector
  Vettore Mult(Vettore &v);
  /// Multiplication by a vector on a vector u
  void Mult(Vettore &v,Vettore &u);
  /// Fill the entries randomly
  void RandomFill(double Max);
  /// Fill the entries of a differential operator
  void FillDiffOperator(SPLINE Wg,int NDim);
  /// Fill the entries for the Canny edge detector
  void FillCanny();
  /// Fill the entries for the Gauss blur
  void FillGaussian(double Sigma,double CutOff);
  /// Fill the entries for the 5x5 Gauss blur
  void FillGaussian5();
  /// Transpose the matrix
  void Transpose();
  /// Normalize the matrix
  void Normalize();
  /// Set all the entries to zero
  void Clear();
  /// Multiply by a scalar
  void Multiply(double Val);
  /// Copy on a matrix
  void CopyOn(Matrice *B);
  /// Print the entries
  void Print();
  /// Returns a value in 1d 
  double Val(int row);
  /// Returns a value in 2d
  double Val(int row,int col);
  /// Returns a value in 3d
  double Val(int row,int col,int zed);
  /// Computes the determinants
  double Det();
  /// Sums two matrices
  Matrice operator+(Matrice&);
  //Matrice operator*(Matrice&) ;
  /// Multiplies two matrices
  Matrice operator*(Matrice &A);
  //  Matrice &operator*(const Matrice&) const;
  /// Multiplies per a scalar
  Matrice operator*(const double&) const;
  /// Multiplies per a vector
  Vettore operator*(const Vettore&) const;
  /// Copy two matrices
  Matrice operator=(Matrice&);
  //  Matrice operator[][](int row,int col);
  /// Tensor product?
  Matrice operator^(Matrice&) const;
  /// Stored entries
  double* data;
  /// Convolute with a matrix
  void ConvoluteMatrix(double *Plot,int NGrid,int NDim,int IfMinImConv);
  /// Convolute with a matrix
  void ConvoluteMatrix1(double *Plot,int NGrid);
  /// Convolute with a matrix
  void ConvoluteMatrix1MinImConv(double *Plot,int NGrid);
  /// Convolute with a matrix
  void ConvoluteMatrix2(double *Plot,int NGrid);
  /// Convolute with a matrix
  void ConvoluteMatrix2MinImConv(double *Plot,int NGrid);
  /// Convolute with a matrix
  void ConvoluteMatrix3(double *Plot,int NGrid);
};
#endif //MATEMATICA_H
