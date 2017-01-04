#ifndef DRIMAGE_H
#define DRIMAGE_H

#include <Matematica.h>
#include <VarData.h>
#ifdef USE_PNG
#include <pngwriter.h>
///Image manipolator
class DrImage{
 private:
  /// Stored information about the piexel colour and position
  double **data;
  /// Stored information about the piexel colour and position
  double **data1;
  /// Width of the image
  int NWidth;
  /// Height of the image
  int NHeight;
  /// Number of files
  int NFile;
  /// File list
  char **FileList;
  /// Loaded initial image
  pngwriter *ImageIn;
 public:
  /// Constructor
  DrImage(int argc,char **argv);
  /// Destructor
  ~DrImage();
  /// Load an image
  void Load(char *FName);
  /// Overwrite the first image
  void ReLoad(char *FName);
  /// Load an image
  void Load2(char *FName);
  /// Write the image
  void Write(char *FName);
  /// Shift
  void BackFold(int Shift);
  /// Gravity
  void Gravity();
  /// Ising 
  void Ising();
  /// Shear 
  void Shear();
  /// Rotor 
  void Rotor();
  /// Gaussian filter
  void Gauss();
  /// Mirror
  void Mirror();
  /// Mirror
  void Difference();
  /// Cut a part of the image
  void Cut();
  /// Convolute the image with a matrix
  void ConvMatrix();
  /// Transpose
  void Transpose();
  /// Rotate the square of 90 degrees
  void EffectOnDataR(double *data2,int l,int w,int h,int ws,int hs,int NSquare);
  /// Transpose the square
  void EffectOnDataT(double *data2,int l,int w,int h,int ws,int hs,int NSquare);
  /// Mirror the square
  void EffectOnDataM(double *data2,int l,int w,int h,int ws,int hs,int NSquare);
  /// Lennard Jones simulation
  void LennardJones();
  /// Fourier
  void Fourier();
  /// Slitscan
  void SlitScan();
};
#endif //USE_PNG
#endif //DRIMAGE_H
