/***********************************************************************
VarData: Header file for the VarData class.
 Copyright (C) 2008 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
***********************************************************************/
#ifndef VARDATA_H
#define VARDATA_H

#include "Matematica.h"
#include <vector> 
#include <list>

#ifdef USE_BOOST
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
using namespace  boost::iostreams;
#endif

;
#define CHAIN_IF_TYPE(c,t) ( ((c)&(t))==(t) )
#define VAR_IF_TYPE(c,t) ( ((c)&(t))==(t) )
#define VAR_ADD_TYPE(c,t) ( (c)|=(t) )
#define VAR_REM_TYPE(c,t) ( (c)&=(~t) )
/// Define the center of the frame
enum ShiftType{
  /// No shift
  SHIFT_NO      =0,
  /// Center of mass of the system
  SHIFT_CM      =1,
  /// Nanoparticle position
  SHIFT_NANO    =2,
  /// xy center of mass z nanoparticle
  SHIFT_CM_NANO =3,
};
/// Type of backfold
enum BackFoldType{
  /// No backfold
  BF_NO       =0,
  /// Back fold every particle
  BF_PART     =1,
  /// Back fold center of mass of the chain
  BF_CHAIN    =2,
  /// Back fold nanoparticle
  BF_NANO     =3,
  /// Do not even allocate
  BF_SKIP     =4,
  /// Plane containing nano axis and normal as reference
  BF_TILT     =5,
};
/// Type of bead
enum BeadType{
  /// Hydrophobic
  BEAD_PHOB   =0,
  /// Hydrophylic
  BEAD_PHIL   =1,
  /// Added chain
  BEAD_OIL    =2,
  /// Nanoparticle
  BEAD_NANO   ,
  /// Every type
  BEAD_EVERY  ,
};
/// Type of chain
enum ChainType{
  /// Every chain
  CHAIN_EVERY   =0x00000,
  /// Polyblock
  CHAIN_POLY    =0x00001,
  /// Homoblock
  CHAIN_HOMO    =0x00002,
  /// Added
  CHAIN_ADDED   =0x00006,
  /// Oriented upwards
  CHAIN_UP      =0x00005,
  /// Oriented downwards
  CHAIN_DOWN    =0x00009,
  /// Inner shell
  CHAIN_INNER   =0x00011,
  /// Outer shell
  CHAIN_OUTER   =0x00021,
  /// Stretched
  CHAIN_STRETCH =0x00200,
  /// Flabby
  CHAIN_FLABBY  =0x00400,
  /// Tilted
  CHAIN_TILTED  =0x00800,
};
/// The different file format used
enum VarFileFormat{
  /// txvl file format
  VAR_SYS_XVT        =  0x0001,
  /// xvt file format
  VAR_SYS_TXVL       =  0x0002,
  /// xyz file format
  VAR_SYS_XYZ        =  0x0004,
  /// xyz file format
  VAR_SYS_XYZT       =  0x0008,
};
/// Flags for the system configuration
enum VarSysPropr{
  /// Presence of a header and box size
  VAR_EDGE           =  0x0001,
  /// To allocate
  VAR_PART_ALLOCATED =  0x0002,
  /// To allocate
  VAR_CH_ALLOCATED   =  0x0004,
  /// Bilayer
  VAR_MEMBRANE       =  0x0008,
  /// Skip particle count/number of particle in the header
  VAR_OPEN_TRUST     =  0x0010,
  /// The structur LINK is allocated
  VAR_LINK_ALLOCATED =  0x0020,
  /// The structure CHAIN is allocated and the chains are defined
  VAR_CHAIN_DEF      =  0x0040,
  /// The structure BLOCK is allocated
  VAR_NBLOCK_ALL     =  0x0080,
};
/// Flags for the creation of an initial system
enum VarCreate{
  /// Absorbing boundary conditions
  VAR_ABSORBING  =  0x0001,
  /// Reflecting boundary conditions
  VAR_REFLECTING =  0x0002,
  /// Added chains
  VAR_ADDED      =  0x0004,
  /// Chains distributed all around
  VAR_DISTRIBUTED=  0x0008,
  /// Two opposed bilayer
  VAR_SOLVENT    =  0x0010,
  /// Planar geometry
  VAR_PLANAR     =  0x0020,
  /// Planar geometry
  VAR_PLANAR_PE  =  0x0040,
  /// Tube geometry
  VAR_TUBE       =  0x0080,
  /// Tube geometry
  VAR_VESICLE    =  0x0100,
  /// Tube geometry
  VAR_COATING    =  0x0200,
  /// Two tails
  VAR_TWOTAILS   =  0x0400,
  /// Quenched obstacle
  VAR_OBSTACLE   =  0x0800,
  /// Quenched obstacle
  VAR_OPPOSED    =  0x0f00,
};
/// Architecture of the block
enum BlockArch{
  /// Points
  ARCH_POINTS   = 0x0001,
  /// Lines
  ARCH_LINES    = 0x0002,
  /// Two tails
  ARCH_TWOTAILS = 0x0004,
  /// Two tails
  ARCH_CLUSTER  = 0x0008,
};
/// Shape of the inclusion
enum NanoShape{
  /// no nanoparticle
  SHAPE_NONE    = 0x00001,
  /// spherical shape
  SHAPE_SPH     = 0x00002,
  /// radius and height
  SHAPE_HEI     = 0x00004,
  /// cylindrical shape
  SHAPE_CYL     = 0x00008,
  /// tilting pill
  SHAPE_TILT    = 0x00010,
  /// cylinder with spherical ends
  SHAPE_PILL    = 0x00020,
  /// planar wall 
  SHAPE_WALL    = 0x00040,
  /// external cylindrical field
  SHAPE_PORE    = 0x00080,
  /// smoothed ends
  SHAPE_SMOOTH  = 0x00100,
  /// smoothed ends
  SHAPE_CLUSTER = 0x00200,
  /// external field
  SHAPE_EXT     = 0x00400,
  /// janus peptide
  SHAPE_JANUS   = 0x00800,
  /// cross links
  SHAPE_CLINKS  = 0x01000,
  /// stalk
  SHAPE_STALK   = 0x02000,
  /// AFM tip
  SHAPE_TIP     = 0x04000,
  /// toroidal shape
  SHAPE_TORUS   = 0x08000,
  /// harmonic external pot
  SHAPE_HARM    = 0x10000,
  /// external potential
  SHAPE_UMBR    = 0x20000,
  /// rigid wall at the boundaries
  SHAPE_BOUND   = 0x40000,
};
/// Information of every particle
typedef struct{
  /// xyz Position of the particle
  double Pos[3];
  /// xyz Backfold distance
  double Bkf[3];
  /// xyzr Velocity of the particle
  double Vel[4];
  /// Particle identifier
  int Idx;
  /// Chain Identifier
  int CId;
  /// Type 
  int Typ;
}PART;
/// Structure with the links of the particles
typedef struct{
  /// with whom is bonded
  int *Link;
  /// How many links per particle
  int NLink;
}LINKS;
/// Information of every chain
typedef struct{
  /// xyzr Postion of the chain
  double Pos[4];
  /// Directive of the chain
  double Dir[3];
  /// Chain velocity
  double Vel[3];
  //  Angle with the normal;
  double Angle;
  /// Type of the chain (see list CHAIN_)
  int Type;
  /// Number of particles per chain
  int NPCh;
  /// Initial bead
  int InitBead;
  /// End bead
  int EndBead;
}CHAIN;
/// Information for every block
typedef struct{
   /// BLock name
  char Name[60];
  /// # particles per chain
  int NPCh;
  /// # particles
  int NPart;
  /// # chains
  int NChain;
  /// Initial particle position
  int InitIdx;
  /// End particle position
  int EndIdx;
  /// Architecture type
  int Arch;
  /// Nano index
  int nNano;
  /// Diblock limit of the chain
  int Asym;
} BLOCK;
/// Information for every soft object
typedef struct{
  /// soft name
  char Name[60];
  /// initial position
  double Pos[3];
  /// bias velocity
  double Vel[3];
  /// dimension xyz/rad height
  double Size[3];
  /// topology
  int Topology;
  /// # particles
  int NPart;
  /// # chains
  int NChain;
  /// # part per chain
  int NPCh;
  /// Initial Position
  int InitIdx;
  /// End position
  int EndIdx;
} SOFT;
/// General information of the system
typedef struct{
  /// Total time
  double Time;
  /// Temperature
  double Temp;
  /// Temperature
  double Beta;
  /// Pot, kinetik, free
  double Energy[3];
  /// xyzr edges of the simulation box
  double Edge[4];
  /// Inverted xyzr edges of the simulation box
  double InvEdge[4];
  /// Center of mass of the system
  double Cm[3];
  /// Velocity of the system
  double Vel[4];
  /// Three components of the pressure
  double Pre[3];
  /// Chemical potential of the water
  double vBB;
  /// Surface tension
  double SurfTens;
  /// Incompatibility
  double chiN;
  /// Incompressibility
  double kappaN;
  /// Density coexistence
  double rho;
  /// Prefactor of the bending potential
  double kappaBend;
  /// Prefactor of the bond spring
  double kappaSpring;
  /// Rest length of the harmonic potential
  double SpringRest;
  /// Convertion unit R_e over CutOff
  double ReOverCutOff;
  /// Timestep
  double Deltat;
  /// Weighting function straight length
  double WFuncStraight2;
  /// Weighting function straight length
  double WFuncStraight3;
  /// Courrent step
  unsigned long Step;
  /// Number of particle
  int NPart;
  /// Number of chain
  int NChain;
  /// Number of particle per chain
  int NPCh;
  /// Number of allocated particles
  int NAllocP;
  /// Number of allocated chains
  int NAllocC;
  /// # of types of the particle
  int NType;
  /// Maximum number of bonds
  int NLink;
  /// Number of blocks
  int NBlock;
  //  Number of nanoparticles
  int NNano;
}GENERAL;
/// Some calculated properties of the system
class Properties{
 public:
  /// Set the values to zero
  Properties(){RePhob=0.;RePhil=0.;VolPhob=0.;VolPhil=0.;FactPhob=0.;FactPhil=0.;GyrPhob=0.;GyrPhil=0.;ChDiff=0.;};
  /// Print the current values
  void Print(){printf("RePhob %lf RePhil %lf\n",RePhob,RePhil);
  printf("RadGyrPhob %lf RadGyrPhil %lf\n",GyrPhob,GyrPhil);
  printf("StructFactPhob %lf StructFactPhil %lf\n",FactPhob,FactPhil);
  printf("VolPhob %lf VolPhil %lf\n",VolPhob,VolPhil);
  printf("ChainDiff %lf\n",ChDiff);};
  /// End2End distance of the 
  double RePhob;//End to end distance
  double RePhil;
  double VolPhob;//volume fraction
  double VolPhil;
  double FactPhob;//structure factor
  double FactPhil;
  double GyrPhob;//radius of gyration
  double GyrPhil;
  double ChDiff;//lateral diffusion of the chains
  Properties operator+(const Properties&) const;
  Properties operator*(const double&) const;
};
/// Energy contribution from the density functional hamiltonian
class MatInt{
 private:
  /// Interaction matrix
  double *IntMatr2;
  /// Interaction matrix
  double *IntMatr3;
  /// Number on spicies
  int NType;
  /// Maximum order
  int NOrd;
 public:
  /// Creates and fill the matrix
  MatInt(int NType,int NOrd);
  /// Fill all the entries
  void FillEntries(double *Matr,int Ord);
  /// Matrix element of the two type interaction
  int IntType(int t1,int t2);
  /// Matrix element of the two type interaction
  int IntType(int t1,int t2,int t3);
  /// Prefactor of the force
  double Coeff(int t1,int t2);
  /// Prefactor of the force
  double Coeff(int t1,int t2,int t3);
  /// Set the prefactor
  void SetCoeff(double Co,int t1,int t2);
  /// Set the prefactor
  void SetCoeff(double Co,int t1,int t2,int t3);
  /// Print the entries
  void Print();
  /// Print n type
  int pNType(){return NType;};
  /// Rescale entries
  void Rescale(double SFactor,int Order);
};
/// Information about the nanoparticle
typedef struct{
  /// Architecture file
  char ArchFile[60];
  /// Position
  double Pos[3];
  /// Backfolded position
  double Bkf[3];
  /// Velocity
  double Vel[3];
  /// Forces
  double Force[3];
  /// Rotation axis
  double Axis[3];
  /// Angular momentum
  double AMom[3];
  /// Temporal angular momentum
  double AMomTemp[3];
  /// Angular velocity
  double AVel[3];
  /// Mass
  double Mass;
  /// Size
  double Rad;
  /// Strength of the interaction
  double Hamaker;
  /// Height of the cylinder
  double Height;
  /// Friction term
  double Gamma;
  /// Stochastic term
  double Zeta;
  /// Viscosity
  double Viscosity;
  /// Reference potential
  double OffSet;
  /// Cut off of the potential
  double CutOff;
  /// Thickness of the LJ well
  double Coating;
  /// Baseline of the potential
  double BaseLine;
  /// Minimum distance threshold
  double DistThr;
  /// Maximum value of the force
  double ForThr;
  /// Maximum value of the potential 
  double PotThr;
  /// Area of a pore or a stalk
  double Area;
  /// 0 none, 1 spherical, 2 cylindrical 3 wall
  int Shape;
  /// Number of links connecting the constituent monomers
  int NLink;
  /// Number of monomers per side
  int NHeight;
  /// Number of monomers per circle
  int NCircle;
  /// In which block is the peptide written
  int nBlock;
}NANO;
/// A Cartesian tern
typedef struct {
  /// Cartesian coordinates
  double x[3];
} XYZ;
/// All the vertex in a cell
typedef struct {
  /// Vertices
  XYZ p[8];
  /// Normals
  XYZ n[8];
  /// Density at the vertices
  double val[8];
} GRIDCELL;
/// Define a triangle
typedef struct {
  /// The three vertices
  XYZ p[3];
  /// Centroid
  XYZ c;
  /// Reference for the vectors
  int v[3];
  /// Normal to the vertices
  XYZ n[3];
} VAR_TRIANGLE;
/// Define a triangle
typedef struct {
  /// The two vertices
  XYZ p[2];
  /// Centroid
  XYZ c;
  /// Reference for the vectors
  int v[2];
  /// Normal to the vertices
  XYZ n[2];
} VAR_LINE;
//-----------------Class-VarData----------------------------------
/// Reads and elaborates a system of chains
class VarData{
 private:
  /// General information on the system
  GENERAL *Gen;
 public:
  /// Implementation of all usefull algorythms
  Matematica *Mat;
  /// Matrix of the prefactor of the interactions
  MatInt *MInt;
  /// Set the constants
  VarData();
  /// Destructor
  ~VarData();
  //-------------------------VarData.cpp-------------------
  /// If enabled call the function position
  void VarMessage(const char * s, ...);
  /// Open the \param InFile and back fold
  bool Open(char *InFile,int BF);
  /// Opens a file without reallocationg
  bool OpenRisk(char *InFile,int BF);
  /// Opens a file checking if the information are correct
  bool OpenTrust(char *InFile,int BF);
  /// Alloc the structures
  void AllocPart();
  /// Alloc the structures
  void AllocChain();
  /// Calculate some basis properties
  Properties SysProperties();
  /// Print a string with the system information
  void SysInfo(char *cSystem);
  /// Print a string with the system definitions
  void SysDef(char *cSystem);
  /// Calculates fundamental quantities
  char *SysState();
  /// Set the virial coefficients from the known values of density coex...
  void SetCoeff();
  /// Set and normalize the virial coefficients from the arrays v2 and v3
  void SetCoeff(double *v2,double *v3);
  /// Return the relative distance between two particles (wrapped)
  double TwoPartDist(int p1,int p2,double *RelDist);
  /// Return the relative distance between two particles (wrapped)
  double TwoPartDist(double *Pos,int p2,double *RelDist);
  /// Return the relative distance between two particles (wrapped)
  double TwoPartDist2(int p1,int p2,double *RelDist);
  /// Return the relative distance between two particles (wrapped)
  double TwoPartDist2(double *Pos,int p2,double *RelDist);
  /// Return the relative distance between two particles (wrapped) if the particles are within the cut off
  int TwoPartDist(int p1,int p2,double *RelDist,double CutOff);
  /// Return the relative distance between two particles (wrapped) if the particles are within the cut off
  int TwoPartDist(double *Pos,int p2,double *RelDist,double CutOff);
  //----------------------VarDataWrite.cpp-----------------
  /// Writes a "system-file" or a "x y z" file" 
  bool Write(char *OutFile);
  /// Writes a "system-file" or a "x y z" file" 
  bool WriteTxvl(char *OutFile);
  /// Writes a "system-file" or a "x y z" file" 
  bool WriteXvt(char *OutFile);
  /// Writes a "system-file" or a "x y z" file" 
  bool WriteXyz(char *OutFile);				       
  /// Header interactions
  void HeaderInteraction(FILE *FileToWrite);
  /// String for the rigid inclusion in the header file
  void StringNano(char *NString,int n);
  /// Write the nano section of the header to the file
  int HeaderNano(FILE *FileToWrite);
  /// Header soft
  int HeaderSoft(char *Line);
  /// Write the positions of the egdes of the rectangles
  void WriteLinkedSurf(FILE *FWrite,double *Plot,int NSample,int NType,double *Bound,int *PId);
  /// Write the particle position as linked edges of squares
  void WriteSurf(FILE *F2Write,double **Plot,int NSample,int OffSet);
  /// Identifier of the shape
  void ShapeId(int iShape,char *Shape);
#ifdef USE_BOOST
  /// Open a file with compression
  void OpenComprWrite(char *FName,filtering_ostream ZipFile);
  /// Open a file with compression
  void OpenComprRead(char *FName,filtering_istream ZipFile);
  /// Close a file with compression
  void CloseCompr(filtering_istream ZipFile);
  /// Close a file with compression
  void CloseCompr(filtering_ostream ZipFile);
#endif
  //------------------------VarDataBackFold.cpp------------
  /// Backfold the particle position
  bool BackFold(int How);
  /// Definition of the chain
  int BfDefChain();
  /// Find the box size if missing
  int BfEdge();
  /// Define the distance form the nanoparticle
  void DistFromNp();
  /// Backfold the system wrt the reference position
  void ShiftRef(int BackFold);
  /// Find the position of the stalk
  int StalkPos(double *OldPos);
  /// Backfold the nano described as a cluster of monomers
  void BfPep();
  /// Describe the backbone of a filament
  void BackBone(double *Line,int NBin);
  /// Describe the line for a linear stalk
  void StalkLineProf(double *Line,int NBin);
  /// Find the position of the stalk second method
  void StalkPos2(double *OldPos,double *CmStalk);
  /// Find the position of the stalk third method
  void StalkPos3(double *OldPos,double *CmStalk);
  /// Find the position of the stalk forth method
  int StalkPos4(double *OldPos,double *CmStalk);
  /// Weight of the neighblorung normal on a vertex
  double NormalWeight(VAR_TRIANGLE *Triang,double *Weight,int NGrid,int NTri);
  /// Connect the lines in a chain
  void ConnectLineChain(VAR_LINE *Triang,int NGrid,int NTri);
  /// Connect the lines in a chain
  void ConnectLineChain2(VAR_LINE *Triang,int NGrid,int NTri);
  /// Connect the lines in a chain
  void ConnectLineChain3(VAR_LINE *Triang,int NGrid,int NTri);
  /// Find the position of the pore
  double PorePos();
  //-----------------VarDataString.cpp------------------
  /// Retrive from a string the information concerning the mask
  int Fetch(char *str,char *mask,char *fmt, ... );
  /// Retrive from a string the position of the brakets
  int BraketPos(char *str,char *mask,int *sPos,int *sLen);
  /// Retrive from a string the information concerning the mask
  int Fetch(char *str,char *mask,int NArg,double *Val);
  /// Copy the value in the @param String to the @param Value 
  bool ReadString(const char *String,char *cLine,double *Value);
  /// Copy the value in the @param String to the @param Value 
  bool ReadString(const char *String,double *Value,char *line);
  /// Copy the value in the @param String to the @param Value 
  bool ReadString(const char *String,char *cLine,int *Value);
  /// Copy the value in the String to the Value referring to the position of pLine
  int  ReadVal(char *pLine,double *Value);
  //-----------------VarDataRead.cpp-------------------
  /// Read a single line in format Xvt
  int ReadLineXvt(char *cLine,double *Pos,int *Type);
  /// Reads a "configuration file"
  bool ReadConf(char *InFile);
  /// Reads a header
  void ReadHeader(FILE *FileToRead);
  /// Reads a header for a txvl file format
  void ReadHeaderTxvl(FILE *FileToRead);
  /// Reads a header of xvl file format
  void ReadHeaderXvt(FILE *FileToRead);
  /// Reads particle type and position
  int ReadPart(FILE *FileToRead);
  /// Reads a type-position-velocity-link file
  int ReadPartTxvl(FILE *FileToRead);
  /// Reads a position-velocity-type file
  int ReadPartXvt(FILE *FileToRead);
  /// Reads a x y z file
  int ReadPartXyz(FILE *FileToRead);
  /// Reads a x y z t file
  int ReadPartXyzt(FILE *FileToRead);
  /// Reads the information to alloc the structure
  int ReadPassThru(FILE *FileToRead);
  /// Reads the specifications about the nano
  int ReadSoft(FILE *ConfFile);
  /// Reads the specifications about the hard object
  void ReadNano(FILE *ConfFile,int NCircle,int NHeight);
  /// Reads and set the specifics of the nano
  int NanoString(char *cLine,int n);
  /// Substitue the nano header
  void SubNanoHeader(char *cFile);
  /// Identifier of the shape
  int ShapeId(char *Shape);
  //--------------------VarDataCreate.cpp------------------------
  /// Define and write the system as described in the conf file
  int DefSoft(char *nome2,char *ConfF);
  /// Creates a trial system
  int TrialSys();
  /// Creates an initial system
  bool CreateSoft(int *arch,double Thickness,int s);
  /// Soft in a tube shape
  void CreateTube(int *arch,double Thickness,int s);
  /// planar membrane
  void CreatePlanar(int *arch,double Thickness,int s);
  /// vesicle
  void CreateVesicle(int *arch,double Thickness,int s);
  /// coating around a cylindrical nanoparticle
  void CreateCoating(int *arch,double Thickness,int s);
  /// Creates obstacles 
  void CreateObstacle(int *arch,double Thickness,int s);
  /// No particle inside the nano
  int CheckNano(double *Pos,int s);
  /// Defines the nanoparticle as a net of monomers
  void AddProtein(int NCircle,int NHeight,int nNano,char *filename);
  /// Defines the nanoparticle as a net of monomers
  void CreateProtein(int nNano,int nStart);
  /// Fill the protein with water
  void AddStuffing(char *filename,int nStuffing,int nNano);
  /// Add phantom solvent at the bottom
  void AddSolvent(char *filename,int nWater);
  /// Add homopolymer chains in the bilayer
  void AddChains(char *filename,double Thickness);
  /// Add cholesterol chains in the bilayer
  void AddCholesterol(char *filename,double Thickness,int s);
  /// Define four different blocks
  void DefBlock(int *NChStep,int How);
  /// set the remaining information
  void DefRest(int *arch,int s);
  /// return the number in the chain of the next particle put 
  int PutPart(int j,int p,int HalfLim,double sigma);
  /** Find the couples of most neighbouring chains */
  void FindNeighbours(char *FileName);
  //---------------------VarDataEl.cpp----------------------------
  /// Swap two chains
  void SwapChain(int c1,int c2,int b);
  /// Swap two cahins
  void SwapChain(int c1,int c2);
  /// Swap two particle
  void SwapPart(int p1,int p2);
  /// Update the new number of chains
  void ChangeNChain(int NChain,int b);
  /// Shift the system accordin to the SHIFT_ definitions
  bool ShiftSys(int How);
  /// Define a normal coordinate for every patch
  void SampleSurface(double *Plot,int NSample,int Type);
  /// Define a normal coordinate for every patch
  MOMENTI SampleSurfacePart(double *Plot,int NSample,int Type);
  /// Define a normal coordinate for every patch
  MOMENTI SampleSurface(Matrice *Plot,int NSample,int Type);
  /// Allocate and fill PlotMem with the particle average position
  MOMENTI SampleSurfaceMem(int NSample);
  /// Load in the array Plot the density of the system 
  void LoadDensFile(double **Plot,int NBin);
  /// Perform a spatial derivative on a surface
  int SpatialDerivative(Matrice *Surface,Matrice *Resp,SPLINE Weight,int NSample);
  /// Shift a block wrt to Shift 
  void ShiftBlock(Vettore *Shift,int b);
  /// Rotate a block wrt to the Axis from the Origin
  void RotateBlock(Vettore *Axis,Vettore *Origin,int b);
  /// Mirror the position wrt to a plane
  void MirrorBlock(Vettore *Px1,Vettore *Px2,Vettore *Px3,int b);
  /// Transform a block
  void Transform(int block);
  /// Point to the shape function
  void Point2Shape(int iShape);
  /// Data type for distance/field functions
  typedef double(VarData::*NANO_DIST)(double *Pos,int n);
  /// Pointer to a distance/field function
  NANO_DIST Nano_Dist;
  /// Pointer to a generic function
  double NanoDist2(double *Pos,int n){return (*this.*Nano_Dist)(Pos,n);}
  /// Distance from the nanoparticle
  double NanoDist2(double x,double y,double z,int n);
  /// No field
  double FieldNo(double *Pos,int n);
  /// Scalar field of a sphere
  double FieldSphere(double *Pos,int n);
  /// Scalar field of a elipsoid
  double FieldElips(double *Pos,int n);
  /// Scalar field of a elipsoid
  double FieldParab(double *Pos,int n);
  /// Scalar field of a cylinder
  double FieldCyl(double *Pos,int n);
  /// Scalar field of a transmembrane protein
  double FieldTransMem(double *Pos,int n);
  /// Scalar field of a janus peptide
  double FieldJanus(double *Pos,int n);
  /// Scalar field of a janus peptide
  double FieldTorus(double *Pos,int n);
  /// Scalar field of a tilted cylinder
  double FieldTilt(double *Pos,int n);
  /// Scalar field of a hard wall at the box edges
  double FieldBound(double *Pos,int n);
  /// Scalar field of a tilted cylinder
  double FieldTiltWall(double *Pos,int n);
  //---------------------VarDataExp.cpp----------------------------
  /// 1-d pair correlation
  int PairCorrelation(double *Point,int NSample,int How,int Type);
  /// Circular 2-d pair correlation
  int PairCorrelationRound(double **Point,int NSample,int Type);
  /// 2-d pair correlation on a square
  int PairCorrelationSquare(double **Point,int NSample,int Type);
  /// 2-d pair correlation on a square fererring to the pep position
  int PairCorrelationPep(double **Point,int NSample,int Type);
  /// 2-d Scattering
  int Scattering2d(double **Point,int NSample,int Type);
  /// 2-d scattering
  int Scattering2D(double **Point,int NSample,int Type);
  /// 1-d spectrum of a surface
  void Spettro2d(double *Points,int NSample,int Type);
  /// 2-d spectrum of a sirface
  void Spettro2d(double *Plot,int NSample);
  //---------------------VarDataContour.cpp------------------------
  /// Calculate the density profile for the x, y, z, r  coordinate 
  int DensityProfile(int coord,int NSample,int NType,double *dDensity);
  /// Sampled three dimentional weighted shape of the system
  int Core(double ***Plot,int NSample,double Border[3][2]);
  /// rzd representation of the system referring to @param How
  int RadDistr(int NSample,double *Plot,double Border[2],int How);
  /// Density profile along a worm like micelle
  int Worm(int Partition,int NSample,double *Border,double *dPoint);
  /** Fill an array of \param NSample values with the volume
   contribution in a rectangular box*/
  void VolumeCircSlab(double *VolContr,int NSample);
  /** Following the contour of a stalk */
  void Stalk(int NSample,int NLevel,double **Plot,double Threshold);
  //---------------------VarDataPos.cpp----------------------------
  /// The naerest @param Vertex -particle close to every chain
  int Arrange(int **Triangle,int Vertex);
  /// Boh...
  int Folding();
  /// A cell list to be fixed
  int OrderPos();
  /// return a univocal index of the chain position
  int CalcnPos(double *Pos);
  /// Boh
  int Neighbour(double *Pos);
  /// Distribution of number of chain per patch
  int NChainPSquare(double *Plot);
  /// Boh
  int LateralFluctuation(double *Plot,int LatValue);
  /// Voronoi tassellation
  int Voronoi();
  /// Return the integer index with respect to the partition NSquare
  int PosVectInt(double *Pos);
  //------------------VarInterp.cpp----------------------------
  /// Discontinous parabolas
  int InterParab(PART *PmIn,PART *PmOut,int NIn,int nOut);
  /// Discontinous parabolas
  int InterParab2(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// Discontinous cubic
  int InterCubica(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// Discontinous forth degree
  int InterForth(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// third order spline
  int InterSpline3(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// forth order spline
  int InterSpline4(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// BSpline
  int InterBSpline(PART *PmIn,PART *PmOut,int NIn,int NOut);
  /// 2-d BSpline
  int InterBSpline2D(double **PlIn,double **PmOut,int NIn,int NOut);
  /// 2-d BSpline
  int InterBSpline2D(double *PlIn,double *PmOut,int NIn,int NOut);
  /// 1-d BSpline
  int InterBSpline1D(double *PlIn,double *PmOut,int NIn,int NOut);
  /// NIn-polynomian
  int InterPoly(PART *PmIn,PART *PmOut,int NIn,int nOut);
  /// Boh
  int InterDerMatrix(PART *Pm, int NMass,SPLINE Weight,double Offset);
  /// Smooth a grid with BSplines
  void SmoothGrid(int NSample,char *FWrite);
  /// Smooth a grid with BSplines and update the particle positions
  void SmoothGrid(int NSample);
  /// Convolute a matrix
  void ConvoluteMatrix(double *Plot,int NGrid,Matrice *Mask,int NDim);
  /// Convolute a matrix 1d 
  void ConvoluteMatrix1(double *Plot,int NGrid,Matrice *Mask);
  /// Convolute a matrix 2d 
  void ConvoluteMatrix2(double *Plot,int NGrid,Matrice *Mask);
  /// Convolute a matrix 3d
  void ConvoluteMatrix3(double *Plot,int NGrid,Matrice *Mask);
  //-------------------VarDataComm.cpp----------------------
  /// Set and reallocate the number of particles
  int SetNPart(int NewNPart);
  /// Set and reallocate the number of chains
  int SetNChain(int NewNCh);
  /// Set and reallocate the number of links
  int SetNLink(int NewNCh);
  /// Set and reallocate the number of particles per chains
  void SetNPCh(int NewNCh);
  /// Set the number of species
  void SetNType(int NewNType);
  /// (re)allocate the links
  int AllocLinks(int NewNCh);
  /// Set NBlock
  int SetNBlock(int Val);
  /// Set NNano
  int SetNNano(int Val);
  /// Copy the part P2 on part P1
  void Copy(PART *P1,PART *P2,int NPartOld);
  /// Copy the chain C2 on chain C1
  void Copy(CHAIN *C1,CHAIN *C2,int NChainOld);
  //-------------------VarDataMarchCubes.cpp----------------
  /// Defines the triangles close to the IsoLevel of the 3d density Plot
  VAR_TRIANGLE *MarchingCubes(double *Plot,int NSample,double IsoLevel,int *NTri);
  /// Defines the triangles close to the IsoLevel of the 3d density Plot
  VAR_LINE *MarchingSquares(double *Plot,int NSample,double IsoLevel,int *NTri);
  //-------------------VarDataCGAL.cpp----------------------
  //#ifdef __CGAL_h__
  /// Calculate the (temporal/radial) area distribution
  void AreaDistr(double *Distr,double *RadDistr,int NSample);
  //#endif //CGAL
  /// Total time
  double pTime(){return Gen->Time;};
  /// Delta t
  double pDeltat(){return Gen->Deltat;};
  /// Temperature
  double pTemp(){return Gen->Temp;};
  /// Beta factor 1/kTB
  double pBeta(){return Gen->Beta;};
  /// Pot, kinetik, free
  double pEnergy(int d){return Gen->Energy[d];};
  /// xyzr edges of the simulation box
  double pEdge(int d){return Gen->Edge[d];};
  /// Inverted xyzr edges of the simulation box
  double pInvEdge(int d){return Gen->InvEdge[d];};
  /// xyzr edges of the simulation box
  double pVol(){return Gen->Edge[0]*Gen->Edge[1]*Gen->Edge[2];};
  /// Center of mass of the system
  double pCm(int d){return Gen->Cm[d];};
  /// Maximum velocity
  double pVelMax(int d){return Gen->Vel[d];};
  /// Incompatibility
  double pchiN(){return Gen->chiN;};
  /// Incompressibility
  double pkappaN(){return Gen->kappaN;};
  /// Bending coupling
  double pkBen(){return Gen->kappaBend;};
  /// Spring coupling
  double pkSpr(){return Gen->kappaSpring;};
  /// Bending coupling
  void SetkBen(double Val){Gen->kappaBend = Val;};
  /// Spring coupling
  void SetkSpr(double Val){Gen->kappaSpring = Val;};
  /// Rest distance of the harmonic potential
  void SetSprRest(double Val){Gen->SpringRest = Val;};
  /// Rest distance of the harmonic potential
  double pSprRest(){return Gen->SpringRest;};
  /// Density coexistence
  double prho(){return Gen->rho;};
  /// Re/CutOff
  double pReOverCutOff(){return Gen->ReOverCutOff;};
  /// Parameter of the second order weighting function
  double pWei2Par(){return Gen->WFuncStraight2;};
  /// Parameter of the third order weighting function
  double pWei3Par(){return Gen->WFuncStraight3;};
  /// Number of steps
  int pStep();
  /// Number of particle
  int pNPart();
  /// Number of chain
  int pNChain();
  /// Number of chain
  int pNChain(int b);
  /// Number of particle per chain
  int pNPCh();
  /// Number of particle per chain
  int pNPCh(int c);
  /// # of types of the particle
  int pNType();
  /// Maximum number of bonds
  int pNLink();
  ///  Number of nanoparticles
  int pNNano();
  /// Number of blocks
  int pNBlock();
  /// Allocated number of particles
  int pNAllocP();
  /// Allocated number of chains
  int pNAllocC();
  /// Set Edge
  void SetEdge(double Val,int d){Gen->Edge[d] = Val;Gen->InvEdge[d] = 1./Val;};
  /// Set Edge
  void SetCNorm(int d){
    if(d<0 || d > 2) return;
    CNorm = d;
    CLat1 = (d+1)%3;
    CLat2 = (d+2)%3;
  };
  /// Set scale factor
  void SetScaleF(double *Scale){
    for(int d=0;d<3;d++){
      ScaleF[d] = Scale[d];
    }
  }
  /// Set reference pos
  void SetShiftPos(double *RefPos){
    for(int d=0;d<3;d++){
      ShiftPos[d] = RefPos[d];
    }
  }
  void SetIfNormalize(int If){
    IfNormalize = If;
  }
  /// Set DeltaT
  void SetDeltat(double Val){Gen->Deltat = Val;};
  /// Set Step
  void SetStep(int Val){Gen->Step = Val;};
  /// Set Temperature
  void SetTemp(double Val){Gen->Temp = Val;Gen->Beta = 1./Val;};
  /// Set Time
  void SetTime(double Val){Gen->Time = Val;};
   /// Increment Step
  void IncrStep(){Gen->Step++;};
  /// Return back folded position
  double pPos(int p,int d);
  /// Return back folded position
  double pChPos(int p,int d);
  /// Return back folded position
  void pPos(int p,double *Pos);
  /// Print the particle position
  double *pPos(int p);
  /// Return the velocity
  double pPosNoBkf(int p,int d);
  /// Return the velocity
  double pVel(int p,int d);
  /// Set the particle position
  void SetPos(int p, double *Pos);
  /// Set the particle position
  void SetPos(int p,int d,double Pos);
  /// Set the particle velocity
  void SetVel(int p, double *Vel);
  /// Set the particle type
  void SetType(int p,int t);
  /// Return the type
  int pType(int p);
  /// Return the chain
  int pChain(int p);
  /// Return back folded nano position
  double pNanoPos(int n,int d);
  /// Set the back folded array for the particle p
  void SetBkf(int p);
  /// Set the back folded array for the nano n
  void SetNanoBkf(int n);
  /// Print a position
  void pPos(double *Pos);
  /// What to draw
  char cWhat2Draw[STRSIZE];
  /// Extra particle
  NANO *Nano;
  /// Particle information of all particle
  PART *Pm;
  /// Array of linking between the particles
  LINKS *Ln;
  /// Information on all chains
  CHAIN *Ch;
  /// Soft bodies
  SOFT *Soft;
  /// Information for every block
  BLOCK *Block;
  /// Particle position/density on the square lattice
  double *PlotMem;
  /// Reference position
  double ShiftPos[3];
  /// Scale factor
  double ScaleF[3];
  ///  Number of soft bodies
  int NSoft;
  /// Number of particle to be considered in the radial density profile
  int NPartNearSphere;
  /// Additional homopolymer chains into the membrane
  int NAddChain;
  /// Additional cholesterol chains into the membrane
  int NAddChol;
  /// Solvent molecules
  int NSolvent;
  /// Stuffing for the cylinder
  int NStuffing;
  /// Normal coordinate
  int CNorm;
  /// lateral coordinate
  int CLat1;
  /// lateral coordinate
  int CLat2;
  /// Type of chain selected
  int NChType;
  /// Type of particle selected
  int NPType;
  /// Number of particles per edge
  int NEdge;
  /// Contains the definition of the system
  int SysType;
  /// Contains the definition of the file format
  int SysFormat;
  /// Contains the information for the creation
  int SysCreate;
  /// If normalize the lateral dimensions to one
  int IfNormalize;
  /// If PlotMem is allocated and filled
  int IfPlotMem;
};
#endif //VARIABILI_H
