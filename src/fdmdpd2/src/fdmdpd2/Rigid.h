/* $Id: Rigid.h 366 2013-08-08 14:19:08Z fuhrmans $ */

#ifndef RIGID_H
#define RIGID_H
#include "common.h"

#define TENS_INTER 0
#define TENS_INTRA 1

#define VAR_IF_TYPE(c,t) ( ((c)&(t))==(t) )
#define VAR_ADD_TYPE(c,t) ( (c)|=(t) )
#define VAR_REM_TYPE(c,t) ( (c)&=(~t) )

#define ABS(A) ((A)>0. ? (A) : -(A))

#ifndef CUBO_H
#define CUBO_H

static int IfTens = 0;

/// The previuos and consecutive particle in the cell linked list
struct DOM_PART{
  /// Next particle in the list
  int Next;
  /// Previous particle in the list
  int Prev;
  /// Coordination number
  int Coord;
};
/** For every cube how many particles and which are the first and the last
 */
struct DOM_CELL{
  /// Number of particles in the cell
  int NPart;
  /// First particle in the cell
  int First;
  /// Last particle in the cell
  int Last;
  /// First pointed particle;
  int Curr1;
  /// Second pointed particle;
  int Curr2;
};
struct DOM_DEC{
  /// Module 10 (nx,ny,nz) -> nc
  int Mod10[3];
  /// Number of cells
  int NCell;
  /// Number of cells per direction
  int nCells[3];
  /// Box sizes
  double Edge[3];
  /// Cell CutOff
  double CutOff;
  /// Number of allocated part
  int NAllocP;
  /// Number of particles and iterators per cell
  struct DOM_CELL *DCell;
  /// List of position of the particles
  struct DOM_PART *Pc;
} ;
typedef struct DIFF_PROF{
  // former positions
  double *OldPos;
  // density normalization
  double *VolSlabInv;
  // density normalization
  double *Count;
  // int particle origin
  int *PartOr;
  // maximum radius
  double BoxRad;
  // initial step
  int InitStep;
  // number of bins
  int NBin;
  // number of steps
  int NStep;
  // enable
  int IfEnable;
};
/// Allocate and fill the cells
int DomDecCutOff(struct DOM_DEC *Dc,double Edge[3],double CutOffCell,VEC2 *Pm,int NPart,double CutOff);
/// Allocate and fill the cells
int DomDec(struct DOM_DEC *Dc,double EdgeExt[3],int *NSect,VEC2 *Pm,int NPart,double CutOff);
/// Rebuild the pairlist
void DCRebuild(struct DOM_DEC *Dc,VEC2 *Pm,int NPart);
/// Add a particle to the cell c
int DCAddPart(struct DOM_DEC *Dc,const int p,double *Pos);
/// Remove a particle form the cell c
int DCRemPart(struct DOM_DEC *Dc,const int p,double *Pos);
/// Add a particle to the cell c
int DCAddPartC(struct DOM_DEC *Dc,const int p,const int c);
/// Remove a particle form the cell c
int DCRemPartC(struct DOM_DEC *Dc,const int p,const int c);
/// Move a particle form the cell c1 to the cell c2
int DCMovePart(struct DOM_DEC *Dc,const int p,double *OldPos,double *NewPos);
/// # part in the cell
int DCNPart(struct DOM_CELL *Cell,const int c);
/// First part in the cell
int DCFirst(struct DOM_CELL *Cell,const int c);
/// Next linked part 
int DCNext(struct DOM_PART *Pc,const int p);
/// Print the number of cells
int DCpNCell(struct DOM_DEC *Dc);
/// Set the coordination number
int DCSetCoorNumb(struct DOM_DEC *Dc,double *Pos,int p);
/// Get the number of the neighbouring cells
int DCGetNeighbours(struct DOM_DEC *Dc,int c,int p,int *NeiList);
/// (x,y,z)->n
int DCnCell(struct DOM_DEC *Dc,const double Pos[3]);
#endif //CUBO_H
enum NanoShape{
  // no nanoparticle
  SHAPE_NONE   = 0x00001,
  // spherical shape
  SHAPE_SPH    = 0x00002,
  // cylindrical shape
  SHAPE_CYL    = 0x00008,
  // tilting pill
  SHAPE_TILT   = 0x00010,
  // cylinder with spherical ends
  SHAPE_PILL   = 0x00020,
  // planar wall 
  SHAPE_WALL   = 0x00040,
  // external cylindrical field
  SHAPE_PORE   = 0x00080,
  // smoothed ends
  SHAPE_SMOOTH = 0x00100,
  // radius and height
  SHAPE_HEI    = 0x00200,
  // tip of a cantilever
  SHAPE_TIP    = 0x00400,
  // external harmonic potential
  SHAPE_EXT    = 0x00800,
  // Janus peptide
  SHAPE_JANUS  = 0x01000,
  // stalk (ref point)
  SHAPE_STALK  = 0x02000,
  // torus
  SHAPE_TORUS  = 0x04000,
  // harmonic
  SHAPE_HARM   = 0x08000,
  // harmonic
  SHAPE_UMBR   = 0x10000,
};
enum TensCalcMode{
  // 1d
  TENS_LINE = 1,
  // 1d projected on axis
  TENS_TILT,
  // 3d
  TENS_3D,
  // 2d 
  TENS_2D,
  // disabled
  TENS_NO,
  // calculate the tension and exit
  TENS_2D_LOOP,
};


typedef struct {double w,x,y,z;} QUAT;
/** structure of the external rigid field */
struct NANO
{
  FILE *Save;// 
  double Pos[3]; // Position
  double OldPos[3]; // Position
  double Vel[3]; // Velocity
  double Axis[3]; // Axis of the cylinder
  double Force[3]; // Sum of all forces
  double Pre[3];// Three component of the pressure
  double AVel[3]; //Angular momentum
  double AMom[3]; // Angolar momentum
  double AMomTemp[3]; // Angolar momentum
  double ForceRMS[4]; // Single components of the forces
  double v2; //Square velocities
  double Pot[2]; //Potential energy/chemical potential
  double Gamma; // Dissipation constant
  double Zeta;// Prefactor for the random noise
  double Mass;// Mass (set to 1)
  double Hamaker;// Strength of the interaction
  double Radius;// Radius
  double Height;// Height (in case of a cylindrical shape)
  double SwitchOff;// Switch off the interactions
  double Viscosity;// Viscosity of the "fluid"
  double Coating;// position of the well
  double Baseline;// value of the potential at the cutoff
  double CutOff;// Cutoff of the interactions
  double DistThr;// Maximum Distance where the potential is over the threshold
  double ForThr;// Threshold for the forces
  double PotThr;// Threshold for the forces
  int n; /* number of nano particles */
  int Inter; // number of interaction
  int Shape; // 0 none 1 spherical 2 cylindrical 3 wall
  int NCircle; //Number of monomers on the basis of the cylinder
  int NSide; //Number of monomers on the side of the cylinder
};
struct PEPTIDE{
  char ArchFile[60];
  int *Link; //Array of the linked particles
  double *Elong;//Rest distance of the linked particles
  double *kSpr; //Spring coupling
  double Pos[3];//Position of the center of mass
  double Vel[3];//Velocity of the peptide
  double Force[3];//Force on the peptide
  double Pot[3];//Potential on the peptide
  double Rad;//Radius
  double Hei;//Height
  double kEl;//Elastic constant
  int Id;//
  int NLink;//Number of links
  int NSide;//Number of monomers per side
  int NCircle;//Number of monomers per circle
  int NInter;//Number of 
};
struct EXTERNAL{
  char ArchFile[60];
  int *Link; //Array of the linked particles
  double *Elong;//Rest distance of the linked particles
  double *kSpr; //Spring coupling
  double kEl;//Elastic constant
  int Id;//
  int NLink;//Number of links
};


/** structure of the tension profile calculation. It writes NComp files which are selected for the direction to which the slabs are perpendicular and the component of the stress tensor */
struct TENS_PROF
{
  /** Array of NLine slabs for the tension profile */
  double **Tension;
  /** Array of summed contributions */
  double **Count;
  /** Slab volume contribution (inverted) */
  double *VolSlabInv;
  /** Rotation matrix to a specific orientation of the cyl */
  double Rot[9];
  /** Reference point for the radial distribution */
  double RefPos[3];
  /** Reference axis for the radial distribution */
  double RefAxis[3];
  /** Which coordinate has to be wrapped */
  int Wrap[3];
  /** Slab volume in one direction */
  double VolSlabInvPerp;
  /** Number of executed sums */
  double FNorma;
  /** Length of the maximal radial distance */
  double BoxRad;
  /** Number of snapshots to average */
  double NAverage;
  /** Coexistance density */
  double RhoCoex;
  /** Force and slab direction */
  int *Coord;
  /** Number of slabs */
  int NSlab;
  /** Number of direction */
  int NComp;
  /** Number of types */
  int NType;
  /** Number of types */
  int NTypeMax;
  /** Total number o NSlab arrays */
  int NLine;
  /** Reference system */
  int Ref;
  /** Step for the temporary file */
  int Step;
  /** Mode of calculation */
  int CalcMode;
  /** Calculate the tension and exit */
  int IfLoop;
};
//In GeomTrans.c
extern double Dist2AxisPerp(double *Pos,double *Axis, double *Perp);
extern void VectProd(double *v,double *u,double *Resp);
double AngleVect(double *u,double *v);
void RotMatr9(double *M,double *Axis,double Angle);
extern void Normalize(double *Axis);
extern double VectAngle(double *u,double *v);
extern double NormVect(double *Axis);
extern QUAT Quaternion(double *Axis,double Angle);
extern void RotMatr16(double *M,double *Axis,double Angle);
struct WIDOM{
  // Distribution of particles in the system
  double *Distr;
  // Variance for the creation of gaussian chains
  double Var;
  // Chemical potential of the chains
  double ChemPot;
  // Edges of the insertion box
  double Edge[3];
  // Maximum step for the insertion
  int NStep;
  // Switch between insertion, deletion, mc ...
  int CalcMode;
  // Number of allocated beads
  int NAlloc;
  // Number of bins for the distribution
  int NBins;
  // Diblock limit of a chain
  int DiblockLim;
  // Domain decomposition
  struct DOM_DEC *Dc;
};


#endif //RIGID_H
