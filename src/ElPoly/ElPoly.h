#ifndef ELPOLY_H
#define ELPOLY_H
#include "../include/VarData.h"
#include "../include/Cubo.h"
//#include "../include/Draw.h"
#ifdef USE_GL
#include <GL/glut.h>
#endif//USE_GL
#include "../include/Draw.h"
#include "../include/SingProc.h"
/// Which Draw function should be called
enum ElVisualize{
  /// Call DrParticle
  EL_PART  =  0x0000,
  /// Call DrQuad
  EL_QUAD    ,
  /// Call DrChain
  EL_CHAIN   ,
  /// Call DrVector
  EL_VECTORS   ,
  /// Call DrDensity
  EL_DENS    ,
  /// Call DrQuad1
  EL_QUAD1   ,
  /// Call DrPotential
  EL_POTENTIAL   ,
  /// Call DrSurface
  EL_SURF    ,
  /// Call DrCross (crosslink)
  EL_CROSS   ,
  /// Call DrIsoipse
  EL_ISOIPSE ,
  /// Call DrSample
  EL_SAMPLE  ,
  /// Call DrSpectrum
  EL_SPECTRUM,
  /// Call DrColor 
  EL_COLOR,
  /// Call DrVoronoi
  EL_VORONOI,
  /// Create a skin
  EL_SKIN,
  /// Stalk view
  EL_STALK,
  /// Call DrTriangulate
  EL_TRIA,
  /// Call DrGenMesh
  EL_MESH,
  /// Call DrMarchingCubes
  EL_ISOLEVEL,
  /// Call DrPolygon
  EL_POLYGON,
  /// Call DrSquareMesh
  EL_SQUAREMESH,
  /// Call OpenGl output
  EL_OPENGL,
  /// Call Povray output
  EL_POVRAY
};
/** Operates and visualizes the information in the VarData class */
class ElPoly : public VarData{
 private:
  /// Lists of file to open
  char **cFile;
 public:
  /// Constructor
  ElPoly(int argc,char **argv,int *FilePos);
  /// Desctructor
  ~ElPoly();
  //-----------------ElPolyProfile.cpp------------------
   /// Calculate the density profile \comm dens
  int DensProf(int NBin,int NSample,int Coord);
   /// Calculate the density profile projecting the particles on the normal density 
  void DensProfNormalSlab(int NBin,int NSample,int Coord);
  /// Density and thickness profile arond the nanoparticle
  int Diff2Files(int NBin,int How);
  /// Pressure difference between the virial and the ideal gas term
  void RestPress(int NBin);
  /// Thickness and density profile along Coord1 
  void SlabProf(int NBin,int nNano,int Coord1);
  /// Print the density profile in the surfaces representation
  void PrintDens(FILE *FileToWrite,double **Plot,double *LatDim,int NBin);
  /// Calculate the thickness from the radial density profile
  void RadDens2Thick(int NBin);
  /// Calculate the thickness from the radial density profile
  void RadDens2Thick2d(int NBin);
  /// Spatial distribution of the temperature 
  int Temperature(int NBin,int Coord);
  /** 3d (rad,norma,dens) density profile in cartesian coordinates*/
  void CartDens(int NBin,int nNano);
  /** 3d (rad,norma,dens) density profile with respect to a initial position \comm  radNano radCm radCmN */
  void RadDistrF(int NBin,int How,int nNano);
  /// Distribution of the bond lengths
  void BondDistr(int NSample);
  /// Distribution of end to end distances
  void E2EDistr(int NSample);
  /// Distribution of end to end distances
  void SplayDistr(int NSample);
  /// Outer shell of a projected (rad,normal) graph
  int RadialShell(int NBin);
  /// Diffusion and density profile around the nanoparticle \comm nano
  int NanoParticle(int NBin);
  /// Density profile of a non straight worm shape membrane \comm worm
  int WormF(int Partition,int NBin);
  /// Separate the leaflets
  void StalkF(int NSample);
  /// Calculate the thickness profile from the density plot
  void ThickFromDens(int NBin);
  /// Distribution of the areas around the protein
  void AreaDistrF(int NBin);
  /// Radial profile of the normal position of the particles
  void RadNormPos(int NBin,int NGrid);
  /// Radial summation of the 3d tension profile
  int SurfTens(int NBin);
  /// Sum more tension profile files
  void SumTens();
  /// Trace of the pressure profile
  int PressTrace();
  /// Contour plot around the inclusion
  int PressRadial();
  /// Change the pressure profile from cartesian to radial
  int Tens2dCartRad();
  /// Write the line describing a linear stalk
  void StalkLineProfF(int NBin);
  /// Dummy function to operate on a sequence of files
  void Prova();
  //-----------------ElPolyMeasure.cpp------------------
  /// 1d pair correlation
  void PairCorr(int NBin,int NDim);
  /// 1d, 2d pair correlation of the chains \comm pairRound pairSquare (obsolete)
  int PairCorrelationF(int NBin,int How);
  /// Diffusivity coefficient of the chains \comm diff
  int Diffusivity();
  /// Diffusivity of particles starting from an initial slab
  void DiffSlab(int NSlab);
  /// 2d Scattering \comm scatt scatt2
  int ScatteringF(int NBin,int How);
  /// Calculates the number of chain per area \comm area
  int NChainPSquareF();
  /// Calculates the number of chain per area \comm area
  void AreaCompr(int NSample);
  /// Calculates the 1d, 2d spectrum \comm spe
  int SpectrumF(int NBin);
  /// Midplanes for a sequence of snapshots
  void Midplane(int NBin);
  /// Calculates the 1d, 2d spectrum \comm spe
  void SpectrumMidplane(int NBin);
  /// Read the header and average the information
  void HeaderAverage(int nNano);
  /// Create a system file with a chain less and the information about the missing chain
  void WidomOut(char *InFile,int NBin);
  /// Create NCh systems with a chain less
  void WidomOut();
  /// Create the histogram of the energy of adding lipids
  void WidomIn(char *InFile,int NBin);
  /// Create a set of files with a chain more
  void WidomIn();
  /// Write the end to end distance of the chains
  void End2EndDistr(char *OutFile);
  /// Write the bond distance between the monomers
  void BondDistr(char *OutFile,int NBin);
  /// Decoupling of the direction/end to end distance of the chains
  void Decoupling(int What);
  /// Measure the elastic coupling between the two sheets
  void ElasticCoupling(int NSample);
  /// Measure the elastic coupling between the two sheets
  void ElasticCouplingNVT();
   /// Measure the elastic coupling between the two sheets
  void BilayerDistance(char *OutFile,int NBin);
  /// Distribution of the end to end distances
  void EndToEndDist();
   //-----------------ElPolyRepr.cpp---------------------
  /// Projection of the system against one direction \comm pro
  int ProjectionF(int NBin,int Coord);
  /// Sampling of the three dimentional space
  int CoreF(int NBin,int How);
  /// Area of a surfuce projected on the coordinate Coord \comm surf
  int Surface(int NBin,int Coord);
  /// Project the velocities of a 3d system on a 2d system wrt a coordinate
  int From3To2d(int Coord,double Param);
  /// Projet a 2d system on one of the two coordinates
  int From2To1d(int Coord);
  /// Project a 3d system on one coordinate
  int From3To1d(int Coord);
  /// Density plot all over different snapshots and calculation of the isolevel surface 
  void IsoSurf(int NSample,double *IsoValue,int NIso);
  /// Calculate the discrete density of the system and the correspondent isolines
  void IsoLine(int NSample,double *IsoValue,int NIso,int How);
  /// Perform the marching cubes on a density plot and print the isolines
  void IsoLine(FILE *F2Write,double *Plot,int NSample,double *IsoLevel,int NIso);
  /// Write the position of the stalk for every snapshot
  void FetchStalk();
  /// Write the position of the pore for every snapshot
  void FetchPore();
  /// Area of hydrophobic in the torus
  void StalkArea();
  /// Average the postion of the lipids over many snapshots
  void AvSnap();
  /// Radial profiles of the slab density between 0⁰ and 90⁰ subdivided in NAngle angles
  void SlabAngleProfs(int NBin,int NAngle,int Coord);
  /// Sample the space in NSample lattice points
  void Sample(int NSample);
  //-----------------ElPolyOutput.cpp-------------------
  /// Prepare a countor plot for tecplot
  void Conv2Tecplot(int NBin,int How);
  /// esport the data in vmd file format
  void Conv2Vmd();  
  /// esport the data in pov file format for rendering
  void Conv2Povray();  
  /// esport the data in radius depth density file format
  void Conv2rzd(int NSample);
  /// esport the data in radius depth density file format
  void Conv2xyzd(int NSample);
  /// Draw a scalar field
  void DrField(int NGrid,double IsoLevel,int nNano,FILE *FWrite);
  /// Draw the bonds
  void DrBondPovRay(double *Pos1,double *Pos2,float *Color);
  /// Print the header for povray
  void HeaderPovRay();
  /// PovRay draw function
  void DrNanoPovRay(int n);
  /// PovRay draw function
  void DrPartPovRay(int p);
  /// Convert into a square lattice
  void ConvLattice(int NSample,char *FName);
  //-----------------ElPolyEl.cpp-----------------------
  /// Remove chains satisfying a condition
  void RemoveChains();
  /// Calculation of the contact angle \comm angle Boh?
  int Angle(int NBin);
  /// Definition of the contact angle function Boh
  double ContactAngle(double x);
  /// Averaged center of mass Boh
  int CenterOfMass(int Coord);
  /// Defines the first and the last file to be elaborated \comm file
  int ChangeFile();
  /// Defines the normal coordinate \comm coord
  int SpecifyCoord();
  /// Opens the f file of the list \comm open
  int OpenFile(int f);
  /// Opens a file
  int OpenFile(char *FileName);
  /// Calculates some properties of the system \comm prop
  int PropertiesF();
  /** Convert the internal definition for the menu of ElPoly in string */
  char *ChooseDraw(int ExtWhat2Draw);
  /** Assing the correct value of What2Draw from \param String */
  void ChooseDraw(char *String);
  /** Reorder the LIPID block in four different layers */
  void DivideLayers(int How);
  /// Set the initial and final number of files 
  void SetBoundFile(int InitFile,int EndFile);
  /// Information on the current file elaborated
  void Processing(int f);
  /// Shift the system to the center
  void Shift2Center();
  /// Set backfold type
  void SetBackFold(int Bf){
    NBackFold = Bf;
  }
  /// Set backfold type
  void SetNVisSkip(int NSkip){
    if(NSkip > 0){
      NVisSkip = NSkip;
    }
  }
#ifdef __glut_h__
  //---------------ElPolyDraw.cpp----------------------
  /// Choose the visualisation (obsolete)
  void RenderPart(void);
  /// Draws in run time (expensive)
  void DrRunTime();
  /// Visualisation of the particle system with velocity color scheme
  void DrColor();
  /// Create a list of all the particles
  void DrPartList();
  /// Draw a vector
  void DrVector(Vettore v,Vettore Origin);
  /// Draw the chains as vectors
  void DrVectors();
  /// Visualize only particles within a certain distance
  int DrIntorno(int p,double Blue);
  /// Visualize the cross links
  void DrCrossLinks();
  /// Boh
  int Graphics(int argc,char **argv);
  /// Sequence of pictures
  void ESlide();
  /// Sequence of pictures
  void ESlide1();
  /// Load and draw the linking of a protein
  void DrProtein(const char *FileName,int nBlock);
  /// Draw the nanoparticle structure
  void DrNano();
  /// Draw the cross linked particles
  void DrCrossLinks(char *FileName);
  /// Data type for distance/field functions
  typedef void(ElPoly::*DRAW_PART)(int p);
  /// Pointer to a distance/field function
  DRAW_PART Draw_Part;
  /// Pointer to a generic function
  void DrawPart(int p){(*this.*Draw_Part)(p);}
  /// OpenGl draw function
  void DrPartOpenGl(int p);
  /// Pointer to a distance
  DRAW_PART Draw_Nano;
  /// Pointer to a generic function
  void DrawNano(int n){(*this.*Draw_Nano)(n);}
  /// OpenGl draw function
  void DrNanoOpenGl(int n);
  /// Assign the pointer to the corrispondent draw function, intilize the visualization
  void DrawFuncHeader();
  /// Dnd the visualization
  void DrawFuncFooter();
  /// Data type for distance/field functions
  typedef void(ElPoly::*DRAW_BOND)(double *Pos1,double *Pos2,float *Color);
  /// Pointer to a distance/field function
  DRAW_BOND Draw_Bond;
  /// Finde the neighbour of the particle p and draw the bond
  void DrBond(int p);
  /// Pointer to a generic function
  void DrawBond(double *Pos1,double *Pos2,float *Color){(*this.*Draw_Bond)(Pos1,Pos2,Color);}
  /// Null function
  void DrBondNo(double *Pos1,double *Pos2,float *Color){};
  /// Draw the bonds
  void DrBondOpenGl(double *Pos1,double *Pos2,float *Color);
  //---------------ElPolyDrawControl.cpp----------------------
  /// Draws all the particles and bonds
  void keyboard(unsigned char key,int x, int y);
  /* /// Definition of the key bindings */
  /* void keyboardDraw(unsigned char key); */
  /// Scaling
  void ElDrawMouse(int button,int state,int x,int y);
  /// Compile the list with some useful primitives
  void CompileList();
  /// Boh
  void DrPosCol();
  /// Creates the menu
  void Menu();
  /// Choose what to visualize
  void ElMenuChoise(int option);
  /// Choose what to visualize  
  void ElMenuVisual(int option);
  /// Creates the menu
  void CreateMenu();
  //---------------ElPolyDrawSurf.cpp----------------------
  /// Draws all the same quotes surfaces
  void DrIsoipse(int NBin,int NIsoipse,int CoordN);
  /// Draws the surface like a sheet using the chains position
  void DrSurface();
  /// Draw the surface in triangles
  void DrSmooth(double *Plot,int NSample,double Min,double Max);
  /// Samples the surface and draws it in triangles
  void DrSample(int NSample);
  /// Calculate the spectrum of the surface and draws it in triangles
  void DrSpectrum();
  /// Samples the surface, applies a 2d matrix derivative and draws the surface in triangles
  void DrDerivative();
  /// Position of every chain in hexagons
  void DrChains();
  /// Every particle is traeted as a vertex
  void DrPolygon();
  /// Call a square in every particle position 
  void DrQuad();
  /// (rad,norm,dens) visualisation
  void DrDensity();
  /// Intensity scheme   
  void DrQuad1();
  /// Potential function
  void DrPotential();
  /// Draw the outer shell
  void DrShell();
  /// Voronoi tassellation
  void DrVoronoi();
  /// Boh
  void DrInterpSurface();
  /// Defines the triangles at the boundaries of the density close to the IsoLevel value
  void DrIsolevel(int NSample,double IsoLevel);
  /// Draw a scalar field
  void DrField(int NGrid,double IsoLevel,int nNano);
  /// Draw the polygons for a square mesh
  void DrSquareMesh();
  /// Draw the normal to a point
  void DrNormalPoint(int p,int NEdge);
  //-----------------------------------------------
  GLuint Quad,Point,*Cylinder,MetalCylinder,Hexagon,Cube,Arrow,GlWall;
  /// Saturation of the color (increase intensity)
  double Saturation;
  /// Name of the block to draw
  char Block2Draw[20];
  /// single triangle defined by three vectors
  void DrTria(Vettore *v00,Vettore *v01,Vettore *v11,Vettore *vN);
  /// Draw the contour of three vectors
  void DrTriaContour(Vettore *v00,Vettore *v01,Vettore *v11);
  /// Double triangle defined by four vectors
  void DrDoubleTria(Vettore *v00,Vettore *v01,Vettore *v11,Vettore *v10,Vettore *vN);
  /// Visualize the surface calculated in Stalk()
  void DrCreateStalk();
  /// Visualize the surface from Stalk.xvl
  void DrStalk();
   /// Cover a regular square grid of points with tiles
  void Tile();
  //----------------------------ElPolyDrawCGAL----------------
  /// Triangulate a surface
  void DrTriangulate();
  /// Build a mesh from the lipid positions
  void DrMesh();
  /// Generate a mesh from a function
  void DrGenMesh();
  /// Construct cells from the lipid positions
  void DrCells();
  /// Find the covering surface for given points
  void DefineSkin(int NSample);
  /// Find the covering surface for given points
  void DefineSurf();
#endif //__glut_h__
  /// First and last file of the list
  int NFile[2];
  /// Total number of file
  int NFileTot;
  /// Boh
  int NPro ;
  /// Current number of the file list
  int quando;
  //  int NBin;
  /// Line size of the gl
  int LineSize;
  /// If it visualizes only the particles/chains within a certain distance
  int IfIntorno;
  /// Boh
  //int IfColour;
  /// Boh
  //int IfChains;
  /// If it draws the bonds
  int IfLine;
  /// Type of the chains to be visualized
  int IfChType;
  /// Type of backfold
  int NBackFold;
  /// Draw the content in the appropriate visulization
  int  What2Draw;
  /// How many lipids are skipped in the visualization
  int NVisSkip;
  /// Boh
  //int *ChainTypes;
  /// Draw output flag
  int DrawOutput;
  /// Boh
  //double Diameter;
  /// Boh
  //double StepDiameter;
  /// Boh
  //double ExtraDiam;
  /// Distance from the nanoparticle
  double Vicinanze;
  /// Define the shrink factor between the box edges
  double InvScaleUn;
  /// External parameter (e.g. for MarchingCubes)
  double ExtParam;
  /// Normal scaling factor (z zooming)
  double ScaleFact;
  /// Output file for drawing
  FILE *DrawOutFile;
#ifdef USE_MPI
  SingProc *Proc;
#endif
};
static int IfImage;

extern void Slide();
#ifdef __glut_h__
#endif//__glut_h__
#endif //ELPOLY_H
