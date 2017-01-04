#ifndef FORCES_H
#define FORCES_H
#include <VarData.h>
#include <Cubo.h>
#include <SingProc.h>
#include <time.h>
#ifdef USE_GL
#include <GL/glut.h>
#endif//USE_GL
/** Describe a different shape for allocating and setting different environment variables */
enum SYS_SHAPE{
  /// Beads allined on a chain
  SYS_1D =     0x0001,
  /// Rectangular square lattice
  SYS_2D =     0x0002,
  /// Cubic square lattice
  SYS_3D =     0x0004,
  /// Two 1d coupled strings
  SYS_LEAVES = 0x0008,
  /// Rotating rigid bodies
  SYS_RIGID  = 0x0010,
  /// 1d representation of a stalk
  SYS_STALK  = 0x0020,
  /// Perfom a md simulation
  SYS_MD     = 0x0040,
  /// Perform a mc simulation
  SYS_MC     = 0x0080,
  /// Perform a Widom insertion and deletion
  SYS_WIDOM  = 0x0100,
  /// Build a trial system
  SYS_TRIAL  = 0x0200,
  /// 1d representation of a pore
  SYS_PORE   = 0x0400,
  /// 1d representation of a pore
  SYS_ROD   = 0x0800,
  /// 1d electric line with branches
  SYS_ELECTRO   = 0x1000,
};
/** Describe a different calculation mode for allocating and setting different environment variables */
enum CALC_MODE{
  /// NVT ensemble
  CALC_NVT     = 0x00001,
  /// mVT ensemble (m chemical potential)
  CALC_mVT     = 0x00002,
  /// Solution of the differential equation
  CALC_SOLVE   = 0x00004,
  /// Pair potential
  CALC_PAIR    = 0x00008,
  /// Harmonic potential
  CALC_HARM    = 0x00010,
  /// Integrated Lennard-Jones potential
  CALC_LJ39    = 0x00020,
  /// Lennard-Jones potential
  CALC_LJ      = 0x00040,
  /// MC density functional hamiltonian
  CALC_DENS_CH = 0x00080,
  /// NcVT ensemble (Nc number of chains)
  CALC_NcVT    = 0x00100,
  /// mcVT ensemble (mc chemical potential of a chain)
  CALC_mcVT    = 0x00200,
  /// Calc 3d pressure profile (cartesian)
  CALC_3d      = 0x00400,
  /// Calc 2d pressure profile (polar)
  CALC_2d      = 0x00800,
  /// MD density functional Hamiltonian
  CALC_DENS    = 0x01000,
  /// Introduce the configurational bias
  CALC_CONF_BIAS = 0x02000,
  /// Introduce the bias on a planar region
  CALC_BIL_BIAS  = 0x04000,
  /// Introduce the bias on a spherical region
  CALC_SPH_BIAS  = 0x08000,
  /// Use a step potential
  CALC_STEP  = 0x10000,
  /// Use a potential for the electrical line
  CALC_ELECTRO  = 0x20000,
};
/** Thermostat flag */
enum THERM_MODE{
  /// no thermostat
  THERM_NO      = 0x0001,
  /// Langevin thermostat
  THERM_LANG    = 0x0002,
  /// Andersen thermostat
  THERM_AND     = 0x0003,
  /// Berendsen thermostat
  THERM_BERE    = 0x0004,
};
/// Interpolation method
enum FIT_TYPE{
  /// Cubic spline
  FIT_SPLINE3 = 0,
  /// Forth order spline
  FIT_SPLINE4 = 1,
  /// Parabolic interpolation
  FIT_PARAB2  = 2,
  /// Cubic interpolation
  FIT_CUBIC   = 3,
  /// 
  FIT_DERIV   = -1,
  /// Finite differences
  FIT_FORTH   = 4,
  /// Polynomial interpolation
  FIT_POLY    = 5,
  /// BSplie interpolation
  FIT_BSPLINE = 6,
};
/// Integration method
enum IntMethod{
  /// MD integration
  INT_MD = 0,
  /// Solution of the differential equation
  INT_DIFF = 1,
  /// MC update
  INT_MC = 2,
};
/// Quantities allocated
enum ALLOCATED{
  /// struct FORCES allocated
  ALL_FORCES = 0x001,
  /// Tabulated force field allocated
  ALL_METR   = 0x002,
  /// Potential force field allocated
  ALL_POT    = 0x004,
  /// MC allocated
  ALL_MC     = 0x008,
  /// MD allocated
  ALL_MD     = 0x010,
  /// Spline interpoints allocated
  ALL_SPLINE = 0x020,
  /// Rigid inclusion allocated
  ALL_CYL    = 0x040,
  /// ?
  ALL_FFIELD = 0x080,
  /// struct TENS allocated
  ALL_TENS   = 0x100,
  /// Local densities allocated
  ALL_DENS   = 0x200,
  /// Energy and position of the trial bonds
  ALL_BIAS   = 0x400,
  /// Energy and position of the trial bonds
  ALL_MATRIX = 0x800,
};
/// Prefactors of the forces
typedef struct{
  /// Prefactor of the Laplacian
  double Lap;
  /// Prefactor of the square laplacian
  double SLap;
  /// Elastic force
  double El[3];
  /// External force
  double Ext;
  /// Lennard-Jones
  double LJ;
  /// Contact/Friction
  double Cont;
  /// Elongation of the springs
  double Elong[3];
  /// Lennard-Jones
  double LJMin;
  /// Lennard-Jones
  double LJCutOff;
  /// CutOff of the lennard jones
  double CutOff2;
  /// Baseline of the potential
  double BaseLine;
  /// Maximum potential allowed
  double PotThr;
  /// Maximum force allowed
  double ForThr;
  /// Correspondent radial distance for this force
  double DistThr;
} KFORCES;
/// Single contribution of the forces
typedef struct{
  /// Direction
  double Dir[3];
  /* /// Helfrich */
  /* double Hel[3]; */
  /// External
  double Ext[3];
  /* /// Elastic */
  /* double El[3]; */
  /* /// Lennard-Jones */
  /* double LJ[3]; */
} FORCES;
/// Spatial tensor profile
typedef struct{
  /// # of components
  int NComp;
  /// # of dimensions
  int NDim;
  /// # of slabs
  int NSlab;
  /// Calculation mode
  int CalcMode;
  /// Pressure array
  double **Pre;
  /// Density array
  double **Dens;
  /// Total pressure
  double PreTot[3];
  /// Reference point
  double RefPos[3];
  /// Edges
  double Edge[3];
  /// Which coordinate to wrap
  int Wrap[3];
  /// Inverse of the edges 
  double EdgeInv[3];
}TENS;
/// Alias to the preferred domain decomposition
//typedef DdDoubleLoop DomDec;
typedef DdLinkedList DomDec;
//typedef DdFixedSize DomDec;
//typedef DdArray DomDec;
/** This class performs the different steps to solve the 
   equations of motion via molecular dynamics simulation,
   simple grancanonical simulation in Monte Carlo and solution
   of 1d differential equations.
*/
class Forces : public VarData{
 private:
  /// Configuration file
  char ConfFile[60];
  /// Initial simulation time
  time_t InitTime;
  /// Current time
  time_t CurrTime;
  /// Density of every particle (second order)
  double *Dens2;
  /// Density of every particle (third order)
  double *Dens3;
  /// Tabulated potential
  double *PTab;
  /// Tabulated force
  double *FTab;
  /// Tabulated metropolis
  double *MTab;
  /// Stored energy of the system
  double OldNrgSys;
  /// Stored energy of the particle
  double *OldNrgBead;
  /// Stored energy of the chain
  double *OldNrgCh;
  /// Stored position of the chain
  double **OldPos;
  /// Stored local density second order
  double *LocDens2;
  /// Stored local density third order
  double *LocDens3;
  /// Distribution of the first bead
  double *FirstBeadDistr;
  /// Position of the bond
  double **BondPosBias;
  /// Cumulative probability of the bias
  double *CumProbBias;
  /// Borders for the linear bias
  double BorderBias[2];
  /// Binning for the distribution
  int NBin;
  /// Binning for the parameter space
  int NGrid;
  /// Number of trial per chain
  int NTrialBias;
  /// Bondary conditions on the box
  int BoundCond[6];
  /// Periodic image convention
  int PeriodicImage[3];
  /// Information for the statistics
  FILE *StatFile1;
  /// Information for the statistics
  FILE *StatFile2 ;
  /// Pointer type to the energy calculation function
  typedef double(Forces::*CALC_NRG)(int p,double *Pot);
  /// Pointer to the energy calculation function
  CALC_NRG NrgBead;
  /// Pointer to the energy calculation function
  CALC_NRG NrgCh;
  /// Choose a calculation mode
  //static CALC_NRG ChooseCalcMode(int Mode);
  void ChooseCalcMode(int Mode);
  /// Pointer type to the potential function
  typedef double(Forces::*CALC_POT)(double Dist,int t1,int t2,double *Pot);
  /// Pointer to the energy calculation function
  CALC_POT CalcPot;
  /// Choose a calculation mode
  //static CALC_POT ChoosePot(int Mode);
  void ChoosePot(int Mode);
  /// Number of points tabulated
  int NTab;
  /// Flags for the system configuration
  int DynFlag;
  /// A generic step in the simulation 
  typedef double(Forces::*STEP_CONF)();
  /// A generic step in the simulation 
  typedef double(Forces::*STEP_SIM)();
  /// Configure the initial system
  STEP_CONF ConfigSys();
  /// Calculate the interactions
  STEP_SIM CalcUpdate();
  /// Matrix with interactions
  Matrice *IntMatrix;
  public:
  //-----------Forces.cpp---------------------
  /// Create an initial system and choose the simulation thecnique
  Forces(int argc,char **argv,int NPart,char *ConfFile);
  /// Read an initial system and choose the simulation thecnique
  Forces(int argc,char **argv,char *ConfFileExt,char *Snapshot);
  /// Internal message
  void Shout(const char * s, ...);
  /// Allocated the structures needed for the corresponding simulation method
  void AllocMethod();
  /// Calculate some initial quantities for the succesive calculations
  void PrepareSys();
  /// Create the grid for the parallelisation
  void PrepareParallel(int argc,char **argv);
  /// Frees the memory
  ~Forces();
  /// System's info
  void Info();
  /// Read the config file
  int ReadConfDinamica(char *File);
  /// Read the config file
  void InitConst();
  /// Realloc the number of particles
  int ReSetNPart(int NewNPart);
  /// Realloc the number of chains
  int ReSetNChain(int NewNChain);
  /// Set a new number of particle per chain
  void ReSetNPCh(int NewNPCh);
  /// Open a new file
  void ReOpen(char *FName,int Bf);
  /// Fill the entries of the interaction matrix
  void FillMatrix();
  //-----------ForcesIntegration.cpp----------
  /** Obtain informations for a better performance in 
      inserting the chains */
  void StudySys();
  /// Minimum of the Helfrich Hamiltonian
  int MinHelfrich();
  /// Boh
  int DynIntegration();
  /// Sum up all the forces and update the positions
  void Dynamics();
  /// Solve a system of four oder differential equation
  void Solve();
  /// Solve a system of four oder differential equation of particles connected by links 
  void SolveLinks();
  /// Solve a system of four oder differential equation of particles connected by links 
  void SolveLinksSparse();
  /// Solve a system of four oder differential equation of particles connected by links 
  void SolveLinksIterative();
  /// Solve a system of four oder differential equation of particles connected by links 
  void SolveRod();
  /// Solve a system of four oder differential equation of particles in a line
  void SolveLeaves();
  /// Boh
  int Update();
  /// Velocity Verlet for a rigid body, first step
  void VelVerletRigid();
  /// Velocity Verlet for a rigid body, second step
  void VelVerletRigid2();
  /// Metropolis acceptance criterium
  int IfMetropolis(double ExpArg,double Weight);
  /// Insert a particle in the box
  int InsertBead(int p);
  /// Ignore a chain in the system (densities, pairlist)
  void IgnoreCh(int c);
  /// Delete a chain in the system
  void RemChFromSys(int c);
  /// Save the chain configuration
  void SaveCh(int c);
  /// Reinsert the chain in the previous position
  void ReInsertCh(int c);
  /// Consider a chain in the system (densities, pairlist)
  void ConsiderCh(int c);
  /// Displace a chain in the box
  double InsertCh(int c);
  /// Build the rest of the chain 
  double InsertRest(int pCurr,int StartPos);
  /// Remove a chain in the box
  double RemoveChBias(int c);
  /// Put a chain in the box
  double InsertChBias(int c);
  /// Weight a set of bonds for the configurational bias
  double WeightSetBond(int p,int t);
  /// Create a set of bonds for the configurational bias
  double CreateSetBond(int p,int t);
  /// Move a particle in the box
  int MoveBead(int p);
  /// Trial insertion
  int TryInsert();
  /// Trial removal
  int TryRemove();
  /// Trial movement
  int TryMove();
  /// Trial desplacement of a chain
  int TryMoveCh();
  /// Trial insertion of a chain
  int TryInsertCh();
  /// Trial removal of a chain
  int TryRemoveCh();
  /// Trial biased insertion of a chain
  int TryInsertChBias();
  /// Trial biased removal of a chain
  int TryRemoveChBias();
  /// Widom insertion
  void WidomInsert(double *NrgDiff);
  /// Widom removal
  void WidomRemove(double *NrgDiff,int p);
  /// Widom insertion
  void WidomInsertCh(double *NrgDiff);
  /// Widom removal
  void WidomRemoveCh(double *NrgDiff,int c);
  /// Widom with Rosenbluth weight
  void WidomBiasChIn(double *Weight);
  /// Widom with Rosenbluth weight
  void WidomBiasChOut(double *Weight,int c);
  /// First step of the velocity Verlet
  void VelVerlet1();
  /// Langevin thermostat
  void LangevinTherm();
  /// Andersen thermostat
  void AndersenTherm();
  /// Berendsen thermostat
  void BerendsenTherm();
  /// No thermostat
  void NoTherm(){};
  /// Second step of the velocity Verlet
  void VelVerlet2();
  /// Choose a calculation mode
  //static CALC_NRG ChooseCalcMode(int Mode);
  void ChooseThermostat(int Mode);
  /// Pointer type to the potential function
  typedef void (Forces::*CALC_THERM)();
  /// Pointer to the energy calculation function
  CALC_THERM CalcTherm;
  /// Pointer to the energy function
  void ApplyTherm(){return (*this.*CalcTherm)();};
  //-----------ForcesForceField.cpp-----------
  /// Helfrich Hamiltonian for a line
  int ForceFieldLine();
  /// Boh
  int SumSomeForces(int HowMany);
  /// Helfrich Hamiltonian with an elastic coupling
  int ForceFieldLeaves();
  /// Armonic potential on a lattice
  int ForceFieldBulk();
  /// Bending potential on a rod
  void ForceFieldRod();
  /// Interaction between rigid bodies
  void ForceFieldRigid();
  /// Calculate the forces for the density functional
  void CalcForcesDensFunc();
  /// Print the force and the potential
  void PrintForce();
  /// Calculate the force field summing every single contribution
  void GetForceField();
  /// Tabulate the values of the force
  void TabForceAlloc(int NTabExt);
  /// Tabulate the values of the potential
  void TabPot();
  /// Calculate all the interaction with the nano
  void NanoInteraction();
  /// Exchange energy with the nano
  double NanoNrg(int p);
  /// Exchange energy with the nano
  double NanoNrg(double *Pos,int t);
  /// Define the parameters for calculating the force
  void DefForceParam();
  /// Define the parameters for calculating the force
  void DefNanoForceParam();
  /// Force between to rigid bodies
  double RigidLJ(int nNano,double Dist,double *Pot,double Sign);
  /// Hamaker potential
  double RigidHamaker(int n,double Dist,double *Pot,double Sign);
  /// Force between to rigid bodies
  double RigidCoulomb(int nNano,double Dist,double *Pot,double Sign);
  /// Distance between two bodies
  double RigidDistanceRad(int n,int nn,double *dr);
  /// Distance between two bodies
  double RigidDistanceAxis(int n,int nn,double *dr);
  // Point to the potential 
  void PointShape(int iShape);
  /* /// Data type for distance/field functions */
  /* typedef double(Forces::*NANO_FORCE)(double Dist2,int n,int p); */
  /* /// Pointer to a distance/field function */
  /* NANO_FORCE Nano_Force; */
  /* /// Pointer to a generic function */
  /* double NanoForce(double Dist2,int n,int p){return (*this.*Nano_Force)(Dist2,n,p);} */
  /// Harmonic potential
  double Harmonic(double Dist2,int t1,int t2,double *Pot);
  /// Step potential
  double StepPot(double Dist2,int t1,int t2,double *Pot);
  ///  Potential for the electrical lines
  double ElectroPot(double Dist2,int t1,int t2,double *Pot);
  /// Integrated Lennard Jones potential
  double LJ39(double Dist2,int t1,int t2,double *Pot);
  /// Classical Lennard Jones potential
  double LJPot(double Dist2,int t1,int t2,double *Pot);
  /// Pointer to a potential
  double Potential(double Dist,int t1,int t2,double *Pot){return (*this.*CalcPot)(Dist,t1,t2,Pot);};  
  //-----------ForcesCalcNrg.cpp-----------
  /// Calculate and sum up the energy of the chains
  double CalcTotNrgCh();
  /// Calculate and sum up the energy of the part
  double CalcTotNrgBead();
  /// Pointer to the energy function
  double CalcNrgBead(int p,double *Pot){return (*this.*NrgBead)(p,Pot);};
  /// Pointer to the chain energy function
  double CalcNrgCh(int c,double *Pot){return (*this.*NrgCh)(c,Pot);};
  /// Calculate the bond and the density functional energies
  double NrgChBondDens(int c,double *Pot);
  /// Calculate the density functional energies
  double NrgChDens(int c,double *Pot);
  /// Calculate the non bonded interaction energy with the neighbouring particles
  double CalcPairwise(int p,double *Pot);
  /// Calculate the non bonded interaction energy with the neighbouring particles
  double CalcPairwiseCh(int c,double *Pot);
  /// Calculate the bonded interaction energy with the neighbouring particles
  double CalcBending(int p);
  /// Calculate the spring interaction energy with the neighbouring particles
  double CalcSpring(int p);
  /// Calculate the spring and the bonded interactions with the other monomers in the chain
  double CalcBonded(int p,double *Pot);
  /// Calculate the bending energy for a ghost particle
  double CalcBendingGhost(double *Pos,int pExt);
  /// Calculate the bonded and spring interaction in a cell
  double CalcBondedCh(int c,double *Pot);
  /* ///Calculate the spring, bending and non bonded interactions */
  /* double CalcNrgCh(int c,double *Pot); */
  /// Calculate the spring, bending and non bonded interactions and write it in OldNrgPm
  void CalcNrgBeadDensFunc();
  /// Calculate the spring, bending and non bonded interactions
  double CalcNrgBeadDensFunc(int p,double *Pot);
  /// Check if all the particles are taken in account
  double CheckDomDec(int p);
  /// Check the pair list
  void CheckPairList();
  /// Iterate all over the particles and calculate the forces
  double SumForcesMD();
  /// Cubic weighting function
  double Wei3(const double r, const double a);
  /// Derivative of the cubic weighting function
  double DerWei3(const double r, const double a);
  /// Quadratic weighting function
  double Wei2(const double r, const double b);
  /// Derivative of the quadratic weighting function
  double DerWei2(const double r, const double b);
  /// Calculation of the energy from the density functional Hamiltonian for the ghost particle
  double DensFuncNrgGhost(double *Pos,int p1,int t1);
  /// Calculation of the energy from the density functional Hamiltonian for the ghost particle
  double DensFuncNrgGhostInternal(double *Pos,int p1,int t1);
  /// Calculation of the energy from the density functional Hamiltonian for the particle p1
  double DensFuncNrgBead(int p1);
   /// Calculate the energy from the density functional Hamiltonian for the chain c
  double DensFuncNrgCh(int c,double *Pot);
   /// Calculate the average energy from the density functional Hamiltonian for the chain c
  double DensFuncNrgChAv(int c);
   /// Calculate the average energy from the density functional Hamiltonian for the chain c
  double DensFuncNrgChInternal(int c);
   /// Calculation of the energy from the density functional Hamiltonian for the system
  double DensFuncNrgSys();
  /// Calculate the local density for the particles between pInit and pEnd
  void CalcDens(int pInit,int pEnd);
  /// Set the local densities to zero
  void ClearDens();
  /// The energy is constant within the cutoff
  double NrgStep(int p);
  /// The energy per chain is constant within the cutoff
  double NrgStepCh(int c,double *Pot);
  /// Energy of an electric line 
  double NrgElectro(int p);
  /// Add the densities connected with the particles between pInit and pEnd
  int AddDens(int pInit,int pEnd);
  /// Substract the densities connected with the particles between pInit and pEnd
  int RemDens(int pInit,int pEnd);
  /** Sum the local density for the particles between pInit and pEnd
      and multiply the factors by the virial coefficients.
   */
  double SumDens(int pInit,int pEnd);
  /// List of the cells close to the particle position
  int ListNeiCell(int p,double *Pos,int *NeiList){
    return Pc->GetNei(Pos,NeiList);
    //return Pc->GetCell(p,NeiList);
  };
  //-----------ForcesBoundary.cpp-------------
  /// Move a particle
  void PullBead();
  /// Move a particle
  void PushBead();
  /// Select a particle
  void SelectBead(int p);
  /// Sinusoidal surface wave
  void Wave();
  /// Add a circle as a boundary condition
  void AddCircle(int nNano);
  /// Add a cylinder as a boundary condition
  void AddCylinder(int nNano);
  /// Add a pore as a boundary condition
  void AddPore(int nNano);
  /// Add all rigid bodies as a boundary condition
  void AddRigid();
  /// Height of the boundary condition depending on the direction
  double HeightBoundary(double *Pos,int dir);
  //-----------ForceCreate.cpp----------------
  /** Create an initial configuration and an appropriate force field */
  void CreateInitial();
  /// Create a plane of connected beads
  void  Create2d();
  /// Create a lattice of connected beads
  void Create3d();
  /// Create two connected sheets and add a protein
  void CreateLeaves();
  /// Create the 1d representation of a stalk
  void CreateStalk();
  /// Create the 1d representation of a pore
  void CreatePore();
  /// Create single line of connected monomers
  void Create1d();
  /// Create rigid bodies
  void CreateRigid();
  /// Create a initial disposition of particle for the MC sim
  void CreateMC();
  /// Create a initial disposition of particle for the MD sim
  void CreateMD();
  /// Create a initial disposition of particle for a stiff rod
  void CreateRod();
  /// Create a initial disposition of houses to collect on a line
  void CreateElectro();
  //---------ElPolyTens.cpp--------------------------------
  /// Alloc the pressure profile 
  void AllocTens();
  /// Calculate the forces for the tension profile
  void CalcTens();
  /// Calculate the densities
  void CalcDens();
  /// Particle positions back folded on the tension reference point (Cartesian coordinates)
  double TensRefCart(double *Pos1,double *Pos2,double *PosP1,double *PosP2);
  ///  Particle positions back folded on the tension reference point (polar coordinates)
  double TensRefPol(double *Pos1,double *Pos2,double *PosP1,double *PosP2);
  /// Data type for distance/field functions
  typedef double(Forces::*TENS_REF)(double *Pos1,double *Pos2,double *PosP1,double *PosP2);
  /// Pointer to a coordinate distance
  TENS_REF Tens_Ref;
  ///  Particle positions back folded on the tension reference point
  double TensRef(double *Pos1,double *Pos2,double *PosP1,double *PosP2){
    return (*this.*Tens_Ref)(Pos1,Pos2,PosP1,PosP2);
  };
  /// Sum the forces on the line joining the points p1 and p2
  void SumTens(int p1,int p2,double Forces,double *DistRel);
  /// Sum the forces on the line joining the points p1 and p2
  void SumTens(int p1,int p2,double *Pre);
  /// Sum the forces on the line joining the points p1 and p2
  void SumTens(double *Pos1,double *Pos2,double *Pre);
  /// Write the pressure and density profile
  void WriteTens(char *TFile,int Comp,double InvNFile);
  /// Write the 2d pressure profile
  void WriteTens2d(FILE *FWrite,int Comp,double InvNFile);
  //---------ForcesLoop.cpp--------------------------------
  /// Find the minimun bilayer thickness for different peptide sizes 
  void ExplorePepSize();
  /// Find the minimun bilayer thickness for different peptide sizes 
  void ExplorePepSize2d();
  /// Find the minimum for different interpeptide distances
  void ExploreDoubleMin();
  /// Total energy of the system
  void CalcTotNrg(char *FName,int nFile);
  /// Exchange energy of the protein
  void CalcNrgPep(char *File2Open,int f);
  /// Run a step further
  void RunDynamics();
  /// Build the widom histograms
  void RunWidom(char *File2Read,int f);
  /// Build the widom histograms
  void RunWidomChIn(char *File2Read,int f);
  /// Build the widom histograms
  void RunWidomChOut(char *File2Read,int f);
  /// Rosenbluth weights for insertion
  void RosenIn(FILE *WidomIn);
  /// Rosenbluth histograms for deletion
  void RosenOut(FILE *WidomIn);
  /// Build the widom histograms with Rosenbluth weight
  void RunWidomBiasChOut(FILE *WidomOut);
  /// Choose the simulation method
  void ChooseSimMode();
  /// Perform a operation every time step
  void Task();
  /// Calculate and sum up the pressure profile from the file list
  void CalcTens(char **argv,int *FilePos,int NFile);
  /// Average of the forces
  void AvForces(char **argv,int *FilePos,int NFile);
  /// Trial loop
  void Trial();
  /// Minmal md
  void MinimalMD();
  /// Minmal nrg
  double MinimalNrg();
  /// Simulation loop for 1d
  void Sim1d(){VelVerlet1();ForceFieldLine();VelVerlet2();};
  /// Simulation loop for 2d
  void Sim2d(){VelVerlet1();Wave();VelVerlet2();};
  /// Simulation loop for 1d
  void Sim3d(){VelVerlet1();ForceFieldBulk();VelVerlet2();};
  /// Simulation loop for leaves
  void SimLeaves(){VelVerlet1();ForceFieldLeaves();AddRigid();VelVerlet2();};
  /// Simulation loop for rigid
  void SimRigid(){VelVerlet1();ForceFieldRigid();VelVerletRigid();VelVerlet2();};
  /// Iterative process to approach to the solution 
  void MinimizeSol();
 //----------------CONSTANTS----------------------------
  /// Step to move a particle
  double IncrDist;
  /// Time step
  double Deltat;
  /// Equilibrium number of particles/chains
  double NChemPotId;
  /// Chemical potential of the particles
  double ChemPotId;
  /// Chemical potential of the particles
  double ChemPotEx;
  /// Standard deviation of the gaussian chain
  double GaussVar;
  /// Cut off of the interactions
  double CutOff;
  /// Spatial separation between particles
  double Dx;
  /// Viscosity of the medium
  double Viscosity;
  /// Total time
  double Time;
  /// Average energy per particle
  double NrgPBead;
  /// Bead to move
  int Bead2Move;
  /// Old part to move
  int Old2Move;
  /// Boh
  int IntMax;
  /// Number of particle per edge
  int nEdge[3];
  /// Shape of system
  int SysShape;
  /// Calculation mode
  int CalcMode;
  /// Thermostat mode
  int ThermMode;
  /// Maximum number of time steps
  int SimLimit;
  /// Which arrays are allocated
  int SysAlloc;
  /// If interpolates with the splines
  int IfInterp;
  /// Boh
  int IfLeaves;
  /// Boh
  int IfMove;
  /// Boh
  int IfNano;
  /// If the matrix has to be changed
  int IfFillMatrix;
  /// Exit from the loop
  int IfExit;
  /// Count accepted moves
  int NInsertion;
  /// Count accepted moves
  int NRemoval;
  /// First and last file of the list
  int NFile[2];
  /// How many timesteps before redrawing
  int NUpdate;
  /// How many timesteps before write the snapshot
  int NWrite;
  /// Total number of points for drawing a spline
  int NSpline;
  /// Prefactor of the forces
  KFORCES Kf;
  /// Array containing the forces for each particle
  FORCES *Fm;
  /// Structure for the pressure calculation
  TENS Tens;
  /// Pair list
  DomDec *Pc;
  /// Splines
  PART *Pl;
#ifdef USE_MPI
  SingProc *Proc;
#endif
#ifdef __glut_h__
  /// Show interpolating lines
  int Interp();
  /// Initialize the scene
  int Graphics(int argc,char **argv);
  /// Additional key bindings
  void keyboard(unsigned char key,int x, int y);
  /// Two dimensional soil
  void DrawSoil();
  /// Two dimensional surface
  void DrawCarpet();
  /// Alternative drawing of the particle position
  void DrawParticles();
  /// Alternative drawing of the particle position
  void DrBondLine(int p);
  /// Draw the scene
  void DrawScene();
  /// Draws cylinder or spheres
  void DrawNano();
  /// Idle function to run the dynamics
  void DynamicsView();
  /// Menu
  void Menu();
  /// List referring the cylinder
  GLuint *Cylinder;
  /// List referring the particle
  GLuint Particles;
  /// Boh
  int NShow;
  /// Produce the images for the video
  int IfMovie;
  /// Current number of frame
  int Frame;
  /// Boh
  int IfExt;
  /// Boh
  int IfSpline;
  /// Boh
  int IfRot;
  /// Visualize spheres or points
  int IfSphere;
  /// Boundary condition for a single particle
  int BeadType;
  /// Menu identifier
  int menu,submenu;
  /// Boh
  int IfLine;
#endif
};
/// Idle function to run the dynamics
extern void DynamicsMotion();
#endif //ELPOLY_H
