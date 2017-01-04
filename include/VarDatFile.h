#ifndef VAR_DAT_FILE_H
#define VAR_DAT_FILE_H
#include "../include/Matematica.h"
#define DIS_IF_TYPE(c,t) ( ((c)&(t))==(t) )
#define DIS_ADD_TYPE(c,t) ( (c)|=(t) )
#define DIS_REM_TYPE(c,t) ( (c)&=(~t) )
/// Visualize log scale
enum VisLog{
  /// No Log
  DIS_NOLOG = 0x001,
  /// Log on x
  DIS_LOGX  = 0x002,
  /// Log on y
  DIS_LOGY  = 0x004,
  /// Log on x and y
  DIS_LOGLOG= 0x006,
};

enum RefType{
  /// sequence of natural numbers
  REF_SEQ  = -1,
  /// is an ascissa his own
  REF_ASC  = -2,
};
/// Reads and stores a data file to be elaborated via Matematica
class VarDatFile:public Matematica{
 private:
  /// Command to execute on x
  char XFormula[60];
  /// Command to execute on y
  char YFormula[60];
  /// Header's string
  char Header[STRSIZE];
  /// Defines some initial values
  void Init();
  /// Allocates
  void Allocate();
  /// Max number of points per array
  int NMax;
  /// Max number of points per array
  int NMaxPunti;
  /// Old Number 
  int NVar;
  /// List of NVar arrays of NMax points
  double **st;
  /// Minima
  double *sMin;
  /// Maxima
  double *sMax;
  /// Names of the arrays
  char **sc;
  /// Number of points per signal
  int *SetNMax;
  /// Line type
  int *sType;
  /// Line color
  int *sColor;
  /// Abscissa column
  int *RefAbsc;
  /// Stored if required for a transformation of the signal
  double *Punti;
  /// Errors of the @see Punti set
  double *PuntiErr;
  /// Stored if required for a additional transformation of the elaborated signal
  double *Punti1;
  /// Values of the empirical distribution
  double *dInter;
  /// Values of the empirical distribution
  double *dError;
  /// Values of the theoretical distribution
  double *dInter1;
  /// Global max
  double GlobMax;
  /// Global Min
  double GlobMin;
  /// x Global max
  double xGlobMax;
  /// x Global Min
  double xGlobMin;
  /// If Punti should be allocated
  int IfPuntiAlloc;
  /// If Punti should be freed
  int IfPuntiFree;
  /// Moltiplication factor for a scaled set (MediaMob)
  //int Incr;
 public:
  /// Allocate the vectors reading from a file
  VarDatFile(char *file,int NBin);
  /// Allocate the vector with a given size
  VarDatFile(int ExtNMax,int ExtNVar,int ExtNBin);
  /// Point to the already allocated vector
  VarDatFile(double **ExtSt,int ExtNMax,int ExtNVar,int ExtNBin);
  /// Load list
  VarDatFile(char **FileList,int *Pos,int NFile,int Values);
  /// Freeing
  ~VarDatFile();
  /// Add the content of the file into the memory
  int Aggiungi(char *file);
  /// Reads information contained in the file
  int FileInfo(FILE *InFIle,int *NCol);
  /// Reads a single line
  int ReadLine(char *cLine,double *Value);
  /// Is it working?
  int  ReadLine(FILE *InFile, const char *Format, ...);
  /// Do a generic elaboration
  int ElabSegnale(int NElMin,int NElMax);
  /// Distribution of a signal
  MOMENTI DistrSegnale(int NElMin,int NElMax,int IfNorm);
  /// Logarithmic distribution of a signal
  MOMENTI DistrLogSegnale(int NElMin,int NElMax,int IfNorm);
  /// Exponential distribution of a signal
  MOMENTI DistrExpSegnale(int NElMin,int NElMax,int IfNorm);
  /// Distribution of a signal
  MOMENTI DistrSegnale(int NElMin,int NElMax,double *Border,int IfNorm);
  /// Distribution of a signal with error
  MOMENTI DistrSignErr(int NElMin,int NElMax,double *Border,int IfNorm);
  /// Tests a gaussian distribution of the signal
  MOMENTI DistrGaussSegnale(int ElMin,int ElMax,int IfNorm);
  /// Weighted distribution
  MOMENTI WeightAverageSet(int CoordY,int ElMin, int ElMax);
  /// Distributions of data organized in different samples
  void DistrSignSample(int ElMin, int ElMax,double **Distr,int NSample,int IfNorm,double *xBound);
  /// Calculate the spectrum
  int SpettroSegnale(int NElMin,int NElMax);
  /// Spectrum of a sequence of lines
  void SpeLine(int ElMin,int ElMax,int NBin,double *Spe);
  /// Normalize the signal
  int NormalizzaSegnale(int NElMin,int NElMax);
  /// Normalize the intervals 
  int NormalizzaInter();
  /// Transforms the signal into its square root
  bool RadiceSegnale();
  /// Autocorrelation of the signal
  int AutocorSegnale(int NElMin,int NElMax);
  /// Sum the values between ElMin and ElMax
  double SumSegnale(int CoordY,int NElMin,int NElMax);
  /// Numerical integration of the signal
  double IntSegnale();
  /// Sums all the yCoordinates
  void SommaSegnali();
  /// Average of all ordinate
  void AverageOrdinate(int ElMin, int ElMax,double *Distr,double *xBound);
  /// Points to a different column
  void Punta(int n);
  /// Points to an external pointer, probably not working
  void Punta(double **ExtSt,int n);
  /// Point an external pointer to the private memory
  void Punta(double *sp,int n);
  /// Point the internal pointer to an external memory
  void PuntaInt(double *sp);
  /// Points the function Elab to the derivative, probably not working
  void DerivataSegnale(){Elab = &Matematica::DerO4;};
  /// Points the function Elab to the square gradient, probably not working
  void VarieSegnale(){Elab = &Matematica::SquareGradient;};
  /// Points the function Elab to the module, probably not working
  void ModuloSegnale(){Elab = &Matematica::Modulo;};
  /// A generic function to be edit
  double VarieSegnale(int NElMin,int NElMax);
  /// Computes a liner fit
  RETTA InterRettSegnale(int CoordY,int NElMin,int NElMax,int LogLog);
  /// Computes a exponential fit
  RETTA InterExpSegnale(int CoordY,int NElMin,int NElMax,int LogLog);
  /// Computes a Gaussian fit
  MOMENTI InterGaussSegnale(int CoordY,int NElMin,int NElMax,int LogLog);
  /// Computes a parabolic fit
  PARABOLA ParabolaSegnale(int CoordY,int NElMin,int NElMax,int LogLog);
  /// Running average
  int MediaMobSegnale(int);
  /// Weighted histogram analysis
  int WeightHistoSign(int NHisto);
  /// To point correlation at a distance dist
  int CorrelaADuePunti(int dist);
  /// Selfsimilarity
  void AutosimilaritaSegnale(int);
  /// Reallocate the pointer dInter dInter1
  void CambiaNBin(int);
  /** Copies Punti in Punti1 to store the information of
      the already elaborated array, to be fixed */
  void CambiaPunti();
  /// Prints the content in memory
  void Print();
  /// Sort with respect to a coordinate
  void Sort();
  /// Writes the content of the elaborated array
  void ScriviPunti(char *file);
  /// Writes the content of the pointed vector
  void ScriviFile(char *file,int CoordY,int LogLog,int NVisMin,int NVisMax);
  /// Writes all the content in memory
  void ScriviTutto(char *file,int LogLog,int NVisMin,int NVisMax);
  /// Rescale to bulk
  void RescaleToBulk(char *FName);
  /// Export in the txvl file format
  void ExportTxvl(char *FName,int NElMin,int NElMax);
  /// Export a contour plot
  void TecPlot(char *FName);
  /// Allocates @see Punti if necessary
  void PuntiAlloc();
  /// Free Punti if it is allocated
  void PuntiFree();
  /// Number of point of the distribution
  int NBin;
  /// Print the number of variables
  int pNVar(){ return NVar;};
  /// Print the maximum number of data in the array
  int pNMax(){return NMax;};
  /// Print the maximum number of data for the column
  int pNRow(int CoordY);
  /// Print global borders
  void pGlobBorder(double *xMin,double *xMax,double *yMin,double *yMax){
    *xMin = xGlobMin;
    *xMax = xGlobMax;
    *yMin = GlobMin;
    *yMax = GlobMax;
  };
  /// Value at the position \param n in the column \param v
  double Val(int CoordY,int n);
  /// Value of the Abscissa at position \param n
  double Abscissa(int CoordY,int n);
  /// Value of the elaborated array at position \param n
  double pPunti(int n);
  /// Value of the error array at postion \param n
  double pPuntiErr(int n);
  /// Find the borders globally
  void pMinMaxGlob(int NVisMin,int NVisMax);
  /// Print the maximum of all sets
  double pMaxGlob(int NVisMin,int NVisMax);//{return GlobMax;};
  /// Print the minimum of all sets
  double pMinGlob(int NVisMin,int NVisMax);//{return GlobMin;};
  /// Print the abscissa minimum of all sets
  double pxMinGlob(int NVisMin,int NVisMax);//{return GlobMin;};
  /// Print the abscissa maximum of all sets
  double pxMaxGlob(int NVisMin,int NVisMax);//{return GlobMax;};
  /// Print the maximum of the column \param n
  double pMax(int CoordY,int NVisMin,int NVisMax);//{return sMax[n];};
  /// Print the minimum of the column \param n
  double pMin(int CoordY,int NVisMin,int NVisMax);//{return sMin[n];};
  /// Print the maximum of all sets
  double pMaxGlobLog(int NVisMin,int NVisMax);//{return GlobMax;};
  /// Print the minimum of all sets
  double pMinGlobLog(int NVisMin,int NVisMax);//{return GlobMin;};
  /// Print the maximum of the column \param n
  double pMaxLog(int CoordY,int NVisMin,int NVisMax);//{return sMax[n];};
  /// Print the minimum of the column \param n
  double pMinLog(int CoordY,int NVisMin,int NVisMax);//{return sMin[n];};
  /// Print the postion \param n of array \param dInter
  double pInter(int n){return dInter[n];};
  /// Print the postion \param n of array \param dError
  double pError(int n){return dError[n];};
  /// Print the postion \param n of array \param dInter1
  double pInter1(int n){return dInter1[n];};
  /// Print the minimum of @see Punti
  double PuntiMin();
  /// Print the maximum of @see Punti
  double PuntiMax();
  /// If the column is an ascissa
  int IsAbscissa(int Col){
    if(RefAbsc[Col] == Col) return 1; return 0;};
  /// If the column is an ascissa
  int pRefAbsc(int Col){return RefAbsc[Col];};
  /// If the column is an ascissa
  int pSetNMax(int Col){return SetNMax[Col];};
  /// If the column is an ascissa
  int IsSequence(int Col){return RefAbsc[Col]==REF_SEQ?1:0;};
  /// Set the natural sequence as abscissa for the column Col
  void ImpSequence(int Col);
  /// Set the natural sequence as abscissa for all sets
  void ImpSequence();
  /// Set the columns for the x array
  void ImpCoordX(int vAbs);
  /// Set the column for the x array
  void ImpCoordX(int vSet,int vAbs);
  /// Set the column for the y array
  void ImpCoordY(int Ext);
  /// Set the value of XFormula
  void setXFormula(char *str);
  /// Set the value of YFormula
  void setYFormula(char *str);
  /// Reverse the sets
  void Reverse();
  /// Smooth the line
  int Smooth(double Fact,int CoordY,int NVisMin,int NVisMax);
  /// Smooth the line
  void SmoothGauss(double Fact,int CoordY,int NVisMin,int NVisMax);
  /// Distribution of distances 
  void DoubleDistFluct();
  /// Execute a formula defined in XFormula and YFormula
  void WriteFormula(char *Exit);
  /// Print the header
  char *PrintHeader();
  /// Points to a column of st
  double *sp;
};
#endif //VARELEMENTIGRAFICI_H
