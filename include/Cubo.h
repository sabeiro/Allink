/***********************************************************************
Cubo: Header file for the Cubo class.
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
#ifndef CUBO_H
#define CUBO_H
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <signal.h>
#include "Matematica.h"
#define MAX_JOINT_VERTEX 16


using namespace std;
/// All the general functions shared among the different domain decomposition structures 
class DomDecBasics{
 private:
 public:
  /// Boundary condition for the cell consider
  double BoundCond[27][3];
  /// Cell where the current particle sits
  int cCurr;
  /// Neighbouring list
  int NeiListCurr[27];
  /// Number of neighbouring cells
  int NNeiCurr;
  /// Current number of neighbours
  int nNeiCurr;
  /// Reference particle
  int p1Curr;
  /// Current particle
  int p2Curr;
  /// 1 if the loop is over
  int IfLoopCurr;
  /// Module 10 (nx,ny,nz) -> nc
  int Mod10[3];
  /// Number of cells
  int NCell;
  /// Number of cells per direction
  int NSect[3];
  /// Box sizes
  double Edge[3];
  /// Inverse box size
  double InvEdge[3];
  /// Cell CutOff
  double CutOff;
  /// Position of the current particle
  double PosCurr[3];
  /// Number of allocated part
  int NAllocP;
  /// Total number of particle
  int NPart;
  /// Number of neighbours
  int NNeighbour;
  /// Allocate 
  DomDecBasics();
  /// Handle errors
  void SigErr(int Condition,const char * s, ...);
  /// Relative position in the cell
  int PosRet(const double Pos[3]);
  /// # part in the class
  int pNPart();
  /// Print the number of cells
  int pNCell();
  /// (x,y,z)->n
  int pCella(const double Pos[3]);
  /// (x,y,z)->(cx,cy,cz,cTot)
  void pCella(const double Pos[3],int c[4]);
  /// Get the coordination number
  int GetCoorNumb(double *Pos);
  /// Get the number of the neighbouring cells for the particles close to a wall
  int GetCellCoord(int c,int Coord,int *NeiList);
  /// Get the number of the neighbouring cells for the particles close to a wall
  int GetCellCoord(double *Pos,int *NeiList);
  /// Get the number of the neighbouring cells
  int GetCellCh(int c,int *NeiList);
  /// Get the number of the neighbouring cells for a ghost particle
  int GetCell(double *Pos,int *NeiList);
  /// Print the content of a cell
  void PrintCell(const int c);
  /// Print the the particles in the cells
  void PrintCells();
  /// Print the the particle list in the cells
  void PrintList(const int c);
  /// Print the the particle list in the cells
  void PrintLists();
  /// Check the list
  void CheckList();
  /// Check the neighbours
  void CheckNei(int p);
  /// Set the CutOff
  int SetCutOff(double CutOffExt);
};
/// The previuos and consecutive particle in the cell linked list
class DomPart {
 private:
  int *Prev;
  int *Next;
  int NPart;
 public:
  DomPart(int NPart);
  /// Returns a entry
  int operator[](int col);
  DomPart& operator++();
};
/// The previuos and consecutive particle in the cell linked list
typedef struct {
  /// Position
  double Pos[3];
  /// Next particle in the list
  int Next;
  /// Previous particle in the list
  int Prev;
  /// Coordination number
  int Coord;
  /// Which cell it's belonging
  int Cell;
  /** Position of the particle, saves time in accessing particle position close in memory */
  //double Pos[3];
} DOMAIN_PART;
/** For every cube how many particles and which are the first and the last
 */
class DomCell{
  /// Increment the current position
  friend void operator+(DomCell &Dc);
 public:
  DomCell(){};
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
  DomCell& operator=(const DomCell &Dc);
  DomCell& operator+=(const int i);
  DomCell& operator++();
};
/// Domain decomposition as pointer to linked particles
class DdLinkedList : public DomDecBasics{
private:
public:
  /// Allocate 
  DdLinkedList(double Edge[3],int NPart,double CutOff);
  /// Erase the pairlist
  void Erase();
  /// Clear the pairlist
  void Clear(){Erase();};
  /// Swap two particles
  int SwapPart(int p1,double *Pos1,int p2,double *Pos2);
  /// Add a particle to the cell c
  void AddPart(const int p,double *Pos);
  /// Remove a particle form the cell c
  void RemPart(const int p,double *Pos);
  /// Remove the particle p
  void RemPart(const int p);
  /// Add a particle to the cell c
  void AddPart(const int p,const int c);
  /// Remove a particle form the cell c
  void RemPart(const int p,const int c);
  /// Move a particle form the cell c1 to the cell c2
  void MovePart(const int p,double *OldPos,double *NewPos);
  /// Move a particle form the cell c1 to the cell c2
  void MovePart(const int p,double *NewPos);
  /// Relative position in the cell
  int PosRet(const double Pos[3]);
  /// # part in the cell
  int pNPart(const int c){return Cella[c].NPart;};
  /// # part in the class
  int pNPart(){return NPart;};
  /// Print the cell which the particle belong
  int pCell(const int p){return Pc[p].Cell;};
  /// First part in the cell
  int First(const int c){return Cella[c].First;};
  /// Next linked part 
  int Next(const int p){return Pc[p].Next;};
  /// Iterate in the cell
  int ItCell(const int c);
  /// Stop the loop and set the counter to zero
  int IfItCell(const int c);
  /// Stop the loop and set the counter to zero
  int IfItCouple(const int c);
  /// Set the coordination umber to the part p
  int SetCoorNumb(double *Pos,int p);
  /// Choose among the different neighbouring lists
  int GetNei(double *Pos,int *NeiList){
    return GetCell(Pos,NeiList);
    // return GetCellCoord(Pos,NeiList);
  }
  /// Set the counters to the initial position
  void SetCounters(int c);
  /// Iterates along all couples
  void Couple(const int c,int *p1,int *p2);
  /// Print the particles in a cell
  void PrintCell(const int c);
  /// Print the the particles in the cells
  void PrintCells();
  /// Gather information of the neighbouring cells
  void SetCurr(int p);
  /// Gather information of the neighbouring cells
  void SetCurrGhost(double *Pos);
  /// Increase the iterator to the next couple
  void NextCurr();
  /// Increase the iterator to the next couple
  void NextCurrGhost();
  /// Tell when the curr loop is over
  int IfCurr();
  /// Tell when the curr loop is over
  int IfCurrGhost();
  /// Retrun the squared current interparticle distance
  void Dist2Curr(double *DistRel);
  /// Retrun the squared current interparticle distance
  void Dist2CurrGhost(double *DistRel);
  /// Print the the particle list in the cells
  void PrintList(const int c);
  /// Print the the particle list in the cells
  void PrintLists();
  /// Check the list
  void CheckList();
  /// Check the neighbours
  void CheckNei(int p);
  /// Increment the current part in the cell
  void IncrCurr(const int c);
  /// Increment the current iterators in the cell
  void IncrCurrList(const int c);
  /// Find the closest neighbour
  int FindClosest(int p1);
  /// Number of particles and iterators per cell
  DomCell *Cella;
  /// List of position of the particles
  DOMAIN_PART *Pc;
};
/// The array of particles in every cell
typedef struct{
  /// Id of the particle
  list <int> Part;
}DdCell;
/// Domain decomposition using arrays
class DdArray : DomDecBasics{
private:
  /// Edges of the system
  double Edge[3];
  /// Number of cells per dimension
  int NSect[3];
  /// Modulus per dimension
  int Mod10[3];
  /// Total number of part
  int NPart;
  /// Total number of cells
  int NCell;
  /// CutOff distance
  double CutOff;
  /// First iterator
  list <int>::iterator NCurr;
  /// Second iterator
  list <int>::iterator NCurr2;
public:
  /// Allocate 
  DdArray(double EdgeExt[3],int NPart,double CutOff);
  /// Cell structure
  DdCell *Cella;
  /// Erase the content of the cells
  void Erase();
  /// Clear the pairlist
  void Clear(){Erase();};
  /// Add a part to the cell
  void AddPart(int p,double *Pos);
  /// Remove a part from the cell
  void RemPart(int p,double *Pos);
  /// Move a part 
  void MovePart(int p,double *OldPos,double *NewPos);
  /// Swap two parts
  void SwapPart(int p1,double *Pos1,int p2,double *Pos2);
  /// Set the counters to the initial pointer
  void SetCounters(int c);
  /// End the iteration in the cell
  int IfItCell(int c);
  /// Increase the counters
  void IncrCurr(int c);
  /// Value of the current particle in the cell c
  int ItCell(int c);
  /// Choose among the different neighbouring lists
  int GetNei(double *Pos,int *NeiList){
    return GetCell(Pos,NeiList);
  }
  /// Returns the pointed couple of particle
  void Couple(const int c,int *p1,int *p2);
  /// Tell when the loop is over
  int IfItCouple(const int c);
  /// Increment the pointers
  void IncrCurrList(const int c);
  /// Print the number of cells
  int pNCell(){return NCell;};
  /// Print the content of one cell
  void PrintCell(const int c);
  /// Print the content of every cell
  void PrintCells();
};
/// The array of particles in every cell
typedef struct{
  /// How many particles
  int NPart;
  /// Id of the particle
  int *Part;
}DdFixCell;
/// Domain decomposition using arrays
class DdFixedSize : DomDecBasics{
private:
  /// Maximum number of particles per cell
  int NPCell;
  /// First iterator
  int NCurr;
  /// Second iterator
  int NCurr2;
public:
  /// Allocate 
  DdFixedSize(double EdgeExt[3],int NPart,double CutOff);
  /// Cell structure
  DdFixCell *Cella;
  /// Erase the content of the cells
  void Erase();
  /// Clear the pairlist
  void Clear(){Erase();};
  /// Add a part to the cell
  void AddPart(int p,double *Pos);
  /// Remove a part from the cell
  void RemPart(int p,double *Pos);
  /// Add a part to the cell
  void AddPart(int p,int c);
  /// Remove a part from the cell
  void RemPart(int p1,int c);
  /// Move a part 
  void MovePart(int p,double *OldPos,double *NewPos);
  /// Swap two parts
  void SwapPart(int p1,double *Pos1,int p2,double *Pos2);
  /// Set the counters to the initial pointer
  void SetCounters(int c);
  /// End the iteration in the cell
  int IfItCell(int c);
  /// Increase the counters
  void IncrCurr(int c);
  /// Value of the current particle in the cell c
  int ItCell(int c);
  /// Choose among the different neighbouring lists
  int GetNei(double *Pos,int *NeiList){
    return GetCell(Pos,NeiList);
  }
  /// Returns the pointed couple of particle
  void Couple(const int c,int *p1,int *p2);
  /// Tell when the loop is over
  int IfItCouple(const int c);
  /// Increment the pointers
  void IncrCurrList(const int c);
  /// Print the content of one cell
  void PrintCell(const int c);
  /// Print the content of every cell
  void PrintCells();
};
/// Stupid double loop to check the other pair loops
class DdDoubleLoop : public DomDecBasics{
private:
public:
  /// Allocate 
  DdDoubleLoop(double Edge[3],int NPart,double CutOff);
  /// Erase the pairlist
  void Erase();
  /// Clear the pairlist
  void Clear(){Erase();};
  /// Swap two particles
  int SwapPart(int p1,double *Pos1,int p2,double *Pos2);
  /// Add a particle to the cell c
  void AddPart(const int p,double *Pos);
  /// Remove a particle form the cell c
  void RemPart(const int p,double *Pos);
  /// Remove the particle p
  void RemPart(const int p);
  /// Add a particle to the cell c
  void AddPart(const int p,const int c);
  /// Remove a particle form the cell c
  void RemPart(const int p,const int c);
  /// Move a particle form the cell c1 to the cell c2
  void MovePart(const int p,double *OldPos,double *NewPos);
  /// Move a particle form the cell c1 to the cell c2
  void MovePart(const int p,double *NewPos);
  /// Relative position in the cell
  int PosRet(const double Pos[3]);
  /// # part in the cell
  int pNPart(const int c){return Cella[c].NPart;};
  /// # part in the class
  int pNPart(){return NPart;};
  /// Print the cell which the particle belong
  int pCell(const int p){return Pc[p].Cell;};
  /// First part in the cell
  int First(const int c){return Cella[c].First;};
  /// Next linked part 
  int Next(const int p){return Pc[p].Next;};
  /// Iterate in the cell
  int ItCell(const int c);
  /// Stop the loop and set the counter to zero
  int IfItCell(const int c);
  /// Stop the loop and set the counter to zero
  int SetCoorNumb(double *Pos,int p);
  /// Choose among the different neighbouring lists
  int GetNei(double *Pos,int *NeiList){
    return GetCell(Pos,NeiList);
    // return GetCellCoord(Pos,NeiList);
  }
  /// Set the counters to the initial position
  void SetCounters(int c);
  /// Print the particles in a cell
  void PrintCell(const int c);
  /// Print the the particles in the cells
  void PrintCells();
  /// Gather information of the neighbouring cells
  void SetCurr(int p);
  /// Increase the iterator to the next couple
  void NextCurr();
  /// Tell when the curr loop is over
  int IfCurr();
  /// Retrun the squared current interparticle distance
  void Dist2Curr(double *DistRel);
  /// Print the the particle list in the cells
  void PrintList(const int c);
  /// Print the the particle list in the cells
  void PrintLists();
  /// Check the list
  void CheckList();
  /// Check the neighbours
  void CheckNei(int p);
  /// Increment the current part in the cell
  void IncrCurr(const int c);
  /// Increment the current iterators in the cell
  void IncrCurrList(const int c);
  /// Number of particles and iterators per cell
  DomCell *Cella;
  /// List of position of the particles
  DOMAIN_PART *Pc;
};
/* typedef struct { */
/*   /// Position */
/*   double Pos[3]; */
/* } PPOS; */
/* typedef struct { */
/*   /// Vertex id */
/*   list <int> v; */
/* } VERTEX; */
typedef struct {
  /// Vertex postion
  double Pos[3];
  /// Vertex id
  int v;
  /// Triangles/lines connected (fixed size!!! should be enough)
  int t[MAX_JOINT_VERTEX];
  /// Number of triangles connected to the vertex
  int NTria;
} VERTEX;
/// Connects the triangles by vertices
class NeiVertex{
 private:
  /// Triangle iterator
  //list <int>::iterator tCurr;
  int tCurr;
  /// Vertex iterator
  int vCurr;
  /// Current memory position of the vertices
  /// List of vertices
  VERTEX *Vertex;
  /* /// List of positions */
  /* list <PPOS> vPos; */
  /// Memory postion of the vertices
  int *vMemPos;
  /// Edge sizes
  double Edge[3];
  /// Number of grid segments
  int NGrid;
  /// Total number of vertex
  int NVert;
  /// Total number of triangles
  int NTria;
  /// Number of vertices per traingle (polygon)
  int NvPt;
 public:
  /// Alloc NVertices
  NeiVertex(int NTriaExt,int NvPtExt,int NGridExt,double *EdgeExt);
  /// Delete the class
  ~NeiVertex();
  /// Add the triangle t at the vertex v
  void Add(int v,int t,double *Pos);
  /// Add the triangle t in the Pos
  void Add(double *Pos,int t);
  /// Remove the triangle t from the vertex v
  void Rem(int v,int t);
  /// Copy two vertices
  void CopyVert2To1(VERTEX Vert1,VERTEX Vert2);
  /// Swap to vertices
  void Swap(int v1,int v2);
  /// Reorder and fill the vertices 
  void Reorder();
  /// Set counters to zero for the point v
  void SetCounters(int v);
  /// Set all the counters to zero
  void SetCounters();
  /// Correspondent vertes for Pos
  int GetVertex(double *Pos);
  /// End of the counters
  int IfItCell(int v);
  /// Increment the counter
  void IncrCurr(int v);
  /// Current vertex iterator for the vertex v
  int VertCurr(int v);
  /// Current triangle for the vertex v
  int TriaCurr(int v);
  /// Position of the vertex v
  void PosVertex(int v,double *Pos);
  /// Print the entire structure
  void Print();
};

#endif //CUBO_H
