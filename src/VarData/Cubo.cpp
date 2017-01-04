#include <Cubo.h>
#define MASK_X 0
#define MASK_Y 4
#define MASK_Z 8
#define GET_MASK_X 0x00f
#define GET_MASK_Y 0x0f0
#define GET_MASK_Z 0xf00
#define DOM_CENTER 0
#define DOM_LEFT   1
#define DOM_RIGHT  2
#define DOM_SET_POS(n,p,m)  ( (n)|=((p)<<(m)) )
#define DOM_GET_POS(n,m,g)    ( ((n)&(g))>>(m) )
//FIXME
#define DOM_IF_POS(n,m,d)   ( ((n)>>(m))==(d) )
//-----------------BASICS-----------------------------
/** Allocator of the general structure */
DomDecBasics::DomDecBasics(){
  //  signal(SIGSTOP,StopProg);
}
// void DomDecBasics::StopProg(int sig){
//   printf("Program stopped\n");
// }
/** Signal an error and exit */
void DomDecBasics::SigErr(int IsWrong,const char * s, ...){
#ifdef DEBUG
  if(IsWrong){
    va_list args;
    va_start(args, s);
    fprintf(stderr, "Dom Dec error: ");
    vfprintf(stderr, s, args);
    fprintf(stderr, "exiting the program\n");
    va_end(args);
    abort();
  }
#else  
  return;
#endif
}
/** Set the cut off of the grid spacing, the cut off should be much smaller than the box size */
int DomDecBasics::SetCutOff(double CutOffExt){
  for(int d=0;d<3;d++){
    if(CutOffExt > Edge[d]/(double)NSect[d]){
      printf("CutOff larger than the cell size\n");
      return 1;
    }
  }
  CutOff = CutOffExt;
  return 0;
}
/** Number of particles */
int DomDecBasics::pNPart(){
  return NPart;
}
/** Number of cells */
int DomDecBasics::pNCell(){
  return NCell;
}
/** Return the list of neighbouring cells*/
int DomDecBasics::GetCell(double *Pos,int *NeiList){
  int c = pCella(Pos);
  return GetCellCh(c,NeiList);
}
/** Neighbouring cells with periodic boundary conditions */
int DomDecBasics::GetCellCh(int c,int *NeiList){
  SigErr(c < 0 || c > NCell,"cell value %d over the ranges\n",c);
  if(NCell == 1){
    return 1;
  }
  int NeiMod[3];
  int NeiCell[3];
  int NNei = 0;
  NeiMod[2] = (int)(c/(double)Mod10[2]);
  NeiMod[1] = (int)((c-NeiMod[2]*Mod10[2])/(double)Mod10[1]);
  NeiMod[0] = (int)(c-NeiMod[1]*Mod10[1]-NeiMod[2]*Mod10[2]);
  for(int n=0;n<27;n++) for(int d=0;d<3;d++) BoundCond[n][d] = 0.;
  for(int cx=-1;cx<=1;cx++){
    NeiCell[0] = NeiMod[0] + cx;
    if(NeiCell[0] >= NSect[0]){ 
      NeiCell[0] -= NSect[0];
      BoundCond[NNei][0] = -1.;
    }
    if(NeiCell[0] < 0){
      NeiCell[0] += NSect[0];
      BoundCond[NNei][0] = 1.;
    }
    for(int cy=-1;cy<=1;cy++){
      NeiCell[1] = NeiMod[1] + cy;
      if(NeiCell[1] >= NSect[1]){
	NeiCell[1] -= NSect[1];
	BoundCond[NNei][1] = -1.;
      }
      if(NeiCell[1] < 0){
	NeiCell[1] += NSect[1];
	BoundCond[NNei][1] = 1.;
      }
     for(int cz=-1;cz<=1;cz++){
	NeiCell[2] = NeiMod[2] + cz;
	if(NeiCell[2] >= NSect[2]){
	  NeiCell[2] -= NSect[2];
	  BoundCond[NNei][2] = -1.;
	}
	if(NeiCell[2] < 0){
	  NeiCell[2] += NSect[2];
	  BoundCond[NNei][2] = 1.;
	}
	NeiList[NNei] = NeiCell[0]*Mod10[0] + NeiCell[1]*Mod10[1] + NeiCell[2]*Mod10[2];
	NNei++;
      }
    }
  }
  //printf("\n %d) ",c);for(int i=0;i<NNei;i++) printf("%d ",NeiList[i]);
  return NNei;
}
/** Coordination number of the particle in the cell */
int DomDecBasics::GetCellCoord(double *Pos,int *NeiList){
  int c = pCella(Pos);
  int Coord = GetCoorNumb(Pos);
  return GetCellCoord(c,Coord,NeiList);
}
/** Coordination number of the particle in the cell, every particle has a flag which tells to which cell border is close to. Saves computational time. */
int DomDecBasics::GetCellCoord(int c,int Coord,int *NeiList){
  SigErr(c < 0,"Probably the particle position is not a number cell %d\n",c);
  int NNei = 0;
  int nNei[3];
  int NeiCell[3];
  int Mask[3]    = {MASK_X,MASK_Y,MASK_Z};
  int GetMask[3] = {GET_MASK_X,GET_MASK_Y,GET_MASK_Z};
  int Move[3] = {0,0,0};
  NeiCell[2] = (int)(c/(double)Mod10[2]);
  NeiCell[1] = (int)((c-NeiCell[2]*Mod10[2])/(double)Mod10[1]);
  NeiCell[0] = (int)(c-NeiCell[1]*Mod10[1]-NeiCell[2]*Mod10[2]);
  for(int n=0;n<27;n++)for(int d=0;d<3;d++) BoundCond[n][d] = 0.;
  for(int d=0;d<3;d++){
    int d1 = (d+1)%3;
    int d2 = (d+2)%3;
    if( DOM_GET_POS(Coord,Mask[d],GetMask[d]) == DOM_RIGHT ){
      Move[d] = DOM_RIGHT;
      int nc = NeiCell[d] + 1.;
      if(nc >= NSect[d]){
	nc = 0;
	BoundCond[NNei][d] = -1.;
      }
      NeiList[NNei++] = nc*Mod10[d] + NeiCell[d1]*Mod10[d1] + NeiCell[d2]*Mod10[d2];
    }
    else if( DOM_GET_POS(Coord,Mask[d],GetMask[d]) == DOM_LEFT ){
      Move[d] = DOM_LEFT;
      int nc = NeiCell[d] - 1.;
      if(nc < 0){
	nc = NSect[d]-1;
	BoundCond[NNei][d] = 1.;
      }
      NeiList[NNei++] = nc*Mod10[d] + NeiCell[d1]*Mod10[d1] + NeiCell[d2]*Mod10[d2];
    }
  }
  if(NNei > 1){
    for(int d=0;d<3;d++){
      nNei[d] = NeiCell[d];
      if(Move[d] == DOM_LEFT){
  	nNei[d] = NeiCell[d] - 1;
  	if(nNei[d] < 0) nNei[d] = NSect[d] - 1;
      }
      if(Move[d] == DOM_RIGHT){
  	nNei[d] = NeiCell[d] + 1;
  	if(nNei[d] >= NSect[d]) nNei[d] = 0;
      }
    }
    NeiList[NNei++] = nNei[0]*Mod10[0] + nNei[1]*Mod10[1] + nNei[2]*Mod10[2];
  }
  if(NNei == 4){
    for(int dd=0;dd<3;dd++){
      for(int d=0;d<3;d++){
  	nNei[d] = NeiCell[d];
  	if(d == dd) continue;
  	if(Move[d] == DOM_LEFT){
  	  nNei[d] = NeiCell[d] - 1;
  	  if(nNei[d] < 0) nNei[d] = NSect[d] - 1;
  	}
  	if(Move[d] == DOM_RIGHT){
  	  nNei[d] = NeiCell[d] + 1;
  	  if(nNei[d] >= NSect[d]) nNei[d] = 0;
  	}
      }
      NeiList[NNei++] = nNei[0]*Mod10[0] + nNei[1]*Mod10[1] + nNei[2]*Mod10[2];
    }
  }
  NeiList[NNei++] = c;
  // for(int n=0;n<NNei;n++) assert(NeiList[n] < NCell);
  // for(int n=0;n<NNei;n++) assert(NeiList[n] >= 0);
  return NNei;
}
/** Return the unique cell identification number for the given position */
int DomDecBasics::pCella(const double Pos[3]){
  int ris[3];
  for(int d=0;d<3;d++){
    double PosBox = Pos[d] - floor(Pos[d]*InvEdge[d])*Edge[d];
    ris[d] = (int)(PosBox/Edge[d]*NSect[d]);
  }
  int risp = ris[0]*Mod10[0]+ris[1]*Mod10[1]+ris[2]*Mod10[2];
  SigErr(risp >= NCell,"Probably the particle position is not a number\n");
  return risp;
}
/** Return the unique cell identification number for the given position */
void DomDecBasics::pCella(const double Pos[3],int c[4]){
  for(int d=0;d<3;d++){
    double PosBox = Pos[d];// - floor(Pos[d]*InvEdge[d])*Edge[d];
    c[d] = (int)(PosBox/Edge[d]*NSect[d]);
  }
  c[4] = c[0]*Mod10[0]+c[1]*Mod10[1]+c[2]*Mod10[2];
  SigErr(c[4] >= NCell,"Probably the particle position is not a number\n");
}
/** Retrun the coordination number for the given position*/
int DomDecBasics::GetCoorNumb(double *Pos){
  double CellLen[3];
  int Border[3];
  double BorderR[3];
  int Mask[3] = {MASK_X,MASK_Y,MASK_Z};
  int Coord = 0;
  for(int d=0;d<3;d++){
    Border[d]  = (int)(Pos[d]/Edge[d]*NSect[d]);
    CellLen[d] = (Edge[d]/(double)NSect[d]);
    if(Pos[d] - Border[d]*CellLen[d] < CutOff){
      DOM_SET_POS(Coord,DOM_LEFT,Mask[d]);
    }
    if(-Pos[d] + (Border[d]+1)*CellLen[d] < CutOff){
      DOM_SET_POS(Coord,DOM_RIGHT,Mask[d]);
    }
  }
  return Coord;
}
//-------------------LINKED-LIST-------------------------------
/** Constructor for the domain decomposition with the linked list*/
DdLinkedList::DdLinkedList(double EdgeExt[3],int NPartExt,double CutOffExt){
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
    InvEdge[d] = 1./Edge[d];
    //NSect[d] = (int)(Edge[d]/(2.*CutOffExt));
    NSect[d] = (int)(Edge[d]/CutOffExt);
    SigErr(CutOffExt > 2.*Edge[d],"Insufficient number of cells \n");
  }
  NPart = 0;
  SetCutOff(CutOffExt);
  Mod10[0] = 1;
  Mod10[1] = NSect[0];
  Mod10[2] = NSect[1]*Mod10[1];
  NCell = NSect[0]*NSect[1]*NSect[2];
  Cella = new DomCell[NCell];
  SigErr(Cella==NULL,"DomCell not allocated\n");
  for(int c=0;c<NCell;c++){
    Cella[c].NPart = 0;
    Cella[c].First = Cella[c].Last = -2;
  }
  NAllocP = NPartExt + 2000;
  Pc = new DOMAIN_PART[NAllocP];
  SigErr(Pc == NULL,"DOMAIN_PART not allocated\n");
  printf("%d\n",NCell);
}
/** Empty the records of the cells*/
void DdLinkedList::Erase(){
  for(int c=0;c<NCell;c++){
    Cella[c].NPart = 0;
    Cella[c].Last = -1;
  }
  for(int p=0;p<NPart;p++){
    Pc[p].Cell = -1;
    Pc[p].Next = -1;
    Pc[p].Prev = -2;
  }
  NPart = 0;
}
/** Coordination number of the particle in the cell*/
int DdLinkedList::SetCoorNumb(double *Pos,int p){
  Pc[p].Coord = GetCoorNumb(Pos);
  return Pc[p].Coord;
}
/** Add a part to the correspondent cell*/
void DdLinkedList::AddPart(const int p,double *Pos){
  // for(int d=0;d<3;d++){
  //   SigErr(Pos[d] < 0. || Pos[d] > Edge[d],"particle %d over the boudaries 0< %lf < %lf\n",p,Pos[d],Edge[d]);}
  int c = pCella(Pos);
  SetCoorNumb(Pos,p);
  for(int d=0;d<3;d++) Pc[p].Pos[d] = Pos[d];
  AddPart(p,c);
}
/** Add a part to the correspondent cell*/
void DdLinkedList::AddPart(const int p,const int c){
  //PrintCell(c);
  int pl = Cella[c].Last;
  SigErr(pl == p,"The particle %d is already the last in the cell %d\n",p,c);
  SigErr(p >= NAllocP,"More particle than allocated %d > %d\n",p,NAllocP);
  SigErr(Cella[c].NPart >= NAllocP,"More particle than allocated\n");
  Cella[c].NPart++;
  NPart++;
  Cella[c].Last = p;
  Pc[p].Cell = c;
  Pc[p].Next = -1;
  if(Cella[c].NPart == 1){
    Cella[c].First = p;
    Pc[p].Prev = -2;
  }
  else{
    if(pl >= 0) Pc[pl].Next = p;
    Pc[p].Prev = pl;
  }
  //printf("%d - %d - %d c) %d last %d coord %x cell %d (%lf %lf %lf)\n",Pc[p].Prev,p,Pc[p].Next,c,Cella[c].Last,Pc[p].Coord,Pc[p].Cell,Pc[p].Pos[0],Pc[p].Pos[1],Pc[p].Pos[2]);
}
/** Remove particle from the cell */
void DdLinkedList::RemPart(const int p,double *Pos){
  int c = pCella(Pos);
  RemPart(p,c);
}
/** Remove a part to the correspondent cell*/
void DdLinkedList::RemPart(const int p,const int c){
  SigErr(p >= NAllocP,"More particles than allocated\n");
  //assert(Cella[c].NPart > 0);
  SigErr(Cella[c].NPart <= 0,"Cannot remove, the cell is already empty\n");
  SigErr(c < 0,"The cell %d does not exist\n",c);
  Cella[c].NPart--;
  NPart--;
  if(Cella[c].NPart == 0){
    Cella[c].First = -2;
    Cella[c].Last  = -1;
  }
  else{
    int pp = Pc[p].Prev;
    int pn = Pc[p].Next;
    if(pp == -2) Cella[c].First = pn;
    else Pc[pp].Next = pn;
    if(pn == -1) Cella[c].Last = pp;
    else Pc[pn].Prev = pp;
    Pc[p].Cell = -1;
  }
  //printf("%d - %d - %d c) %d last %d coord %x cell %d\n",pp,p,pn,c,Cella[c].Last,Pc[p].Coord,Pc[p].Cell);
  //  assert(pp != -2 && Pc[pp].Next != pp);
}
/** Remove a part to the correspondent cell*/
void DdLinkedList::RemPart(const int p){
  int c = Pc[p].Cell;
  SigErr(c < 0,"The cell %d or the particle %d does not exist\n",c,p);
  RemPart(p,c);
}
/** Swap two particles*/
int DdLinkedList::SwapPart(int p1,double *Pos1,int p2,double *Pos2){
  if(p1==p2) return 0;
  DOMAIN_PART Pn;
  int c1 = Pc[p1].Cell;
  int c2 = Pc[p2].Cell;
  for(int d=0;d<3;d++){
    Pc[p1].Pos[d] = Pos2[d];
    Pc[p2].Pos[d] = Pos1[d];
  }
  if(c1!=c2){
    if(Cella[c2].First == p2)
      Cella[c2].First  = p1;
    if(Cella[c2].Last == p2)
      Cella[c2].Last  = p1;
    if(Cella[c1].First == p1)
      Cella[c1].First  = p2;
    if(Cella[c1].Last == p1)
      Cella[c1].Last  = p2;
  }
  if(c1==c2){
    if(Cella[c1].First == p2)
      Cella[c1].First = p1;
    else if(Cella[c1].First == p1)
      Cella[c1].First = p2;
    if(Cella[c1].Last == p2)
      Cella[c1].Last = p1;
    else if(Cella[c1].Last == p1)
      Cella[c1].Last = p2;
  }
  int p1p = Pc[p1].Prev;
  int p2p = Pc[p2].Prev;
  int p1n = Pc[p1].Next;
  int p2n = Pc[p2].Next;
  if(p1p >= 0) Pc[p1p].Next = p2;
  if(p1n >= 0) Pc[p1n].Prev = p2;
  if(p2p >= 0) Pc[p2p].Next = p1;
  if(p2n >= 0) Pc[p2n].Prev = p1;
  Pn.Prev  = Pc[p1].Prev;
  Pn.Next  = Pc[p1].Next;
  Pn.Coord = Pc[p1].Coord;
  Pn.Cell  = Pc[p1].Cell;
  Pc[p1].Prev  = Pc[p2].Prev;
  Pc[p1].Next  = Pc[p2].Next;
  Pc[p1].Coord = Pc[p2].Coord;
  Pc[p1].Cell  = Pc[p2].Cell;
  Pc[p2].Prev  = Pn.Prev;
  Pc[p2].Next  = Pn.Next;
  Pc[p2].Coord = Pn.Coord;
  Pc[p2].Cell  = Pn.Cell;
}
/** Shift a particle from one position to its new*/
void DdLinkedList::MovePart(const int p,double *NewPos){
  int cn = pCella(NewPos);
  int co = pCell(p);
  SetCoorNumb(NewPos,p);
  for(int d=0;d<3;d++)Pc[p].Pos[d] = NewPos[d];
  if(cn == co) return;
  RemPart(p,co);
  AddPart(p,cn);
}
/** Shift a particle from one position to its new*/
void DdLinkedList::MovePart(const int p,double *OldPos,double *NewPos){
  int cn = pCella(NewPos);
  int co = pCella(OldPos);
  SetCoorNumb(NewPos,p);
  for(int d=0;d<3;d++)Pc[p].Pos[d] = NewPos[d];
  if(cn == co) return ;
  RemPart(p,co);
  AddPart(p,cn);
}
/** Set the counters to the first particle of the cell c1 for the first loop*/
void DdLinkedList::SetCounters(int c1){
  for(int c=0;c<NCell;c++){
    Cella[c].Curr1 = Cella[c].First;
    Cella[c].Curr2 = Cella[c].First;
    if(Cella[c].NPart > 1)
      Cella[c].Curr2 = Pc[Cella[c].Curr1].Next;
  }
}
int DdLinkedList::ItCell(const int c){
  SigErr(c >= NCell,"The cell %d does not exist\n",c);
  SigErr(Cella[c].Curr1 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr1,NAllocP);
  return Cella[c].Curr1;
}
/** Return 0 when the loop inside the cell is over*/
int DdLinkedList::IfItCell(const int c){
  if(Cella[c].Curr1 < 0 || Cella[c].NPart == 0){
    Cella[c].Curr1 = Cella[c].First;
    return 0;
  }
  return 1;
}
/** Increment the iterator to the next particle*/
void DdLinkedList::IncrCurr(const int c){
  Cella[c].Curr1 = Pc[Cella[c].Curr1].Next;
}
/** Print the content of the cell */
void DdLinkedList::PrintCell(const int c){
  for(SetCounters(c);IfItCell(c);IncrCurr(c)){
    int p = ItCell(c);
    printf("%d) # %d %d_%d_%d %x\n",c,Cella[c].NPart,Pc[p].Prev,ItCell(c),Pc[p].Next,Pc[p].Coord);
  }
}
/** Print the content of all cells*/
void DdLinkedList::PrintCells(){
  for(int c=0;c<NCell;c++)
    PrintCell(c);
}
/** 
    Set the iterators for the current particle and build the list of neighbouring cells.
*/
void DdLinkedList::SetCurr(int p){
  cCurr = Pc[p].Cell;
  SigErr(cCurr < 0,"The cell %d or the particle %d does not exists\n",cCurr,p);
  NNeiCurr = GetNei(Pc[p].Pos,NeiListCurr);
  nNeiCurr = 0;
  IfLoopCurr = 1;
  cCurr = NeiListCurr[nNeiCurr];
  p1Curr = p;
  p2Curr = Cella[cCurr].First;
  if(p2Curr == p1Curr) Pc[p2Curr].Next;
  while(p2Curr < 0){
    nNeiCurr++;
    //printf(" %d)  %d/%d %d %d\n",p,nNeiCurr,NNeiCurr,cCurr,p2Curr);
    if(nNeiCurr==NNeiCurr){
      IfLoopCurr = 0;
      break;
    }
    cCurr = NeiListCurr[nNeiCurr];
    p2Curr = Cella[cCurr].First;
  }
  //printf("%d) %d %d/%d %d\n",p,p2Curr,nNeiCurr,NNeiCurr,cCurr);
  //for(int n=0;n<27;n++)printf("  %d] %0.f %0.f %0.f\n",n,BoundCond[n][0],BoundCond[n][1],BoundCond[n][2]);
  //PrintCell(cCurr);
  Cella[cCurr].Curr1 = p1Curr;
  Cella[cCurr].Curr2 = p2Curr;
}
void DdLinkedList::NextCurr(){
  p2Curr = Pc[p2Curr].Next;
  if(p2Curr == p1Curr) Pc[p2Curr].Next;
  while(p2Curr < 0){
    nNeiCurr++;
    if(nNeiCurr==NNeiCurr){
      IfLoopCurr = 0;
      break;
    }
    cCurr = NeiListCurr[nNeiCurr];
    p2Curr = Cella[cCurr].First;
  }
  //printf("   %d) %d %d\n",p1Curr,p2Curr,nNeiCurr);
  Cella[cCurr].Curr1 = p1Curr;
  Cella[cCurr].Curr2 = p2Curr;
}
/** 
    Set the iterators for the current ghost particle and build the list of neighbouring cells.
*/
int DdLinkedList::IfCurr(){
  return IfLoopCurr;
}
/** Iterate one step and return the position */
void DdLinkedList::Dist2Curr(double *DistRel){
  for(int d=0;d<3;d++){
    DistRel[d] = Pc[p1Curr].Pos[d] - Pc[p2Curr].Pos[d];
    //TOFIX
    //printf("%d-%d) %.0f %.0f  %d %d-%d\n",p1Curr,p2Curr,BoundCond[nNeiCurr][d]*Edge[d],-floor(DistRel[d]*InvEdge[d] + .5)*Edge[d],nNeiCurr,Pc[p1Curr].Cell,Pc[p2Curr].Cell);
    //DistRel[d] += BoundCond[nNeiCurr][d]*Edge[d];
    DistRel[d] -= floor(DistRel[d]*InvEdge[d] + .5)*Edge[d];
  }
  DistRel[3] = SQR(DistRel[0]) + SQR(DistRel[1]) + SQR(DistRel[2]);
}
void DdLinkedList::SetCurrGhost(double *Pos){
  cCurr = pCella(Pos);
  SigErr(cCurr < 0,"The cell %d does not exists\n",cCurr);
  NNeiCurr = GetNei(Pos,NeiListCurr);
  nNeiCurr = 0;
  IfLoopCurr = 1;
  cCurr = NeiListCurr[nNeiCurr];
  for(int d=0;d<3;d++) PosCurr[d] = Pos[d];
  p2Curr = Cella[cCurr].First;
  while(p2Curr < 0){
    nNeiCurr++;
    //printf(" %d)  %d/%d %d %d\n",p,nNeiCurr,NNeiCurr,cCurr,p2Curr);
    if(nNeiCurr==NNeiCurr){
      IfLoopCurr = 0;
      break;
    }
    cCurr = NeiListCurr[nNeiCurr];
    p2Curr = Cella[cCurr].First;
  }
  //printf("%d) %d %d/%d %d\n",p,p2Curr,nNeiCurr,NNeiCurr,cCurr);
  //for(int n=0;n<27;n++)printf("  %d] %0.f %0.f %0.f\n",n,BoundCond[n][0],BoundCond[n][1],BoundCond[n][2]);
  //PrintCell(cCurr);
  Cella[cCurr].Curr2 = p2Curr;
}
void DdLinkedList::NextCurrGhost(){
  p2Curr = Pc[p2Curr].Next;
  while(p2Curr < 0){
    nNeiCurr++;
    if(nNeiCurr==NNeiCurr){
      IfLoopCurr = 0;
      break;
    }
    cCurr = NeiListCurr[nNeiCurr];
    p2Curr = Cella[cCurr].First;
  }
  //printf("   %d) %d %d\n",p1Curr,p2Curr,nNeiCurr);
  Cella[cCurr].Curr2 = p2Curr;
}
int DdLinkedList::IfCurrGhost(){
  return IfLoopCurr;
}
/** Iterate one step and return the position */
void DdLinkedList::Dist2CurrGhost(double *DistRel){
  for(int d=0;d<3;d++){
    DistRel[d] = PosCurr[d] - Pc[p2Curr].Pos[d];
    //TOFIX
    //printf("%d-%d) %.0f %.0f  %d %d-%d\n",p1Curr,p2Curr,BoundCond[nNeiCurr][d]*Edge[d],-floor(DistRel[d]*InvEdge[d] + .5)*Edge[d],nNeiCurr,Pc[p1Curr].Cell,Pc[p2Curr].Cell);
    //DistRel[d] += BoundCond[nNeiCurr][d]*Edge[d];
    DistRel[d] -= floor(DistRel[d]*InvEdge[d] + .5)*Edge[d];
  }
  DistRel[3] = SQR(DistRel[0]) + SQR(DistRel[1]) + SQR(DistRel[2]);
}
/** Associate the two iterator of the cell to the particles p1 and p2*/
void DdLinkedList::Couple(const int c,int *p1,int *p2){
  SigErr(Cella[c].Curr1 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr1,NAllocP);
  SigErr(Cella[c].Curr2 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr2,NAllocP);
  *p1 = Cella[c].Curr1;
  *p2 = Cella[c].Curr2;
}
/** Return zero when both iterators are over the loop */
int DdLinkedList::IfItCouple(const int c){
  if(Cella[c].NPart < 2) return 0;
  if(Cella[c].Curr1 == -1){
    Cella[c].Curr1 = Cella[c].First;
    Cella[c].Curr2 = Pc[Cella[c].First].Next;
    return 0;
  }
  return 1;
}
/** Increment the second iterator and the first when the second is over the loop*/
void DdLinkedList::IncrCurrList(const int c){
  int p2 = Cella[c].Curr2;
  Cella[c].Curr2 = Pc[p2].Next;
  if(Pc[p2].Next == -1){
    int p1 = Cella[c].Curr1;
    Cella[c].Curr1 = Pc[p1].Next;
    Cella[c].Curr2 = Pc[Pc[p1].Next].Next;
    if(Cella[c].Curr2 == -1){
      Cella[c].Curr1 = -1;
    }
  }
  SigErr(Cella[c].Curr2 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr2,NAllocP);
}
int DdLinkedList::FindClosest(int p1){
  double DistRel[4];
  int pCloser = -1;
  double DistCurr = 1000000.;
  for(SetCurr(p1);IfCurr();NextCurr()){
    if(p2Curr <= p1Curr) continue;
    Dist2Curr(DistRel);
    if(DistRel[3] < DistCurr){
      DistCurr = DistRel[3];
      pCloser = p2Curr;
    }
  }
  return pCloser;
}
//--------------------DomPart-----DomCell---------------------
/** Increase the iterator, deactivated*/
DomCell& DomCell::operator++(){
  //It++;
  return *this;
}
/** Increase the iterator, deactivated*/
DomPart& DomPart::operator++(){
}
/** Return the cell class*/
DomCell& DomCell::operator=(const DomCell &Dc){
  return *this;
}
DomPart::DomPart(int ExtNPart){
  NPart = ExtNPart;
  Next = new int[NPart];
  Prev = new int[NPart];
}
int DomPart::operator[](int col){
  return Next[col];
} 
//-------------------------------DdArray------------------------
/**
   Every cell contains a stl list of particle.
 */
DdArray::DdArray(double EdgeExt[3],int NPartExt,double CutOffExt){
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
    InvEdge[d] = 1./Edge[d];
    NSect[d] = (int)(Edge[d]/(double)(CutOffExt));
  }
  NPart = 0;
  CutOff = CutOffExt;
  Mod10[0] = 1;
  Mod10[1] = NSect[0];
  Mod10[2] = NSect[1]*Mod10[1];
  NCell = NSect[0]*NSect[1]*NSect[2];
  Cella = new DdCell[NCell];
  if(Cella == NULL) {printf("DomCell not allocated\n");}
}
void DdArray::Erase(){
  for(int c=0;c<NCell;c++){
    Cella[c].Part.clear();
  }
}
void DdArray::AddPart(int p,double *Pos){
  int c = pCella(Pos);
  Cella[c].Part.push_back(p);
  NPart++;
}
void DdArray::RemPart(int p,double *Pos){
  int c = pCella(Pos);
  Cella[c].Part.remove(p);
  NPart--;
}
void DdArray::MovePart(int p,double *OldPos,double *NewPos){
  int c1 = pCella(OldPos);
  int c2 = pCella(NewPos);
  if(c1 == c2) return;
  Cella[c1].Part.remove(p);
  Cella[c2].Part.push_back(p);
}
void DdArray::SwapPart(int p1,double *Pos1,int p2,double *Pos2){
  int c1 = pCella(Pos1);
  int c2 = pCella(Pos2);
  if(c1 == c2) return;
  Cella[c1].Part.remove(p1);
  Cella[c1].Part.push_back(p2);
  Cella[c2].Part.remove(p2);
  Cella[c2].Part.push_back(p2);
}
void DdArray::SetCounters(int c){
  //if(Cella[c].Part.size() == 0) return;
  NCurr = Cella[c].Part.begin();
  NCurr2 = Cella[c].Part.begin();
  if(Cella[c].Part.size() > 1)
    ++NCurr2;
};
int DdArray::IfItCell(int c){
  if(NCurr == Cella[c].Part.end()) return 0;
  return 1;
};
void DdArray::IncrCurr(int c){
  ++NCurr;
}
int DdArray::ItCell(int c){
  return *NCurr;
}
void DdArray::Couple(const int c,int *p1,int *p2){
  *p1 = *NCurr;
  *p2 = *NCurr2;
}
int DdArray::IfItCouple(const int c){
  if(NCurr == Cella[c].Part.end())
    return 0;
  return 1;
}
void DdArray::IncrCurrList(const int c){
  ++NCurr2;
  if(NCurr2 == Cella[c].Part.end()){
    NCurr2 = NCurr;
    ++NCurr2;
    ++NCurr;
  }
}
void DdArray::PrintCell(const int c){
  for(SetCounters(c);IfItCell(c);IncrCurr(c)){
    int p = ItCell(c);
    printf("%d) # %d %d\n",c,Cella[c].Part.size(),ItCell(c));
  }
}
void DdArray::PrintCells(){
  for(int c=0;c<NCell;c++)
    PrintCell(c);
}
//--------------------------------DdFixedSize----------------------------
/**
   Every cell contains a stl list of particle.
 */
DdFixedSize::DdFixedSize(double EdgeExt[3],int NPartExt,double CutOffExt){
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
    InvEdge[d] = 1./Edge[d];
    NSect[d] = (int)(Edge[d]/(double)(CutOffExt));
  }
  NPart = 0;
  NPCell = 30;
  CutOff = CutOffExt;
  Mod10[0] = 1;
  Mod10[1] = NSect[0];
  Mod10[2] = NSect[1]*Mod10[1];
  NCell = NSect[0]*NSect[1]*NSect[2];
  Cella = new DdFixCell[NCell];
  for(int c=0;c<NCell;c++){
    Cella[c].Part = new int[NPCell];
  }
  if(Cella == NULL) {printf("DomCell not allocated\n");}
}
void DdFixedSize::Erase(){
  for(int c=0;c<NCell;c++){
    for(int p=0;p<Cella[c].NPart;p++){
      Cella[c].Part[p] = -1;
    }
    Cella[c].NPart = 0;
  }
}
void DdFixedSize::AddPart(int p,double *Pos){
  int c = pCella(Pos);
  AddPart(p,c);
}
void DdFixedSize::AddPart(int p,int c){
  //printf("Adding %d in %d\n",p,c);
  if(Cella[c].NPart > NPCell){
    printf("Maximum allocated size in the domain decomposition (%d) reached, terminating\n",NPCell);
    assert(Cella[c].NPart < NPCell);
  }
  Cella[c].Part[Cella[c].NPart++] = p;
  NPart++;
}
void DdFixedSize::RemPart(int p1,double *Pos){
  int c = pCella(Pos);
  RemPart(p1,c);
}
void DdFixedSize::RemPart(int p1,int c){
  //PrintCell(c);
  //printf("removing %d from %d\n",p1,c);
  int FoundPartInCell = 0;
  int Temp = 0;
  for(int p=0;p<Cella[c].NPart;p++){
    if(p1==Cella[c].Part[p]){
      for(int pp=p;pp<Cella[c].NPart-1;pp++){
	Temp = Cella[c].Part[pp];
	Cella[c].Part[pp] = Cella[c].Part[pp+1];
	Cella[c].Part[pp+1] = Temp;
      }
      Cella[c].NPart--;
      NPart--;
      FoundPartInCell = 1;
      break;
    }
  }
  //PrintCell(c);
  assert(FoundPartInCell);
}
void DdFixedSize::MovePart(int p,double *OldPos,double *NewPos){
  int c1 = pCella(OldPos);
  int c2 = pCella(NewPos);
  if(c1 == c2) return;
  RemPart(p,c1);
  AddPart(p,c2);
}
void DdFixedSize::SwapPart(int p1,double *Pos1,int p2,double *Pos2){
  int c1 = pCella(Pos1);
  int c2 = pCella(Pos2);
  if(c1 == c2) return;
  RemPart(p1,c1);
  RemPart(p2,c2);
  AddPart(p1,c2);
  AddPart(p2,c1);
}
void DdFixedSize::SetCounters(int c){
  //if(Cella[c].Part.size() == 0) return;
  NCurr = 0;
  NCurr2 = 0;
  if(Cella[c].NPart > 0) NCurr2++;
};
int DdFixedSize::IfItCell(int c){
  if(NCurr == Cella[c].NPart) return 0;
  return 1;
};
void DdFixedSize::IncrCurr(int c){
  NCurr++;
}
int DdFixedSize::ItCell(int c){
  return Cella[c].Part[NCurr];
}
void DdFixedSize::Couple(const int c,int *p1,int *p2){
  //  printf("%d %d %d\n",NCurr,NCurr2,Cella[c].NPart);
  *p1 = Cella[c].Part[NCurr];
  *p2 = Cella[c].Part[NCurr2];
}
int DdFixedSize::IfItCouple(const int c){
  if(NCurr >= Cella[c].NPart)
    return 0;
  return 1;
}
void DdFixedSize::IncrCurrList(const int c){
  NCurr2++;
  if(NCurr2 >= Cella[c].NPart){
    NCurr2 = NCurr;
    NCurr2++;
    NCurr++;
  }
}
void DdFixedSize::PrintCell(const int c){
  for(SetCounters(c);IfItCell(c);IncrCurr(c)){
    int p = ItCell(c);
    printf("%d) # %d %d\n",c,Cella[c].NPart,ItCell(c));
  }
}
void DdFixedSize::PrintCells(){
  for(int c=0;c<NCell;c++)
    PrintCell(c);
}
//-------------------DOUBLE-LOOP-------------------------------
/** Constructor for the domain decomposition with the linked list*/
DdDoubleLoop::DdDoubleLoop(double EdgeExt[3],int NPartExt,double CutOffExt){
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
    InvEdge[d] = 1./Edge[d];
    NSect[d] = 1;
  }
  NPart = 0;
  SetCutOff(CutOffExt);
  Mod10[0] = 1;
  Mod10[1] = 1;
  Mod10[2] = 1;
  NCell = 0;
  Cella = new DomCell[NCell];
  SigErr(Cella==NULL,"DomCell not allocated\n");
  for(int c=0;c<NCell;c++){
    Cella[c].NPart = 0;
    Cella[c].First = Cella[c].Last = -2;
  }
  NAllocP = NPartExt + 2000;
  Pc = new DOMAIN_PART[NAllocP];
  SigErr(Pc == NULL,"DOMAIN_PART not allocated\n");
}
/** Empty the records of the cells*/
void DdDoubleLoop::Erase(){
  for(int c=0;c<NCell;c++){
    Cella[c].NPart = 0;
    Cella[c].Last = -1;
  }
  for(int p=0;p<NPart;p++){
    Pc[p].Cell = -1;
    Pc[p].Next = -1;
    Pc[p].Prev = -2;
  }
  NPart = 0;
}
/** Coordination number of the particle in the cell*/
int DdDoubleLoop::SetCoorNumb(double *Pos,int p){
  Pc[p].Coord = GetCoorNumb(Pos);
  return Pc[p].Coord;
}
/** Add a part to the correspondent cell*/
void DdDoubleLoop::AddPart(const int p,double *Pos){
  for(int d=0;d<3;d++) Pc[p].Pos[d] = Pos[d];
  AddPart(p,0);
}
/** Add a part to the correspondent cell*/
void DdDoubleLoop::AddPart(const int p,const int c){
  Cella[c].NPart++;
  NPart++;
  Cella[c].Last = p;
  Pc[p].Cell = c;
  Pc[p].Next = -1;
}
/** Remove particle from the cell */
void DdDoubleLoop::RemPart(const int p,double *Pos){
  RemPart(p,0);
}
/** Remove a part to the correspondent cell*/
void DdDoubleLoop::RemPart(const int p,const int c){
  Cella[c].NPart--;
  NPart--;
  if(Cella[c].NPart == 0){
    Cella[c].First = -2;
    Cella[c].Last  = -1;
  }
  else{
    int pp = Pc[p].Prev;
    int pn = Pc[p].Next;
    if(pp == -2) Cella[c].First = pn;
    else Pc[pp].Next = pn;
    if(pn == -1) Cella[c].Last = pp;
    else Pc[pn].Prev = pp;
    Pc[p].Cell = -1;
  }
}
/** Remove a part to the correspondent cell*/
void DdDoubleLoop::RemPart(const int p){
  RemPart(p,0);
}
/** Swap two particles*/
int DdDoubleLoop::SwapPart(int p1,double *Pos1,int p2,double *Pos2){
  printf("Swap part disabled\n");
}
/** Shift a particle from one position to its new*/
void DdDoubleLoop::MovePart(const int p,double *NewPos){
  printf("Move part disabled\n");
}
/** Shift a particle from one position to its new*/
void DdDoubleLoop::MovePart(const int p,double *OldPos,double *NewPos){
  printf("Move part disabled\n");
}
/** Set the counters to the first particle of the cell c1 for the first loop*/
void DdDoubleLoop::SetCounters(int c){
  p1Curr = 0;
  p2Curr = 1;
  cCurr = 0;
}
int DdDoubleLoop::ItCell(const int c){
  return 0;
}
/** Return 0 when the loop inside the cell is over*/
int DdDoubleLoop::IfItCell(const int c){
  if(p1Curr == NPart) return 0;
  return 1;
}
/** Increment the iterator to the next particle*/
void DdDoubleLoop::IncrCurr(const int c){
  p2Curr++;
  if(p2Curr <= NPart){
      p1Curr++;
      p2Curr = p1Curr+1;
  }
}
/** Print the content of the cell */
void DdDoubleLoop::PrintCell(const int c){
  for(SetCounters(c);IfItCell(c);IncrCurr(c)){
    int p = ItCell(c);
    printf("%d) # %d %d_%d_%d %x\n",c,Cella[c].NPart,Pc[p].Prev,ItCell(c),Pc[p].Next,Pc[p].Coord);
  }
}
/** Print the content of all cells*/
void DdDoubleLoop::PrintCells(){
  for(int c=0;c<NCell;c++)
    PrintCell(c);
}
/** */
void DdDoubleLoop::SetCurr(int p){
  cCurr = 0;
  p1Curr = p;
  p2Curr = p+1;
}
void DdDoubleLoop::NextCurr(){
  p2Curr++;
}
int DdDoubleLoop::IfCurr(){
  if(p2Curr == NPart) return 0;
  return 1;
}
/** Iterate one step and return the position */
void DdDoubleLoop::Dist2Curr(double *DistRel){
  for(int d=0;d<3;d++){
    DistRel[d] = Pc[p1Curr].Pos[d] - Pc[p2Curr].Pos[d];
    //TOFIX
    //DistRel[d] += BoundCond[nNeiCurr][d]*Edge[d];
    DistRel[d] -= floor(DistRel[d]*InvEdge[d] + .5)*Edge[d];
  }
  DistRel[3] = SQR(DistRel[0]) + SQR(DistRel[1]) + SQR(DistRel[2]);
}
//-----------------------------NeiVertex------------------------
/**
   
 */
NeiVertex::NeiVertex(int NTriaExt,int NvPtExt,int NGridExt,double *EdgeExt){
  NTria = NTriaExt;
  NvPt = NvPtExt;
  NVert = NTria*NvPt;
  NGrid = 2*NGridExt-1;
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
  }
  Vertex = new VERTEX[NVert];
  vMemPos = new int[NVert];
  for(int v=0;v<NVert;v++){
    Vertex[v].NTria = 0;
  }
  NVert = 0;
}
NeiVertex::~NeiVertex(){
  delete [] Vertex;
}
int NeiVertex::GetVertex(double *Pos){
  int s[3] = {0,0,0};
  for(int d=0;d<3;d++){
    s[d] = (int)(Pos[d]/Edge[d]*NGrid);
  }
  return (s[0]*NGrid+s[1])*NGrid+s[2];
}
void NeiVertex::Add(double *Pos,int t){
  int v = GetVertex(Pos);
  Add(v,t,Pos);
}
void NeiVertex::Add(int v,int t,double *Pos){
  Vertex[NVert].v = v;
  for(int d=0;d<3;d++)
    Vertex[NVert].Pos[d] = Pos[d];
  Vertex[NVert].t[0] = t;
  //printf("Add %d %d %d\n",NVert,Vertex[NVert].v,Vertex[NVert].t[0]);
  NVert++;
}
void NeiVertex::CopyVert2To1(VERTEX Vert1,VERTEX Vert2){
  Vert1.v = Vert2.v;
  Vert1.NTria = Vert2.NTria;
  for(int d=0;d<3;d++) Vert1.Pos[d] = Vert2.Pos[d];
  for(int t=0;t<Vert1.NTria;t++) Vert1.t[t] = Vert2.t[t];
}
void NeiVertex::Swap(int v1,int v2){
  VERTEX Vert = Vertex[v1];
  Vertex[v1] = Vertex[v2];
  Vertex[v2] = Vert;
  CopyVert2To1(Vert,Vertex[v1]);
  CopyVert2To1(Vertex[v1],Vertex[v2]);
  CopyVert2To1(Vertex[v2],Vert);
}
void NeiVertex::Reorder(){
  //sorting
  for(int i=0;i<NVert;i++){
    for(int j=i;j>0;j--){
      if(Vertex[j].v < Vertex[j-1].v){
  	Swap(j-1,j);
      }
      else
  	break;
    }
  }
  for(int v=0;v<NVert;v++){
    int NDel = 0;
    Vertex[v].NTria = 1;
    for(int vv=v+1;Vertex[vv].v == Vertex[vv-1].v;vv++){
      if(vv > NVert) break;
      int t = vv-v;
      Vertex[v].t[t] = Vertex[vv].t[0];
      Vertex[v].NTria = t+1;
      NDel++;
      SigErr(NDel>MAX_JOINT_VERTEX,"Neighbour Vertices: the vertex has more conjunction points than expected\n");
    }
    for(int d=1;d<NDel+1;d++){
      for(int vv=v+1;vv<NVert-1;vv++){
        Swap(vv,vv+1);
      }
    }
    NVert -= NDel;
  }
}
void NeiVertex::SetCounters(){
  vCurr = 0;
  tCurr = 0;
}
/** Initialize the */
void NeiVertex::SetCounters(int v){
  for(int vv=vCurr;vv<NVert;vv++){
    if(Vertex[vv].v == v){
      vCurr = vv;
      break;
    }
  }
  tCurr = 0;
}
int NeiVertex::IfItCell(int v){
  if(tCurr == Vertex[vCurr].NTria) return 0;
  return 1;
}
void NeiVertex::IncrCurr(int v){
  tCurr++;
}
int NeiVertex::VertCurr(int v){
  return vCurr;
}
int NeiVertex::TriaCurr(int v){
  return Vertex[vCurr].t[tCurr];
}
void NeiVertex::Print(){
  SetCounters();
  for(int v=0;v<NVert;v++){
    printf("---- (%d) ----\n",Vertex[v].v);
    for(SetCounters(v);IfItCell(v);IncrCurr(v)){
      printf("%d - ",VertCurr(v));
    }
    printf("\n");
  }
  SetCounters();
}
