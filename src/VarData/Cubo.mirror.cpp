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
void DomDecBasics::SigErr(int IsWrong,const char * s, ...){
#ifdef DEBUG
  if(IsWrong){
    va_list args;
    va_start(args, s);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, s, args);
    fprintf(stderr, "exiting the program\n");
    va_end(args);
    abort();
  }
#else  
  return;
#endif
}
DomDecBasics::DomDecBasics(double EdgeExt[3],int NPartExt,double CutOffExt){
  IfMirror = 1;
  IfBorder = 0;
  double CellSize = IfBorder ? CutOffExt : 2.*CutOffExt;
  for(int d=0;d<3;d++){
    Edge[d] = EdgeExt[d];
    InvEdge[d] = 1./Edge[d];
    NSect[d] = (int)(Edge[d]/CellSize) + 2*IfMirror;
    SigErr(CutOffExt > 2.*Edge[d],"Insufficient number of cells \n");
  }
  NPart = 0;
  SetCutOff(CutOffExt);
  Mod10[0] = 1;
  Mod10[1] = NSect[0];
  Mod10[2] = NSect[1]*Mod10[1];
  NCell = NSect[0]*NSect[1]*NSect[2];
}
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
int DomDecBasics::pNPart(){
  return NPart;
}
int DomDecBasics::pNCell(){
  return NCell;
}
int DomDecBasics::GetCell(double *Pos,int *NeiList){
  int c = pCella(Pos);
  return GetCellCh(c,NeiList);
}
int DomDecBasics::pCella(const double Pos[3]){
  int ris[3];
  for(int d=0;d<3;d++){
    ris[d] = (int)(Pos[d]*InvEdge[d]*NSect[d]);
  }
  return pCella(ris);
}
int DomDecBasics::pCella(const double *Pos,int *cd){
  for(int d=0;d<3;d++){
    cd[d] = (int)(Pos[d]*InvEdge[d]*NSect[d]);
  }
  return pCella(cd);  
}
int DomDecBasics::pCella(const int c[3]){
  //+1 for the mirror notation
  int risp = (c[0]+IfMirror)*Mod10[0]+(c[1]+IfMirror)*Mod10[1]+(c[2]+IfMirror)*Mod10[2];
  SigErr(risp >= NCell,"Probably the particle position is not a number or not in the box\n");
  return risp;
}
int DomDecBasics::GetCellChMirror(int c[3],int *NeiList){
  int NeiMod[3];
  int NeiCell[3];
  int NNei = 0;
  for(int cx=-1;cx<=1;cx++){
    NeiCell[0] = NeiMod[0] + cx;
    for(int cy=-1;cy<=1;cy++){
      NeiCell[1] = NeiMod[1] + cy;
      for(int cz=-1;cz<=1;cz++){
	NeiList[NNei++] = NeiCell[0]*Mod10[0] + NeiCell[1]*Mod10[1] + NeiCell[2]*Mod10[2];
      }
    }
  }
  return NNei;
}
int DomDecBasics::GetCellCh(int c,int *NeiList){
  int cd[3];
  cd[2] = (int)(c/(double)Mod10[2]);
  cd[1] = (int)((c-cd[2]*Mod10[2])/(double)Mod10[1]);
  cd[0] = (int)(c-cd[1]*Mod10[1]-cd[2]*Mod10[2]);
  GetCellCh(cd,NeiList);
}
int DomDecBasics::GetCellCh(int *cd,int *NeiList){
  int NeiMod[3] = {cd[0],cd[1],cd[2]};
  int NeiCell[3];
  int NNei = 0;
  for(int cx=-1;cx<=1;cx++){
    NeiCell[0] = NeiMod[0] + cx;
    if(NeiCell[0] >= NSect[0]) NeiCell[0] -= NSect[0];
    if(NeiCell[0] < 0) NeiCell[0] += NSect[0];
    for(int cy=-1;cy<=1;cy++){
      NeiCell[1] = NeiMod[1] + cy;
      if(NeiCell[1] >= NSect[1]) NeiCell[1] -= NSect[1];
      if(NeiCell[1] < 0) NeiCell[1] += NSect[1];
     for(int cz=-1;cz<=1;cz++){
	NeiCell[2] = NeiMod[2] + cz;
	if(NeiCell[2] >= NSect[2]) NeiCell[2] -= NSect[2];
	if(NeiCell[2] < 0) NeiCell[2] += NSect[2];
	NeiList[NNei++] = NeiCell[0]*Mod10[0] + NeiCell[1]*Mod10[1] + NeiCell[2]*Mod10[2];
      }
    }
  }
  return NNei;
}
int DomDecBasics::GetCellCoord(double *Pos,int *NeiList){
  int c = pCella(Pos);
  int Coord = GetCoorNumb(Pos);
  return GetCellCoord(c,Coord,NeiList);
}
int DomDecBasics::GetCellCoord(int c,int Coord,int *NeiList){
  SigErr(c < 0,"Probably the particle position is not a number\n");
  int NNei = 0;
  int nNei[3];
  int NeiCell[3];
  int Mask[3]    = {MASK_X,MASK_Y,MASK_Z};
  int GetMask[3] = {GET_MASK_X,GET_MASK_Y,GET_MASK_Z};
  int Move[3] = {0,0,0};
  NeiCell[2] = (int)(c/(double)Mod10[2]);
  NeiCell[1] = (int)((c-NeiCell[2]*Mod10[2])/(double)Mod10[1]);
  NeiCell[0] = (int)(c-NeiCell[1]*Mod10[1]-NeiCell[2]*Mod10[2]);
  for(int d=0;d<3;d++){
    int d1 = (d+1)%3;
    int d2 = (d+2)%3;
    if( DOM_GET_POS(Coord,Mask[d],GetMask[d]) == DOM_LEFT ){
      Move[d] = DOM_LEFT;
      int nc = NeiCell[d] - 1;
      if(nc < 0) nc = NSect[d]-1;
      NeiList[NNei++] = nc*Mod10[d] + NeiCell[d1]*Mod10[d1] + NeiCell[d2]*Mod10[d2];
    }
    else if( DOM_GET_POS(Coord,Mask[d],GetMask[d]) == DOM_RIGHT ){
      Move[d] = DOM_RIGHT;
      int nc = NeiCell[d] + 1;
      if(nc >= NSect[d]) nc = 0;
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
void DdLinkedList::Erase(){
  for(int c=0;c<NCell;c++){
    Cella[c].NPart = 0;
    Cella[c].Last = -1;
  }
  // for(int p=0;p<NPart;p++){
  //   Pc[p].Cell = 0;
  //   Pc[p].Next = -1;
  //   Pc[p].Prev = -2;
  // }
  NPart = 0;
}
int DdLinkedList::SetCoorNumb(double *Pos,int p){
  Pc[p].Coord = GetCoorNumb(Pos);
  return Pc[p].Coord;
}
void DdLinkedList::SetCounters(int c1){
  for(int c=0;c<NCell;c++){
    Cella[c].Curr1 = Cella[c].First;
    Cella[c].Curr2 = Cella[c].First;
    if(Cella[c].NPart > 1)
      Cella[c].Curr2 = Pc[Cella[c].Curr1].Next;
  }
}
void DdLinkedList::AddPart(const int p,double *Pos){
  int cd[3];
  pCella(Pos,cd);
  SetCoorNumb(Pos,p);
  int c = pCella(cd);
  int pl = Cella[c].Last;
  SigErr(pl == p,"The particle %d is already the last in the cell %d\n",p,c);
  SigErr(p >= NAllocP,"More particle than allocated\n");
  SigErr(Cella[c].NPart >= NAllocP,"More particle than allocated\n");
  Cella[c].NPart++;
  NPart++;
  SigErr(pl >= NAllocP,"More particle than allocated\n");
  Cella[c].Last = p;
  for(int d=0;d<3;d++){
    Pc[p].Cell[d] = cd[d];
    Pc[p].Pos[d] = Pos[d];
  }
  AddPartMirror(p);
  Pc[p].Next = -1;
  if(Cella[c].NPart == 1){
    Cella[c].First = p;
    Pc[p].Prev = -2;
    return ;
  }
  if(pl >= 0) Pc[pl].Next = p;
  Pc[p].Prev = pl;
  //printf("%d - %d - %d c) %d last %d coord %x cell %d\n",pp,p,pn,c,Cella[c].Last,Pc[p].Coord,Pc[p].Cell);
}
void DdLinkedList::CheckMirror(int p){
  int NCellMirror = 0;
  for(int d=0;d<3;d++){
    double Pos[3] = {Pc[p].Pos[0],Pc[p].Pos[1],Pc[p].Pos[2]};
    int cd[3] = {Pc[p].Cell[0],Pc[p].Cell[1],Pc[p].Cell[2]};
    if(Pc[p].Pos[d] - CutOff < 0.){
      Pos[d] = Pc[p].Pos[d] - CutOff;
      cd[d] -= 1;
      AddPart(p,cd,Pos);
    }
    else if(Edge[d] - Pc[p].Pos[d] < CutOff){
      Pos[d] = Pc[p].Pos[d] + CutOff;
      cd[d] += 1;
      AddPart(p,cd,Pos);
    }
  }
}
void DdLinkedList::AddMirrorPart(const int p,int *cd,double *Pos){
  SetCoorNumb(Pos,p);
  int c = pCella(cd);
  int pl = Cella[c].Last;
  SigErr(pl == p,"The particle %d is already the last in the cell %d\n",p,c);
  SigErr(p >= NAllocP,"More particle than allocated\n");
  SigErr(Cella[c].NPart >= NAllocP,"More particle than allocated\n");
  Cella[c].NPart++;
  NPart++;
  SigErr(pl >= NAllocP,"More particle than allocated\n");
  Cella[c].Last = p;
  for(int d=0;d<3;d++){
    Pc[p].Cell[d] = cd[d];
    Pc[p].Pos[d] = Pos[d];
  }
  Pc[p].Next = -1;
  if(Cella[c].NPart == 1){
    Cella[c].First = p;
    Pc[p].Prev = -2;
    return ;
  }
  if(pl >= 0) Pc[pl].Next = p;
  Pc[p].Prev = pl;
}
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
    return ;
  }
  int pp = Pc[p].Prev;
  int pn = Pc[p].Next;
  if(pp == -2) Cella[c].First = pn;
  else Pc[pp].Next = pn;
  if(pn == -1) Cella[c].Last = pp;
  else Pc[pn].Prev = pp;
  for(int d=0;d<3;d++)
    Pc[p].Cell[d] = -1;
  //printf("%d - %d - %d c) %d last %d coord %x cell %d\n",pp,p,pn,c,Cella[c].Last,Pc[p].Coord,Pc[p].Cell);
  //  assert(pp != -2 && Pc[pp].Next != pp);
}
void DdLinkedList::SwapPart(int p1,double *Pos1,int p2,double *Pos2){
  if(p1==p2) return ;
  DOMAIN_PART Pn;
  int c1 = pCella(Pc[p1].Cell);
  int c2 = pCella(Pc[p2].Cell);
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
  for(int d=0;d<3;d++){
    Pn.Cell[d]  = Pc[p1].Cell[d];
    Pc[p1].Cell[d]  = Pc[p2].Cell[d];
    Pc[p2].Cell[d]  = Pn.Cell[d];
  }
  Pc[p1].Prev  = Pc[p2].Prev;
  Pc[p1].Next  = Pc[p2].Next;
  Pc[p1].Coord = Pc[p2].Coord;
  Pc[p2].Prev  = Pn.Prev;
  Pc[p2].Next  = Pn.Next;
  Pc[p2].Coord = Pn.Coord;
}
void DdLinkedList::RemPart(const int p,double *Pos){
  int c = pCella(Pos);
  RemPart(p,c);
}
void DdLinkedList::RemPart(const int p){
  int c = pCella(Pc[p].Cell);
  SigErr(c < 0,"The cell %d does not exist\n",c);
  return RemPart(p,c);
}
void DdLinkedList::MovePart(const int p,double *NewPos){
  int cd[3];
  int cn = pCella(NewPos,cd);
  int co = pCell(p);
  SetCoorNumb(NewPos,p);
  if(cn == co) return;
  RemPart(p,co);
  AddPart(p,cd);
}
void DdLinkedList::MovePart(const int p,double *OldPos,double *NewPos){
  int cd[3];
  int cn = pCella(NewPos,cd);
  int co = pCella(OldPos);
  SetCoorNumb(NewPos,p);
  if(cn == co) return ;
  RemPart(p,co);
  AddPart(p,cd);
}
int DdLinkedList::ItCell(const int c){
  SigErr(c >= NCell,"The cell %d does not exist\n",c);
  SigErr(Cella[c].Curr1 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr1,NAllocP);
  return Cella[c].Curr1;
}
int DdLinkedList::IfItCell(const int c){
  if(Cella[c].Curr1 < 0 || Cella[c].NPart == 0){
    Cella[c].Curr1 = Cella[c].First;
    return 0;
  }
  return 1;
}
void DdLinkedList::IncrCurr(const int c){
  Cella[c].Curr1 = Pc[Cella[c].Curr1].Next;
}
void DdLinkedList::PrintCell(const int c){
  for(SetCounters(c);IfItCell(c);IncrCurr(c)){
    int p = ItCell(c);
    printf("%d) # %d %d_%d_%d %x\n",c,Cella[c].NPart,Pc[p].Prev,ItCell(c),Pc[p].Next,Pc[p].Coord);
  }
}
void DdLinkedList::PrintCells(){
  for(int c=0;c<NCell;c++)
    PrintCell(c);
}
DomCell& DomCell::operator++(){
  //It++;
  return *this;
}
DomPart& DomPart::operator++(){
  

}
void DdLinkedList::Couple(const int c,int *p1,int *p2){
  SigErr(Cella[c].Curr1 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr1,NAllocP);
  SigErr(Cella[c].Curr2 >= NAllocP,"Poiting to a particle %d over the number of allocated particles\n",Cella[c].Curr2,NAllocP);
  *p1 = Cella[c].Curr1;
  *p2 = Cella[c].Curr2;
}
int DdLinkedList::IfItCouple(const int c){
  if(Cella[c].NPart < 2) return 0;
  if(Cella[c].Curr1 == -1){
    Cella[c].Curr1 = Cella[c].First;
    Cella[c].Curr2 = Pc[Cella[c].First].Next;
    return 0;
  }
  return 1;
}
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
DomCell& DomCell::operator=(const DomCell &Dc){
  return *this;
}
void operator+(DomCell &Dc){
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
//-----------------------------NeiVertex------------------------
/**
   
 */
NeiVertex::NeiVertex(int NGridExt){
  NGrid = 2*NGridExt-1;
  int NVertex = NGrid*NGrid*NGrid;
  Vertex = new VERTEX[NVertex];
  vPos   = new PPOS[NVertex];
}
NeiVertex::~NeiVertex(){
  delete [] Vertex;
  delete [] vPos;
}
int NeiVertex::GetVertex(double *Pos){
  int s[3];
  for(int d=0;d<3;d++){
    s[d] = (int)(Pos[d]*NGrid);
  }
  return (s[0]*NGrid+s[1])*NGrid+s[2];
}
void NeiVertex::Add(double *Pos,int t){
  int v = GetVertex(Pos);
  Vertex[v].v.push_back(t);
  for(int d=0;d<3;d++)
    vPos[v].Pos[d] = Pos[d];  
}
void NeiVertex::Add(int v,int t,double *Pos){
  Vertex[v].v.push_back(t);
  for(int d=0;d<3;d++)
    vPos[v].Pos[d] = Pos[d];
}
void NeiVertex::Rem(int v,int t){
  Vertex[v].v.remove(t);
}
void NeiVertex::SetCounters(int v){
  tCurr = Vertex[v].v.begin();
}
int NeiVertex::IfItCell(int v){
  if(tCurr == Vertex[v].v.end()) return 0;
  return 1;
};
void NeiVertex::IncrCurr(int v){
  ++tCurr;
}
int NeiVertex::VertCurr(int v){
  return *tCurr;
}
void NeiVertex::PosVertex(int v,double *Pos){
  for(int d=0;d<3;d++)
    Pos[d] = vPos[v].Pos[d];
}


//-----------------------------
