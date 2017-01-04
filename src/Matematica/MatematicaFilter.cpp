#include "../include/Matematica.h"

int Matematica::PermuteRandomAll(PERMUTE *Perm,int NMass){
  int *Sequence = (int *)calloc(NMass,sizeof(int));
  for(int n=0;n<NMass;n++){
    int nSquare = n;
    int Ran = (int)(NMass*Casuale());
    int IfGood=0;
    int IfAlready=1;
    for(int gg= nSquare-1;gg>=0;gg--){
      //printf("%d %d) %d\n",nSquare,gg,Sequence[gg]);
      if(nSquare == Sequence[gg]){
	IfAlready = 0;
	break;
      }
    }
    if(!IfAlready) continue;
    while(Ran==nSquare || !IfGood){// Ranw != Squarew && Ranh != Squareh){
      IfGood = 1;
      Ran = (int)(NMass*Casuale());
      for(int gg= nSquare-1;gg>=0;gg--){
	//printf("%d %d) %d %d\n",nSquare,gg,Ran,Sequence[gg]);
	if(Ran == Sequence[gg]){
	  IfGood = 0;
	  break;
	}
      }
    }
    Sequence[nSquare] = Ran;
    Sequence[Ran] = nSquare;
    //	for(int g=0;g<NGrid*NGrid;g++)printf("%d %d\n",g,Sequence[g]);
  }
  for(int n=0;n<NMass;n++){
    Perm[n].n = n;
    Perm[n].m = Sequence[n];
  }
  int nTemp = 0;
  for(int n=0;n<NMass;n++){
    if(Perm[n-nTemp].m < n ){
      for(int nn=n-nTemp;nn<NMass-1-nTemp;nn++){
	Perm[nn].n = nn+1+nTemp;
	Perm[nn].m = Sequence[nn+1+nTemp];
      }
      nTemp++;
    }
  }
  //for(int n=0;n<NMass;n++)printf("%d %d - %d %d \n",Perm[n].n,Perm[n].m,n,Sequence[n]);
  return 1;
}
int Matematica::PermuteRandomAll(int *Sequence,int NMass){
  for(int n=0;n<NMass;n++){
    int nSquare = n;
    int Ran = (int)(NMass*Casuale());
    int IfGood=0;
    int IfAlready=1;
    for(int gg= nSquare-1;gg>=0;gg--){
      //printf("%d %d) %d\n",nSquare,gg,Sequence[gg]);
      if(nSquare == Sequence[gg]){
	IfAlready = 0;
	break;
      }
    }
    if(!IfAlready) continue;
    while(Ran==nSquare || !IfGood){// Ranw != Squarew && Ranh != Squareh){
      IfGood = 1;
      Ran = (int)(NMass*Casuale());
      for(int gg= nSquare-1;gg>=0;gg--){
	//printf("%d %d) %d %d\n",nSquare,gg,Ran,Sequence[gg]);
	if(Ran == Sequence[gg]){
	  IfGood = 0;
	  break;
	}
      }
    }
    Sequence[nSquare] = Ran;
    Sequence[Ran] = nSquare;
    //if(Sequence[nSquare] < nSquare) Sequence[nSquare] = -1;
    //	for(int g=0;g<NGrid*NGrid;g++)printf("%d %d\n",g,Sequence[g]);
  }
  return 1;
}
int Matematica::ApplyFilter(Matrice *Point,Matrice *Res,Matrice *Mask){
  if(Point->Size() != Res->Size()){
    printf("Matrices differ! %d %d \n",Point->Size(),Res->Size());
    return 0;
  }
  int NMaskh = Mask->Size();
  int NMaskw = Mask->Size();
  int height = Point->Size();
  int width  = Point->Size();
  int NHalf = Mask->Size()/2;
  for(int h=0;h<height; h++) {
    for(int w=0;w<width; w++) {
      //Res->setvalue(h,w,0.);
      for(int lh=0;lh<NMaskh;lh++){
	for(int lw=0;lw<NMaskw;lw++){
	  if(h+lh-NHalf <0 ) continue;
	  if(w+lw-NHalf <0 ) continue;
	  if(h+lh-NHalf > height-1) continue;
	  if(w+lw-NHalf > width-1) continue;
	  Res->Add(h,w,Mask->Val(lh,lw)*Point->Val(h+lh-NHalf,w+lw-NHalf));
	}
      }
    }
  }
  return 1;
}
int Matematica::ApplyFilter(Matrice *Res,Matrice *Mask){
  int NMaskh = Mask->pNRow();
  int NMaskw = Mask->pNCol();
  //  double Temp[NMaskh][NMaskw];
  int width = Res->pNRow();
  int height = Res->pNCol();
  Matrice *Temp = new Matrice(width,height);
  Res->CopyOn(Temp);
  //int NHalf = Mask->Size()/2;
  int NHalf = (int)Mask->pNCol()/2;
  for(int h=0;h<height; h++) {
    for(int w=0;w<width; w++) {
      //Res->setvalue(h,w,0.);
      double dTemp=0.;
      for(int lh=0;lh<NMaskh;lh++){
	int l1h = h + lh - NHalf;
	if(l1h >= height) continue;//g1x -= NGrid;
	if(l1h < 0) continue;//g1x + NGrid;
	for(int lw=0;lw<NMaskw;lw++){
	  int l1w = w + lw - NHalf;
	  if(l1w >= width) continue;//g1x -= NGrid;
	  if(l1w < 0) continue;//g1x + NGrid;
	  dTemp += Mask->Val(lw,lh)*Temp->Val(l1w,l1h);
	}
      }
      Res->Set(w,h,dTemp);
    }
  }
  delete Temp;
  return 0;
}
int Matematica::Transform(int *Out,int *In,int NEdge,int operation){
  int NMin = NEdge-1;
  if( OP_IF(operation,OP_INVERT) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[r*NEdge+c] == 1 ? 0 : 1;
  }
  if( OP_IF(operation,OP_ROT_90) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[(NMin-c)*NEdge+r];
  }
  if( OP_IF(operation,OP_ROT_180) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[(NMin-r)*NEdge+(NMin-c)];
  }
  if( OP_IF(operation,OP_ROT_270) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[c*NEdge+(NMin-r)];
  }
  if( OP_IF(operation,OP_MIRROR) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[r*NEdge+(NMin-c)];
  }
  if( OP_IF(operation,OP_TRANSPOSE) ){
    for(int r=0;r<NEdge;r++)
      for(int c=0;c<NEdge;c++)
	Out[r*NEdge+c] = In[c*NEdge+r];
  }
  return 1;
}
void Matematica::BackFold(Matrice *In,Matrice *Out,int NShift){
  for(int r=0;r<In->pNRow();r++){
    int r1 = r + NShift;
    if(r1 >= In->pNRow()) r1 -= In->pNRow();
    if(r1 < 0) r1 += In->pNRow();
    for(int c=0;c<In->pNCol();c++){
      Out->Set(r1,c,In->Val(r,c));
    }
  }
}
