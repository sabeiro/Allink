#include <VarData.h>
int VarData::SetNPart(int NewNPart){
  if(Gen->NPart == NewNPart) return 1;
  int OldNPart = pNPart();
  Gen->NPart = NewNPart;
  if( !VAR_IF_TYPE(SysType,VAR_PART_ALLOCATED) ){
    Gen->NAllocP = Gen->NPart;
    if(pNPCh() > 0 && pNPCh() < 200 ) Gen->NAllocP += 2000*Gen->NPCh;
    else Gen->NAllocP += 2000;
    Gen->NAllocP += 2000;
    Pm = (PART *)calloc(Gen->NAllocP,sizeof(PART));
    if(Pm == NULL){
      printf("Couldn't alloc the particles\n");
      exit(1);
    }
    VAR_ADD_TYPE(SysType,VAR_PART_ALLOCATED);
    if(pNLink() > 0){
      AllocLinks(NewNPart);
    }
  }
  if(Gen->NAllocP < NewNPart){
    Gen->NAllocP = NewNPart;
    if(pNPCh() > 0 && pNPCh() < 200 ) Gen->NAllocP += 2000*Gen->NPCh;
    else Gen->NAllocP += 2000;
    PART *Pn = (PART *)calloc(Gen->NPart,sizeof(PART));
    Copy(Pn,Pm,OldNPart);
    Pm = (PART *)realloc(Pm,Gen->NAllocP*sizeof(PART));
    if(Pm == NULL){
      printf("Couldn't realloc the particles\n");
      exit(1);
    }
    if(pNLink() > 0){
      AllocLinks(NewNPart);
    }
    Copy(Pm,Pn,OldNPart);
    free(Pn);
  }
  else return 1;
  return 0;
}
int VarData::SetNLink(int NewNLink){
  if(Gen->NLink == NewNLink) return 1;
  Gen->NLink = NewNLink;
  return 0;
}
int VarData::AllocLinks(int NewNPart){
  if( !VAR_IF_TYPE(SysType,VAR_LINK_ALLOCATED) ){
    Ln = (LINKS *)calloc(Gen->NAllocP,sizeof(LINKS));
    if(Ln == NULL){
      printf("Couldn't alloc the links\n");
      exit(1);
    }
    for(int p=0;p<Gen->NAllocP;p++){
      Ln[p].Link = (int *)calloc(Gen->NLink,sizeof(int));
      if(Ln[p].Link == NULL){
	printf("Couldn't alloc the links\n");
	exit(1);
      }
    }
    VAR_ADD_TYPE(SysType,VAR_LINK_ALLOCATED);
  }
  else if(Gen->NAllocP < NewNPart){
    Ln = (LINKS *)realloc(Ln,Gen->NAllocP*sizeof(LINKS));
    if(Ln == NULL){
      printf("Couldn't realloc the links\n");
      exit(1);
    }
    for(int p=0;p<Gen->NAllocP;p++){
      Ln[p].Link = (int *)realloc(Ln[p].Link,Gen->NLink*sizeof(int));
      if(Ln[p].Link == NULL){
	printf("Couldn't realloc the links\n");
	exit(1);
      }
    }
  }
  return 0;
}
int VarData::SetNChain(int NewNChain){
  int OldNChain = pNChain();
  Gen->NChain = NewNChain;
  if( !VAR_IF_TYPE(SysType,VAR_CH_ALLOCATED) ){
    Gen->NAllocC = Gen->NChain+200;
    Ch = (CHAIN *)calloc(Gen->NAllocC,sizeof(CHAIN));
    VAR_ADD_TYPE(SysType,VAR_CH_ALLOCATED);
  }
  if(Gen->NAllocC < NewNChain){
    Gen->NAllocC = NewNChain + 200;
    CHAIN *Cn = (CHAIN *)calloc(Gen->NChain,sizeof(CHAIN));
    Copy(Cn,Ch,OldNChain);
    Ch = (CHAIN *)realloc(Ch,Gen->NAllocC*sizeof(CHAIN));
    if(Ch == NULL){
      printf("Couldn't realloc the chains\n");
      exit(1);
    }
    Copy(Ch,Cn,OldNChain);
    free(Cn);
  }
  else return 1;
  return 0;
}
//obsolete?
void VarData::AllocPart(){
  if(Gen->NPart == 0){
    printf("No particles or no particle tag: {} \n");
    return ;
  }
  if( !VAR_IF_TYPE(SysType,VAR_PART_ALLOCATED) ){
    Gen->NAllocC = Gen->NChain+100;
    Gen->NAllocP = Gen->NPart+100*Gen->NPCh;
    Pm = (PART *)calloc(Gen->NAllocP,sizeof(PART));
    Ln = (LINKS *)calloc(Gen->NAllocP,sizeof(LINKS));
    if(Ln == NULL || Pm == NULL){
      return ;
    }
    for(int p=0;p<Gen->NAllocP;p++){
      Ln[p].Link = (int *)calloc(Gen->NLink,sizeof(int));
    }
    Ch = (CHAIN *)calloc(Gen->NAllocC,sizeof(CHAIN));
    VAR_ADD_TYPE(SysType,VAR_PART_ALLOCATED);
    VAR_ADD_TYPE(SysType,VAR_CH_ALLOCATED);
  }
  //FIXME: doesn't realloc
  else {
    //free(Ch);
    //for(int p=0;p<Gen->NPart;p++)  free(Pm[p].Link);
    //free(Pm);
    // Pm = (PART *)realloc(Pm,Gen->NPart*sizeof(PART));
    //Pm = (PART *)calloc(Gen->NPart,sizeof(PART));
    // for(int p=0;p<Gen->NPart;p++){
    //   Pm[p].Link = (int *)realloc(Pm[p].Link,Gen->NLink*sizeof(int));
    //   //Pm[p].Link = (int *)calloc(Gen->NLink,sizeof(int));
    // }
    // Ch = (CHAIN *)realloc(Ch,Gen->NChain*sizeof(CHAIN));
    //Ch = (CHAIN *)calloc(Gen->NChain,sizeof(CHAIN));
  }
  if(!Pm){ printf("Non alloca Pm\n");return ;}
  if(!Ch){ printf("Non alloca Ch\n");return ;}
}
void VarData::Copy(PART *P1,PART *P2,int NPartOld){
  for(int p=0;p<pNPart();p++){
    if(p >= NPartOld) break;
    for(int d=0;d<3;d++){
      P1[p].Pos[d] = P2[p].Pos[d];
      P1[p].Vel[d] = P2[p].Vel[d];
    }
    P1[p].Typ = P2[p].Typ;
  }
}
void VarData::Copy(CHAIN *C1,CHAIN *C2,int NChainOld){
  for(int c=0;c<pNChain();c++){
    if(c >= NChainOld) break;
    for(int d=0;d<3;d++){
      C1[c].Pos[d] = C2[c].Pos[d];
      C1[c].Vel[d] = C2[c].Vel[d];
    }
    C1[c].Type = C2[c].Type;
  }
}
void VarData::SetNPCh(int NewNPCh){
  //blocks?
  Gen->NPCh = NewNPCh;
}
int VarData::SetNBlock(int Val){
  if(Val == Gen->NBlock) return 1;
  if(!VAR_IF_TYPE(SysType,VAR_NBLOCK_ALL)){
    Gen->NBlock = Val;
    Block = (BLOCK *)calloc(Gen->NBlock,sizeof(BLOCK));
    VAR_ADD_TYPE(SysType,VAR_NBLOCK_ALL);
  }
  else{
    Gen->NBlock = Val;
    Block = (BLOCK *) realloc(Block,Gen->NBlock*sizeof(BLOCK));
  }
  return 0;
};
int VarData::SetNNano(int Val){
  if(Val == Gen->NNano) return 1;
  //if(!VAR_IF_TYPE(SysType,VAR_NNANO_ALL)){
  //  VAR_ADD_TYPE(SysType,VAR_NBLOCK_ALL);
  Gen->NNano = Val;
  if(Val == 0) Val = 1;
  Nano = (NANO *) realloc(Nano,Val*sizeof(NANO));
  return 0;
};
void VarData::SetNType(int NewNType){
  //blocks?
  Gen->NType = NewNType;
}
double VarData::pChPos(int c,int d){
  return Ch[c].Pos[d];
};
double VarData::pPos(int p,int d){
  return Pm[p].Pos[d] + Pm[p].Bkf[d];
};
void VarData::pPos(int p,double *Pos){
  for(int d=0;d<3;d++) 
    Pos[d] = Pm[p].Pos[d] + Pm[p].Bkf[d];
};
double VarData::pPosNoBkf(int p,int d){
  return Pm[p].Pos[d];
};
double VarData::pVel(int p,int d){
  return Pm[p].Vel[d];
};
void VarData::SetPos(int p, double *Pos){
  for(int d=0;d<3;d++)
    Pm[p].Pos[d] = Pos[d];
};
void VarData::SetPos(int p,int d,double Pos){
  Pm[p].Pos[d] = Pos;
};
void VarData::SetVel(int p, double *Vel){
  for(int d=0;d<3;d++)
    Pm[p].Vel[d] = Vel[d];
};
void VarData::SetType(int p,int t){
  Pm[p].Typ = t;
};
int VarData::pType(int p){
  return Pm[p].Typ;
};
int VarData::pChain(int p){
  return Pm[p].CId;
};
double VarData::pNanoPos(int n,int d){
  return Nano[n].Pos[d] + Nano[n].Bkf[d];
};
void VarData::SetBkf(int p){
  for(int d=0;d<3;d++) 
    Pm[p].Bkf[d] = -floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
};
void VarData::SetNanoBkf(int n){
  for(int d=0;d<3;d++) 
    Nano[n].Bkf[d] = -floor(Nano[n].Pos[d]*pInvEdge(d))*pEdge(d);
};
// void VarData::pPos(int p){
//   printf("%lf %lf %lf\n",Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2]);
// };
double *VarData::pPos(int p){
  return Pm[p].Pos;
}
void VarData::pPos(double *Pos){
  printf("%lf %lf %lf\n",Pos[0],Pos[1],Pos[2]);
};
int VarData::pStep(){return Gen->Step;};
int VarData::pNPart(){return Gen->NPart;};
int VarData::pNChain(){return Gen->NChain;};
int VarData::pNChain(int b){return Block[b].NChain;};
int VarData::pNPCh(){return Gen->NPCh;};
int VarData::pNPCh(int b){return Block[b].NPCh;};
int VarData::pNType(){return Gen->NType;};
int VarData::pNLink(){return Gen->NLink;};
int VarData::pNNano(){return Gen->NNano;};
int VarData::pNBlock(){return Gen->NBlock;};
int VarData::pNAllocP(){return Gen->NAllocP;};
int VarData::pNAllocC(){return Gen->NAllocC;};
