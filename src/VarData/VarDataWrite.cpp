/***********************************************************************
VarData: This  Program reads and writes a specific file format 
storing all the information relative to a set of equal structure
polymers to the CHAIN, PART and GENERAL structures. It provides 
two different ways to backfold the coordinates, a fanction that 
creates an initial system with different option and some function
for the data analisys. The first calculate the distribution of the
monomer in the box, the second the distribution of the bonds.
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
#include "../include/VarData.h"
;
#ifdef USE_BOOST
void VarData::OpenComprWrite(char *FName,filtering_ostream ZipFile){
   ofstream OpenZipFile(FName, ios_base::binary);
   // create filter:
   ZipFile.push(bzip2_compressor());
   ZipFile.push(OpenZipFile);
}
void VarData::OpenComprRead(char *FName,filtering_istream ZipFile){
   ifstream OpenZipFile(FName, ios_base::binary);
   // create filter:
   ZipFile.push(bzip2_compressor());
   ZipFile.push(OpenZipFile);
}
void VarData::CloseCompr(filtering_istream ZipFile){
  ZipFile.reset();
}
void VarData::CloseCompr(filtering_ostream ZipFile){
  ZipFile.reset();
}
#endif
bool VarData::Write(char *OutFile){
  if( SysFormat == VAR_SYS_TXVL ){
    if(WriteTxvl(OutFile))
      return 1;
  }
  else if( SysFormat == VAR_SYS_XVT ){
    if(WriteXvt(OutFile))
      return 1;
  }
  else 
    WriteXyz(OutFile);
  return 0;
}
bool VarData::WriteTxvl(char *OutFile){
  VarMessage("Write");
  FILE *FileToWrite = fopen(OutFile,"w");
  SigErr(FileToWrite==NULL,"The output file %s could not be opened\n",OutFile);
  fprintf(FileToWrite,"# ");
  if( VAR_IF_TYPE(SysType,VAR_EDGE) ){
    fprintf(FileToWrite,"l( %.2g %.2g %.2g) ",
	    Gen->Edge[0],Gen->Edge[1],Gen->Edge[2]);
  }
  fprintf(FileToWrite,"c[%d] ",Gen->NChain);
  fprintf(FileToWrite,"s[%d] ",Gen->Step);
  fprintf(FileToWrite,"n[%d] ",Gen->NPart);
  fprintf(FileToWrite,"L[%d] ",Gen->NLink);
  fprintf(FileToWrite,"v[%d] ",NEdge);
  fprintf(FileToWrite,"d[%s] ",cWhat2Draw);
  fprintf(FileToWrite,"D(%lf) ",Gen->Deltat);
  fprintf(FileToWrite,"i(%.2f %.2f %.2f) ",Gen->rho,Gen->chiN,Gen->kappaN);
  fprintf(FileToWrite,"\n");
  HeaderNano(FileToWrite);
  for(int p=0;p<Gen->NPart;p++){
    fprintf(FileToWrite,"{ t[%d %d %d] ",
	    Pm[p].Idx,Pm[p].CId,Pm[p].Typ);
    fprintf(FileToWrite,"x(%lf %lf %lf) ",
	    Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2]);
    fprintf(FileToWrite,"v(%lf %lf %lf) ",
	    Pm[p].Vel[0],Pm[p].Vel[1],Pm[p].Vel[2]);
    //	printf("Pm[%d].NLink %d\n",p,Pm[p].NLink);
    if(pNLink() > 0){
      for(int i=0;i<Ln[p].NLink;i++)
	fprintf(FileToWrite,"l[%d] ",Ln[p].Link[i]);
    }
    fprintf(FileToWrite,"}\n");//End of particle info
  }
  // for(int n=0;n<Gen->NNano;n++){
  //   fprintf(FileToWrite,"{t[%d %d %d] ",Gen->NPart-1+n,-1,2);
  //   fprintf(FileToWrite,"x(%lf %lf %lf) ",Nano[n].Pos[0],Nano[n].Pos[1],Nano[n].Pos[2]);
  //   fprintf(FileToWrite,"v(%lf %lf %lf)}\n",0.,0.,0.);
  // }
  fclose(FileToWrite);
  return 0;
}
void VarData::HeaderInteraction(FILE *FileToWrite){
  double Norm2 = CUBE(pReOverCutOff()) / SQR(pNPCh());
  double Norm3 = CUBE(SQR(pReOverCutOff()) / pNPCh());
  MInt->Rescale(1./Norm2,2);
  MInt->Rescale(1./Norm3,3);
  fprintf(FileToWrite,"# v=%.2f %.2f %.2f %.2f %.2f %.2f ",MInt->Coeff(0,0),MInt->Coeff(0,1),MInt->Coeff(0,2),MInt->Coeff(1,1),MInt->Coeff(1,2),MInt->Coeff(2,2));
  fprintf(FileToWrite,"w=%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",MInt->Coeff(0,0,0),MInt->Coeff(0,0,1),MInt->Coeff(0,0,2),MInt->Coeff(0,1,1),MInt->Coeff(0,1,2),MInt->Coeff(0,2,2),MInt->Coeff(1,1,1),MInt->Coeff(1,1,2),MInt->Coeff(1,2,2),MInt->Coeff(2,2,2));
  fprintf(FileToWrite,"# a2=%.1f a3=%.1f Re=%.1f N=%d ks=%.3f kb=%.3f l0=0\n",Gen->WFuncStraight2,Gen->WFuncStraight3,Gen->ReOverCutOff,Gen->NPCh,Gen->kappaSpring,Gen->kappaBend);
  MInt->Rescale(Norm2,2);
  MInt->Rescale(Norm3,3);
}
void VarData::ShapeId(int iShape,char *Shape){
  sprintf(Shape,"no");
  if(VAR_IF_TYPE(iShape,SHAPE_SPH)) sprintf(Shape,"sph");
  else if(VAR_IF_TYPE(iShape,SHAPE_TIP)) sprintf(Shape,"tip");
  else if(VAR_IF_TYPE(iShape,SHAPE_CYL)) sprintf(Shape,"cyl");
  else if(VAR_IF_TYPE(iShape,SHAPE_PILL)) sprintf(Shape,"pill");
  else if(VAR_IF_TYPE(iShape,SHAPE_TILT)) sprintf(Shape,"tilt");
  else if(VAR_IF_TYPE(iShape,SHAPE_WALL)) sprintf(Shape,"wall"); 
  else if(VAR_IF_TYPE(iShape,SHAPE_PORE)) sprintf(Shape,"pore"); 
  else if(VAR_IF_TYPE(iShape,SHAPE_EXT)) sprintf(Shape,"ext"); 
  else if(VAR_IF_TYPE(iShape,SHAPE_JANUS)) sprintf(Shape,"janus");
  else if(VAR_IF_TYPE(iShape,SHAPE_STALK)) sprintf(Shape,"stalk");
  else if(VAR_IF_TYPE(iShape,SHAPE_TORUS)) sprintf(Shape,"torus");
  else if(VAR_IF_TYPE(iShape,SHAPE_HARM)) sprintf(Shape,"harm");
  else if(VAR_IF_TYPE(iShape,SHAPE_UMBR)) sprintf(Shape,"umbr");
  else if(VAR_IF_TYPE(iShape,SHAPE_BOUND)) sprintf(Shape,"bound");
}
void VarData::StringNano(char *String,int n){
  char Shape[60];
  ShapeId(Nano[n].Shape,Shape);
  if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CLUSTER)){ 
    sprintf(Shape,"Architecture%d.dat",n);
    sprintf(String,"# Pep x(%.2f %.2f %.2f) a(%.2f %.2f %.2f) g(%.2f %.2f %.2f) i[%d] d[%d %d] fn{%s}",Nano[n].Pos[0],Nano[n].Pos[1],Nano[n].Pos[2],Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2],Nano[n].Rad,Nano[n].Height,Nano[n].Hamaker,n,Nano[n].NCircle,Nano[n].NHeight,Shape);
  }
  else if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CLINKS)){ 
    sprintf(Shape,"CrossLinks.dat");
    sprintf(String,"# Ext i[%d] fn{%s}",0,Shape);
  }
  else 
    sprintf(String,"# Rigid x(%.2f %.2f %.2f) a(%.2f %.2f %.2f) c(%.2f %.2f %.2f %.4f) s{%s}",Nano[n].Pos[0],Nano[n].Pos[1],Nano[n].Pos[2],Nano[n].Axis[0],Nano[n].Axis[1],Nano[n].Axis[2],Nano[n].Rad,Nano[n].Hamaker,Nano[n].Height,Nano[n].OffSet,Shape);
}
int VarData::HeaderNano(FILE *FileToWrite){
  if(Gen->NNano == 0) return 1;
  for(int n=0;n<Gen->NNano;n++){
    if(Nano->Shape == SHAPE_NONE)continue;
    char String[120];
    StringNano(String,n);
    fprintf(FileToWrite,"%s\n",String);
  }
  return 0;
}
int VarData::HeaderSoft(char *Line){
  if(NSoft == 0) return 1;
  for(int n=0;n<NSoft;n++){
    char Shape[60];
    sprintf(Shape,"no");
    if(Soft[n].Topology == VAR_PLANAR) sprintf(Shape,"planar");
    if(Soft[n].Topology == VAR_VESICLE) sprintf(Shape,"vesicle");
    if(Soft[n].Topology == VAR_COATING) sprintf(Shape,"coating");
    if(Soft[n].Topology == VAR_DISTRIBUTED) sprintf(Shape,"distributed");
    sprintf(Line,"# Soft x(%.2f %.2f %.2f) c(%.2f %.2f) s{%s} n{%s}\n",Soft[n].Pos[0],Soft[n].Pos[1],Soft[n].Pos[2],Soft[n].Size[0],Soft[n].Size[1],Shape,Soft[n].Name);
  }
  return 0;
}
bool VarData::WriteXvt(char *OutFile){
  VarMessage("Write");
  FILE *FileToWrite = NULL;
  char Line2Put[512];
  FileToWrite = fopen(OutFile,"w");
  SigErr(FileToWrite==NULL,"The output file %s could not be opened\n",OutFile);
  fprintf(FileToWrite,"# L=%lf %lf %lf t=%lf blocks=%d\n",Gen->Edge[0],Gen->Edge[1],Gen->Edge[2],pTime(),Gen->NBlock);
  HeaderInteraction(FileToWrite);
  HeaderNano(FileToWrite);
  for(int b=0,cOff=0,pCurr=0;b<Gen->NBlock;b++,cOff+=Block[b].NChain){
    if(Block[b].NChain == 0) continue;
    fprintf(FileToWrite,"# n=%d N=%d name=%s\n",Block[b].NChain,Block[b].NPCh,Block[b].Name);
    for(int c=cOff;c<Block[b].NChain+cOff;c++){
      for(int p=0;p<Block[b].NPCh;p++,pCurr++){
	// if(p > 0){
	//   for(int d=0;d<3;d++){
	//     if(Pm[pCurr].Pos[d] - Pm[pCurr-1].Pos[d] > .5*pEdge(d))
	//       Pm[pCurr].Pos[d] -= pEdge(d);
	//     else if(Pm[pCurr].Pos[d] - Pm[pCurr-1].Pos[d] < -.5*pEdge(d))
	//       Pm[pCurr].Pos[d] += pEdge(d);
	//   }
	// }
	fprintf(FileToWrite,"%lf %lf %lf ",
		pPos(pCurr,0),pPos(pCurr,1),pPos(pCurr,2));
	fprintf(FileToWrite,"%lf %lf %lf ",
		Pm[pCurr].Vel[0],Pm[pCurr].Vel[1],Pm[pCurr].Vel[2]);
	fprintf(FileToWrite," %d ",
		Pm[pCurr].Typ);
	fprintf(FileToWrite,"\n");
      }
    }
  }
  fclose(FileToWrite);
  return 0;
}
bool VarData::WriteXyz(char *OutFile){
  VarMessage("Write");
  FILE *FileToWrite = NULL;
  if( (FileToWrite = fopen(OutFile,"w")) == 0){
    printf("The output file %s could not be opened\n",OutFile);
    return 0;
  }
  for(int p=0;p<Gen->NPart;p++){
    fprintf(FileToWrite,"%lf %lf %lf\n",
	    Pm[p].Pos[0],Pm[p].Pos[1],Pm[p].Pos[2]);
  }
  fclose(FileToWrite);
  return 0;
}
void VarData::WriteLinkedSurf(FILE *FWrite,double *Plot,int Values,int NType,double *Bound,int* PId){
  int p = PId[0];
  int c = PId[1];
  int t = PId[2];
  double InvValues = 1./(double)Values;
  int link[4] = {0,0,0,0};
  for(int vv=0;vv<Values-1;vv++){
    for(int v=0;v<Values-1;v++){
      link[0] = (v+0)*Values + (vv+0);
      if(Plot[(link[0])*NType+t] < .1) continue;
      //------------defines-the-squares---------------------
      link[1] = (v+0)*Values + (vv+1);
      link[2] = (v+1)*Values + (vv+0);
      link[3] = (v+1)*Values + (vv+1);
      int pRef = p;
      for(int lx=0;lx<2;lx++){
	for(int ly=0;ly<2;ly++){
	  int l = 2*lx+ly;
	  double NanoAdded = Plot[(link[l])*NType+2]/Bound[2]+Plot[(link[l])*NType+3]/Bound[3];
	  double Phob = t == 0 ? Plot[(link[l])*NType+0]/Bound[0] : 0.;
	  double Phil = t == 1 ? Plot[(link[l])*NType+1]/Bound[1] : 0.;
	  double r = (v +lx)*InvValues*pEdge(3);
	  double z = (vv+ly)*InvValues*(pEdge(CNorm));
	  double dens = (Plot[(link[l])*NType+t]);//+Plot[(v*Values+vv)*NType+1]+Plot[(v*Values+vv)*NType+2]);
	  int l1 = pRef + (p+1)%4;
	  int l2 = pRef + (p+2)%4;
	  int l3 = pRef + (p+3)%4;
	  fprintf(FWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p,c,t,r,z,dens,NanoAdded,Phob,Phil,l1,l2,l3);
	  p++;
	}
      }
      c++;
    }
  }
  PId[0] = p;
  PId[1] = c;
  PId[2] = t;
  return ;
  double Threshold = .1;
  int *IdSerial = (int *)calloc(Values*Values*2,sizeof(int));
  for(int vv=0,p=0;vv<Values-1;vv++){
    for(int v=0;v<Values-1;v++){
      int l = (vv+0)*Values + (v+0);
      IdSerial[l] = p;
      printf("%d %d %d %d\n",l,p,v,vv);
      if(Plot[l*NType+t] < Threshold) continue;
      p++;
    }
  }
  for(int vv=0,p=0,c=0;vv<Values-1;vv++){
    for(int v=0;v<Values-1;v++){
      link[0] = (v+0)*Values + (vv+0);
      if(Plot[(link[0])*NType+t] < Threshold) continue;
      //------------defines-the-squares---------------------
      link[1] = (vv+0)*Values + (v+1);
      link[2] = (vv+1)*Values + (v+0);
      link[3] = (vv+1)*Values + (v+1);
      printf("%d %d %d %d\n",link[0],link[1],link[2],link[3]);
      for(int l=0;l<4;l++){
	link[l] = IdSerial[link[l]];
      }
      printf("%d %d %d %d\n",link[0],link[1],link[2],link[3]);
      double NanoAdded = Plot[(link[0])*NType+2]/Bound[2]+Plot[(link[0])*NType+3]/Bound[3];
      double Phob = t == 0 ? Plot[(link[0])*NType+0]/Bound[0] : 0.;
      double Phil = t == 1 ? Plot[(link[0])*NType+1]/Bound[1] : 0.;
      double r = (v)*InvValues*pEdge(3);
      double z = (vv)*InvValues*(pEdge(CNorm));
      double dens = (Plot[(link[0])*NType+t]);//+Plot[(v*Values+vv)*NType+1]+Plot[(v*Values+vv)*NType+2]);
      fprintf(FWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",link[0],c,t,r,z,dens,NanoAdded,Phob,Phil,link[1],link[2],link[3]);
      p++;
      c++;
   }
  }
  free(IdSerial);
}
void VarData::WriteSurf(FILE *F2Write,double **Plot,int NSample,int OffSet){
  double InvNSample = 1./(double)NSample;
  for(int sx=0;sx<NSample;sx++){
    for(int sy=0;sy<NSample;sy++){
      double x = (sx+.5)*pEdge(CLat1)*InvNSample;
      double y = (sy+.5)*pEdge(CLat2)*InvNSample;
      int l01 = (sx+0)*NSample+(sy+1)+OffSet;
      if(sy == NSample - 1)
	l01 = (sx+0)*NSample+(0)+OffSet;
      int l10 = (sx+1)*NSample+(sy+0)+OffSet;
      if(sx == NSample - 1)
	l10 = (0)*NSample+(sy+0)+OffSet;
      int lm0 = (sx-1)*NSample+(sy+0)+OffSet;
      if(sx == 0)
	lm0 = (NSample - 1)*NSample+(sy+0)+OffSet;
      int l0m = (sx+0)*NSample+(sy-1)+OffSet;
      if(sy == 0)
	l0m = (sx)*NSample+(NSample-1)+OffSet;
      fprintf(F2Write,"{x(%lf %lf %lf) l[%d] l[%d] l[%d] l[%d]}\n",x,y,Plot[sx][sy],l01,l10,lm0,l0m);
    }
  }
}
