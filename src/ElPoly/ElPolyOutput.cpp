#include "ElPoly.h"

void ElPoly::Conv2Tecplot(int NBin,int How){
  int NType = 3;
  double Round = 0.001;
  int Nx = MIN(NBin,NEdge);
  int Ny = MIN(NBin,NEdge);
  double **Plot = (double **)calloc(NType,sizeof(double));
  double **Count = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot[t] = (double *)calloc(Nx*Ny,sizeof(double));
    Count[t] = (double *)calloc(Nx*Ny,sizeof(double));    
  }
  FILE *TecPlot = fopen("TecPlot.dat","w");
  if(How == 0){//Dens
    fprintf(TecPlot,"VARIABLES = \"R\", \"Z\", \"oil density\",\"density phob\",\"density phil\"\n");
    fprintf(TecPlot,"ZONE J=%d, K=%d, F=POINT\n",Nx,Ny);
  }
  else if(How == 1){
    fprintf(TecPlot,"VARIABLES = \"R\", \"Z\", \"pressure\",\"density phob\",\"density phil\"\n");
    fprintf(TecPlot,"ZONE J=%d, K=%d, F=POINT\n",Nx,Ny);
  }
  else if(How == 2){
    fprintf(TecPlot,"VARIABLES = \"R\", \"Z\", \"height\",\"density phob\",\"density phil\"\n");
    fprintf(TecPlot,"ZONE J=%d, K=%d, F=POINT\n",Nx,Ny);
  }
  for(int p=0;p<pNPart();p++){
    int vx = (int)((pPos(p,0)+Round)*pInvEdge(0)*Nx);
    if(vx < 0 || vx >= Nx) continue;
    int vy = (int)((pPos(p,1)+Round)*pInvEdge(1)*Ny);
    if(vy < 0 || vy >= Ny) continue;
    if(How == 0){//Dens
      int t = Pm[p].Typ;
      t = (t+1)%3;
      // if(t==2){t=0;Plot[1][vx*Nx+vy] += -1.;}
      // else 
	Plot[t][vx*Nx+vy] += pPos(p,2);
      Count[t][vx*Nx+vy] += 1.;
    }
    else if(How == 1){
      for(int t=0;t<3;t++){
	Plot[t][vx*Nx+vy] += pVel(p,t);
	Count[t][vx*Nx+vy] += 1.;
      }
    }
    else if(How == 2){
      int t = 0;//Pm[p].Typ;
      Plot[t][vx*Nx+vy] += pPos(p,2);
      Count[t][vx*Nx+vy] += 1.;
    }
  }
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<3;t++){
    //Mask.ConvoluteMatrix(Plot[t],Nx,2);
    Mask.ConvoluteMatrix(Plot[t],Nx,NDim,IfMinImConv);
  }
  for(int vx=0;vx<Nx;vx++){
    for(int vy=0;vy<Ny;vy++){
      double r = vx*pEdge(0)/(double)Nx;
      double z = vy*pEdge(1)/(double)Ny - pEdge(1)*.5;
      double Norm0 = Count[0][vx*Nx+vy] > 0. ? 1./Count[0][vx*Nx+vy] : 1.;
      double Norm1 = Count[1][vx*Nx+vy] > 0. ? 1./Count[1][vx*Nx+vy] : 1.;
      double Norm2 = Count[2][vx*Nx+vy] > 0. ? 1./Count[2][vx*Nx+vy] : 1.;
      fprintf(TecPlot,"%lf %lf %lf %lf %lf\n",r,z,
	      Plot[0][vx*Nx+vy]*Norm0,Plot[1][vx*Nx+vy]*Norm1,Plot[2][vx*Nx+vy]*Norm2);
    }
  }
  // if(Pm[p].Vel[0]+Pm[p].Vel[1]+Pm[p].Vel[2]<.1){
  //   Dens1 = -5.;
  //   Dens2 = -5.;
  //   Dens3 = -5.;
  // }
  // if(pPos(p,0) > 0 && pPos(p,0) < 5.)
  //   if(pPos(p,1)>13. && pPos(p,1)<19.)
  // 	if(Pm[p].Vel[0] < 0.01){
  // 	  Dens1 = -3.;
  // 	  Dens2 = -3.;
  // 	  Dens3 = -3.;
  // 	}
  fclose(TecPlot);
  for(int t=0;t<NType;t++){
    free(Plot[t]);
    free(Count[t]);
  }
  free(Plot);
  free(Count);
}
void ElPoly::Conv2Vmd(){
  FILE *OutVmd = fopen("Sim.vtf","w");
  int cOff = 0;
  int pOff = 0;
  char Shape[30];
  int NanoPoint = 10;
  for(int b=0;b<pNBlock();b++){
    int bType = 0;
    if(!strncmp(Block[b].Name,"PEP",3)){
      bType = 2;
    }    
    for(int c=cOff;c<cOff+pNChain(b);c++){
      for(int p=pOff;p<pOff+pNPCh(b);p++){
	int Type = pType(p) + bType;
	if(Type == 0)
	  fprintf(OutVmd,"a %d r 0.8 n A resid %d res %s\n",p,pChain(p),Block[b].Name);
	else if(Type == 1)
	  fprintf(OutVmd,"a %d r 0.8 n B resid %d res %s\n",p,pChain(p),Block[b].Name);
	else if(Type == 2)
	  fprintf(OutVmd,"a %d r 0.8 n D resid %d res %s\n",p,pChain(p),Block[b].Name);
	else if(Type == 3)
	  fprintf(OutVmd,"a %d r 0.8 n E resid %d res %s\n",p,pChain(p),Block[b].Name);
      }
      pOff += pNPCh(b);
      fprintf(OutVmd,"b %d::%d\n",c*pNPCh(b),(c+1)*pNPCh(b)-1);
    }
    cOff += pNChain(b);
  }
  for(int n=0;n<pNNano();n++){
    ShapeId(Nano[n].Shape,Shape);
    for(int i=0;i<NanoPoint;i++){
      fprintf(OutVmd,"a %d r %lf n Nano resid %d res %s\n",n*NanoPoint+pNPart()+i,Nano[n].Rad,cOff+n,Shape);
    }
    fprintf(OutVmd,"b %d::%d\n",n*NanoPoint+pNPart(),n*NanoPoint+pNPart()+NanoPoint-1);
  }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_CHAIN)) return ;
    fprintf(OutVmd,"\ntimestep\npbc %lf %lf %lf\n",pEdge(0),pEdge(1),pEdge(2));
    for(int p=0;p<pNPart();p++){
      fprintf(OutVmd,"%lf %lf %lf\n",pPos(p,0),pPos(p,1),pPos(p,2));
    }
    for(int n=0;n<pNNano();n++){
      double Pos[3];
      for(int d=0;d<3;d++){
	Pos[d] = pNanoPos(n,d) - .5*Nano[n].Height*Nano[n].Axis[d];
      }
      for(int i=0;i<NanoPoint;i++){
	for(int d=0;d<3;d++){
	  Pos[d] += Nano[n].Height/(double)NanoPoint*Nano[n].Axis[d];
	}
	fprintf(OutVmd,"%lf %lf %lf\n",Pos[0],Pos[1],Pos[2]);
      }
    }
  }
}
//-----------------------POVRAY--------------------------
void ElPoly::DrPartPovRay(int p){
  double Pos[3];
  for(int d=0;d<3;d++) Pos[d] = pPos(p,d)*InvScaleUn;
  int Typ = pType(p) < 6 ? pType(p) : 5;
  int Chc = Ch[Pm[p].CId].Type;
  if(pType(p) == 0 && CHAIN_IF_TYPE(Chc,CHAIN_ADDED) )Typ = 3;
  fprintf(DrawOutFile,"ElSphere(<%.4f, %.4f, %.4f>,",Pos[0],Pos[1],Pos[2]);
  fprintf(DrawOutFile,"SphRad,rgbt<%.4f, %.4f, %.4f, %.4f>)\n",ColorType[Typ][0],ColorType[Typ][1],ColorType[Typ][2],ColorType[Typ][3]-1.);
}
void ElPoly::DrNanoPovRay(int n){
  if(VAR_IF_TYPE(Nano[n].Shape,SHAPE_CYL)){
    double PosP[3];
    double PosN[3];
    int Typ = 2;
    for(int d=0;d<3;d++){
      PosP[d] = pNanoPos(n,d) + .5*Nano[n].Height*Nano[n].Axis[d];
      PosN[d] = pNanoPos(n,d) - .5*Nano[n].Height*Nano[n].Axis[d];
      PosP[d] *= InvScaleUn;
      PosN[d] *= InvScaleUn;
    }
    fprintf(DrawOutFile,"ElCylinder(<%.2f, %.2f, %.2f>,",PosP[0],PosP[1],PosP[2]);
    fprintf(DrawOutFile,"<%.2f, %.2f, %.2f>%.2f,",PosN[0],PosN[1],PosN[2],Nano[n].Rad*InvScaleUn);
    fprintf(DrawOutFile,"rgbt<%.4f, %.4f, %.4f, %.4f>,1)\n",ColorType[Typ][0],ColorType[Typ][1],ColorType[Typ][2],ColorType[Typ][3]-1.);
  }
  else{
    // Point2Shape(Nano[n].Shape);
    // DrField(128,SQR(Nano[n].Rad),n,DrawOutFile);
  }
}
void ElPoly::DrBondPovRay(double *Pos1,double *Pos2,float *Color){
  fprintf(DrawOutFile,"ElCylinder(<%lf, %lf, %lf>,",Pos1[0],Pos1[1],Pos1[2]);
  fprintf(DrawOutFile,"<%lf, %lf, %lf>",Pos2[0],Pos2[1],Pos2[2]);
  fprintf(DrawOutFile,"CylRad,rgbt<%.4f, %.4f, %.4f, %.4f>,1)\n",Color[0],Color[1],Color[2],Color[3]-1.);
}
void ElPoly::Conv2Povray(){
  //particles, lines
  HeaderPovRay();
  for(int f=NFile[0];f<NFile[1];f++){
    SigErr(DrawOutFile != NULL,"DrawOutFile already allocated, can't use the file");
    char FName[60];
    sprintf(FName,"PovSnap%05d.pov",pStep());
    DrawOutFile = fopen(FName,"w");
    int ImSize[2] = {1000,1000};
    fprintf(DrawOutFile,"// POV 3.x input script : plot.pov\n// command: povray +W%d +H%d -I%s -O%s.tga +P +X +A +FT +C\n",ImSize[0],ImSize[1],FName,FName);
    fprintf(DrawOutFile,"#if (version < 3.5)\n#error \"POV3DisplayDevice has been compiled for POV-Ray 3.5 or above.\\nPlease upgrade POV-Ray.\"\n#end\n");
    fprintf(DrawOutFile,"#include \"PovHeader.inc\"\n");
    fprintf(DrawOutFile,"// System \n");
    for(int b=0,NPep=0;b<pNBlock();b++){
      for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx;p++){
	DrPartPovRay(p);
	int Typ = Pm[p].Typ;
	double Pos1[3];
	double Pos2[3];
	for(int l=0;l<Ln[p].NLink;l++){
	  int link = Ln[p].Link[l];
	  if(p == link) continue;
	  for(int d=0;d<3;d++){
	    Pos1[d] = pPos(p,d)*InvScaleUn;
	    Pos2[d] = pPos(p,link)*InvScaleUn;
	    if(Pos1[d] - Pos2[d] > .5*InvScaleUn)
	      Pos2[d] += pEdge(d)*InvScaleUn;
	    else if(Pos1[d] - Pos2[d] < -.5*InvScaleUn)
	      Pos2[d] -= pEdge(d)*InvScaleUn;
	  }
	  DrBondPovRay(Pos1,Pos2,ColorType[Typ]);
	}
      }
    }
    for(int n=0;n<pNNano();n++){
      DrNanoPovRay(n);
    }
    fclose(DrawOutFile);
  }
}
#include "ElPolyDrawSurf.h"
void ElPoly::DrField(int NGrid,double IsoLevel,int nNano,FILE *FWrite){
  int Typ = 2;
  double CubeDist[8];
  double EdgeVertex[12][3];
  double EdgeNormal[12][3];
  double InvNGrid = 1./(double)NGrid;
  double Pos[3];
  double Pos1[3];
  double Cm[3];
  for(int d=0;d<3;d++){
    Cm[d] = .5*pEdge(d)*InvScaleUn;
  }
  for(int gx=0;gx<NGrid;gx++){
    Pos[0] = gx*InvNGrid*pEdge(0);
    for(int gy=0;gy<NGrid;gy++){
      Pos[1] = gy*InvNGrid*pEdge(1);
      for(int gz=0;gz<NGrid;gz++){
	Pos[2] = gz*InvNGrid*pEdge(2);
	for(int v=0;v<8;v++){
	  for(int d=0;d<3;d++){
	    Pos1[d] = Pos[d] + VertCube[v][d]*InvNGrid*pEdge(d);
	  }
	  CubeDist[v] = NanoDist2(Pos1,nNano);
	}
	int Flag = 0;
	for(int v=0;v<8;v++){
	  if(CubeDist[v] <= IsoLevel)
	    Flag |= 1<<v;
	}
	int CubeShape = CubeTop[Flag];
	if(CubeShape==0) continue;
	for(int e=0;e<12;e++){
	  if(CubeShape & (1<<e)){
	    double Delta = CubeDist[EdgeConn[e][1]] - CubeDist[EdgeConn[e][0]];
	    double OffSet = (IsoLevel-CubeDist[EdgeConn[e][0]])/Delta;
	    if(Delta == 0.0){
	      OffSet = .5;
	    }
	    EdgeNormal[e][0] = NanoDist2(Pos[0]-0.01,Pos[1],Pos[2],nNano) 
	      - NanoDist2(Pos[0]+0.01,Pos[1],Pos[2],nNano);
	    EdgeNormal[e][1] = NanoDist2(Pos[0],Pos[1]-0.01,Pos[2],nNano) 
	      - NanoDist2(Pos[0],Pos[1]+0.01,Pos[2],nNano);
	    EdgeNormal[e][2] = NanoDist2(Pos[0],Pos[1],Pos[2]-0.01,nNano) 
	      - NanoDist2(Pos[0],Pos[1],Pos[2]+0.01,nNano);
	    double Norm = sqrt(SQR(EdgeNormal[e][0])+SQR(EdgeNormal[e][1])+SQR(EdgeNormal[e][2]));
	    for(int d=0;d<3;d++){
	      EdgeVertex[e][d] = Pos[d] + (VertCube[EdgeConn[e][0]][d]+OffSet*EdgeDir[e][d])*InvNGrid*pEdge(d);
	      EdgeVertex[e][d] *= InvScaleUn;
	      EdgeNormal[e][d] *= 1./Norm;
	    }
	  }
	}
	for(int t=0;t<5;t++){
	  if(TrConnTable[Flag][3*t] < 0.) break;
	  fprintf(FWrite,"ElTriangle (");
	  for(int d=0;d<3;d++){
	    int v = TrConnTable[Flag][3*t+d];
	    fprintf(FWrite,"<%.2f, %.2f, %.2f>,",EdgeVertex[v][0]-Cm[0],EdgeVertex[v][1]-Cm[1],EdgeVertex[v][2]-Cm[2]);
	    fprintf(FWrite,"<%.2f, %.2f, %.2f>,",EdgeNormal[v][0],EdgeNormal[v][1],EdgeNormal[v][2]);
	  }
	  fprintf(FWrite,"rgbt<%.3f, %.3f, %.3f, %.3f>)\n",ColorType[Typ][0],ColorType[Typ][1],ColorType[Typ][2],ColorType[Typ][3]-1.);
	}
      }
    }
  }
}
void ElPoly::Conv2rzd(int NBin){
  // FILE *FOut = fopen("data.dat","w");
  // for(int p=0;p<pNPart();p++){
  //   fprintf(FOut,"%lf %lf %lf %lf\n",pPos(p,0),pPos(p,1),pPos(p,2),pVel(p,0));
  // }
  // fclose(FOut);
  // return;
  int NType = 3;
  double Round = 0.001;
  double **Plot = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++)
    Plot[t] = (double *)calloc(NBin*NBin,sizeof(double));
  double *Count = (double *)calloc(NBin*NBin*NType,sizeof(double));
  FILE *TecPlot = fopen("data.dat","w");
  for(int p=0;p<pNPart();p++){
    int vx = (int)((pPos(p,0)+Round)*pInvEdge(0)*NBin);
    if(vx < 0 || vx >= NBin) continue;
    int vy = (int)((pPos(p,1)+Round)*pInvEdge(1)*NBin);
    if(vy < 0 || vy >= NBin) continue;
    for(int t=0;t<3;t++){
      Plot[t][vx*NBin+vy] += pVel(p,t);
      Count[(vx*NBin+vy)*NType+t] += 1.;
    }
  }
  //smooth and write dens
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<3;t++){
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NBin,NDim,IfMinImConv);
  }
  for(int vx=0;vx<NBin;vx++){
    for(int vy=0;vy<NBin;vy++){
      double r = vx*pEdge(0)/(double)NBin;
      double z = vy*pEdge(1)/(double)NBin - pEdge(1)*.5;
      double Norm0 = Count[(vx*NBin+vy)*NType+0] > 0. ? 1./Count[(vx*NBin+vy)*NType+0] : 1.;
      double Norm1 = Count[(vx*NBin+vy)*NType+1] > 0. ? 1./Count[(vx*NBin+vy)*NType+1] : 1.;
      double Norm2 = Count[(vx*NBin+vy)*NType+2] > 0. ? 1./Count[(vx*NBin+vy)*NType+2] : 1.;
      fprintf(TecPlot,"%lf %lf %lf %lf %lf\n",r,z,
	      Plot[0][vx*NBin+vy]*Norm0,Plot[1][vx*NBin+vy]*Norm1,Plot[2][vx*NBin+vy]*Norm2);
      // fprintf(TecPlot,"%lf %lf %lf \n",r,z,Plot[1][vx*NBin+vy]*Norm0);
    }
  }
  fclose(TecPlot);
  for(int t=0;t<NType;t++)
    free(Plot[t]);
  free(Plot);
  free(Count);
}
void ElPoly::Conv2xyzd(int NBin){
  // FILE *FOut = fopen("data.dat","w");
  // for(int p=0;p<pNPart();p++){
  //   fprintf(FOut,"%lf %lf %lf %lf\n",pPos(p,0),pPos(p,1),pPos(p,2),pVel(p,0));
  // }
  // fclose(FOut);
  // return;
  int NType = 3;
  double Round = 0.001;
  double **Plot = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++)
    Plot[t] = (double *)calloc(NBin*NBin*NBin,sizeof(double));
  double *Count = (double *)calloc(NBin*NBin*NBin*NType,sizeof(double));
  FILE *TecPlot = fopen("data.dat","w");
  for(int p=0;p<pNPart();p++){
    int vx = (int)((pPos(p,0)+Round)*pInvEdge(0)*NBin);
    if(vx < 0 || vx >= NBin) continue;
    int vy = (int)((pPos(p,1)+Round)*pInvEdge(1)*NBin);
    if(vy < 0 || vy >= NBin) continue;
    int vz = (int)((pPos(p,2)+Round)*pInvEdge(2)*NBin);
    if(vz < 0 || vz >= NBin) continue;
    for(int t=0;t<3;t++){
      Plot[t][(vx*NBin+vy)*NBin+vz] += pVel(p,t);
      Count[((vx*NBin+vy)*NBin+vz)*NType+t] += 1.;
    }
  }
  //smooth and write dens
  for(int vx=0;vx<NBin;vx++){
    for(int vy=0;vy<NBin;vy++){
      for(int vz=0;vz<NBin;vz++){
	double x = vx*pEdge(0)/(double)NBin;
	double y = vy*pEdge(1)/(double)NBin;
	double z = vz*pEdge(1)/(double)NBin;
	double Norm0 = Count[((vx*NBin+vy)*NBin+vz)*NType+0] > 0. ? 1./Count[((vx*NBin+vy)*NBin+vz)*NType+0] : 1.;
	double Norm1 = Count[((vx*NBin+vy)*NBin+vz)*NType+1] > 0. ? 1./Count[((vx*NBin+vy)*NBin+vz)*NType+1] : 1.;
	double Norm2 = Count[((vx*NBin+vy)*NBin+vz)*NType+2] > 0. ? 1./Count[((vx*NBin+vy)*NBin+vz)*NType+2] : 1.;
	fprintf(TecPlot,"%lf %lf %lf %lf %lf %lf\n",x,y,z,
		Plot[0][(vx*NBin+vy)*NBin+vz]*Norm0,Plot[1][(vx*NBin+vy)*NBin+vz]*Norm1,Plot[2][(vx*NBin+vy)*NBin+vz]*Norm2);
	// fprintf(TecPlot,"%lf %lf %lf %lf\n",x,y,z,Plot[0][vx*NBin+vy]*Norm0);
      }
    }
  }
  fclose(TecPlot);
  for(int t=0;t<NType;t++)
    free(Plot[t]);
  free(Plot);
  free(Count);
}
void ElPoly::HeaderPovRay(){
  FILE *HeaderPov = fopen("PovHeader.inc","w");
  double SphRad = 0.005;
  double CylRad = 0.002;
  double PerspAngle = 15;
  int ImSize[2] = {1000,1000};
  fprintf(HeaderPov,"// POV 3.x input script : plot.pov\n// command: povray +W%d +H%d -ISnap.pov -OSnap.pov.tga +P +X +A +FT +C\n",ImSize[0],ImSize[1]);
  fprintf(HeaderPov,"#if (version < 3.5)\n#error \"POV3DisplayDevice has been compiled for POV-Ray 3.5 or above.\\nPlease upgrade POV-Ray.\"\n#end\n");
  //include
    fprintf(HeaderPov,"// #include \"colors.inc\"\n	\
//  #include \"stones.inc\"\n\
//  #include \"textures.inc\"\n\
//  #include \"shapes.inc\"\n\
//  #include \"glass.inc\"\n\
//  #include \"metals.inc\"\n\
//  #include \"woods.inc\"\n");
  fprintf(HeaderPov,"\
#declare clip_on=array[3] {0, 0, 0};\n\
#declare clip=array[3];\n\
#declare scaledclip=array[3];\n\
#declare SphRad=%lf;\n\
#declare CylRad=%lf;\n\
#declare line_width=0.0020;\n",SphRad,CylRad);
  // //materials
  // fprintf(HeaderPov,"\
  // #declare Fotoni  = photons {\n\
  //       target\n\
  //       refraction on\n\
  //       reflection on\n\
  //   }\n");
  //define a point
  fprintf(HeaderPov,"\
#macro ElPoint (P1, R1, C1)\n\
  #local T = texture { finish { ambient 1.0 diffuse 0.0 phong 0.0 specular 0.0 } pigment { C1 } }\n \
  #if(clip_on[2])\n\
  intersection {\n\
    sphere {P1, R1 texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
    clip[2]\n\
  }\n\
  #else\n\
  sphere {P1, R1 texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
  #end\n\
#end\n");
  //define a line
  fprintf(HeaderPov,"\
#macro ElLine (P1, P2, C1)\n\
  #local T = texture { finish { ambient 1.0 diffuse 0.0 phong 0.0 specular 0.0 } pigment { C1 } }\n\
  #if(clip_on[2])\n\
  intersection {\n\
    cylinder {P1, P2, line_width texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
    clip[2]\n\
  }\n\
  #else\n\
  cylinder {P1, P2, line_width texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
  #end\n\
#end\n");
  //define a sphere
  fprintf(HeaderPov,"#macro ElSphere(P1, R1, C1)\n\
  #local T = texture { pigment { C1 } }\n\
  #local M = material{\n		 \
    texture{\n\
      pigment{gradient y color_map{[0.4 C1][0.4 C1]}}\n\
      finish{ambient 0 diffuse 0.4 specular 1 roughness 0.0001 reflection 0.25}\n\
    }\n\
    interior{ior 1.33}\n				\
  }\n\
  #if(clip_on[2])\n\
  intersection {\n\
    sphere {P1, R1 material {M} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
    clip[2]\n\
  }\n\
  #else\n\
  sphere {P1, R1 material {M} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
  #end\n\
#end\n");
  //define a cylinder
  fprintf(HeaderPov,"\
#macro ElCylinder (P1, P2, R1, C1, O1)\n\
  #local T = texture { pigment { C1 } }\n\
  #if(clip_on[2])\n\
  intersection {\n\
    cylinder {P1, P2, R1 #if(O1) open #end texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow\n\
}\n\
    clip[2]\n\
  }\n\
  #else\n\
  cylinder {P1, P2, R1 #if(O1) open #end texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
  #end\n\
#end\n");
  //define a triangle
  fprintf(HeaderPov,"\
#macro ElTriangle (P1, N1, P2, N2, P3, N3, C1)\n\
  #local T = texture { pigment { C1 } }\n\
  smooth_triangle {P1, N1, P2, N2, P3, N3 texture {T} #if(clip_on[1]) clipped_by {clip[1]} #end no_shadow}\n\
#end\n");
  //set the camera
  double Loc[3] = {0.5,4.0,4.};
  double LAt[3] = {0.5,0.5,0.5};
  double Up[3] = {0.,6.,0.};
  double Right[3] = {4.8,0.,0.};
  double Dir[3] = {0.,0.,4.};
  double Light0[3] = {0.5,4.0,4.0};
  double Light1[3] = {0.5,4.0,2.0};
  double CBack[3] = {1.,1.,1.};
  #ifdef __glut_h__
  // Draw *Dr;
  // Loc[0] = Dr->xp; Loc[1] = Dr->yp; Loc[2] = Dr->zp+Dr->zw;
  // Light0[0] = Dr->xl0; Light0[1] = Dr->yl0; Light0[2] = Dr->zl0;
  // Light1[0] = Dr->xl1; Light1[1] = Dr->yl1; Light1[2] = Dr->zl1;
  // CBack[0] = Dr->Rback; CBack[1] = Dr->Gback; CBack[2] = Dr->Bback;
  #endif //__glut_h__
  fprintf(HeaderPov,"camera {\n\
//perspective orthographic fisheye ultra_wide_angle...\n\
  location <%.3f, %.3f, %.3f>\n\
  look_at <%.3f, %.3f, %.3f>\n\
  sky <0. 0. 1.>\n\
//  up y\n\
//  right 1.33*x\n\
//  direction <%.3f, %.3f, %.3f>\n\
//  translate <0.0, 1.0, 1.0>\n\
//  rotate <0.0, 1.0, 1.0>\n\
    angle %f\n		    \
}\n",
	  Loc[0],Loc[1],Loc[2],
	  LAt[0],LAt[1],LAt[2],
	  Dir[0],Dir[1],Dir[2],
	  PerspAngle);
  //set the light0
  fprintf(HeaderPov,"\
light_source {\n\
  <%lf, %lf, %lf>\n\
  color rgb<1.000, 1.000, 1.000>\n\
  parallel\n\
  point_at <.5, .5, .5>\n\
}\n",Light0[0],Light0[1],Light0[2]);
  //set the light1
  fprintf(HeaderPov,"\
light_source {\n\
  <%lf, %lf, %lf>\n\
  color rgb<1.000, 1.000, 1.000>\n\
  parallel\n\
  point_at <.5, .5, .5>\n\
}\n",Light1[0],Light1[1],Light1[2]);
  //set the fog
  fprintf(HeaderPov,"\
  fog { fog_type   2\n\
      distance   50\n\
      color       rgb<1.000, 1.000, 1.000>\n\
      fog_offset 0.1\n\
      fog_alt    1.5\n\
      turbulence 1.8\n\
    }\n");
  //background
  fprintf(HeaderPov,"\
background {\n\
  color rgb<%lf, %lf, %lf>\n\
}\n",CBack[0],CBack[1],CBack[2]);
  //texture
  fprintf(HeaderPov,"\
#default { texture {\n\
 finish { ambient 0.000 diffuse 0.650 phong 0.1 phong_size 40.000 specular 0.500 }\n\
} }\n");
}
void ElPoly::ConvLattice(int NSample,char *FName){
  int NType = 3;
  double **Plot = (double **) calloc(3,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot[t] = (double *)calloc(NSample*NSample,sizeof(double));
  }
  LoadDensFile(Plot,NSample);
  FILE *FWrite = fopen(FName,"w");
  fprintf(FWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(0),pEdge(1),pEdge(2),NSample,ChooseDraw(EL_PART));
  double InvNSample = 1./(double)NSample;
  for(int t=0;t<3;t++){
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
	//if(fabs(Plot[t][sx*NSample+sy]) < 0.01) continue;
	double x = sx*InvNSample*pEdge(CLat1);
	double y = sy*InvNSample*pEdge(CLat2);
	fprintf(FWrite,"{t[%d %d %d] x(%lf %lf %lf)}\n",t*pNPart()+sx*NSample+sy,t,t,x,y,Plot[t][sx*NSample+sy]);
      }
    }
  }
  for(int t=0;t<3;t++) free(Plot[t]);
  free(Plot);
}
