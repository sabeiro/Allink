#include "VarData.h"

//
// Marching Cubes Example Program 
// by Cory Bloyd (corysama@yahoo.com)
//
// A simple, portable and complete implementation of the Marching Cubes
// and Marching Tetrahedrons algorithms in a single source file.
// There are many ways that this code could be made faster, but the 
// intent is for the code to be easy to understand.
//
// For a description of the algorithm go to
// http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
//
// This code is public domain.
//
XYZ VertexInterp(double isolevel,XYZ p1,XYZ p2,double valp1,double valp2);
int PolygoniseCube(GRIDCELL g,double iso,VAR_TRIANGLE *tri);
int PolygoniseSquare(GRIDCELL g,double iso,VAR_LINE *tri);
XYZ TriangNorm(XYZ p1,XYZ p2,XYZ p3);
VAR_TRIANGLE *VarData::MarchingCubes(double *Plot,int NSample,double IsoLevel,int *NTriang){
  int NTri = 0;
  double InvNSample = 1./(double)NSample;
  int NVertex = 2*NSample-1;
  GRIDCELL Grid;
  VAR_TRIANGLE triangles[10];
  VAR_TRIANGLE *Triang = NULL;
  for(int sx = 0;sx<NSample-1;sx++){
    for(int sy = 0;sy<NSample-1;sy++){
      for(int sz = 0;sz<NSample-1;sz++){
	Grid.p[0].x[0] = (sx)*InvNSample*pEdge(0);
	Grid.p[0].x[1] = (sy)*InvNSample*pEdge(1);
	Grid.p[0].x[2] = (sz)*InvNSample*pEdge(2);
	Grid.val[0] = Plot[((sx)*NSample+(sy))*NSample+sz];
	Grid.p[1].x[0] = (sx+1)*InvNSample*pEdge(0);
	Grid.p[1].x[1] = (sy)*InvNSample*pEdge(1);
	Grid.p[1].x[2] = (sz)*InvNSample*pEdge(2);
	Grid.val[1] = Plot[((sx+1)*NSample+(sy))*NSample+sz];
	Grid.p[2].x[0] = (sx+1)*InvNSample*pEdge(0);
	Grid.p[2].x[1] = (sy+1)*InvNSample*pEdge(1);
	Grid.p[2].x[2] = (sz)*InvNSample*pEdge(2);
	Grid.val[2] = Plot[((sx+1)*NSample+(sy+1))*NSample+sz];
	Grid.p[3].x[0] = (sx)*InvNSample*pEdge(0);
	Grid.p[3].x[1] = (sy+1)*InvNSample*pEdge(1);
	Grid.p[3].x[2] = (sz)*InvNSample*pEdge(2);
	Grid.val[3] = Plot[((sx)*NSample+(sy+1))*NSample+sz];
	Grid.p[4].x[0] = (sx)*InvNSample*pEdge(0);
	Grid.p[4].x[1] = (sy)*InvNSample*pEdge(1);
	Grid.p[4].x[2] = (sz+1)*InvNSample*pEdge(2);
	Grid.val[4] = Plot[((sx)*NSample+(sy))*NSample+sz+1];
	Grid.p[5].x[0] = (sx+1)*InvNSample*pEdge(0);
	Grid.p[5].x[1] = (sy)*InvNSample*pEdge(1);
	Grid.p[5].x[2] = (sz+1)*InvNSample*pEdge(2);
	Grid.val[5] = Plot[((sx+1)*NSample+(sy))*NSample+sz+1];
	Grid.p[6].x[0] = (sx+1)*InvNSample*pEdge(0);
	Grid.p[6].x[1] = (sy+1)*InvNSample*pEdge(1);
	Grid.p[6].x[2] = (sz+1)*InvNSample*pEdge(2);
	Grid.val[6] = Plot[((sx+1)*NSample+(sy+1))*NSample+sz+1];
	Grid.p[7].x[0] = (sx)*InvNSample*pEdge(0);
	Grid.p[7].x[1] = (sy+1)*InvNSample*pEdge(1);
	Grid.p[7].x[2] = (sz+1)*InvNSample*pEdge(2);
	Grid.val[7] = Plot[((sx)*NSample+(sy+1))*NSample+sz+1];
	int n = PolygoniseCube(Grid,IsoLevel,triangles);
	Triang = (VAR_TRIANGLE *)realloc(Triang,(NTri+n)*sizeof(VAR_TRIANGLE));
	int s[3];
	for(int l=0;l<n;l++){
	  for(int v=0;v<3;v++){
	    for(int d=0;d<3;d++){
	      s[d] = (int)(triangles[l].p[v].x[d]*pInvEdge(d)*NVertex);
	      if(s[d] < 0 || s[d] >= NVertex){
		printf("Marching: 0 <= sx %d < NVertex %d  0. < %lf < %lf\n",sx,NVertex,triangles[l].p[v].x[d],pEdge(d));
		continue;
	      }
	    }
	    triangles[l].v[v] = (s[0]*NVertex+s[1])*NVertex+s[2];
	  }
	  Triang[NTri+l] = triangles[l];
	}
	NTri += n;
      }
    }
  }
  *NTriang = NTri;
  return Triang;
}
/**--------------------------------------------------------------
  Given a grid cell and an isolevel, calculate the triangular
  facets requied to represent the isosurface through the cell.
  Return the number of triangular facets, the array "triangles"
  will be loaded up with the vertices at most 5 triangular facets.
  0 will be returned if the grid cell is either totally above
  of totally below the isolevel.
*/
int PolygoniseCube(GRIDCELL g,double iso,VAR_TRIANGLE *tri){
  extern int edgeTable[256];
  extern int triTable[256][16];
  int i,ntri = 0;
  int cubeindex;
  XYZ vertlist[12];
  double OneThird = 1./3.;
  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  cubeindex = 0;
  if (g.val[0] < iso) cubeindex |= 1;
  if (g.val[1] < iso) cubeindex |= 2;
  if (g.val[2] < iso) cubeindex |= 4;
  if (g.val[3] < iso) cubeindex |= 8;
  if (g.val[4] < iso) cubeindex |= 16;
  if (g.val[5] < iso) cubeindex |= 32;
  if (g.val[6] < iso) cubeindex |= 64;
  if (g.val[7] < iso) cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return(0);

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1) {
    vertlist[0] = VertexInterp(iso,g.p[0],g.p[1],g.val[0],g.val[1]);
  }
  if (edgeTable[cubeindex] & 2) {
    vertlist[1] = VertexInterp(iso,g.p[1],g.p[2],g.val[1],g.val[2]);
  }
  if (edgeTable[cubeindex] & 4) {
    vertlist[2] = VertexInterp(iso,g.p[2],g.p[3],g.val[2],g.val[3]);
  }
  if (edgeTable[cubeindex] & 8) {
    vertlist[3] = VertexInterp(iso,g.p[3],g.p[0],g.val[3],g.val[0]);
  }
  if (edgeTable[cubeindex] & 16) {
    vertlist[4] = VertexInterp(iso,g.p[4],g.p[5],g.val[4],g.val[5]);
  }
  if (edgeTable[cubeindex] & 32) {
    vertlist[5] = VertexInterp(iso,g.p[5],g.p[6],g.val[5],g.val[6]);
  }
  if (edgeTable[cubeindex] & 64) {
    vertlist[6] = VertexInterp(iso,g.p[6],g.p[7],g.val[6],g.val[7]);
  }
  if (edgeTable[cubeindex] & 128) {
    vertlist[7] = VertexInterp(iso,g.p[7],g.p[4],g.val[7],g.val[4]);
  }
  if (edgeTable[cubeindex] & 256) {
    vertlist[8] = VertexInterp(iso,g.p[0],g.p[4],g.val[0],g.val[4]);
  }
  if (edgeTable[cubeindex] & 512) {
    vertlist[9] = VertexInterp(iso,g.p[1],g.p[5],g.val[1],g.val[5]);
  }
  if (edgeTable[cubeindex] & 1024) {
    vertlist[10] = VertexInterp(iso,g.p[2],g.p[6],g.val[2],g.val[6]);
  }
  if (edgeTable[cubeindex] & 2048) {
    vertlist[11] = VertexInterp(iso,g.p[3],g.p[7],g.val[3],g.val[7]);
  }
  /* Create the triangles */
  for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
    tri[ntri].p[0] = vertlist[triTable[cubeindex][i  ]];
    tri[ntri].p[1] = vertlist[triTable[cubeindex][i+1]];
    tri[ntri].p[2] = vertlist[triTable[cubeindex][i+2]];
    for(int d=0;d<3;d++){
      tri[ntri].c.x[d] = (tri[ntri].p[0].x[d]+tri[ntri].p[1].x[d]+tri[ntri].p[2].x[d])*OneThird;
    }
    //TODO: Find the normal to the vertices, 
    // i.e. how the triangles are connected
    tri[ntri].n[0] = TriangNorm(tri[ntri].p[0],tri[ntri].p[1],tri[ntri].p[2]);
    tri[ntri].n[1] = TriangNorm(tri[ntri].p[0],tri[ntri].p[1],tri[ntri].p[2]);
    tri[ntri].n[2] = TriangNorm(tri[ntri].p[0],tri[ntri].p[1],tri[ntri].p[2]);
    //if(fabs(tri[ntri].n[0].x+tri[ntri].n[0].y+tri[ntri].n[0].z) <= 0.) continue;
    ntri++;
  }
  return(ntri);
}
VAR_LINE *VarData::MarchingSquares(double *Plot,int NSample,double IsoLevel,int *NTriang){
  int NTri = 0;
  double InvNSample = 1./(double)NSample;
  int NVertex = 2*NSample-1;
  GRIDCELL Grid;
  VAR_LINE triangles[4];
  VAR_LINE *Triang = NULL;
  for(int sx = 0;sx<NSample-1;sx++){
    for(int sy = 0;sy<NSample-1;sy++){
      Grid.p[0].x[0] = (sx)*InvNSample*pEdge(0);
      Grid.p[0].x[1] = (sy)*InvNSample*pEdge(1);
      Grid.p[0].x[2] = IsoLevel;
      Grid.val[0] = Plot[(sx)*NSample+(sy)];
      Grid.p[1].x[0] = (sx+1)*InvNSample*pEdge(0);
      Grid.p[1].x[1] = (sy)*InvNSample*pEdge(1);
      Grid.p[1].x[2] = IsoLevel;
      Grid.val[1] = Plot[(sx+1)*NSample+(sy)];
      Grid.p[2].x[0] = (sx+1)*InvNSample*pEdge(0);
      Grid.p[2].x[1] = (sy+1)*InvNSample*pEdge(1);
      Grid.p[2].x[2] = IsoLevel;
      Grid.val[2] = Plot[(sx+1)*NSample+(sy+1)];
      Grid.p[3].x[0] = (sx)*InvNSample*pEdge(0);
      Grid.p[3].x[1] = (sy+1)*InvNSample*pEdge(1);
      Grid.p[3].x[2] = IsoLevel;
      Grid.val[3] = Plot[(sx)*NSample+(sy+1)];
      int n = PolygoniseSquare(Grid,IsoLevel,triangles);
      Triang = (VAR_LINE *)realloc(Triang,(NTri+n)*sizeof(VAR_LINE));
      int s[3];
      for(int l=0;l<n;l++){
	for(int v=0;v<2;v++){
	  for(int d=0;d<2;d++){
	    s[d] = (int)(triangles[l].p[v].x[d]*pInvEdge(d)*NVertex);
	    if(s[d] < 0 || s[d] >= NVertex){
	      printf("Marching squares: 0 <= sx %d < NVertex %d  0. < %lf < %lf\n",s[d],NVertex,triangles[l].p[v].x[d],pEdge(d));
	      continue;
	    }
	  }
	  triangles[l].v[v] = s[0]*NVertex+s[1];
	}
	Triang[NTri+l] = triangles[l];
      }
      NTri += n;
    }
  }
  *NTriang = NTri;
  return Triang;
}
/**--------------------------------------------------------------
  Given a grid cell and an isolevel, calculate the triangular
  facets requied to represent the isosurface through the cell.
  Return the number of triangular facets, the array "triangles"
  will be loaded up with the vertices at most 5 triangular facets.
  0 will be returned if the grid cell is either totally above
  of totally below the isolevel.
*/
int PolygoniseSquare(GRIDCELL g,double iso,VAR_LINE *tri){
  extern int edgeTable2d[16];
  extern int triTable2d[16][4];
  int i,ntri = 0;
  int cubeindex;
  XYZ vertlist[4];
  XYZ Normal;
  Normal.x[0] = 0.;
  Normal.x[1] = 0.;
  Normal.x[2] = 1.;
  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  cubeindex = 0;
  if (g.val[0] < iso) cubeindex |= 1;
  if (g.val[1] < iso) cubeindex |= 2;
  if (g.val[2] < iso) cubeindex |= 4;
  if (g.val[3] < iso) cubeindex |= 8;
  /* Square is entirely in/out of the surface */
  if (edgeTable2d[cubeindex] == 0)
    return(0);
  /* Find the vertices where the surface intersects the cube */
  if (edgeTable2d[cubeindex] & 1) {
    vertlist[0] = VertexInterp(iso,g.p[0],g.p[1],g.val[0],g.val[1]);
  }
  if (edgeTable2d[cubeindex] & 2) {
    vertlist[1] = VertexInterp(iso,g.p[1],g.p[2],g.val[1],g.val[2]);
  }
  if (edgeTable2d[cubeindex] & 4) {
    vertlist[2] = VertexInterp(iso,g.p[2],g.p[3],g.val[2],g.val[3]);
  }
  if (edgeTable2d[cubeindex] & 8) {
    vertlist[3] = VertexInterp(iso,g.p[3],g.p[0],g.val[3],g.val[0]);
  }
  /* Create the triangles */
  for(i=0;triTable2d[cubeindex][i]!=-1;i+=2){
    tri[ntri].p[0] = vertlist[triTable2d[cubeindex][i  ]];
    tri[ntri].p[1] = vertlist[triTable2d[cubeindex][i+1]];
    for(int d=0;d<2;d++){
      tri[ntri].c.x[d] = (tri[ntri].p[0].x[d]+tri[ntri].p[1].x[d])*.5;
    }
    //TODO: Find the normal to the vertices, 
    // i.e. how the triangles are connected
    tri[ntri].n[0] = Normal;
    tri[ntri].n[1] = Normal;
    //if(fabs(tri[ntri].n[0].x+tri[ntri].n[0].y+tri[ntri].n[0].z) <= 0.) continue;
    ntri++;
  }
  return(ntri);
}
XYZ TriangNorm(XYZ p1,XYZ p2,XYZ p3){
  XYZ n;
  double Norm = 0.;
  n.x[0] = (p2.x[1]-p1.x[1])*(p2.x[2]-p3.x[2]) 
    - (p2.x[1]-p3.x[1])*(p2.x[2]-p1.x[2]);
  n.x[1] = (p2.x[2]-p1.x[2])*(p2.x[0]-p3.x[0]) 
    - (p2.x[2]-p3.x[2])*(p2.x[0]-p1.x[0]);
  n.x[2] = (p2.x[0]-p1.x[0])*(p2.x[1]-p3.x[1]) 
    - (p2.x[0]-p3.x[0])*(p2.x[1]-p1.x[1]);
  Norm = sqrt(SQR(n.x[0])+SQR(n.x[1])+SQR(n.x[2]));
  if(Norm > 0.){
    Norm = 1./Norm;
    n.x[0] *= Norm;
    n.x[1] *= Norm;
    n.x[2] *= Norm;
  }
  return n;
}
/**-------------------------------------------------------------------------
  Return the point between two points in the same ratio as
  isolevel is between valp1 and valp2
*/
XYZ VertexInterp(double isolevel,XYZ p1,XYZ p2,double valp1,double valp2){
  double mu;
  XYZ p;
  if (ABS(isolevel-valp1) < 0.00001)
    return(p1);
  if (ABS(isolevel-valp2) < 0.00001)
    return(p2);
  if (ABS(valp1-valp2) < 0.00001)
    return(p1);
  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x[0] = p1.x[0] + mu * (p2.x[0] - p1.x[0]);
  p.x[1] = p1.x[1] + mu * (p2.x[1] - p1.x[1]);
  p.x[2] = p1.x[2] + mu * (p2.x[2] - p1.x[2]);
  return(p);
}
/**
  int edgeTable[256].  It corresponds to the 2^8 possible combinations of
  of the eight (n) vertices either existing inside or outside (2^n) of the
  surface.  A vertex is inside of a surface if the value at that vertex is
  less than that of the surface you are scanning for.  The table index is
  constructed bitwise with bit 0 corresponding to vertex 0, bit 1 to vert
  1.. bit 7 to vert 7.  The value in the table tells you which edges of
  the table are intersected by the surface.  Once again bit 0 corresponds
  to edge 0 and so on, up to edge 12.
  Constructing the table simply consisted of having a program run thru
  the 256 cases and setting the edge bit if the vertices at either end of
  the edge had different values (one is inside while the other is out).
  The purpose of the table is to speed up the scanning process.  Only the
  edges whose bit's are set contain vertices of the surface.
  Vertex 0 is on the bottom face, back edge, left side.
  The progression of vertices is clockwise around the bottom face
  and then clockwise around the top face of the cube.  Edge 0 goes from
  vertex 0 to vertex 1, Edge 1 is from 2->3 and so on around clockwise to
  vertex 0 again. Then Edge 4 to 7 make up the top face, 4->5, 5->6, 6->7
  and 7->4.  Edge 8 thru 11 are the vertical edges from vert 0->4, 1->5,
  2->6, and 3->7.
      4--------5     *---4----*
     /|       /|    /|       /|
    / |      / |   7 |      5 |
   /  |     /  |  /  8     /  9
  7--------6   | *----6---*   |
  |   |    |   | |   |    |   |
  |   0----|---1 |   *---0|---*
  |  /     |  /  11 /     10 /
  | /      | /   | 3      | 1
  |/       |/    |/       |/
  3--------2     *---2----*
*/
int edgeTable[256]={
  0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   
};
/**
  int triTable[256][16] also corresponds to the 256 possible combinations
  of vertices.
  The [16] dimension of the table is again the list of edges of the cube
  which are intersected by the surface.  This time however, the edges are
  enumerated in the order of the vertices making up the triangle mesh of
  the surface.  Each edge contains one vertex that is on the surface.
  Each triple of edges listed in the table contains the vertices of one
  triangle on the mesh.  The are 16 entries because it has been shown that
  there are at most 5 triangles in a cube and each "edge triple" list is
  terminated with the value -1.
  For example triTable[3] contains
  {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
  This corresponds to the case of a cube whose vertex 0 and 1 are inside
  of the surface and the rest of the verts are outside (00000001 bitwise
  OR'ed with 00000010 makes 00000011 == 3).  Therefore, this cube is
  intersected by the surface roughly in the form of a plane which cuts
  edges 8,9,1 and 3.  This quadrilateral can be constructed from two
  triangles: one which is made of the intersection vertices found on edges
  1,8, and 3; the other is formed from the vertices on edges 9,8, and 1.
  Remember, each intersected edge contains only one surface vertex.  The
  vertex triples are listed in counter clockwise order for proper facing.
*/
int triTable[256][16] =
  {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
   {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
   {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
   {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
   {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
   {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
   {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
   {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
   {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
   {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
   {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
   {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
   {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
   {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
   {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
   {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
   {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
   {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
   {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
   {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
   {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
   {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
   {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
   {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
   {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
   {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
   {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
   {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
   {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
   {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
   {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
   {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
   {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
   {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
   {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
   {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
   {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
   {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
   {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
   {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
   {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
   {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
   {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
   {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
   {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
   {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
   {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
   {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
   {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
   {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
   {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
   {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
   {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
   {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
   {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
   {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
   {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
   {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
   {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
   {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
   {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
   {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
   {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
   {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
   {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
   {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
   {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
   {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
   {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
   {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
   {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
   {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
   {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
   {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
   {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
   {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
   {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
   {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
   {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
   {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
   {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
   {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
   {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
   {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
   {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
   {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
   {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
   {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
   {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
   {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
   {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
   {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
   {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
   {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
   {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
   {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
   {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
   {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
   {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
   {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
   {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
   {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
   {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
   {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
   {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
   {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
   {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
   {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
   {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
   {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
   {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
   {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
   {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
   {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
   {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
   {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
   {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
   {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
   {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
   {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
   {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
   {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
   {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
   {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
   {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
   {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
   {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
   {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
   {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
   {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
   {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
   {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
   {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
   {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
   {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
   {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
   {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
   {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
   {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
   {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
   {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
   {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
   {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
   {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
   {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
   {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
   {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
   {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
   {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
   {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
   {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
   {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
   {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
   {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
   {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
   {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
   {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
   {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
   {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
   {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
   {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
   {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
   {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
   {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
   {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
   {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
   {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
   {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
   {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
   {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
   {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
   {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
   {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
   {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
   {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
   {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
   {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
   {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
   {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
   {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
   {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
/**
  int edgeTable[16].  It corresponds to the 2^4 possible combinations of
  of the eight (n) vertices either existing inside or outside (2^n) of the
  surface.  A vertex is inside of a surface if the value at that vertex is
  less than that of the surface you are scanning for.  The table index is
  constructed bitwise with bit 0 corresponding to vertex 0, bit 1 to vert
  1.. bit 3 to vert 3.  The value in the table tells you which edges of
  the table are intersected by the surface.  Once again bit 0 corresponds
  to edge 0 and so on, up to edge 3.
  Constructing the table simply consisted of having a program run thru
  the 256 cases and setting the edge bit if the vertices at either end of
  the edge had different values (one is inside while the other is out).
  The purpose of the table is to speed up the scanning process.  Only the
  edges whose bit's are set contain vertices of the surface.
  Vertex 0 is on the bottom face, back edge, left side.
  The progression of vertices is clockwise around the bottom face
  and then clockwise around the top face of the cube.  Edge 0 goes from
  vertex 0 to vertex 1, Edge 1 is from 2->3 and so on around clockwise to
  vertex 0 again. 

  3---2---2
  |       |
  |       |
  3       1
  |       |
  |       |
  0---0---1
*/
int edgeTable2d[16]={
  0x0 , 0x09, 0x03, 0x0a, 0x06, 0x0f, 0x05, 0x0c, 
  0x0c, 0x05, 0x0f, 0x06, 0x0a, 0x03, 0x09, 0x0
};
int triTable2d[16][4] = {
  {-1,-1,-1,-1},
  {0,3,-1,-1},
  {0,1,-1,-1},
  {3,1,-1,-1},
  {1,2,-1,-1},
  {0,1,2,3},
  {0,2,-1,-1},
  {2,3,-1,-1},
  {2,3,-1,-1},
  {2,0,-1,-1},
  {3,0,1,2},
  {1,2,-1,-1},
  {1,3,-1,-1},
  {0,1,-1,-1},
  {3,0,-1,-1},
  {-1,-1,-1,-1}
};
  