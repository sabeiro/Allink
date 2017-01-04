/***********************************************************************
ElPoly:This progam provide a graphical visualisation of the data 
opend by VarData using openGL glut. The most important option are 
the possibility of changing the backfold of the polymers with 'c', 
see the subsequent file in the list with '>', see the bond with 'b'. 
Copyright (C) 2008-2010 by Giovanni Marelli <sabeiro@virgilio.it>


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include "ElPoly.h"
			 //#include <stdafx.h>

#ifdef __glut_h__
extern Draw *Dr;
  //#ifndef __CGAL_h__
#ifdef USE_CGAL
//#include <CGAL/basic.h>
//#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
			 //#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
//using namespace CGAL;

//default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

//c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;

typedef GT::Sphere_3 Sphere_3;
//typedef GT::Point_3 Point_3;
typedef GT::FT FT;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_euclidean_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
// typedef Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
// typedef Triangulation_2<Gt,Tds> Triangulation;
template < class Gt >
class My_face_base : public CGAL::Triangulation_face_base_2<Gt>
{
public:
  CGAL::Color color;
  My_face_base() :
    CGAL::Triangulation_face_base_2<Gt>() {}
  My_face_base(void* v0, void* v1, void* v2) : 
    CGAL::Triangulation_face_base_2<Gt>(v0,v1,v2) {}
  My_face_base(void* v0, void* v1, void* v2, void* n0, void* n1, void* n2) : 
    CGAL::Triangulation_face_base_2<Gt>(v0,v1,v2,n0,n1,n2) {}
};
typedef struct {double r;double g;double b;double a;} COLOR;

void ElPoly::DrTriangulate(){
  if(Dr->IfMaterial){
    glEnable(GL_LIGHTING);
    glEnable( GL_LIGHT0 );
  }
  else 
    glDisable(GL_LIGHTING);
  Delaunay dtUp;
  Delaunay dtDown;
  //vector <COLOR> HueUp;
  //vector <COLOR> HueDown;
  COLOR ColChain;
  GLfloat Color[4];
  ColChain.a = 1.;
  double AreaMean = pEdge(CLat1)*pEdge(CLat2)*.5/(double)pNChain();
  BfDefChain();
  double Pos[3];
  Delaunay::Finite_vertices_iterator vUp = dtUp.finite_vertices_begin();
  for(int c=0;c<pNChain();c++){
    int Chc = Ch[c].Type;
    //if(Ch[c].Pos[CNorm] > .7*pEdge(CNorm)) continue;
    //if(Ch[c].Pos[CNorm] < .3*pEdge(CNorm)) continue;    
    if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
    ColChain.r = 0.;
    if(CHAIN_IF_TYPE(Chc,CHAIN_UP))
      ColChain.r = 1.;
    ColChain.b = 0.;
    if(CHAIN_IF_TYPE(Chc,CHAIN_FLABBY))
      ColChain.b = 1.;
    ColChain.g=.7;
    if(CHAIN_IF_TYPE(Chc,CHAIN_TILTED))
      ColChain.g = .5;
    for(int d=0;d<3;d++)
      Pos[d] = Ch[c].Pos[d] - floor(Ch[c].Pos[d]*pInvEdge(d))*pEdge(d);
    if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP) ){
      K::Point_3 ChPos(Pos[0],Pos[1],Pos[2]);
      dtUp.insert(ChPos);
      //HueUp.push_back(ColChain);
    }
    else{
      K::Point_3 ChPos(Pos[0],Pos[1],Pos[2]);
      dtDown.insert(ChPos);
      //HueDown.push_back(ColChain);
    }
  }
  Delaunay::Face_iterator fcTrUp   = dtUp.finite_faces_begin();
  Delaunay::Face_iterator fcTrDown = dtDown.finite_faces_begin();
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  for(int p=0;fcTrUp != dtUp.faces_end(); ++fcTrUp,p++){
    Delaunay::Vertex_handle vf1 = fcTrUp->vertex(0),
      vf2 = fcTrUp->vertex(1),vf3 = fcTrUp->vertex(2);
    K::Point_3 pf1 = vf1->point();
    K::Point_3 pf2 = vf2->point();
    K::Point_3 pf3 = vf3->point();
    Vettore v1(pf1.x(),pf1.y(),pf1.z());
    Vettore v2(pf2.x(),pf2.y(),pf2.z());
    Vettore v3(pf3.x(),pf3.y(),pf3.z());
    Vettore vN(3);
    v1.Mult(InvScaleUn);
    v2.Mult(InvScaleUn);
    v3.Mult(InvScaleUn);
    vN = (v1-v2) ^ (v3-v2);
    glEnable(GL_LIGHTING);
    double Depth = pf1.z()*pInvEdge(CNorm)*Saturation+ExtParam;
    Dr->DepthMap(Depth,Color);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color); 
    glColor4fv(Color);
    //if(vN.Norm() > 2.*AreaMean) continue;
    DrTria(&v1,&v2,&v3,&vN);
    glDisable(GL_LIGHTING);
    glColor4f(0.,0.,1.,1.);
    DrTriaContour(&v1,&v2,&v3);
  }
  for(int p=0;fcTrDown != dtDown.faces_end(); ++fcTrDown,p++){
    Delaunay::Vertex_handle vf1 = fcTrDown->vertex(0),
      vf2 = fcTrDown->vertex(1),vf3 = fcTrDown->vertex(2);
    K::Point_3 pf1 = vf1->point();
    K::Point_3 pf2 = vf2->point();
    K::Point_3 pf3 = vf3->point();
    Vettore v1(pf1.x(),pf1.y(),pf1.z());
    Vettore v2(pf2.x(),pf2.y(),pf2.z());
    Vettore v3(pf3.x(),pf3.y(),pf3.z());
    Vettore vN(3);
    v1.Mult(InvScaleUn);
    v2.Mult(InvScaleUn);
    v3.Mult(InvScaleUn);
    vN = (v1-v2) ^ (v3-v2);
    glEnable(GL_LIGHTING);
    double Depth = pf1.z()*pInvEdge(CNorm)*Saturation+ExtParam;
    Dr->DepthMap(Depth,Color);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,Color); 
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,Color); 
    glColor4fv(Color);
    //if(vN.Norm() > 2.*AreaMean) continue;
    //    glColor4f(HueDown[p].r,HueDown[p].g,HueDown[p].b,HueDown[p].a);
    DrTria(&v1,&v2,&v3,&vN);
    glDisable(GL_LIGHTING);
    glColor4f(1.,0.,0.,1.);
    DrTriaContour(&v1,&v2,&v3);
  }
  glDisable(GL_LIGHTING);
  for(int n=0;n<pNNano();n++) DrawNano(n);
  glEndList();
}
#include <CGAL/Triangulation_3.h>
typedef CGAL::Triangulation_3<K> Triang;
void ElPoly::DrCells(){
  Delaunay dtUp;
  COLOR ColChain;
  ColChain.a = 1.;
  double AreaMean = pEdge(CLat1)*pEdge(CLat2)*.5/(double)pNChain();
  std::list<Triang::Point> LipPos;
  Delaunay::Finite_vertices_iterator vUp = dtUp.finite_vertices_begin();
  for(int b=0,cOff=0;b<pNBlock();b++,cOff+=Block[b].NChain){
    for(int c=cOff;c<pNChain(b)+cOff;c+=20){
      int Chc = Ch[c].Type;
      if(!CHAIN_IF_TYPE(Chc,NChType)) continue;
      ColChain.r = 0.;
      if(CHAIN_IF_TYPE(Chc,CHAIN_UP))
	ColChain.r = 1.;
      ColChain.b = 0.;
      if(CHAIN_IF_TYPE(Chc,CHAIN_FLABBY))
	ColChain.b = 1.;
      ColChain.g=.7;
      if(CHAIN_IF_TYPE(Chc,CHAIN_TILTED))
	ColChain.g = .5;
      Triang::Point LipPoint(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2],Ch[c].Pos[CNorm]);
      LipPos.push_front(LipPoint);
      //HueUp.push_back(ColChain);
    }
  }
  Triang Tr(LipPos.begin(),LipPos.end());
  Triang::size_type NVert = Tr.number_of_vertices();
  std::vector<Triang::Point> V(3);
  assert(Tr.is_valid());
  // Triang::Locate_type lt;
  // int li,lj;
  // Triang::Point p(0,0,0);
  // Triang::Cell_handle c = Tr.locate(p,lt,li,lj);
  // assert( lt == Triang::VERTEX);
  // assert(c->vertex(li)->point() == p);
  // Triang::Vertex_handle v = c->vertex( (li+1)&3);
  // Triang::Cell_handle nc = c->neighbor(li);
  // int nli;
  // assert(nc->has_vertex(v,nli));
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  //Triang::Cell_iterator CIt = Tr.finite_cells_begin();
  Triang::Finite_cells_iterator CIt = Tr.finite_cells_begin();
  for(int p=0;CIt != Tr.finite_cells_end();++CIt,p++){
    Triang::Vertex_handle vf1 = CIt->vertex(0),vf2 = CIt->vertex(1),vf3 = CIt->vertex(2);
    K::Point_3 pf1 = vf1->point();
    K::Point_3 pf2 = vf2->point();
    K::Point_3 pf3 = vf3->point();
    Vettore v1(pf1.x(),pf1.y(),pf1.z());
    Vettore v2(pf2.x(),pf2.y(),pf2.z());
    Vettore v3(pf3.x(),pf3.y(),pf3.z());
    Vettore vN(3);
    v1.Mult(InvScaleUn);
    v2.Mult(InvScaleUn);
    v3.Mult(InvScaleUn);
    vN = (v1-v2) ^ (v3-v2);
    //if(vN.Norm() > 2.*AreaMean) continue;
    double Sfumatura = .3*Mat->Casuale();
    glColor4f(0.1,.4+Sfumatura,0.2,1.);
    //glColor4f(HueUp[p].r,HueUp[p].g,HueUp[p].b,HueUp[p].a);
    DrTria(&v1,&v2,&v3,&vN);
    glColor4f(0.,.0,0.,1.);
    DrTriaContour(&v1,&v2,&v3);
  }
  for(int n=0;n<pNNano();n++) DrawNano(n);
  glEndList();
}
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

FT sphere_function (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}
void ElPoly::DrGenMesh(){
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  Surface_3 surface(sphere_function,Sphere_3(CGAL::ORIGIN, 2.));
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.1,  // radius bound
                                                     0.1); // distance bound
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  Tr::Finite_cells_iterator CIt = tr.finite_cells_begin();
  for(int p=0;CIt != tr.finite_cells_end();++CIt,p++){
    Tr::Vertex_handle vf1 = CIt->vertex(0),vf2 = CIt->vertex(1),vf3 = CIt->vertex(2);
    K::Point_3 pf1 = vf1->point();
    K::Point_3 pf2 = vf2->point();
    K::Point_3 pf3 = vf3->point();
    Vettore v1(pf1.x(),pf1.y(),pf1.z());
    Vettore v2(pf2.x(),pf2.y(),pf2.z());
    Vettore v3(pf3.x(),pf3.y(),pf3.z());
    Vettore vN(3);
    v1.Mult(InvScaleUn);
    v2.Mult(InvScaleUn);
    v3.Mult(InvScaleUn);
    vN = (v1-v2) ^ (v3-v2);
    //if(vN.Norm() > 2.*AreaMean) continue;
    double Sfumatura = .3*Mat->Casuale();
    glColor4f(0.1,.4+Sfumatura,0.2,1.);
    //glColor4f(HueUp[p].r,HueUp[p].g,HueUp[p].b,HueUp[p].a);
    DrTria(&v1,&v2,&v3,&vN);
    glColor4f(1.,.0,0.,1.);
    DrTriaContour(&v1,&v2,&v3);
  }
  for(int n=0;n<pNNano();n++) DrawNano(n);
  glEndList();
}
#else
void ElPoly::DrTriangulate(){
  printf("DefSurface without CGAL not implemented\n");
}
void ElPoly::DrMesh(){
  printf("DefSurface without CGAL not implemented\n");
}
void ElPoly::DrCells(){
  printf("DefSurface without CGAL not implemented\n");
}
#endif // __CGAL_h__
#endif //__glut_h__
