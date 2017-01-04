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
  glLineWidth(.5);
  Delaunay dtUp;
  Delaunay dtDown;
  //vector <COLOR> HueUp;
  //vector <COLOR> HueDown;
  COLOR ColChain;
  ColChain.a = 1.;
  double AreaMean = Gen->Edge[CLat1]*Gen->Edge[CLat2]*.5/(double)Gen->NChain;
  Delaunay::Finite_vertices_iterator vUp = dtUp.finite_vertices_begin();
  for(int c=0;c<Gen->NChain;c++){
    int Chc = Ch[c].Type;
    if(Ch[c].Pos[CNorm] > .7*Gen->Edge[CNorm]) continue;
    if(Ch[c].Pos[CNorm] < .3*Gen->Edge[CNorm]) continue;    
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
    if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP) ){
      K::Point_3 ChPos(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2],Ch[c].Pos[CNorm]);
      dtUp.insert(ChPos);
      //HueUp.push_back(ColChain);
    }
    else{
      K::Point_3 ChPos(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2],Ch[c].Pos[CNorm]);
      dtDown.insert(ChPos);
      //HueDown.push_back(ColChain);
    }
  }
  Delaunay::Face_iterator fcTrUp = dtUp.finite_faces_begin();
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
    //if(vN.Norm() > 2.*AreaMean) continue;
    glColor4f(0.,.7,0.,1.);
    //glColor4f(HueUp[p].r,HueUp[p].g,HueUp[p].b,HueUp[p].a);
    DrTria(&v1,&v2,&v3,&vN);
    glColor4f(0.,.0,0.,1.);
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
    //if(vN.Norm() > 2.*AreaMean) continue;
    glColor4f(0.,.7,0.,1.);
    //    glColor4f(HueDown[p].r,HueDown[p].g,HueDown[p].b,HueDown[p].a);
    DrTria(&v1,&v2,&v3,&vN);
    glColor4f(0.,.0,0.,1.);
    DrTriaContour(&v1,&v2,&v3);
  }
  DrNano();
  glEndList();
}

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/make_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_3.h>
#include <list>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
//#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <fstream>
#include "skin_surface_writer.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_refinement_policy_3.h>
//#include <CGAL/mesh_skin_surface_3.h>
//#include <CGAL/subdivide_skin_surface_mesh_3.h>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;
typedef Polyhedron::Traits::Vector_3                 Vector;
//CGAL::Skin_surface_refinement_policy_3<SkinSurface, Polyhedron> policy(skin);
#include <CGAL/Skin_surface_refinement_policy_3.h>
// typedef K::Point_3                                          Bare_point;
// typedef CGAL::Weighted_point<Bare_point,K::RT>              Weighted_point;
// typedef CGAL::Polyhedron_3<K>                               Polyhedron;


void ElPoly::DefineSkin(){
  list <Weighted_point> WeiPoint;
  double shrinkfactor = 1.0;

  for(int c=0;c<4;c++){//Gen->NChain;c++){
    Bare_point ChPoint(Ch[c].Pos[0]*InvScaleUn,Ch[c].Pos[1]*InvScaleUn,Ch[c].Pos[2]*InvScaleUn);
    WeiPoint.push_front(Weighted_point(ChPoint, 0.5));
  }
  Polyhedron Polyhe;
  //CGAL::make_skin_surface_mesh_3(Polyhe, WeiPoint.begin(), WeiPoint.end(), shrinkfactor);
  Skin_surface_3 skin_surface(WeiPoint.begin(),WeiPoint.end(),shrinkfactor);
  CGAL::mesh_skin_surface_3(skin_surface,Polyhe);
  CGAL::subdivide_skin_surface_mesh_3(skin_surface,Polyhe);
  std::ofstream out("mesh.off");
  out << Polyhe;
  // Polyhedron::Facet_iterator fcUp = Polyhe.facet_begin();
  // for(;fcUp != Polyhe.faces_end(); ++fcUp){
  //   Halfedge_around_facet_circulator heUp = fcUp.halfedge();
  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);
  glColor4f(0.0,0.0,0.0,1.0);
  double Pos[3];
  //  for (Polyhedron::Vertex_iterator vit = Polyhe.vertices_begin();vit != Polyhe.vertices_end(); vit ++) {
    //Vector n = policy.normal(vit);
    //n = n/sqrt(n*n);
    //    out << vit->point() << " " << n << std::endl;
  // Polyhedron::Halfedge_iterator heUp = Polyhe.halfedges_begin();
  // for(;heUp != Polyhe.halfedges_end(); ++heUp){
  //   //Polyhedron::Halfedge_handle Half = *heUp;
  //   Polyhedron::Vertex_handle veUp = heUp->vertex();
//K::Point_3 pf1 = vit->point();


  // for(Polyhedron::Facet_iterator fi = Polyhe.facets_begin();fi != Polyhe.facets_end(); ++fi) {
  //   Skin_surface_3::HFC hc = fi->facet_begin();
  //   Polyhedron::HFC hc_end = hc;
  //   glPushMatrix();//Particle
  //   glBegin(GL_LINES);
  //   do {
  //     Polyhedron::Vertex_handle vh = (*hc).vertex();
  //     K::Point_3 pf1 = vh->point();
  //     Pos[0] = pf1.x();Pos[1] = pf1.y();Pos[2] = pf1.z();
  //     glVertex3d(Pos[0],Pos[1],Pos[2]);
  //   } while (++hc != hc_end);
  //   glEnd();
  //   glPopMatrix();//Particle
    
  // }
  glEndList();
  
  // Tr tr;     // 3D-Delaunay triangulation
  // C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // Surface_3 surface(sphere_function,Sphere_3(CGAL::ORIGIN, 2.)); 
  // Surface_mesh_default_criteria_3<Tr> criteria(30., 0.1,0.1); 
  // // meshing surface
  // make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  // std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

  // DT dt;
  // for(int c=0;c<Gen->NChain;c++){
  //   if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP) )continue;
  //   Point_3<K> ChPos(Ch[c].Pos[CLat1],Ch[c].Pos[CLat2],Ch[c].Pos[CNorm]);
  //   dt.insert(ChPos);
  // }
  // Face_iterator fcTr = dt.finite_faces_begin();
  // glDeleteLists(Dr->Particles,1);
  // Dr->Particles = glGenLists(1);
  // glNewList(Dr->Particles,GL_COMPILE);
  // for(;fcTr != dt.faces_end(); ++fcTr){
  //   Vertex_handle vf1 = fcTr->vertex(0),
  //     vf2 = fcTr->vertex(1),vf3 = fcTr->vertex(2);
  //   Point pf1 = vf1->point();
  //   Point pf2 = vf2->point();
  //   Point pf3 = vf3->point();
  //   Vettore v1(pf1.x() - pf2.x(),pf1.y() - pf2.y(),pf1.z()-pf2.z());
  //   Vettore v2(pf3.x() - pf2.x(),pf3.y() - pf2.y(),pf1.z()-pf2.z());
  //   Vettore vN(3);
  //   vN = v1 ^ v2;
  //   DrTira(v1,v2,v3,vN);
    
  // }
  // glEndList();
}
#else
void ElPoly::DrTriangulate(){
  printf("DefSurface without CGAL not implemented\n");
}
#endif // __CGAL_h__
#endif //__glut_h__
