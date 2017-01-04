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
#ifdef USE_CGAL
#ifndef __CGAL_h__
#include <fstream>
#include <list>
// #include <CGAL/basic.h>
// #include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/make_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
//#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include "skin_surface_writer.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_refinement_policy_3.h>
//#include <CGAL/mesh_skin_surface_3.h>
//#include <CGAL/subdivide_skin_surface_mesh_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K> Traits;
typedef CGAL::Skin_surface_3<Traits>   Skin_surface_3;
typedef Skin_surface_3::FT             FT;
typedef Skin_surface_3::Weighted_point Weighted_point;
//typedef CGAL::Weighted_point<Bare_point,K::RT>              Weighted_point;
typedef Weighted_point::Point          Bare_point;
// typedef K::Point_3                  Bare_point;
typedef CGAL::Polyhedron_3<K,
			   CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;
//typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef Polyhedron::Traits::Vector_3   Vector;
typedef Polyhedron::Vertex_iterator    Vertex_iterator;
typedef Polyhedron::Halfedge_iterator  Halfedge_iterator;
typedef Polyhedron::Facet_iterator     Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator HFC;
typedef Polyhedron::Vertex_handle      Vertex_handle;
// CGAL::Skin_surface_refinement_policy_3<SkinSurface, Polyhedron> policy(skin);



void ElPoly::DefineSkin(int NSample){
  std::list<Weighted_point> l;
  FT shrinkfactor = 0.5;
  double *Plot  = new double[pNType()*CUBE(NSample)];
  double *Count = new double[CUBE(NSample)];
  double Thre = 10.;
  double Radius = pEdge(0)/(double)NSample;
  for(int p=0;p<pNPart();p++){
    int t = pType(p);
    int vx = (int)(pPos(p,0)/pEdge(0)*NSample);
    int vy = (int)(pPos(p,1)/pEdge(1)*NSample);
    int vz = (int)(pPos(p,2)/pEdge(2)*NSample);
    int vTot = (vz*NSample+vy)*NSample+vx;
    Plot[vTot*pNType()+t] += 1.;
  }
  double *Norm = (double *)calloc(pNType(),sizeof(double));
  for(int t=0;t<pNType();t++){
    for(int v=0;v<CUBE(NSample);v++){
      if(Norm[t] < Plot[v*pNType()+t])
	Norm[t] = Plot[v*pNType()+t];
    }
    Norm[t] = Norm[t] <= 0. ? 1. : Norm[t];
  }
  for(int vx=0;vx<NSample;vx++){
    double x = vx*pEdge(0)/(double)NSample;
    for(int vy=0;vy<NSample;vy++){
      double y = vy*pEdge(1)/(double)NSample;
      for(int vz=0;vz<NSample;vz++){
	double z = vz*pEdge(2)/(double)NSample;
	int vTot = (vz*NSample+vy)*NSample+vx;
	if(Plot[vTot*pNType()] > Thre){
	  l.push_front(Weighted_point(Bare_point(x,y,z),Radius));
	}
      }
    }
  }

  Polyhedron Polyhe;

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);
  CGAL::mesh_skin_surface_3(skin_surface, Polyhe);

  //  CGAL::subdivide_skin_surface_mesh_3(skin_surface, Polyhe);

  // std::ofstream out("mesh.off");
  // out << Polyhe;

  glDeleteLists(Dr->Particles,1);
  Dr->Particles = glGenLists(1);
  glNewList(Dr->Particles,GL_COMPILE);

  // Polyhedron::Facet_iterator fcUp = Polyhe.facets_begin();
  // for(;fcUp != Polyhe.facets_end(); ++fcUp){
  //   Polyhedron::Supports_facet_halfedge = fcUp.halfedge();
  //   //Halfedge_around_facet_circulator heUp = fcUp.halfedge();
  // }

  // for (Vertex_iterator vit = Polyhe.vertices_begin();vit != Polyhe.vertices_end(); vit++){
  //   // Vector n = policy.normal(vit);
  //   // n = n/sqrt(n*n);
  //   cout << vit->point() << std::endl;
  //   Halfedge_iterator heUp = Polyhe.halfedges_begin();
  //   for(;heUp != Polyhe.halfedges_end(); ++heUp){
  //     //Polyhedron::Halfedge_handle Half = *heUp;
  //     Vertex_handle veUp = heUp->vertex();
  //     K::Point_3 pf1 = vit->point();
  //   }
  // }
  CGAL::Inverse_index<Vertex_handle> index(Polyhe.vertices_begin(),
  					   Polyhe.vertices_end());

  for(Facet_iterator fi = Polyhe.facets_begin();fi != Polyhe.facets_end(); ++fi) {
    HFC hc = fi->facet_begin();
    HFC hc_end = hc;
    Polyhedron::Vertex_handle vf1 = (*hc).vertex();
    hc++;
    Polyhedron::Vertex_handle vf2 = (*hc).vertex();
    hc++;
    Polyhedron::Vertex_handle vf3 = (*hc).vertex();
    hc++;
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


    // glPushMatrix();//Particle
    // glBegin(GL_LINES);
    // do {
    //   Polyhedron::Vertex_handle vh = (*hc).vertex();
    //   K::Point_3 pf1 = vh->point();
    //   glVertex3d(pf1.x(),pf1.y(),pf1.z());
    // } while (++hc != hc_end);
    // glEnd();
    // glPopMatrix();//Particle
  }
  
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
  //   Point_3<K> ChPos(pPos(c,CLat1),pPos(c,CLat2),pPos(c,CNorm));
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
void ElPoly::DefineSkin(){
  printf("DefSurface without CGAL not implemented\n");
}
#endif // USE_CGAL
#endif // __CGAL_h__
#endif //__glut_h__
