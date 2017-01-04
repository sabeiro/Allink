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
#include <CGAL/trace.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

void ElPoly::DefineSurf(){
  // }
  // int main(void)
  // {
  // Poisson options
  FT sm_angle = 20.0; // Min triangle angle in degrees.
  FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
  FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.

  // Reads the point set file in points[].
  // Note: read_xyz_points_and_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
  // The position property map can be omitted here as we use iterators over Point_3 elements.
  PointList points;
  std::ifstream stream("data/kitten.xyz");
  if (!stream ||
      !CGAL::read_xyz_points_and_normals(
					 stream,
					 std::back_inserter(points),
					 CGAL::make_normal_of_point_with_normal_pmap(std::back_inserter(points))))
    {
      std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
      return;
    }

  // Creates implicit function from the read points using the default solver (TAUCS).
  // Note: this method requires an iterator over points
  // + property maps to access each point's position and normal.
  // The position property map can be omitted here as we use iterators over Point_3 elements.
  Poisson_reconstruction_function function(
					   points.begin(), points.end(),
					   CGAL::make_normal_of_point_with_normal_pmap(points.begin()));

  // Computes the Poisson indicator function f()
  // at each vertex of the triangulation.
  if ( ! function.compute_implicit_function() )
    return;

  // Computes average spacing
  FT average_spacing = CGAL::compute_average_spacing(points.begin(), points.end(),
						     6 /* knn = 1 ring */);
    
  // Gets one point inside the implicit surface
  // and computes implicit function bounding sphere radius.
  Point inner_point = function.get_inner_point();
  Sphere bsphere = function.bounding_sphere();
  FT radius = std::sqrt(bsphere.squared_radius());

  // Defines the implicit surface: requires defining a
  // conservative bounding sphere centered at inner point.
  FT sm_sphere_radius = 5.0 * radius;
  FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
  Surface_3 surface(function,
		    Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
		    sm_dichotomy_error/sm_sphere_radius);

  // Defines surface mesh generation criteria
  CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
						      sm_radius*average_spacing,  // Max triangle size
						      sm_distance*average_spacing); // Approximation error

  // Generates surface mesh with manifold option
  STr tr; // 3D Delaunay triangulation for surface mesh generation
  C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
  CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
			  surface,                              // implicit surface
			  criteria,                             // meshing criteria
			  CGAL::Manifold_with_boundary_tag());  // require manifold mesh

  if(tr.number_of_vertices() == 0)
    return ;

  // saves reconstructed surface mesh
  std::ofstream out("kitten_poisson-20-30-0.375.off");
  Polyhedron output_mesh;
  CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
  out << output_mesh;

  return ;
}
#else
void ElPoly::DefineSurf(){
  printf("DefSurface without CGAL not implemented\n");
}
#endif // __CGAL_h__
#endif //__glut_h__
