/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
/* H0 X                                                                            */
/* H0 X   libAtoms+QUIP: atomistic simulation library                              */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by                                    */
/* H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,      */
/* H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.          */
/* H0 X                                                                            */
/* H0 X   Copyright 2006-2010.                                                     */
/* H0 X                                                                            */
/* H0 X   These portions of the source code are released under the GNU General     */
/* H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html          */
/* H0 X                                                                            */
/* H0 X   If you would like to license the source code under different terms,      */
/* H0 X   please contact Gabor Csanyi, gabor@csanyi.net                            */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by Noam Bernstein as part of          */
/* H0 X   his employment for the U.S. Government, and are not subject              */
/* H0 X   to copyright in the USA.                                                 */
/* H0 X                                                                            */
/* H0 X                                                                            */
/* H0 X   When using this software, please cite the following reference:           */
/* H0 X                                                                            */
/* H0 X   http://www.libatoms.org                                                  */
/* H0 X                                                                            */
/* H0 X  Additional contributions by                                               */
/* H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras                 */
/* H0 X                                                                            */
/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

/* This file contains routines which use the Computational Geometry
   Algorithms Library (CGAL, http://www.cgal.org) to compute alpha
   shapes. It is used to determine a crack front given the set of 
   crack surface atoms. */

#define _USE_MATH_DEFINES
#include <math.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/bounding_box.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <boost/iterator/transform_iterator.hpp>

extern "C" {
#include "libatoms.h"
}

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef K::FT FT;

typedef K::Point_2  Point;
typedef K::Segment_2  Segment;

typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;

typedef Alpha_shape_2::Face  Face;
typedef Alpha_shape_2::Vertex Vertex;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Face_handle  Face_handle;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef Alpha_shape_2::Face_circulator  Face_circulator;
typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape_2::Locate_type Locate_type;

typedef Alpha_shape_2::Face_iterator  Face_iterator;
typedef Alpha_shape_2::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator  Edge_iterator;
typedef Alpha_shape_2::Edge_circulator  Edge_circulator;

typedef Alpha_shape_2::Alpha_iterator Alpha_iterator;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_shape_vertices_iterator;

typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

std::map<Point,int>::key_type get_key(std::map<Point,int>::value_type aPair) {
  return aPair.first;
}

typedef std::map<Point,int>::key_type (*get_key_t)(std::map<Point,int>::value_type);
typedef std::map<Point,int>::iterator map_iterator;
typedef boost::transform_iterator<get_key_t, map_iterator> mapkey_iterator;


extern "C" void c_alpha_shape_2(int *n, double x[], double y[],
			      double *alpha, int *shape_n, int shape_list[],
			      int *error)
{
  INIT_ERROR;

  // Map from points to original indices in input data
  std::map<Point,int> point_map;

  for (int i=0; i<*n; i++) {
    point_map[Point(x[i],y[i])] = i;
  }

  mapkey_iterator points_begin(point_map.begin(), get_key);
  mapkey_iterator points_end(point_map.end(), get_key);
  
  Alpha_shape_2 A(points_begin, points_end, FT(*alpha), Alpha_shape_2::GENERAL);

  Arrangement_2 arr;

  // Traverse alpha shape edges and add to arrangement
  for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
      it != A.alpha_shape_edges_end();
      ++it){
    CGAL::insert(arr, A.segment(*it));
  }

  // Check we have one bounded face
  if (arr.number_of_faces() != 2) {
    RAISE_ERROR("c_alpha_shape_2: number of bounded faces != 1, try increasing alpha");
  }

  Arrangement_2::Face_const_iterator fi = arr.faces_begin();
  ++fi; // skip over the unbounded face which all arrangements have
 
  Arrangement_2::Ccb_halfedge_const_circulator ccb = fi->outer_ccb();

  // Check shape_list is big enough to store alpha shape
  int tmp_shape_length = 0;
  do {
    tmp_shape_length++;
    ++ccb;
  } while (ccb != fi->outer_ccb());
  if (tmp_shape_length > *shape_n) {
    RAISE_ERROR("c_alpha_shape_2: size of alpha shape exceeds length of shape_list");
  }

  // traverse outer CCB again, adding results to shape_list
  ccb = fi->outer_ccb();
  *shape_n = 0;
  do {
    shape_list[(*shape_n)++] = point_map[ccb->source()->point()];
    ++ccb;
  } while (ccb != fi->outer_ccb());
}

extern "C" void c_crack_front_alpha_shape(int *n, double x[], double y[], 
					  double *alpha, double *angle_threshold, 
					  int *front_n, int front_list[], int *error)
{
  INIT_ERROR;

  // Map from points to original indices in input data
  std::map<Point,int> point_map;

  for (int i=0; i<*n; i++) {
    point_map[Point(x[i],y[i])] = i;
  }

  mapkey_iterator points_begin(point_map.begin(), get_key);
  mapkey_iterator points_end(point_map.end(), get_key);
  
  Alpha_shape_2 A(points_begin, points_end, FT(*alpha), Alpha_shape_2::GENERAL);

  Arrangement_2 arr;

  // Traverse alpha shape edges and add to arrangement
  for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
      it != A.alpha_shape_edges_end();
      ++it){
    CGAL::insert(arr, A.segment(*it));
  }

  // Check we have one bounded face
  if (arr.number_of_faces() != 2) {
    RAISE_ERROR("c_crack_front_alpha_shape: number of bounded faces != 1, try increasing alpha");
  }

  // Compute bounding box of original point set and use it
  // to perform a vertical ray shooting query
  K::Iso_rectangle_2 bbox = bounding_box(points_begin, points_end);
  Point p = Point(bbox.xmax(), bbox.ymin());
  Walk_pl walk_pl(arr);
  CGAL::Object    obj = walk_pl.ray_shoot_up (p);
  Arrangement_2::Vertex_const_handle xmax_vertex;
  if (!CGAL::assign (xmax_vertex, obj)) {
    RAISE_ERROR("c_crack_front_alpha_shape: vertical line search did not hit a vertex");
  }
    
  Arrangement_2::Face_const_iterator fi = arr.faces_begin();
  ++fi; // skip over the unbounded face which all arrangements have

  // find half-line with source at xmax_vertex
  Arrangement_2::Ccb_halfedge_const_circulator curr, xmax_edge;
  for(xmax_edge = fi->outer_ccb(); xmax_edge->source() != xmax_vertex; ++xmax_edge);

  double bearing;
  Point p1, p2;
  std::list<int> front;

  // start at xmax_edge and traverse outer CCB until angle along line 
  // exceeds angle_threshold
  curr = xmax_edge;
  do {
    p1 = curr->source()->point();
    p2 = curr->target()->point();
    
    bearing = 180./M_PI*std::atan2(CGAL::to_double(p2.y())-CGAL::to_double(p1.y()),
				   CGAL::to_double(p2.x())-CGAL::to_double(p1.x()));

    front.push_back(point_map[curr->source()->point()]);
  } while (++curr != xmax_edge && bearing < *angle_threshold);

  // now go back to xmax_edge, and go around backwards until
  // angle is less than -angle_theshold
  curr = xmax_edge;
  --curr; // don't want xmax_edge vertex to appear twice in output list
  do {
    p1 = curr->twin()->source()->point();
    p2 = curr->twin()->target()->point();
    
    bearing = 180./M_PI*std::atan2(CGAL::to_double(p2.y())-CGAL::to_double(p1.y()),
				   CGAL::to_double(p2.x())-CGAL::to_double(p1.x()));

    front.push_front(point_map[curr->source()->point()]);
  } while (--curr != xmax_edge && bearing > -*angle_threshold);

  // Copy output into array front_list
  if (front.size() > *front_n) {
    RAISE_ERROR("c_crack_front_alpha_shape: Not enough space in front_list");
  }
  *front_n=0;
  for (std::list<int>::iterator it = front.begin(); it != front.end(); ++it) {
    front_list[(*front_n)++] = *it;
  }
}

#ifdef MAIN_PROGRAM
int main()
{
  FILE *stream = fopen("pos_xz","r");
  char buffer[256];
  int i,n;

  fgets(buffer,256,stream);
  sscanf(buffer,"%d",&n);

  double *x = new double[n];
  double *y = new double[n];

  for (i=0; i<n; i++) {
    fgets(buffer,256,stream);
    sscanf(buffer,"%lf %lf",&(x[i]),&(y[i]));
  }
  fclose(stream);

  double alpha = 100.0;
  double angle_threshold = 160.0;
  
  int front_n = 100;
  int *front_list = new int[front_n];

  int error = 0;
  crack_front_alpha_shape(&n, x, y, &alpha, &angle_threshold, &front_n, front_list, &error);

  if (error == 0) {
    for (i = 0; i < front_n; i++)
      std::cout << front_list[i] << std::endl;
  }			       
}
#endif
