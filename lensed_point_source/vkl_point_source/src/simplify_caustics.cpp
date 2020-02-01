#include <vector>

#include "contour-tracing.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace PS = CGAL::Polyline_simplification_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polygon_2<K>                                   Polygon_2;
typedef K::Point_2                                           Point;
typedef Polygon_2::Vertex_iterator                           VertexIterator;
typedef PS::Squared_distance_cost                            Cost1;
typedef PS::Scaled_squared_distance_cost                     Cost2;
typedef PS::Stop_below_count_ratio_threshold                 Stop;
typedef PS::Stop_above_cost_threshold                        Stop2;


void simplifyPolygon(std::vector<Contour*>& caustics,std::vector<Contour*>& simplified){

  for(int i=0;i<caustics.size();i++){
    Polygon_2 p;
    for(int j=0;j<caustics[i]->x.size();j++){
      p.push_back(Point(caustics[i]->x[j],caustics[i]->y[j]));
    }
    //    CGAL::set_pretty_mode(std::cout);
    //    std::cout << "created the polygon p:" << std::endl;
    //    std::cout << p << std::endl;
    //    std::cout << std::endl;
  
    Cost1 cost;
    p = PS::simplify(p,cost,Stop(0.1));
    //    CGAL::set_pretty_mode(std::cout);
    //    std::cout << "created the polygon p:" << std::endl;
    //    std::cout << p << std::endl;
    //    std::cout << std::endl;
    
    for(VertexIterator vi=p.vertices_begin();vi!=p.vertices_end();vi++){
      simplified[i]->x.push_back(vi->x());
      simplified[i]->y.push_back(vi->y());
    }
    simplified[i]->x.push_back(p.vertices_begin()->x());
    simplified[i]->y.push_back(p.vertices_begin()->y());
    
    
  }

}
