#ifndef SOURCE_PLANE_HPP
#define SOURCE_PLANE_HPP

#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <iostream>

/*
#include <boost/geometry.hpp> 
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/algorithms/append.hpp> 
BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
*/

#include "covKernels.hpp"
#include "tableDefinition.hpp"

class Nlpar;
class ImagePlane;
class CollectionMassModels;


class BaseSourcePlane {
public:
  std::string type;                // adaptive or regular
  int Si;                          // pixels in x, if adaptive then Si is set to 0
  int Sj;                          // pixels in y, if adaptive then Sj is set to 0
  int Sm;                          // total number of pixels
  double* src;                     // values of the pixels
  double* x;                       // pixel x-coordinates
  double* y;                       // pixel y-coordinates
  double* s_dx;                    // source derivative in x
  double* s_dy;                    // source derivative in y
  int* mask_vertices;              // source masked pixels indices (0 masked, 1, not masked)
  std::map<int,int> in_total;      // map matching the index of the pixel inside the mask with the one in the total source plane.
  std::map<int,int> total_in;      // the opposite from above
  double* lambda_out;              // regularization parameters lambda_out for the source pixels outside the mask (0 inside the mask)
  int Smask;                       // number of source pixels in the image mask
  std::vector<double> bound_x;     // embedding polygon vertex x-coordinates
  std::vector<double> bound_y;     // embedding polygon vertex y-coordinates
  std::string reg;                 // name of the regularization scheme
  int eigenSparseMemoryAllocForH;  // estimate of the non-zero elements per row of the regularization matrix H
  bool sample_reg = false;         // sampling regularization matrix related parameters
  double lambda_out_sum;           // the sum of the log (ln) entries of lambda_out
  BaseCovKernel* kernel;           // pointer to kernel class
  mytable L;
  mytable H;
  mytable Ds;
  

  BaseSourcePlane(){};
  BaseSourcePlane(const BaseSourcePlane& other);
  ~BaseSourcePlane(){
    free(src);
    free(x);
    free(y);
    free(s_dx);
    free(s_dy);
    free(mask_vertices);
    free(lambda_out);
    if( this->reg == "covariance_kernel" ){
      delete this->kernel;
    }
  };


  //virtual members
  virtual BaseSourcePlane* clone() = 0;
  virtual void constructH() = 0;
  virtual void outputSource(const std::string path) = 0;
  virtual void outputSourceErrors(double* errors,const std::string path) = 0;
  virtual void constructDerivatives() = 0;             // constructs the derivatives of the source on the source grid
  virtual void constructDs(ImagePlane* image) = 0;     // interpolates the derivatives of the source on the deflected image grid
  virtual void createInterpolationWeights(ImagePlane* image) = 0;
  virtual void inMask(ImagePlane* image) = 0;
  virtual void outputMask(const std::string path) = 0;

  //non-virtual members
  void normalize();
  void constructL(ImagePlane* image);
  void setSourceCovariance(std::vector<Nlpar> reg_pars){};
};


class FixedSource: public BaseSourcePlane {
public:
  FixedSource(int source_i,int source_j,double size,std::string reg_scheme);
  FixedSource(int i,int j,double width,double height,std::string reg_scheme);
  FixedSource(int i,int j,double xmin,double xmax,double ymin,double ymax,std::string reg_scheme);
  FixedSource(const FixedSource& source) : BaseSourcePlane(source) {};
  virtual FixedSource* clone(){
    return new FixedSource(*this);
  };

  virtual void createInterpolationWeights(ImagePlane* image);
  virtual void constructH();
  virtual void constructDerivatives(){};
  virtual void constructDs(ImagePlane* image){};
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path);
  virtual void inMask(ImagePlane* image);
  virtual void outputMask(const std::string path);

  void setGridRect(double width,double height);
  void setGridRect(double xmin,double xmax,double ymin,double ymax);
  void boundPolygon();
  bool pointInPolygon(double x,double y);

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double width;
  double height;
};

/*
class FloatingSource: public BaseSourcePlane {
public:

  //virtual members
  FloatingSource(int source_i,int source_j,double size,double x0,double y0,std::string reg_scheme);
  void createInterpolationWeights(ImagePlane* image);
  void constructH();
  void constructDs(ImagePlane* image){};
  void outputSource(const std::string path);
  void outputSourceErrors(double* errors,const std::string path);

  //non-virtual members
  void setGrid(std::map<std::string,std::string> pars);
  void boundPolygon();
  bool pointInPolygon(double x,double y);

  double x0;
  double y0;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};
*/

class AdaptiveSource: public BaseSourcePlane {
public:
  AdaptiveSource(int Sm,std::string reg_scheme);
  AdaptiveSource(std::string mode,int Sm,int spacing,std::string reg_scheme);
  AdaptiveSource(const AdaptiveSource& source);
  virtual AdaptiveSource* clone(){
    return new AdaptiveSource(*this);
  };
  ~AdaptiveSource();

  virtual void createInterpolationWeights(ImagePlane* image);
  virtual void constructH();
  virtual void constructDerivatives();
  virtual void constructDs(ImagePlane* image);
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path);
  virtual void inMask(ImagePlane* image);
  virtual void outputMask(const std::string path){};

  void createAdaGrid(ImagePlane* image,CollectionMassModels* mycollection);
  void createDelaunay();
  void writeTriangles();
  void setGrid(std::map<std::string,std::string> pars){};
  bool pointInPolygon(double x,double y){};
  void writeVertexFaces(int index,std::string output);
  //  void boundPolygon();

private:

  struct xypoint {
    double x;
    double y;
  };
  
  struct a_triangle {
    int a;
    int b;
    int c;
  };

  xypoint intersection_point_x(xypoint p0,xypoint p1,xypoint p2);
  xypoint intersection_point_y(xypoint p0,xypoint p1,xypoint p2);
  double triangleArea(a_triangle triangle);

  std::string mode;
  int spacing;
  int n_triangles;
  std::vector<a_triangle> triangles;
  std::vector< std::vector<int> > opposite_edges_per_vertex;

  std::vector<int> mask_triangles;
};



class FactorySourcePlane{//This is a singleton class.
public:
  FactorySourcePlane(FactorySourcePlane const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactorySourcePlane const&) = delete;

  static FactorySourcePlane* getInstance(){
    static FactorySourcePlane dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseSourcePlane* createSourcePlane(std::map<std::string,std::string> source){
    if( source["type"] == "fixed" ){
      return new FixedSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),source["reg_s"]);
      //    } else if( source["type"] == "floating" ){
      //      return new FloatingSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),stof(source["x0"]),stof(source["y0"]),source["reg"]);
    } else if( source["type"] == "adaptive" ){
      if( source["mode"] == "random" ){
	return new AdaptiveSource(stoi(source["sm"]),source["reg_s"]);
      } else if( source["mode"] == "image" || source["mode"] == "grid" ){
	return new AdaptiveSource(source["mode"],stoi(source["sm"]),stoi(source["spacing"]),source["reg_s"]);
      }
    } else {
      return NULL;
    }
  }

private:
  FactorySourcePlane(){};
};


#endif /* SOURCE_PLANE_HPP */
