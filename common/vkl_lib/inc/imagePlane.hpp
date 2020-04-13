#ifndef IMAGE_PLANE_HPP
#define IMAGE_PLANE_HPP

#include <valarray>
#include <string>
#include <vector>
#include <map>

#include "tableDefinition.hpp"

class Cross {
public:
  int i0;
  int j0;
  double* coeff_x;
  double* coeff_y;

  Cross(int a,int b){
    this->i0 = a;
    this->j0 = b;
    int n = 8;
    this->coeff_x = (double*) malloc(n*sizeof(double));
    this->coeff_y = (double*) malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
      this->coeff_x[i] = 0.0;
      this->coeff_y[i] = 0.0;
    }
  }
  ~Cross(){
    free(coeff_x);
    free(coeff_y);
  }
};


class ImagePlane {
public:
  int Ni;                    //pixels in x direction
  int Nj;                    //pixels in y direction
  int Nm;                    //total pixels in the image data
  //  std::map<int,int> lookup;  //matching the indices of the data pixels to the indices of the mask pixels
  int Nmask;                 //pixels in the mask
  double width;              //in arcsec
  double height;             //in arcsec
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double* img;               //values of the pixels
  double* x;                 //pixel x-coordinates in arcsec
  double* y;                 //pixel y-coordinates in arcsec
  double* defl_x;            //deflected x coordinates
  double* defl_y;            //deflected y coordinates
  int* active;               //active image pixels used in the construction of the adaptive grid
  InterpolationCell** cells = 0;    //for each image pixel, the corresponding source pixels and interpolation weights
  Cross** crosses = 0;       //array of cross structures for creating the D_s(s_p)*D_psi for reconstructing perturbations
  InterpolationCell** dpsi_cells = 0; //for each image pixel, the corresponding dpsi pixels and interpolation weights, if applicable
  mytable B;
  mytable C;
  mytable S;
  std::string noise_flag;

  ImagePlane(){};
  ImagePlane(const std::string filepath,int i,int j,double w,double h);                      //used to read images
  ImagePlane(const std::string filepath,int i,int j,double w,double h,double x0,double y0);  //used to read images
  ImagePlane(int i,int j,double w,double h);                                                 //used to create images
  ImagePlane(int i,int j,double w,double h,double x0,double y0);                             //used to create images
  ImagePlane(int i,double w,double h);                                                       //used to create masked images
  ImagePlane(const ImagePlane& image);
  ~ImagePlane();

  void writeImage(const std::string filename);
  void writeBin(const std::string filename);
  void readFits(const std::string filename,std::valarray<float>& contents);
  void readB(const std::string filepath,int i,int j,int ci,int cj);
  void readC(const std::string flag,const std::string filepath);
  void readS(const std::string filepath);
  //  void setMaskedC(mytable* Cout,mytable* S,mytable* C);
  void maskData(std::map<int,int> lookup,ImagePlane* masked);
  //  void printCross(int k,mytable Ds);
  void printCross(int k);
  void lowerResRebinAdditive(ImagePlane* newImage);
  void lowerResRebinIntegrate(ImagePlane* newImage);

  
private:
  void setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
  void setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
};

#endif /* IMAGE_PLANE_HPP */
