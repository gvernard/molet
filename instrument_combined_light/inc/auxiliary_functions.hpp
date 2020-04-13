#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>
#include "json/json.h"

class ImagePlane;

class offsetPSF {
public:
  int offset_image;
  int offset_cropped;
  int nj;
  int ni;
};


class PSF {
public:
  ImagePlane* original_psf = NULL;
  ImagePlane* scaled_psf   = NULL;
  ImagePlane* cropped_psf  = NULL;
  double* kernel           = NULL;
  
  PSF(std::string fname,int pix_x,int pix_y,double width,double height,ImagePlane* mydata){
    this->original_psf = new ImagePlane(fname,pix_x,pix_y,width,height);
    interpolatePSF(mydata);
  }

  ~PSF(){
    delete(original_psf);
    delete(scaled_psf);
    delete(cropped_psf);
    free(kernel);
  }

  void interpolatePSF(ImagePlane* image);
  void cropPSF(double threshold);
  void createKernel(int Ni,int Nj);
  void convolve(ImagePlane* image);
  offsetPSF offsetPSFtoPosition(double x,double y,ImagePlane* image);
};


class TransformPSF {
public:
  double x0; // in arcsec
  double y0; // in arcsec
  double rot; // degrees normal cartesian
  bool flip_x;
  bool flip_y;

  TransformPSF(double a,double b,double c,bool d,bool e);
  
  void applyTransform(double xin,double yin,double& xout,double& yout);
  double interpolateValue(double x,double y,PSF* mypsf);

private:
  double cosrot;
  double sinrot;
};


class LightCurve {
public:
  std::vector<double> time;
  std::vector<double> signal;
  
  LightCurve(){};
  LightCurve(const Json::Value lc);
  LightCurve(std::vector<double> time);
  LightCurve(std::vector<double> time,std::vector<double> signal);
  
  void interpolate(std::vector<double> obs_time,double delay,double* interpolated);
  void interpolate(LightCurve* int_lc,double delay);
  Json::Value jsonOut();
  Json::Value jsonOutMag();
};


void outputLightCurvesJson(std::vector<LightCurve*> lcs,std::string filename);
  

#endif /* AUXILIARY_HPP */
