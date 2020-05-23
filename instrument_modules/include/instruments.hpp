#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#include <string>
#include "json/json.h"

class ImagePlane;
class BaseNoise;

class offsetPSF {
public:
  int offset_image;
  int offset_cropped;
  int nj;
  int ni;

  void print();
  void printFrame(FILE* fh,int Ni,int Nj,double w,double h);
};

class Instrument {
public:
  static std::string path;
  std::string name;
  double lambda_min; // in nm
  double lambda_max; // in nm
  double resolution; // in arcsec
  ImagePlane* original_psf = NULL;
  ImagePlane* scaled_psf   = NULL;
  ImagePlane* cropped_psf  = NULL;
  double* kernel           = NULL;
  BaseNoise* noise         = NULL;
  
  Instrument(std::string name,Json::Value noise_pars);
  ~Instrument();

  static double getResolution(std::string name);
  std::string getName();
  void interpolatePSF(ImagePlane* image);
  void cropPSF(double threshold);
  void createKernel(int Ni,int Nj);
  void convolve(ImagePlane* image);
  offsetPSF offsetPSFtoPosition(double x,double y,ImagePlane* image);
};

#endif /* INSTRUMENT_HPP */
