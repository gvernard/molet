#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>
#include <vector>

#include "json/json.h"
#include "vkllib.hpp"

class PSF;
class Instrument;
class offsetPSF;

vkl::RectGrid createObsStatic(Instrument* mycam,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,int res_x,int res_y,double& F_conv_extended,double& F_conv_lens,bool convolve_lens);
vkl::RectGrid createObsPS(vkl::RectGrid* supersim,std::vector<double> image_signal,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y);
void writeCutout(vkl::RectGrid* obs,std::string fname);
void writeAllCutouts(std::vector<double> tobs,Json::Value images,Json::Value samp_LC,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,Instrument* mycam,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y,std::string mock,bool convolve_lens,std::string out_path);
void convertGridFromFlux(vkl::RectGrid* obs,double ZP);
double convertFromFlux(double flux,double ZP);

// user provided function for getting the file names of a time dependent PSF
std::vector<std::string> getFileNames(std::vector<double> tobs,std::string path);


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


#endif /* AUXILIARY_HPP */
