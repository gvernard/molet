#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#include <string>
#include "json/json.h"
#include "vkllib.hpp"

class BaseNoise;

class offsetPSF {
public:
  int offset_image;
  int offset_cropped;
  int nj;
  int ni;

  void print();
  void printFrame(FILE* fh,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
};

class Instrument {
public:
  static std::string path;
  std::string name;
  double lambda_min; // in nm
  double lambda_max; // in nm
  double resolution; // in arcsec
  double readout;    // in electrons
  double ZP;
  
  vkl::RectGrid* original_psf = NULL;
  vkl::RectGrid* scaled_psf   = NULL;
  vkl::RectGrid* cropped_psf  = NULL;
  double* kernel         = NULL;
  BaseNoise* noise       = NULL;

  Instrument(std::string name,Json::Value noise_pars);
  Instrument(std::string name,double ZP,Json::Value noise_pars);
  ~Instrument();

  void common_constructor(Json::Value noise_pars);
  static bool checkInstrumentExists(std::string name);
  static void createNewInstrument(Json::Value pars,std::string path_to_psf);
  static double getResolution(std::string name);
  std::string getName();
  void interpolatePSF(vkl::RectGrid* grid);
  void cropPSF(double threshold);
  void createKernel(int Nx,int Ny);
  void convolve(vkl::RectGrid* grid);
  void convolve(vkl::RectGrid* grid_in,vkl::RectGrid* grid_out);
  offsetPSF offsetPSFtoPosition(double x,double y,vkl::RectGrid* grid);

  void replacePSF(std::string path_to_file);
  void preparePSF(vkl::RectGrid* grid,double ratio);
  double sumPSF(offsetPSF* psf_offset);
};

#endif /* INSTRUMENT_HPP */
