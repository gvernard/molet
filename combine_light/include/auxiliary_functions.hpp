#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>
#include <vector>

#include "json/json.h"
#include "vkllib.hpp"

class PSF;
class Instrument;
class offsetPSF;

class LightCurve {
public:
  std::vector<double> time;
  std::vector<double> signal;
  std::vector<double> dsignal;
  
  LightCurve(){};
  LightCurve(const Json::Value lc);
  LightCurve(std::vector<double> time);
  LightCurve(std::vector<double> time,std::vector<double> signal,std::vector<double> dsignal);
  
  void interpolate(std::vector<double> obs_time,double delay,double* interpolated);
  void interpolate(LightCurve* int_lc,double delay);
  Json::Value jsonOut();
  Json::Value jsonOutMag();
};


Json::Value readLightCurvesJson(std::string lc_type,std::string type,std::string instrument_name,std::string in_path,std::string out_path);
std::vector<LightCurve*> conversions(Json::Value lcs_json,double zs,double scale,double ZP);
void extendIntrinsicLightCurve(LightCurve* lc_in,LightCurve* lc_out,double td,double time_start,double time_end);
void extendExtrinsicLightCurve(LightCurve* lc_ex,LightCurve* lc_out,double td,double time_start,double time_end,double time_start_in);
void combineSupernovaInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target);
void justSupernovaInSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* target);
void combineInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target);
void combineInExUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* LC_unmicro,LightCurve* target);
void combineInUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_unmicro,LightCurve* target);
void justOneSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC,LightCurve* target);
void outputLightCurvesJson(std::vector<LightCurve*> lcs,std::string filename);

vkl::RectGrid createObsBase(Instrument* mycam,vkl::RectGrid* supersim,int res_x,int res_y,std::string out_path);
vkl::RectGrid createPointSourceLight(vkl::RectGrid* supersim,std::vector<double> image_signal,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y);
void writeCutout(std::string cutout_scale,vkl::RectGrid* obs_pp_light,vkl::RectGrid* obs_base,std::string fname);
void writeAllCutouts(std::vector<double> tobs,Json::Value images,std::vector<LightCurve*> samp_LC,vkl::RectGrid* supersim,Instrument* mycam,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y,std::string mock,std::string instrument_name,std::string cutout_scale,std::string out_path);

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
