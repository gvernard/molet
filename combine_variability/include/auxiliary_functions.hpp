#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <string>
#include <vector>

#include "json/json.h"


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
  Json::Value jsonOutMag(double ZP);
};


Json::Value readLightCurvesJson(std::string lc_type,std::string type,std::string instrument_name,std::string in_path,std::string out_path);
std::vector<LightCurve*> conversions(Json::Value lcs_json,double zs,double scale,double ZP);
std::vector<LightCurve*> conversions_SN(Json::Value lcs_json,double zs,double scale,double ZP,double tobs_min,double td_max);
void tailorLightCurve(LightCurve* lc_input,LightCurve* lc_output,double td,double tstart,double tend,double floor_val);
void combineSupernovaInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target);
void justSupernovaInSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* target);
void combineInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target);
void combineInExUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* LC_unmicro,LightCurve* target,double unmicro_ratio);
void combineInUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_unmicro,LightCurve* target,double unmicro_ratio);
void justOneSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC,LightCurve* target);
void outputLightCurvesJson(std::vector<LightCurve*> lcs,double ZP,std::string filename);



class BaseLagKernel {
public:
  virtual void getKernel(std::vector<double> time,std::vector<double>& kernel) = 0;
};


class DeltaKernel : public BaseLagKernel {
public:
  double t_peak; // in days
  DeltaKernel(double t_peak) : t_peak(t_peak) {};
  void getKernel(std::vector<double> time,std::vector<double>& kernel);
};

class TopHatKernel : public BaseLagKernel {
public:
  double t_limit; // in days
  TopHatKernel(double radius); // radius must be in 10^14 cm
  void getKernel(std::vector<double> time,std::vector<double>& kernel);  
private:
  double speed_of_light = 25.9; // in 10^14 cm /day
};


#endif /* AUXILIARY_HPP */
