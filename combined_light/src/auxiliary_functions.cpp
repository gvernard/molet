#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "auxiliary_functions.hpp"


// START:LIGHTCURVE ========================================================================================
LightCurve::LightCurve(const Json::Value lc){
  for(int t=0;t<lc["time"].size();t++){
    this->time.push_back( lc["time"][t].asDouble() );
    this->signal.push_back( lc["signal"][t].asDouble() );
    this->dsignal.push_back( lc["dsignal"][t].asDouble() );
  }
}

LightCurve::LightCurve(std::vector<double> t){
  this->time = t;
  this->signal.resize(t.size());
  this->dsignal.resize(t.size());
}

LightCurve::LightCurve(std::vector<double> t,std::vector<double> s,std::vector<double> ds){
  this->time   = t;
  this->signal = s;
  this->signal = ds;
}

// I think the following function is not used somewhere anymore
void LightCurve::interpolate(std::vector<double> obs_time,double delay,double* interpolated){
  // Assumption 1: obs_time and interpolated have the same length
  // Assumption 2: the light curve is at least as long as obs_time (in actual time values)
  int i = 0;
  int j  = 1;
  while( i < obs_time.size() ){
    double t = obs_time[i] + delay;
    if( this->time[j] > t ){
      //interpolate
      interpolated[i] = this->signal[j-1]+(t-this->time[j-1])*(this->signal[j]-this->signal[j-1])/(this->time[j]-this->time[j-1]);
      i++;
    } else {
      j++;
    }
  }
}

void LightCurve::interpolate(LightCurve* int_lc,double delay){
  // Assumption 1: obs_time and interpolated have the same length
  // Assumption 2: the light curve is at least as long as obs_time (in actual time values)
  int i = 0;
  int j = 1;
  while( i < int_lc->time.size() ){
    double t = int_lc->time[i] + delay;
    if( this->time[j] > t ){
      //interpolate
      int_lc->signal[i]  = this->signal[j-1]+(t-this->time[j-1])*(this->signal[j]-this->signal[j-1])/(this->time[j]-this->time[j-1]);
      int_lc->dsignal[i] = this->dsignal[j-1]+(t-this->time[j-1])*(this->dsignal[j]-this->dsignal[j-1])/(this->time[j]-this->time[j-1]);
      i++;
    } else {
      j++;
    }
  }
}

Json::Value LightCurve::jsonOut(){
  Json::Value json_lc,json_time,json_signal,json_dsignal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(this->signal[i]);
    json_dsignal.append(this->dsignal[i]);
  }
  json_lc["time"]    = json_time;
  json_lc["signal"]  = json_signal;
  json_lc["dsignal"] = json_dsignal;
  return json_lc;
}

Json::Value LightCurve::jsonOutMag(){
  Json::Value json_lc,json_time,json_signal,json_dsignal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(-2.5*log10(this->signal[i]));
    json_dsignal.append( -2.5*log10(sqrt(1.0-pow(this->dsignal[i]/this->signal[i],2))) );
    // This is a symmetric error derived from the asymmetric logarithmic errors as: dm = [(m+)+(m-)]/2
  }
  json_lc["time"]    = json_time;
  json_lc["signal"]  = json_signal;
  json_lc["dsignal"] = json_dsignal;
  return json_lc;
}
// END:LIGHTCURVE ========================================================================================



Json::Value readLightCurvesJson(std::string lc_type,std::string type,std::string instrument_name,std::string in_path,std::string out_path){
  std::ifstream fin;
  Json::Value lcs_json;
  if( type == "custom" ){
    fin.open(in_path+"/input_files/"+instrument_name+"_LC_"+lc_type+".json",std::ifstream::in);
  } else {
    fin.open(out_path+"output/"+instrument_name+"_LC_"+lc_type+".json",std::ifstream::in);
  }
  fin >> lcs_json;
  fin.close();
  return lcs_json;
}

std::vector<LightCurve*> conversions(Json::Value lcs_json,double zs,double scale){
  int N = lcs_json.size();
  double fac = 1.0 + zs;
  std::vector<LightCurve*> lcs(N);
  for(int i=0;i<N;i++){
    lcs[i] = new LightCurve(lcs_json[i]);    
    for(int j=0;j<lcs[i]->signal.size();j++){
      lcs[i]->signal[j] = scale*pow(10.0,-0.4*lcs[i]->signal[j]); // Convert from magnitudes to intensities and scale by a factor if necessary
      lcs[i]->time[j] *= fac; // Convert time to the observer's frame
    }
  }	
  return lcs;
}

void combineInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target){
  // === Combining two signals: intrinsic and extrinsic
  LightCurve* base = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  LC_extrinsic->interpolate(target,0.0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = target->signal[t]*macro_mag*base->signal[t];
    target->dsignal[t] = target->dsignal[t]*macro_mag*base->signal[t]; // no error from macro-mag, intrinsic, and unmicrolensed flux, only from the microlensing mag
  }
  delete(base);
}

void combineInExUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* LC_unmicro,LightCurve* target){
  // === Combining three signals: intrinsic, intrinsic unmicrolensed, and extrinsic
  LightCurve* base = new LightCurve(time);
  LightCurve* base_unmicro = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  LC_unmicro->interpolate(base_unmicro,td);
  LC_extrinsic->interpolate(target,0.0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*(target->signal[t]*base->signal[t] + base_unmicro->signal[t]);
    target->dsignal[t] = macro_mag*target->dsignal[t]*base->signal[t]; // no error from macro-mag, intrinsic, and unmicrolensed flux, only from the microlensing mag    
  }
  delete(base);
  delete(base_unmicro);
}

void combineInUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_unmicro,LightCurve* target){
  // === Combining two signals: intrinsic and unmicrolensed
  LightCurve* base = new LightCurve(time);
  LightCurve* base_unmicro = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  LC_unmicro->interpolate(base_unmicro,0.0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*(base->signal[t] + base_unmicro->signal[t]);
    target->dsignal[t] = 0.0;
  }
  delete(base);
  delete(base_unmicro);
}

void justInSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* target){
  // === Just one signal: intrinsic
  LightCurve* base = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*base->signal[t]; // this line includes only the intrinsic signal and excludes microlensing
    target->dsignal[t] = 0.0;
  }
  delete(base);
}

void outputLightCurvesJson(std::vector<LightCurve*> lcs,std::string filename){
  Json::Value lcs_json;
  for(int q=0;q<lcs.size();q++){
    lcs_json.append( lcs[q]->jsonOutMag() );
  }
  std::ofstream lcs_file(filename);
  lcs_file << lcs_json;
  lcs_file.close();									      
}





// START:TRANSFORM PSF =====================================================================================
TransformPSF::TransformPSF(double a,double b,double c,bool d,bool e){
  this->x0 = a;
  this->y0 = b;
  this->rot = c * 0.01745329251;// converting to rad
  this->flip_x = d;
  this->flip_y = e;

  this->cosrot = cos(this->rot);
  this->sinrot = sin(this->rot);
}

void TransformPSF::applyTransform(double xin,double yin,double& xout,double& yout){
  xout =  (xin-this->x0)*this->cosrot + (yin-this->y0)*this->sinrot;
  yout = -(xin-this->x0)*this->sinrot + (yin-this->y0)*this->cosrot;
  if( this->flip_x ){
    xout *= -1.0;
  }
  if( this->flip_y ){
    yout *= -1.0;
  }
}

double TransformPSF::interpolateValue(double x,double y,PSF* mypsf){
  return 1.0;
}
// END:TRANSFORM PSF =====================================================================================
