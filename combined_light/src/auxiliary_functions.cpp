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
  }
}

LightCurve::LightCurve(std::vector<double> t){
  this->time = t;
  this->signal.resize(t.size());
}

LightCurve::LightCurve(std::vector<double> t,std::vector<double> s){
  this->time   = t;
  this->signal = s;
}

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
  int j  = 1;
  while( i < int_lc->time.size() ){
    double t = int_lc->time[i] + delay;
    if( this->time[j] > t ){
      //interpolate
      int_lc->signal[i] = this->signal[j-1]+(t-this->time[j-1])*(this->signal[j]-this->signal[j-1])/(this->time[j]-this->time[j-1]);
      i++;
    } else {
      j++;
    }
  }
}

Json::Value LightCurve::jsonOut(){
  Json::Value json_lc,json_time,json_signal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(this->signal[i]);
  }
  json_lc["time"]   = json_time;
  json_lc["signal"] = json_signal;
  return json_lc;
}

Json::Value LightCurve::jsonOutMag(){
  Json::Value json_lc,json_time,json_signal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(-2.5*log10(this->signal[i]));
  }
  json_lc["time"]   = json_time;
  json_lc["signal"] = json_signal;
  return json_lc;
}
// END:LIGHTCURVE ========================================================================================




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
