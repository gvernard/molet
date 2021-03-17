#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "json/json.h"

#include "instruments.hpp"
#include "auxiliary_functions.hpp"

int main(int argc,char* argv[]){

  /*
    Requires:
    - multiple_images.json
  */

  //=============== BEGIN:INITIALIZE =======================================================
  std::ifstream fin;  

  // Read the main parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();
  const Json::Value intrinsic_pars = root["point_source"]["variability"]["intrinsic"];
  int N_in = intrinsic_pars["N_in"].asInt();

  // Read the multiple images' parameters from JSON
  std::string out_path = argv[2];
  Json::Value images;
  fin.open(out_path+"output/multiple_images.json",std::ifstream::in);
  fin >> images;
  fin.close();

  double zs = root["source"]["redshift"].asDouble();
  double Mi = intrinsic_pars["absolute_i_mag"].asDouble();
  //========================================================================================



  //=============== SORT INSTRUMENTS BY WAVELENGTH AND GET THE SMALLEST ONE =================
  // Get the rest wavelength for each instrument
  std::map<std::string,double> tmp;
  for(int b=0;b<root["instruments"].size();b++){
    const Json::Value instrument = root["instruments"][b];
    std::string instrument_name = root["instruments"][b]["name"].asString();
    Instrument mycam(instrument_name,root["instruments"][b]["noise"]);
    double lobs  = (mycam.lambda_min + mycam.lambda_max)/2.0;
    tmp.insert( std::pair<std::string,double>(instrument_name,lobs/(1.0+zs)) );
  }  
  std::map<double,std::string> lambda_rest = flip_map(tmp);
  std::string min_instrument_name = lambda_rest.begin()->second;
  double lambda_min = lambda_rest.begin()->first;
  // Now lambda_rest is sorted in ascending wavelength

  // Get the mean magnitude for each instrument
  std::map<std::string,double> mean_mags;
  for(Json::Value::const_iterator it=intrinsic_pars["mean_mag"].begin();it!=intrinsic_pars["mean_mag"].end();it++){
    mean_mags.insert( std::pair<std::string,double>(it.key().asString(),it->asDouble()) );
  }
  double mean_lambda_min = mean_mags[min_instrument_name];
  //=========================================================================================


  
  //=============== CREATE THE TIME VECTOR IN THE SOURCE REFERENCE FRAME ====================
  // Get maximum image time delay
  double td_max = 0.0;
  for(int q=0;q<images.size();q++){
    double td = images[q]["dt"].asDouble();
    if( td > td_max ){
      td_max = td;
    }
  }
  
  // Get time duration for the intrinsic light curve
  double t;
  double tmin = root["instruments"][0]["time"][0].asDouble();
  double tmax = root["instruments"][0]["time"][0].asDouble();
  for(int i=0;i<root["instruments"].size();i++){
    t = root["instruments"][i]["time"][0].asDouble();
    if( tmin > t ){
      tmin = t;
    }
    int last_t = root["instruments"][i]["time"].size() - 1;
    t = root["instruments"][i]["time"][last_t].asDouble();
    if( tmax < t ){
      tmax = t;
    }
  }

  // Extend time window by the maximum time delay and 10% of the provided duration in both directions
  double duration = tmax - tmin;
  tmax = tmax + td_max + 0.1*duration;
  tmin = tmin - 0.1*duration;
  if( tmin < 0 ){
    tmin = 0;
  }
  //std::cout << tmin << " " << tmax << std::endl;

  // Get time vector in the source reference frame
  // I add a single day to tmin every time in the observer's frame and then convert to the source reference frame.
  // The time vector in the source reference frame is sampled in such a way that when dilated to the observer's frame it is sampled every day.
  int N = (int) ceil(tmax -tmin);
  std::vector<double> time(N);
  for(int t=0;t<N;t++){
    time[t] = (tmin + t)/(1+zs);
  }
  //=========================================================================================



  //=============== GET THE LIGHT CURVES IN THE SMALLEST WAVELENGTH ==========================
  std::vector< std::vector<double> > all_signals = getDRWLightCurve(time,zs,Mi,lambda_min,mean_lambda_min,N_in,123);  
  writeLightCurves(time,all_signals,out_path+"output/"+min_instrument_name+"_LC_intrinsic.json");
  //==========================================================================================



  //=============== GET THE LIGHT CURVES IN THE OTHER WAVELENGTHS ============================
  for(std::map<double,std::string>::iterator it=lambda_rest.begin();it!=lambda_rest.end();it++){
    if( it->second == min_instrument_name ){
      continue;
    } else {
      // calculate the light curves in the other filters here
      writeLightCurves(time,all_signals,out_path+"output/"+it->second+"_LC_intrinsic.json");
    }
  }
  //==========================================================================================


  
  return 0;
}
