// Checks that the duration of the following time vectors is consistent (everything is in the observer's time frame):
// - intrinsic light curve time
// - unmicrolensed light curve time
// - observing time vector for each instrument
// To compare these, we also need the maximum time delay.

// This code outputs the final minimum and maximum observational time encompassing all instruments

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "json/json.h"



int main(int argc,char* argv[]){
  std::ifstream fin;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];

  bool unmicro = false;
  if( root["point_source"]["variability"].isMember("unmicro") ){
    unmicro = true;
  }

  bool supernova = false;
  if( root["point_source"]["source_type"].asString() == "supernova" ){
    supernova = true;
  }


  // Read the multiple images' parameters from JSON to get the maximum time delay
  Json::Value images;
  fin.open(out_path+"output/multiple_images.json",std::ifstream::in);
  fin >> images;
  fin.close();
  // Get maximum image time delay
  double td_max = 0.0;
  for(int q=0;q<images.size();q++){
    double td = images[q]["dt"].asDouble();
    if( td > td_max ){
      td_max = td;
    }
  }


  // Find the tmin and tmax across all instrument times.
  std::vector<double> tmins;
  std::vector<double> tmaxs;
  std::vector<std::string> names;
  for(int b=0;b<root["instruments"].size();b++){
    int Ntime = root["instruments"][b]["time"].size();
    tmins.push_back( root["instruments"][b]["time"][0].asDouble() );
    tmaxs.push_back( root["instruments"][b]["time"][Ntime-1].asDouble() );
    names.push_back( root["instruments"][b]["name"].asString() );
  }
  double tobs_max = tmaxs[0];
  double tobs_min = tmins[0];
  for(int k=0;k<tmins.size();k++){
    if( tobs_max < tmaxs[k] ){
      tobs_max = tmaxs[k];
    }
    if( tobs_min > tmins[k] ){
      tobs_min = tmins[k];
    }
  }


  

  bool check = false;
  


  
  //=============================================== Perform various checks ================================================      
  if( supernova ){

    std::string intrinsic_path;
    if( root["point_source"]["variability"]["intrinsic"]["type"].asString() == "custom" ){
      intrinsic_path = in_path+"input_files/";
    } else {
      intrinsic_path = out_path+"output/";
    }

    // Check intrinsic light curves
    std::vector<double> duration(names.size());
    for(int n=0;n<names.size();n++){
      std::string iname = names[n];      
      Json::Value json_intrinsic;
      fin.open(intrinsic_path+iname+"_LC_intrinsic.json",std::ifstream::in);
      fin >> json_intrinsic;
      fin.close();

      // Loop for all different realizations of intrinsic light curves
      std::vector<double> dur;
      for(int i=0;i<json_intrinsic.size();i++){
	int Ntime = json_intrinsic[i]["time"].size();
	double tin_min = json_intrinsic[i]["time"][0].asDouble();
	double tin_max = json_intrinsic[i]["time"][Ntime-1].asDouble();
	dur.push_back( tin_max-tin_min );
	if( tobs_min > (tin_min-td_max) ){
	  fprintf(stderr,"To include intrinsic light curve %d, an earlier starting time by at least %f days is required for instrument %s!\n",i,tobs_min-(tin_min-td_max),iname.c_str());
	  check = true;
	}
	if( tobs_max < (tin_max+1) ){
	  fprintf(stderr,"To include intrinsic light curve %d, a later ending observing time by at least %f day is required in instrument %s!\n",i,(tin_max+1)-tobs_max,iname.c_str());
	  check = true;
	}
      }

      if( dur.size()>1 ){
	for(int i=0;i<dur.size();i++){
	  if( dur[i-1] != dur[i] ){
	    fprintf(stderr,"The length of all the intrinsic light curves for the %s instrument MUST be the same. Currently the length of light curve %d is different!\n",iname.c_str(),i);
	    check = true;
	  }
	}
      }
      duration[n] = dur[0];
    }


    // Check custom extrinsic light curves for SN
    std::string extrinsic_path;
    if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "custom" ){
      extrinsic_path = in_path+"input_files/";

      for(int n=0;n<names.size();n++){
	std::string iname = names[n];      
	Json::Value json_extrinsic;
	fin.open(extrinsic_path+iname+"_LC_extrinsic.json",std::ifstream::in);
	fin >> json_extrinsic;
	fin.close();


	for(int q=0;q<json_extrinsic.size();q++){
	  if( json_extrinsic[q].size() > 0 ){
	    for(int i=0;i<json_extrinsic[q].size();i++){
	      if( json_extrinsic[q][i]["time"][0].asDouble() != 1.0 ){
		fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must begin at time=1, i.e. on the first day of the explosion.\n",i,iname.c_str());
		check = true;
	      }
	    }
	  }
	  //if( json_extrinsic[i]["signal"][0] != 1.0 ){
	  //  fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a magnification of 1 at time=0.\n",i,iname.c_str());
	  //  check = true;
	  //}
	  //int Ntime = json_extrinsic[i]["time"].size();
	  //if( json_extrinsic[i]["time"][Ntime-1] != duration[n] ){
	  //  fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have the same duration as the corresponding intrinsic light curve(s), viz. %f days.\n",i,iname.c_str(),duration[n]);
	  //  check = true;
	  //}
	}
      }

    }

    // Check maximum time doesn't correspond to too big profiles (combined check with size_cutoff)

    // Print convolution information for the standard GERLUMPH map resolution.

  } else {
  //=======================================================================================================================    
    std::string intrinsic_path;
    if( root["point_source"]["variability"]["intrinsic"]["type"].asString() == "custom" ){
      intrinsic_path = in_path+"input_files/";
    } else {
      intrinsic_path = out_path+"output/";
    }
    std::string extrinsic_path;
    if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "custom" ){
      extrinsic_path = in_path+"input_files/";
    } else {
      extrinsic_path = out_path+"output/";
    }
    std::string unmicro_path;
    if( unmicro ){
      if( root["point_source"]["variability"]["unmicro"]["type"].asString() == "custom" ){
	unmicro_path = in_path+"input_files/";
      } else {
	unmicro_path = out_path+"output/";
      }
    }

    // Quick loop to check if the accretion dize does not become too large for some wavelength, e.g. above several Rein
    for(int n=0;n<names.size();n++){

    }
      
    for(int n=0;n<names.size();n++){
      std::string iname = names[n];
      
      // Check if all intrinsic light curves begin at t < tobs_min - td_max
      // and that the maximum observational time is at least equal to the maximum intrinsic time (in the observer's frame).
      Json::Value json_intrinsic;
      fin.open(intrinsic_path+iname+"_LC_intrinsic.json",std::ifstream::in);
      fin >> json_intrinsic;
      fin.close();
      
      for(int i=0;i<json_intrinsic.size();i++){
	if( json_intrinsic[i]["time"][0].asDouble() > (tobs_min-td_max) ){
	  fprintf(stderr,"Instrument %s: Intrinsic light curve %d requires an earlier starting time by at least %f days!\n",iname.c_str(),i,json_intrinsic[i]["time"][0].asDouble() - (tobs_min-td_max));
	  check = true;
	}
	int Ntime = json_intrinsic[i]["time"].size();
	if( json_intrinsic[i]["time"][Ntime-1].asDouble() < tobs_max ){
	  fprintf(stderr,"Instrument %s: Intrinsic light curve %d requires a later ending time by at least %f days!\n",iname.c_str(),i, tobs_max - json_intrinsic[i]["time"][Ntime-1].asDouble());
	  check = true;
	}
      }

      // Check custom extrinsic light curves
      // They must extend from 0 to above tobs_max

      
      // Check unmicrolensed light curves
      if( unmicro ){
	Json::Value json_unmicro;
	fin.open(unmicro_path+iname+"_LC_unmicro.json",std::ifstream::in);
	fin >> json_unmicro;
	fin.close();
	
	if( json_unmicro.size() != json_intrinsic.size() ){
	  fprintf(stderr,"Instrument %s: Number of intrinsic (%d) and unmicrolensed (%d) light curves should be the same!\n",iname.c_str(),json_intrinsic.size(),json_unmicro.size());
	  check = true;
	}
	
	// Same check as above for unmicrolensed light curves
	for(int i=0;i<json_unmicro.size();i++){
	  if( json_unmicro[i]["time"][0].asDouble() > (tobs_min-td_max) ){
	    fprintf(stderr,"Instrument %s: Unmicrolensed light curve %d requires an earlier starting time by at least %f days!\n",iname.c_str(),i,json_unmicro[i]["time"][0].asDouble() - (tobs_min-td_max));
	    check = true;
	  }
	  int Ntime = json_unmicro[i]["time"].size();
	  if( json_unmicro[i]["time"][Ntime-1].asDouble() < tobs_max ){
	    fprintf(stderr,"Instrument %s: Unmicrolensed light curve %d requires a later ending time by at least %f days!\n",iname.c_str(),i,tobs_max - json_unmicro[i]["time"][Ntime-1].asDouble());
	    check = true;
	  }
	}    
      }
    }
    
  }
  //=======================================================================================================================    

  
  if( check ){
    return 1;
  }

  



  
  // Output tmin and tmax for the observerational frame in all instruments
  Json::Value out;
  out["tobs_min"] = tobs_min;
  out["tobs_max"] = tobs_max;
  std::ofstream file_tobs(out_path+"output/tobs.json");
  file_tobs << out;
  file_tobs.close();
  
  return 0;
}
