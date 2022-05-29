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

#include "gerlumph.hpp"
#include "instruments.hpp"
#include "noise.hpp"


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
  std::vector<double> lrest;
  std::vector<std::string> names;
  double zs = root["source"]["redshift"].asDouble();
  for(int b=0;b<root["instruments"].size();b++){
    int Ntime = root["instruments"][b]["time"].size();
    tmins.push_back( root["instruments"][b]["time"][0].asDouble() );
    tmaxs.push_back( root["instruments"][b]["time"][Ntime-1].asDouble() );
    names.push_back( root["instruments"][b]["name"].asString() );
    Instrument mycam(root["instruments"][b]["name"].asString(),root["instruments"][b]["noise"]);
    double lobs  = (mycam.lambda_min + mycam.lambda_max)/2.0;
    lrest.push_back( lobs/(1.0+zs) );
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
    std::string ex_type = root["point_source"]["variability"]["extrinsic"]["type"].asString();
    if( ex_type == "custom" ){
      std::string extrinsic_path = in_path+"input_files/";

      for(int n=0;n<names.size();n++){
	std::string iname = names[n];      
	Json::Value json_extrinsic;
	fin.open(extrinsic_path+iname+"_LC_extrinsic.json",std::ifstream::in);
	fin >> json_extrinsic;
	fin.close();
	std::cout << duration[n] << std::endl;

	for(int q=0;q<json_extrinsic.size();q++){
	  if( json_extrinsic[q].size() > 0 ){
	    for(int i=0;i<json_extrinsic[q].size();i++){
	      if( json_extrinsic[q][i]["time"][0].asDouble() != 1.0 ){
		fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must begin at time=1, i.e. on the first day of the explosion.\n",i,iname.c_str());
		check = true;
	      }
	      int Ntime = json_extrinsic[q][i]["time"].size();
	      if( json_extrinsic[q][i]["time"][Ntime-1].asDouble() < duration[n] ){
		fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a longer duration than the corresponding intrinsic light curve(s), viz. %f days.\n",i,iname.c_str(),duration[n]);
		check = true;
	      }
	    }
	  }
	  //if( json_extrinsic[i]["signal"][0] != 1.0 ){
	  //  fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a magnification of 1 at time=0.\n",i,iname.c_str());
	  //  check = true;
	  //}
	}
      }

    } else if( ex_type == "expanding_supernova" ){
      double cutoff_fac = root["point_source"]["variability"]["extrinsic"]["size_cutoff"].asDouble();  // in Rein
      if( cutoff_fac > 7.0 ){
	fprintf(stderr,"Cutoff size for the largest supernova profile too big. Consider reducing it below 7 Einstein radii.\n");
	check = true;
      }

      double v_expand  = root["point_source"]["variability"]["extrinsic"]["v_expand"].asDouble();     // in 10^5 km/s

      Json::Value cosmo;
      fin.open(out_path+"output/angular_diameter_distances.json",std::ifstream::in);
      fin >> cosmo;
      fin.close();
      double Dl  = cosmo[0]["Dl"].asDouble();
      double Ds  = cosmo[0]["Ds"].asDouble();
      double Dls = cosmo[0]["Dls"].asDouble();
      double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
      double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm

      // Loop over the intrinsic light curve duration in each filter
      for(int i=0;i<duration.size();i++){
	double R_max = duration[i]*v_expand*8.64/Rein; // the numerator is in 10^14 cm
	if( R_max > cutoff_fac ){
	  fprintf(stderr,"Maximum physical size of the supernova in instrument %s is above the cutoff size of %f Einstein radii.\n",names[i].c_str(),cutoff_fac);
	  check = true;
	}
      }
    }


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
    Json::Value::Members json_members = root["point_source"]["variability"]["extrinsic"]["profiles"].getMemberNames();
    std::map<std::string,std::string> main_map; // size should be json_members.size()+2
    for(int i=0;i<json_members.size();i++){
      main_map.insert( std::pair<std::string,std::string>(json_members[i],root["point_source"]["variability"]["extrinsic"]["profiles"][json_members[i]].asString()) );
    }
    main_map.insert( std::pair<std::string,std::string>("rhalf","") );
    
    std::vector<double> rhalfs(names.size());
    for(int i=0;i<names.size();i++){
      if( main_map["type"] == "vector" ){
	rhalfs[i] = root["point_source"]["variability"]["extrinsic"]["profiles"]["rhalf"][i].asDouble();
      } else {
	rhalfs[i] = BaseProfile::getSize(main_map,lrest[i]);
      }
    }

    Json::Value cosmo;
    fin.open(out_path+"output/angular_diameter_distances.json",std::ifstream::in);
    fin >> cosmo;
    fin.close();
    double Dl  = cosmo[0]["Dl"].asDouble();
    double Ds  = cosmo[0]["Ds"].asDouble();
    double Dls = cosmo[0]["Dls"].asDouble();
    double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
    double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm
    for(int i=0;i<names.size();i++){
      if( rhalfs[i]/Rein > 7 ){
	fprintf(stderr,"Accretion disc size for instrument %s is too big. Consider reducing rest wavelength %f so that the disc becomes smaller than 7 Einstein radii.\n",names[i].c_str(),lrest[i]);
	check = true;	
      }
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

      // Check custom extrinsic light curves, they must extend from below tobs_min to above tobs_max
      if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "custom" ){
	Json::Value json_extrinsic;
	fin.open(extrinsic_path+iname+"_LC_extrinsic.json",std::ifstream::in);
	fin >> json_extrinsic;
	fin.close();
	
	for(int q=0;q<json_extrinsic.size();q++){
	  if( json_extrinsic[q].size() > 0 ){
	    for(int i=0;i<json_extrinsic[q].size();i++){
	      if( json_extrinsic[q][i]["time"][0] > tobs_min ){
		fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a starting time earlier than the minimum observing time, viz. <%f days.\n",i,iname.c_str(),tobs_min);
		check = true;
	      }
	      int Ntime = json_extrinsic[q][i]["time"].size();
	      if( json_extrinsic[q][i]["time"][Ntime-1] < tobs_max ){
		fprintf(stderr,"Custom extrinsic light curve %d for instrument %s must have a later ending time than the maximum observing time, viz. >%f days.\n",i,iname.c_str(),tobs_max);
		check = true;
	      }
	    }
	  }
	}
      }
      
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
