#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "json/json.h"

#include "gerlumph.hpp"


int main(int argc,char* argv[]){

  /*
    Requires:
    - angular_diameter_distances.json
    - gerlumph_maps.json
    - multiple_images.json
  */
  
  //=============== BEGIN:INITIALIZE =======================
  std::ifstream fin;
  Json::Value::Members jmembers;
  
  // Read the main parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string out_path = argv[2];
  std::string output = out_path + "output/";

  
  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(output+"angular_diameter_distances.json",std::ifstream::in);
  fin >> cosmo;
  fin.close();

  // Read matching gerlumph map parameters
  Json::Value maps;
  fin.open(output+"gerlumph_maps.json",std::ifstream::in);
  fin >> maps;
  fin.close();

  // Read multiple image parameters
  Json::Value multiple_images;
  fin.open(output+"multiple_images.json",std::ifstream::in);
  fin >> multiple_images;
  fin.close();

  int Nfixed   = root["point_source"]["variability"]["extrinsic"]["Nex"].asInt();
  int Nfilters = root["instruments"].size();

  
  // Calculate the Einstein radius of the microlenses on the source plane
  double Dl  = cosmo[0]["Dl"].asDouble();
  double Ds  = cosmo[0]["Ds"].asDouble();
  double Dls = cosmo[0]["Dls"].asDouble();
  double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
  double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm

  // Maximum sized profile and consequent maxOffset
  int k_max = 0;
  double rhalf_max = 0.0;
  for(int k=0;k<Nfilters;k++){
    double v     = root["point_source"]["variability"]["extrinsic"]["profiles"][k]["v_exp"].asDouble();     // velocity in 10^14 cm/day
    int Ntime    = root["instruments"][k]["time"].size();
    double tmax  = root["instruments"][k]["time"][Ntime-1].asDouble();
    double rhalf = v*tmax;
    if(  rhalf > rhalf_max ){
      rhalf_max = rhalf;
      k_max = k;
    }
  }
  double incl_max   = root["point_source"]["variability"]["extrinsic"]["profiles"][k_max]["incl"].asDouble();    // inclination in degrees
  double orient_max = root["point_source"]["variability"]["extrinsic"]["profiles"][k_max]["orient"].asDouble();  // orientation in degrees
  //================= END:INITIALIZE =======================



  






  //=============== BEGIN:MAP LOOP =======================
  Json::Value images;

  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() == "none" ){

      Json::Value image;
      for(int k=0;k<Nfilters;k++){
 	std::string instrument_name = root["instruments"][k]["name"].asString();
	image[instrument_name] = Json::Value(Json::arrayValue);
      }
      images.append(image);
      
    } else {
      MagnificationMap map(maps[m]["id"].asString(),Rein);

      UniformDisc* maxprofile = new UniformDisc(map.pixSizePhys,rhalf_max,incl_max,orient_max);
      int maxOffset = (int) ceil(maxprofile->Nx/2);
      delete(maxprofile);
      
      FixedLocationCollection fixed(Nfixed,maxOffset,maxOffset);
      fixed.createGridLocations();

      Json::Value image;
      for(int k=0;k<Nfilters;k++){
	double v      = root["point_source"]["variability"]["extrinsic"]["profiles"][k]["v_exp"].asDouble();     // velocity in 10^14 cm/day
	double incl   = root["point_source"]["variability"]["extrinsic"]["profiles"][k]["incl"].asDouble();    // inclination in degrees
	double orient = root["point_source"]["variability"]["extrinsic"]["profiles"][k]["orient"].asDouble();  // orientation in degrees

	int Ntime = root["instruments"][k]["time"].size();
	std::vector<double> time(Ntime);
	for(int t=0;t<Ntime;t++){
	  time[t] = root["instruments"][k]["time"][t].asDouble();
	}

	double** lcs_raw = (double**) malloc(Nfixed*sizeof(double*));
	for(int f=0;f<Nfixed;f++){
	  lcs_raw[f] = (double*) malloc(Ntime*sizeof(double));
	}

	for(int t=0;t<Ntime;t++){
	  double rhalf = v*time[t]; // half light radius of a Uniform disc in 10^14cm
	  UniformDisc profile(map.pixSizePhys,rhalf,incl,orient);

	  EffectiveMap emap(maxOffset,&map);
	  Kernel kernel(map.Nx,map.Ny);
	  kernel.setKernel(&profile);
	  map.convolve(&kernel,&emap);

	  fixed.setEmap(&emap);
	  fixed.extract();
	  for(int f=0;f<Nfixed;f++){
	    lcs_raw[f][t] = fixed.m[f];
	  }
	}

	// Output light curve for filter for image
	std::string instrument_name = root["instruments"][k]["name"].asString();
	Json::Value lcs;
	for(int f=0;f<Nfixed;f++){
	  Json::Value lc;
	  Json::Value jtime;
	  Json::Value jsignal;
	  for(int t=0;t<Ntime;t++){
	    jtime.append( time[t] );
	    jsignal.append( lcs_raw[f][t] );
	  }
	  lc["time"]   = jtime;
	  lc["signal"] = jsignal;
	  lcs.append(lc);
	}
	image[instrument_name] = lcs;

	
	for(int f=0;f<Nfixed;f++){
	  delete(lcs_raw[f]);
	}
	delete(lcs_raw);
	
      }
      images.append(image);

    }
  }
  //================= END:MAP LOOP =======================




  // Write light curves
  for(int k=0;k<Nfilters;k++){
    std::string instrument_name = root["instruments"][k]["name"].asString();

    Json::Value filter;
    for(int m=0;m<maps.size();m++){
      Json::Value lcs;
      if( maps[m]["id"].asString() == "none" ){
	lcs = Json::Value(Json::arrayValue);
      } else {
	lcs = images[m][instrument_name];
      }
      filter.append(lcs);      
    }
    
    std::ofstream file_filter(output+instrument_name+"_LC_extrinsic.json");
    file_filter << filter;
    file_filter.close();
  }


  

  
  return 0;
}
