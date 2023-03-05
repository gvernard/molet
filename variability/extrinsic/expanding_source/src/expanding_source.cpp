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
    - rhalfs_per_map.json
    - intrinsic_times.json
  */
  
  //=============== BEGIN:INITIALIZE =======================
  std::ifstream fin;
  
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

  std::vector<double> phig;
  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() != "none" ){
      phig.push_back( (multiple_images[m]["phig"].asDouble() + 90) ); // convert phig from EoN to normal cartesian
    } else {
      phig.push_back( 0.0 ); // I need this line
    }
  }
  
  int Nfixed   = root["point_source"]["variability"]["extrinsic"]["Nex"].asInt();
  int Nfilters = root["instruments"].size();

  
  // Calculate the Einstein radius of the microlenses on the source plane
  double Dl  = cosmo[0]["Dl"].asDouble();
  double Ds  = cosmo[0]["Ds"].asDouble();
  double Dls = cosmo[0]["Dls"].asDouble();
  double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
  double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm

  // Read half light radii for each map
  Json::Value rhalfs;
  fin.open(output+"rhalfs_per_map.json",std::ifstream::in); // in 10^14 cm
  fin >> rhalfs;
  fin.close();

  // Expansion velocities
  std::vector<double> v_expand(Nfilters);
  for(int k=0;k<Nfilters;k++){
    std::string name = root["instruments"][k]["name"].asString();
    
  }

  // Read start and end times per filter
  Json::Value tins_json;
  fin.open(output+"intrinsic_times.json",std::ifstream::in);
  fin >> tins_json;
  fin.close();
  double tins_min = tins_json["tins_min"].asDouble();
  double tins_max = tins_json["tins_max"].asDouble();
  double duration = tins_max - tins_min;
	
  
  // Maximum sized profile and consequent maxOffset
  double rhalf_max = 0.0;
  for(int m=0;m<maps.size();m++){
    int size = rhalfs[m].size();
    if( size > 0 ){
      if( rhalfs[m][size-1] > rhalf_max ){
	rhalf_max = rhalfs[m][size-1].asDouble();
      }
    }
  }
  double incl      = root["point_source"]["variability"]["extrinsic"]["incl"].asDouble();    // inclination in degrees
  double orient    = root["point_source"]["variability"]["extrinsic"]["orient"].asDouble();  // orientation in degrees
  //================= END:INITIALIZE =======================

  

  //=============== BEGIN:MAP LOOP =======================
  Json::Value images;
  Json::Value maps_locs;
  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() == "none" ){

      Json::Value image;
      for(int k=0;k<Nfilters;k++){
 	std::string instrument_name = root["instruments"][k]["name"].asString();
	image[instrument_name] = Json::Value(Json::arrayValue);
      }
      images.append(image);
      maps_locs.append(Json::Value(Json::arrayValue));
      
    } else {
      gerlumph::MagnificationMap map(maps[m]["id"].asString(),Rein);
      
      // Get maximum profile offset
      gerlumph::UniformDisc* maxprofile = new gerlumph::UniformDisc(map.pixSizePhys,rhalf_max,incl,orient);
      int maxOffset = (int) ceil(maxprofile->Nx/2);
      delete(maxprofile);
      
      // Create fixed locations. This, in fact, is map independent but we need to do it here because it requires maxOffset to be known
      gerlumph::FixedLocationCollection fixed(Nfixed,map.Nx - 2*maxOffset,map.Ny - 2*maxOffset);
      fixed.createGridLocations();

      Json::Value image;
      int Nrhalf = rhalfs[m].size();

      double** raw_signal  = (double**) malloc(Nfixed*sizeof(double*));
      double** raw_dsignal = (double**) malloc(Nfixed*sizeof(double*));
      for(int f=0;f<Nfixed;f++){
	raw_signal[f]  = (double*) malloc(Nrhalf*sizeof(double));
	raw_dsignal[f] = (double*) malloc(Nrhalf*sizeof(double));
      }

      for(int r=0;r<Nrhalf;r++){
	double rhalf = rhalfs[m][r].asDouble(); // half light radius in 10^14cm
	gerlumph::UniformDisc profile(map.pixSizePhys,rhalf,incl,orient); // shape of the brightness profile

	gerlumph::EffectiveMap emap(maxOffset,&map);
	gerlumph::Kernel kernel(map.Nx,map.Ny);
	kernel.setKernel(&profile);
	map.convolve(&kernel,&emap);
	
	fixed.setEmap(&emap);
	fixed.extract();
	for(int f=0;f<Nfixed;f++){
	  raw_signal[f][r]  = fixed.m[f];
	  raw_dsignal[f][r] = fixed.dm[f];
	}
      }
      
      for(int k=0;k<Nfilters;k++){
	// Store light curve for filter for image
	std::string instrument_name = root["instruments"][k]["name"].asString();
	Json::Value lcs;

	// Convert rhalfs to time using the filter's velocity
	Json::Value jtime;
	for(int r=0;r<Nrhalf;r++){
	  double rhalf = rhalfs[m][r].asDouble(); // half light radius in 10^14cm
	  double v_expand = 8.64*root["point_source"]["variability"]["extrinsic"]["v_expand"][instrument_name].asDouble();   // velocity in 10^14 cm / day
	  double t = rhalf/v_expand; // in days
	  jtime.append(t); // Keeping time normalized to 0.
	  //jtime.append(t + tins_min); // Converting to the intrinsic light curve's time
	  if( t > duration ){
	    break; // Check after assignment to make sure the full duration is included in the light curve
	  }
	}

	for(int f=0;f<Nfixed;f++){
	  Json::Value lc;
	  Json::Value jsignal;
	  Json::Value jdsignal;
	  for(int t=0;t<jtime.size();t++){
	    jsignal.append( raw_signal[f][t] );
	    jdsignal.append( raw_dsignal[f][t] );
	  }
	  lc["time"]    = jtime;
	  lc["signal"]  = jsignal;
	  lc["dsignal"] = jdsignal;
	  lcs.append(lc);
	}
	image[instrument_name] = lcs;
      }
      printf("Map %d done.\n",m);
      images.append(image);
	
      // Some clean up
      for(int f=0;f<Nfixed;f++){
	delete(raw_signal[f]);
	delete(raw_dsignal[f]);
      }
      delete(raw_signal);
      delete(raw_dsignal);
      
      // Fixed locations (different for each map, but always the same orientation - minus the shear angle),
      // in normalized units (0 to 1)
      Json::Value locs;
      for(int i=0;i<Nfixed;i++){
	Json::Value point;
	point["x"] = (maxOffset + fixed.A[i].x)/map.Nx;
	point["y"] = (maxOffset + fixed.A[i].y)/map.Ny;
	locs.append(point);
      }
      maps_locs.append(locs);
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

  
  // Write locations in normalized coordinates
  std::ofstream file_maps_locs(output+"fixed_xy_locations.json");
  file_maps_locs << maps_locs;
  file_maps_locs.close();

  // Write supernova parameters
  Json::Value out;
  out["Rein"] = Rein;
  std::ofstream file_sn_prop(output+"supernova_profile_properties.json");
  file_sn_prop << out;
  file_sn_prop.close();
  
  return 0;
}
