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

  // Monitoring time duration in each filter
  // Read tobs_min and tobs_max
  Json::Value tobs_json;
  fin.open(output+"tobs.json",std::ifstream::in);
  fin >> tobs_json;
  fin.close();
  double tobs_max = tobs_json["tobs_max"].asDouble();
  double tobs_min = tobs_json["tobs_min"].asDouble();

  
  // Maximum sized profile and consequent maxOffset
  double v_expand  = root["point_source"]["variability"]["extrinsic"]["v_expand"].asDouble();   // velocity in 10^5 km/s
  double R_max     = (tobs_max-tobs_min)*v_expand*8.64; // this is in 10^14 cm
  double rhalf_max = R_max/sqrt(2.0); // this is the corresponding half-light radius for a uniform disc of radius R_max
  double incl      = root["point_source"]["variability"]["extrinsic"]["incl"].asDouble();    // inclination in degrees
  double orient    = root["point_source"]["variability"]["extrinsic"]["orient"].asDouble();  // orientation in degrees
  double ff        = root["point_source"]["variability"]["extrinsic"]["fractional_increase"].asDouble();  // fractional increase in flux of a uniform disc at each timestep
  double cutoff    = root["point_source"]["variability"]["extrinsic"]["size_cutoff"].asDouble();  // Size of the SN profile above which we assume it is too big to be microlensed
  //================= END:INITIALIZE =======================



  






  //=============== BEGIN:MAP LOOP =======================
  Json::Value images;
  Json::Value maps_locs;
  std::vector< std::vector<double> > times;
  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() == "none" ){

      Json::Value image;
      for(int k=0;k<Nfilters;k++){
 	std::string instrument_name = root["instruments"][k]["name"].asString();
	image[instrument_name] = Json::Value(Json::arrayValue);
      }
      images.append(image);
      maps_locs.append(Json::Value(Json::arrayValue));
      std::vector<double> time; // just an empty vector
      times.push_back(time);
      
    } else {
      MagnificationMap map(maps[m]["id"].asString(),Rein);

      // Get maximum profile offset
      UniformDisc* maxprofile = new UniformDisc(map.pixSizePhys,rhalf_max,incl,orient);
      int maxOffset = (int) ceil(maxprofile->Nx/2);
      delete(maxprofile);
      
      // Create fixed locations. This, in fact, is map independent but we need to do it here because it requires maxOffset to be known
      FixedLocationCollection fixed(Nfixed,map.Nx - 2*maxOffset,map.Ny - 2*maxOffset);
      fixed.createGridLocations();

      // Create sampling time (map dependent)
      // Here we compute the time vector so that at each timestep the area of a circle with the given r1/2 (=v_expand*t) is ff% larger than before.
      std::vector<double> time{1}; // in days
      int i = 0;
      double t  = time.back();
      int cutoff_pix = (int) ceil(cutoff*Rein/map.pixSizePhys);
      do{
	double dt = t*(sqrt(1.0 + ff) - 1.0);
	if( dt < 1.0 ){
	  dt = 1.0; // Set minimum dt to 1 day
	}
	int dr_pix = (int) ceil( v_expand*dt*8.64/map.pixSizePhys ); // dr converted to 10^14 cm and then to map pixels.
	t += dt;
	int r_pix  = (int) ceil( v_expand*t*8.64/map.pixSizePhys ); // r at t+dt converted to 10^14 cm and then to map pixels.
	if( dr_pix >= 1 && r_pix < cutoff_pix ){
	  //printf("%f %f %d %d\n",t,dt,dr_pix,r_pix);
	  time.push_back(t); // consider this time step only if it corresponds to +1 or more pixels in the map, otherwise it would be like performing the same convolution twice.
	}
      } while( t<tobs_max ); // This is an 'until' loop (condition checked at the end) to make sure tobs_max will exist within the final light curve
      times.push_back(time);
      int Ntime = time.size();


      Json::Value image;
      for(int k=0;k<Nfilters;k++){
	double** raw_signal  = (double**) malloc(Nfixed*sizeof(double*));
	double** raw_dsignal = (double**) malloc(Nfixed*sizeof(double*));
	for(int f=0;f<Nfixed;f++){
	  raw_signal[f]  = (double*) malloc(Ntime*sizeof(double));
	  raw_dsignal[f] = (double*) malloc(Ntime*sizeof(double));
	}


	for(int t=0;t<Ntime;t++){
	  double rhalf = v_expand*time[t]; // half light radius of a Uniform disc in 10^14cm
	  UniformDisc profile(map.pixSizePhys,rhalf,incl,orient);

	  EffectiveMap emap(maxOffset,&map);
	  Kernel kernel(map.Nx,map.Ny);
	  kernel.setKernel(&profile);
	  map.convolve(&kernel,&emap);

	  fixed.setEmap(&emap);
	  fixed.extract();
	  for(int f=0;f<Nfixed;f++){
	    raw_signal[f][t]  = fixed.m[f];
	    raw_dsignal[f][t] = fixed.dm[f];
	  }
	}

	
	// Output light curve for filter for image
	std::string instrument_name = root["instruments"][k]["name"].asString();
	Json::Value lcs;
	for(int f=0;f<Nfixed;f++){
	  Json::Value lc;
	  Json::Value jtime;
	  Json::Value jsignal;
	  Json::Value jdsignal;
	  for(int t=0;t<Ntime;t++){
	    jtime.append( time[t] ); // must NOT convert to absolute time
	    jsignal.append( raw_signal[f][t] );
	    jdsignal.append( raw_dsignal[f][t] );
	  }
	  lc["time"]    = jtime;
	  lc["signal"]  = jsignal;
	  lc["dsignal"] = jsignal;
	  lcs.append(lc);
	}
	image[instrument_name] = lcs;
	
	for(int f=0;f<Nfixed;f++){
	  delete(raw_signal[f]);
	  delete(raw_dsignal[f]);
	}
	delete(raw_signal);
	delete(raw_dsignal);
	
      }
      images.append(image);

      
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

  // Write accretion disc parameters
  Json::Value json_rhalfs;
  for(int m=0;m<maps.size();m++){
    Json::Value rhalf;
    std::cout << times[m].size() << std::endl;
    for(int t=0;t<times[m].size();t++){
      double uni_rhalf = v_expand*times[m][t]*8.64/sqrt(2.0);
      //double uni_rhalf = times[m][t];
      rhalf.append(uni_rhalf); // in 10^14 cm
    }
    json_rhalfs[maps[m]["id"].asString()] = rhalf;
  }
  Json::Value out;
  out["Rein"] = Rein;
  out["half-light"] = json_rhalfs;
  std::ofstream file_disc_prop(output+"supernova_profile_properties.json");
  file_disc_prop << out;
  file_disc_prop.close();
  
  return 0;
}
