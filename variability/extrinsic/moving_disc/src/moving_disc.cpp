#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "json/json.h"

#include "gerlumph.hpp"
#include "auxiliary_functions.hpp"
#include "instruments.hpp"
#include "noise.hpp"

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
	

  int Nlc = root["point_source"]["variability"]["extrinsic"]["Nex"].asInt();
  int Nfilters = root["instruments"].size();

  
  // Calculate the Einstein radius of the microlenses on the source plane
  double Dl  = cosmo[0]["Dl"].asDouble();
  double Ds  = cosmo[0]["Ds"].asDouble();
  double Dls = cosmo[0]["Dls"].asDouble();
  double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
  double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm

  // Create total velocity vectors
  velocityComponents vel(Nlc);
  std::vector<double> vtot(Nlc);
  std::vector<double> phi_vtot(Nlc);
  double ra  = root["point_source"]["variability"]["extrinsic"]["pars"]["ra"].asDouble();
  double dec = root["point_source"]["variability"]["extrinsic"]["pars"]["dec"].asDouble();
  double sigma_pec_l = root["point_source"]["variability"]["extrinsic"]["pars"]["sigma_pec_l"].asDouble();
  double sigma_pec_s = root["point_source"]["variability"]["extrinsic"]["pars"]["sigma_pec_s"].asDouble();
  double sigma_disp  = root["point_source"]["variability"]["extrinsic"]["pars"]["sigma_disp"].asDouble();
  double zl = root["lenses"][0]["redshift"].asDouble();
  double zs = root["source"]["redshift"].asDouble();
  vel.createVelocitiesK04(123,ra,dec,sigma_pec_l,sigma_pec_s,sigma_disp,1.0,zl,zs,Dl,Ds,Dls);
  for(int i=0;i<Nlc;i++){
    vtot[i]     = vel.tot[i].v;
    phi_vtot[i] = vel.tot[i].phi;
  }

  // Monitoring time duration in each filter
  // Read tobs_min and tobs_max
  Json::Value tobs_json;
  fin.open(output+"tobs.json",std::ifstream::in);
  fin >> tobs_json;
  fin.close();
  double tobs_max = tobs_json["tobs_max"].asDouble();
  double tobs_min = tobs_json["tobs_min"].asDouble();
  
  // Create light curve collection
  LightCurveCollection mother(Nlc);

  // Get rest wavelength at the middle of each instrument's range
  std::vector<double> lrest(Nfilters);
  for(int k=0;k<Nfilters;k++){
    Instrument mycam(root["instruments"][k]["name"].asString(),root["instruments"][k]["noise"]);
    double lobs  = (mycam.lambda_min + mycam.lambda_max)/2.0;
    lrest[k] = lobs/(1.0+zs);
  }

  // Get half-light radii
  std::vector<double> rhalfs(Nfilters);
  if( root["point_source"]["variability"]["extrinsic"]["profiles"]["type"].asString() == "parametric" ){
    double r0 = root["point_source"]["variability"]["extrinsic"]["profiles"]["r0"].asDouble();
    double l0 = root["point_source"]["variability"]["extrinsic"]["profiles"]["l0"].asDouble();
    double nu = root["point_source"]["variability"]["extrinsic"]["profiles"]["nu"].asDouble();
    for(int k=0;k<Nfilters;k++){
      rhalfs[k] = BaseProfile::sizeParametric(r0,l0,nu,lrest[k]);
    }
  } else if( root["point_source"]["variability"]["extrinsic"]["profiles"]["type"].asString() == "ss_disc" ){
    double mbh  = root["point_source"]["variability"]["extrinsic"]["profiles"]["mbh"].asDouble();
    double fedd = root["point_source"]["variability"]["extrinsic"]["profiles"]["fedd"].asDouble();
    double eta  = root["point_source"]["variability"]["extrinsic"]["profiles"]["eta"].asDouble();
    for(int k=0;k<Nfilters;k++){
      rhalfs[k] = BaseProfile::sizeSS(mbh,fedd,eta,lrest[k]);
    }
  } else if( root["point_source"]["variability"]["extrinsic"]["profiles"]["type"].asString() == "vector" ){
    for(int k=0;k<Nfilters;k++){
      rhalfs[k] = root["point_source"]["variability"]["extrinsic"]["profiles"]["rhalf"][k].asDouble();
    }    
  }
  
  // Create profile parameter maps
  Json::Value::Members json_members = root["point_source"]["variability"]["extrinsic"]["profiles"].getMemberNames();
  std::map<std::string,std::string> main_map; // size should be json_members.size()+2
  for(int k=0;k<json_members.size();k++){
    main_map.insert( std::pair<std::string,std::string>(json_members[k],root["point_source"]["variability"]["extrinsic"]["profiles"][json_members[k]].asString()) );
  }
  main_map.insert( std::pair<std::string,std::string>("pixSizePhys","") );
  main_map.insert( std::pair<std::string,std::string>("rhalf","") );
  std::vector< std::map<std::string,std::string> > profile_parameter_map(Nfilters);
  for(int k=0;k<Nfilters;k++){
    main_map["rhalf"] = std::to_string(rhalfs[k]);
    profile_parameter_map[k] = main_map;
  }
  //================= END:INITIALIZE =======================



  
  //=============== BEGIN:MAP LOOP =======================
  int profMaxOffset;
  int map_res;
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
      MagnificationMap map(maps[m]["id"].asString(),Rein);
      map_res = map.Nx;

      // Read profiles. Has to be done separately for each map because of the possibly different map.pixSizePhys
      std::vector<BaseProfile*> profiles(Nfilters);
      for(int k=0;k<Nfilters;k++){
	profile_parameter_map[k]["pixSizePhys"] = std::to_string(map.pixSizePhys);
	profiles[k] = FactoryProfile::getInstance()->createProfileFromHalfRadius(profile_parameter_map[k]);
      }
      
      Json::Value image;
      for(int k=0;k<Nfilters;k++){
	// set convolution kernel
	profMaxOffset = (int) ceil(profiles[Nfilters-1]->Nx/2);
	EffectiveMap emap(profMaxOffset,&map);
	Kernel kernel(map.Nx,map.Ny);
	
	// Set light curves
	mother.setEmap(&emap);
	mother.createVelocityLocations(254,tobs_max-tobs_min,vtot,phi_vtot,phig[m]); // Same in all filters. Will change only if duration_max is replaced by duration[k]
	
	kernel.setKernel(profiles[k]);
	map.convolve(&kernel,&emap);
	
	mother.extractFull();
	// Filter light curves
	//	int lc_index = filterMaxVelTot(vtot);
	

	
	// Output light curve
	std::string instrument_name = root["instruments"][k]["name"].asString();
	Json::Value lcs;
	for(int i=0;i<Nlc;i++){
	  Json::Value lc;
	  Json::Value time;
	  Json::Value signal;
	  Json::Value dsignal;
	  double t_interval = 11574/vtot[i]; // 11574 = 1/86400 * 10^9, first term from [day] in [s], second from 10^14 cm pixel size
	  for(int j=0;j<mother.lightCurves[i]->Nsamples;j++){
	    time.append(mother.lightCurves[i]->t[j]*t_interval + tobs_min); // light curve is originally in units of 10^14 cm
	    signal.append(mother.lightCurves[i]->m[j]);
	    dsignal.append(mother.lightCurves[i]->dm[j]);
	  }
	  lc["time"]    = time;
	  lc["signal"]  = signal;
	  lc["dsignal"] = dsignal;
	  lcs.append(lc);
	}
	image[instrument_name] = lcs;
      }
      images.append(image);


      // Light curve start and end points (different for each map, but always the same orientation - minus the shear angle),
      // in normalized units (0 to 1)
      Json::Value locs;
      for(int i=0;i<Nlc;i++){
	Json::Value lc;
	lc["Ax"] = mother.A[i].x/map_res;
	lc["Ay"] = mother.A[i].y/map_res;
	lc["Bx"] = mother.B[i].x/map_res;
	lc["By"] = mother.B[i].y/map_res;
	locs.append(lc);
      }
      maps_locs.append(locs);
      
      for(int k=0;k<Nfilters;k++){
	delete(profiles[k]);
      }
      
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

  // Write velocities
  vel.writeVelocities(output+"lc_velocities.dat");

  // Write start and end points in normalized coordinates
  std::ofstream file_maps_locs(output+"lc_xy_start_end.json");
  file_maps_locs << maps_locs;
  file_maps_locs.close();

  // Write accretion disc parameters
  Json::Value json_lobs;
  Json::Value json_lrest;
  Json::Value json_rhalf;
  for(int k=0;k<Nfilters;k++){
    json_lobs.append((1.0+zs)*lrest[k]);
    json_lrest.append(lrest[k]);
    json_rhalf.append(rhalfs[k]);
  }
  Json::Value out;
  out["Rein"] = Rein;
  out["lambda_obs"] = json_lobs;
  out["lambda_rest"] = json_lrest;
  out["half-light"] = json_rhalf;
  std::ofstream file_disc_prop(output+"accretion_disc_properties.json");
  file_disc_prop << out;
  file_disc_prop.close();

  
  return 0;
}
