#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "json/json.h"

#include "vkllib.hpp"
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

  std::string in_path = argv[2];
  std::string out_path = argv[3];
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
  gerlumph::velocityComponents vel(Nlc);
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
  gerlumph::LightCurveCollection mother(Nlc);




  // Get effective wavelength, i.e. the average wavelength, taking into account the instruments transmission function
  std::vector<double> lrest(Nfilters);
  std::vector<Instrument*> instruments(Nfilters);
  for(int k=0;k<Nfilters;k++){
    instruments[k] = new Instrument(root["instruments"][k]["name"].asString(),root["instruments"][k]["noise"]);
    double leff = 0.0;
    double norm = 0.0;
    for(int i=1;i<instruments[k]->wavelength.size();i++){
      norm += (instruments[k]->wavelength[i] - instruments[k]->wavelength[i-1])*(instruments[k]->throughput[i] + instruments[k]->throughput[i-1])/2.0;
      leff += (instruments[k]->wavelength[i] - instruments[k]->wavelength[i-1])*(instruments[k]->throughput[i]*instruments[k]->wavelength[i] + instruments[k]->throughput[i-1]*instruments[k]->wavelength[i-1])/2.0;
    }
    lrest[k] = (leff/norm)/(1.0+zs);
    //lrest[k] = (instruments[k]->lambda_min + instruments[k]->lambda_max)/2.0/(1.0+zs);
  }


  std::vector<double> rhalfs(Nfilters);
  std::vector< std::map<std::string,std::string> > profile_parameter_map(Nfilters);
  if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "moving_fixed_source_custom" ){

    std::map<std::string,std::string> main_map;
    main_map.insert( std::pair<std::string,std::string>("type","custom") );
    main_map.insert( std::pair<std::string,std::string>("shape","custom") );
    main_map.insert( std::pair<std::string,std::string>("incl","0") );
    main_map.insert( std::pair<std::string,std::string>("orient","0") );
    main_map.insert( std::pair<std::string,std::string>("pixSizePhys","1") ); // This is dummy here, but it is set properly in the loop over the maps
    main_map.insert( std::pair<std::string,std::string>("filename","") );
    main_map.insert( std::pair<std::string,std::string>("profPixSizePhys","") );
    main_map.insert( std::pair<std::string,std::string>("rhalf","") );
    
    for(int k=0;k<Nfilters;k++){
      std::string name = root["instruments"][k]["name"].asString();
      main_map["filename"] = in_path+"input_files/cs_"+name+".fits";
      main_map["profPixSizePhys"] = std::to_string( root["point_source"]["variability"]["extrinsic"][name]["pixSize"].asDouble() );
      gerlumph::BaseProfile* profile = gerlumph::FactoryProfile::getInstance()->createProfileFromPars(main_map);
      rhalfs[k] = profile->getHalfRadius();
      delete(profile);
      main_map["rhalf"] = std::to_string( rhalfs[k] );
      profile_parameter_map[k] = main_map;
    }
    
  } else {
    
    // Create profile parameter maps
    Json::Value::Members json_members = root["point_source"]["variability"]["extrinsic"]["profiles"].getMemberNames();
    std::map<std::string,std::string> main_map; // size should be json_members.size()+2
    for(int k=0;k<json_members.size();k++){
      main_map.insert( std::pair<std::string,std::string>(json_members[k],root["point_source"]["variability"]["extrinsic"]["profiles"][json_members[k]].asString()) );
    }
    main_map.insert( std::pair<std::string,std::string>("pixSizePhys","") );
    main_map.insert( std::pair<std::string,std::string>("rhalf","") );
    
    for(int k=0;k<Nfilters;k++){
      if( main_map["type"] == "vector" ){
	rhalfs[k] = root["point_source"]["variability"]["extrinsic"]["profiles"]["rhalf"][k].asDouble();
      } else {
	//rhalfs[k] = gerlumph::BaseProfile::getSize(main_map,lrest[k]);
	rhalfs[k] = gerlumph::BaseProfile::getSize(main_map,instruments[k]->wavelength,instruments[k]->throughput);
	std::cout << rhalfs[k] << std::endl;
      }
      main_map["rhalf"] = std::to_string(rhalfs[k]);
      profile_parameter_map[k] = main_map;
    }
    
  }

  for(int k=0;k<Nfilters;k++){
    delete(instruments[k]);
  }
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

      // Read profiles. Has to be done separately for each map because of the possibly different map.pixSizePhys
      std::vector<gerlumph::BaseProfile*> profiles(Nfilters);
      for(int k=0;k<Nfilters;k++){
	profile_parameter_map[k]["pixSizePhys"] = std::to_string(map.pixSizePhys);
	profiles[k] = gerlumph::FactoryProfile::getInstance()->createProfileFromHalfRadius(profile_parameter_map[k]);
	//profiles[k]->writeImageFITS(output+"cs_"+root["instruments"][k]["name"].asString()+".fits",1);
      }

      // Get maximum profile offset
      int maxOffset = (int) ceil(profiles[Nfilters-1]->Nx/2);
      
      Json::Value image;
      for(int k=0;k<Nfilters;k++){
	// set convolution kernel
	gerlumph::EffectiveMap emap(maxOffset,&map);
	gerlumph::Kernel kernel(map.Nx,map.Ny);
	
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
	lc["Ax"] = (maxOffset+mother.A[i].x)/map.Nx;
	lc["Ay"] = (maxOffset+mother.A[i].y)/map.Ny;
	lc["Bx"] = (maxOffset+mother.B[i].x)/map.Nx;
	lc["By"] = (maxOffset+mother.B[i].y)/map.Ny;
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
