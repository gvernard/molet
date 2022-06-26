#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "json/json.h"

#include "gerlumph.hpp"
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



  // Check the maximum pixel separation between any timesteps for each velocity and report back.

  
  
  // Create light curve collection
  LightCurveCollection mother(Nlc);

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
      MagnificationMap map(maps[m]["id"].asString(),Rein);

      // Get maximum profile offset by reading the first profile in each filter
      int maxOffset = 0;
      for(int k=0;k<Nfilters;k++){
	std::string instrument_name = root["instruments"][k]["name"].asString();
	int Nx = root["point_source"]["variability"]["extrinsic"][instrument_name]["Nx"].asInt();
	double profPixSizePhys = root["point_source"]["variability"]["extrinsic"][instrument_name]["pixSize"].asDouble();
	double N_map_pixels = Nx*profPixSizePhys/map.pixSizePhys;
	int offset = (int) ceil( N_map_pixels/2.0 );
	if( offset > maxOffset ){
	  maxOffset = offset;
	}
      }
      //std::cout << maxOffset << std::endl;
      
      Json::Value image;
      std::map<std::string,std::string> profile_parameter_map;
      profile_parameter_map["shape"] = "custom";
      profile_parameter_map["incl"] = std::to_string(0.0);
      profile_parameter_map["orient"] = std::to_string(0.0);
      profile_parameter_map["pixSizePhys"] = std::to_string(map.pixSizePhys);
      for(int k=0;k<Nfilters;k++){
	std::string instrument_name = root["instruments"][k]["name"].asString();

	// set convolution kernel
	EffectiveMap emap(maxOffset,&map);
	Kernel kernel(map.Nx,map.Ny);
	
	// Set light curves
	mother.setEmap(&emap);
	double tmax = tobs_max-tobs_min;
	mother.createVelocityLocations(254,tmax,vtot,phi_vtot,phig[m]); // Same in all filters. Will change only if duration_max is replaced by duration[k]

	// Get time vector and set light curve number of sample points
	const Json::Value jtime = root["point_source"]["variability"]["extrinsic"][instrument_name]["time"];
	int Nsnapshots = jtime.size();
	std::vector<double> len_frac(Nsnapshots);
	std::vector<double> length(mother.Ncurves);
	std::vector<double> cosphi(mother.Ncurves);
	std::vector<double> sinphi(mother.Ncurves);
	for (int jj=0;jj<Nsnapshots;jj++){
	  len_frac[jj] = jtime[jj].asDouble()/tmax;
	}
	mother.setLightCurves(len_frac,length,cosphi,sinphi);

        // Start loop over images
	profile_parameter_map["profPixSizePhys"] = root["point_source"]["variability"]["extrinsic"][instrument_name]["pixSize"].asString();
	for(int jj=0;jj<Nsnapshots;jj++){
	  std::cout << "Snap: " << jj << std::endl;
	  // Read profiles. Has to be done separately for each map because of the possibly different map.pixSizePhys. And obviously separately for each filter.
	  char timestep[11];
	  int dum = sprintf(timestep,"%04d",jj);
	  profile_parameter_map["filename"] = in_path+"input_files/vs_"+instrument_name+"/"+timestep+".fits";
	  
	  BaseProfile* profile = FactoryProfile::getInstance()->createProfileFromPars(profile_parameter_map);
	  kernel.setKernel(profile);
	  map.convolve(&kernel,&emap);
	  mother.sampleLightCurveTimestep(jj,len_frac,length,cosphi,sinphi);
	  delete(profile);
	}
	
	// Output light curve
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
  Json::Value out;
  out["Rein"] = Rein;
  std::ofstream file_disc_prop(output+"accretion_disc_properties.json");
  file_disc_prop << out;
  file_disc_prop.close();

  
  return 0;
}
