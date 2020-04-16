#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "json/json.h"

#include "gerlumph.hpp"
#include "auxiliary_functions.hpp"

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



  int Nlc = root["point_source"]["variability"]["extrinsic"]["Nex"].asInt();
  int Nfilters = root["instrument"]["bands"].size();

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
  vel.createVelocitiesK04(321,ra,dec,sigma_pec_l,sigma_pec_s,sigma_disp,1.0,zl,zs,Dl,Ds,Dls);
  for(int i=0;i<Nlc;i++){
    vtot[i]     = vel.tot[i].v;
    phi_vtot[i] = vel.tot[i].phi;
  }

  // Monitoring time duration in each filter
  std::vector<double> duration;
  for(int k=0;k<Nfilters;k++){
    //  double duration = dateDifference(root["instrument"]["start"].asString(),root["instrument"]["end"].asString()); // in days
    int Ntime   = root["instrument"]["bands"][k]["time"].size();
    double t0   = root["instrument"]["bands"][k]["time"][0].asDouble();
    double tmax = root["instrument"]["bands"][k]["time"][Ntime-1].asDouble();
    duration.push_back( ceil(tmax) - floor(t0) );
  }	
  //================= END:INITIALIZE =======================



  
  //=============== BEGIN:MAP LOOP =======================
  Json::Value images;
  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() == "none" ){

      Json::Value image;
      for(int k=0;k<Nfilters;k++){
 	std::string band_name = root["instrument"]["bands"][k]["name"].asString();
	image[band_name] = Json::Value(Json::arrayValue);
      }
      images.append(image);

    } else {
      MagnificationMap map(maps[m]["id"].asString(),Rein);
      
      std::vector<BaseProfile*> profiles(Nfilters);
      for(int k=0;k<Nfilters;k++){
	Json::Value json_profile = root["point_source"]["variability"]["extrinsic"]["profiles"][k];
	BaseProfile* profile = createProfileFromJson(json_profile,map.pixSizePhys);
	profiles[k] = profile;
      }
      
      Json::Value image;
      for(int k=0;k<Nfilters;k++){
	// set convolution kernel
	int profMaxOffset = (int) ceil(profiles[Nfilters-1]->Nx/2);
	EffectiveMap emap(profMaxOffset,&map);
	Kernel kernel(map.Nx,map.Ny);
	
	// Rotate map (i.e. rotate the velocity vectors in an opposite way)
	for(int i=0;i<Nlc;i++){
	  vel.tot[i].phi = vel.tot[i].phi - (multiple_images[m]["phig"].asDouble() - 90);
	  phi_vtot[i]    = phi_vtot[i] - (multiple_images[m]["phig"].asDouble() - 90);
	}
	
	// Set light curves
	LightCurveCollection mother(Nlc,&emap);
	mother.createVelocityLocations(213,duration[k],vtot,phi_vtot);
	
	kernel.setKernel(profiles[k]);
	map.convolve(&kernel,&emap);
	
	mother.extractFull();
	// Filter light curves
	//	int lc_index = filterMaxVelTot(vtot);
	

	
	// Output light curve
	std::string band_name = root["instrument"]["bands"][k]["name"].asString();
	Json::Value lcs;
	for(int i=0;i<Nlc;i++){
	  Json::Value lc;
	  Json::Value time;
	  Json::Value signal;
	  double t_interval = 11574*map.pixSizePhys/vtot[i]; // 11574 = 1/86400 * 10^9, first term from [day] in [s], second from 10^14 cm pixel size
	  for(int j=0;j<mother.lightCurves[i]->Nsamples;j++){
	    time.append(mother.lightCurves[i]->t[j]*t_interval);
	    signal.append(mother.lightCurves[i]->m[j]);
	  }
	  lc["time"] = time;
	  lc["signal"] = signal;
	  lcs.append(lc);
	}
	image[band_name] = lcs;
      }
      images.append(image);

      
      for(int k=0;k<Nfilters;k++){
	delete(profiles[k]);
      }
      
    }
  }
  //================= END:MAP LOOP =======================
  
  std::ofstream file_images(output+"extrinsic_light_curves.json");
  file_images << images;
  file_images.close();

  return 0;
}
