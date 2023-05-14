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
    - tobs.json
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

  // Calculate the Einstein radius of the microlenses on the source plane
  double Dl  = cosmo[0]["Dl"].asDouble();
  double Ds  = cosmo[0]["Ds"].asDouble();
  double Dls = cosmo[0]["Dls"].asDouble();
  double M   = root["point_source"]["variability"]["extrinsic"]["microlens_mass"].asDouble();
  double Rein = 13.5*sqrt(M*Dls*Ds/Dl); // in 10^14 cm

  // Read tobs_min and tobs_max
  Json::Value tobs_json;
  fin.open(output+"tobs.json",std::ifstream::in);
  fin >> tobs_json;
  fin.close();
  double tobs_max = tobs_json["tobs_max"].asDouble();
  double tobs_min = tobs_json["tobs_min"].asDouble();
  double duration = tobs_max - tobs_min;
  
  // Number of filters
  int Nfilters = root["instruments"].size();

  // Quantities to be used in determining the half light radii
  std::vector<double> vels(Nfilters);
  for(int k=0;k<Nfilters;k++){
    std::string name = root["instruments"][k]["name"].asString();
    vels[k] = root["point_source"]["variability"]["extrinsic"]["v_expand"][name].asDouble();   // velocity in 10^5 km/s
  }
  double v_expand  = *std::max_element(std::begin(vels),std::end(vels));
  double ff        = root["point_source"]["variability"]["extrinsic"]["fractional_increase"].asDouble();  // Fractional increase in flux of a uniform disc at each timestep
  double start     = root["point_source"]["variability"]["extrinsic"]["start"].asDouble();  // Days after explosion to start from
  double cutoff    = root["point_source"]["variability"]["extrinsic"]["size_cutoff"].asDouble();  // Size of the SN profile above which we assume it is too big to be microlensed
  //================= END:INITIALIZE =======================


  
  //=============== BEGIN:CREATE TIME AND SIZE VECTORS ===============
  // This part is for the fastest filter!!!
  std::vector< std::vector< double > > rhalfs_per_map(maps.size());

  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() == "none" ){
      std::vector<double> rhalfs; // just an empty vector
      rhalfs_per_map[m] = rhalfs; // can be the same empty vector
    } else {
      // Create sampling time (map dependent)
      // Here we compute the time vector so that at each timestep the area of a circle with the given r1/2 (=v_expand*t) is ff% larger than before.
      double pixSizePhys = gerlumph::MagnificationMap::getPixSizePhys(maps[m]["id"].asString(),Rein);
      std::vector<double> time{start}; // in days
      std::vector<double> rhalfs{start*v_expand*8.64}; // in 10^14 cm
      int i = 0;
      double t  = time.back();
      int cutoff_pix = (int) ceil(cutoff*Rein/pixSizePhys);
      do{
	double dt = t*(sqrt(1.0 + ff) - start);
	if( dt < 1.0 ){
	  dt = 1.0; // Set minimum dt to 1 day
	}
	int dr_pix = (int) ceil( v_expand*dt*8.64/pixSizePhys ); // dr converted to 10^14 cm and then to map pixels.
	t += dt;
	int r_pix  = (int) ceil( v_expand*t*8.64/pixSizePhys ); // r at t+dt converted to 10^14 cm and then to map pixels.
	if( dr_pix >= 1 && r_pix < cutoff_pix ){
	  //printf("%f %f %d %d\n",t,dt,dr_pix,r_pix);
	  time.push_back(t); // consider this time step only if it corresponds to +1 or more pixels in the map, otherwise it would be like performing the same convolution twice.
	  rhalfs.push_back(t*v_expand*8.64); // in 10^14 cm
	}
      } while( t<duration ); // This is an 'until' loop (condition checked at the end) to make sure tobs_max will exist within the final light curve
      rhalfs_per_map[m] = rhalfs;
    }
  }
  //================= END:CREATE TIME AND SIZE VECTORS ===============



  //=============== BEGIN:REPORT NUMBER OF CONVOLUTIONS AND ASK FOR CONFIRMATION ===============
  printf(">>>>>> Please pay ATTENTION to the following: <<<<<<\n");
  printf("   The number of convolutions per map is:\n");
  int convs_tot = 0;
  for(int m=0;m<maps.size();m++){
    if( maps[m]["id"].asString() != "none" ){
      printf(" map %d (%s,%.2f,%.2f,%.2f): %lu\n",m,maps[m]["id"].asString().c_str(),maps[m]["k"].asDouble(),maps[m]["g"].asDouble(),maps[m]["s"].asDouble(),rhalfs_per_map[m].size());
      convs_tot += rhalfs_per_map[m].size();
    }
  }
  printf("   The total number of convolutions will be %d.\n",convs_tot);
  double conv_fac = convs_tot/3600.0;
  printf("   It will approximately take %.2f hours (%d mins) on the CPU and %.2f hours (%d mins) on the GPU.\n",conv_fac*13.0,(int)std::ceil(conv_fac*13.0*60),conv_fac*3.0,(int)std::ceil(conv_fac*3.0*60)); // convolution on a CPU lasts 13s and on the GPU 3s (approximate numbers)
  printf("   (if you want to compille gerlumphpp for GPUs see here: https://github.com/gvernard/gerlumphpp)\n");
  printf("   If that is too long, consider shortening the observing time for each instrument, or increasing the fractional_increase parameter, or decreasing the size_cutoff parameter.\n");
  char choice;
  bool run = true;
  do{
    std::cout << "-> Would you like to proceed with the calculations? (answer 'y' or 'n')" << std::endl;
    std::cin >> choice;
    choice = std::tolower(choice);//Put your letter to its lower case
  } while( choice != 'n' && choice != 'y' );
  if( choice =='n' ){
    return 1;
  }
  //================= END:REPORT NUMBER OF CONVOLUTIONS AND ASK FOR CONFIRMATION ===============



  // Write hald-light radii
  Json::Value rhalf;
  for(int m=0;m<maps.size();m++){
    Json::Value jobj;
    if( maps[m]["id"].asString() == "none" ){
      jobj = Json::arrayValue;
    } else {
      jobj = Json::arrayValue;
      for(int i=0;i<rhalfs_per_map[m].size();i++){
	jobj.append( rhalfs_per_map[m][i] );
      }
    }
    rhalf.append(jobj);
  }
    
  std::ofstream file_rhalf(output+"rhalfs_per_map.json");
  file_rhalf << rhalf;
  file_rhalf.close();


  
  return 0;
}
