#include <cmath>
#include <fstream>
#include <iostream>

#include "json/json.h"

#include "auxiliary_functions.hpp"

// Based on the python version of Ned Wright's javascript cosmology calculator by James Schombert,
// which can be found here: http://www.astro.ucla.edu/~wright/CC.python

int main(int argc,char* argv[]){


  // Read input_molet.json
  std::ifstream fin;
  Json::Value::Members jmembers;

  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string out_path = argv[2];
  std::string output = out_path+"output/";
  

  
  // Assume Benchmark model
  double H0 = root["cosmology"]["H0"].asDouble();
  double Wm = root["cosmology"]["Wm0"].asDouble(); // Omega matter at t=t0 (now)
  double Wv = 1.0 - Wm - 0.4165/(H0*H0);           // Omega vacuum, or Lambda
  // initialize constants
  double c = 299792.458;  // velocity of light in km/sec
  double h  = H0/100.0;
  double Wr = 4.165E-5/(h*h);   // includes 3 massless neutrino species, T0 = 2.72528
  double Wk = 1-Wm-Wr-Wv;
  double DH = c/H0;
  


  std::vector<double> zs;
  std::vector<double> zl;
  
  zl.push_back( root["lenses"][0]["redshift"].asDouble() );
  zs.push_back( root["source"]["redshift"].asDouble() );
  
  if( root["lenses"].size() > 1 ){
    for(int k=1;k<root["lenses"].size()+1;k++){
      zl.push_back( root["lenses"][-1-k]["redshift"].asDouble() );
      zs.push_back( root["lenses"][-1-(k-1)]["redshift"].asDouble() );
    }
  }


  Json::Value distances;
  for(int k=0;k<zl.size();k++){
    // Distance to the lens
    double DCMT_l = transverse_comoving_distance(Wm,Wr,Wk,Wv,zl[k]);
    double DA_l   = DH*DCMT_l/(1.0+zl[k]); // in Mpc
    
    // Distance to the source
    double DCMT_s = transverse_comoving_distance(Wm,Wr,Wk,Wv,zs[k]);
    double DA_s   = DH*DCMT_s/(1.0+zs[k]); // in Mpc
    
    // Distance between lens and source
    double DA_ls = (DCMT_s*sqrt(1.0 + Wk*pow(DCMT_l/DH,2)) - DCMT_l*sqrt(1.0 + Wk*pow(DCMT_s/DH,2)))/(1.0+zs[k]);
    DA_ls        = DH*DA_ls; // in Mpc

    Json::Value dum;
    dum["Dl"]  = DA_l;
    dum["Ds"]  = DA_s;
    dum["Dls"] = DA_ls;
    distances.append(dum);
  }
  

  std::ofstream file_distances(output+"angular_diameter_distances.json");
  file_distances << distances;
  file_distances.close();
 

  return 0;
}
