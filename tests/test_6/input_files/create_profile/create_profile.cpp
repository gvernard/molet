#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "json/json.h"
#include "gerlumph.hpp"

int main(int argc,char* argv[]){
  // Generic options
  double rhalf = 17.16877;        // in 10^14 cm
  double pixSizePhys = rhalf/10.0;
  double ZP = 20.0;  // This has to be defined by hand
  
  std::cout << "Pixel physical size: " << pixSizePhys << std::endl;
  
  // Convert all profile parameters to string and add 'pixSizePhys' and 'rhalf'
  std::map<std::string,std::string> main_map; // size should be json_members.size()+2
  main_map["type"] = "value";
  main_map["rhalf"] = std::to_string(rhalf);
  main_map["shape"] = "gaussian";
  main_map["incl"] = "0.0";
  main_map["orient"] = "0.0";
  main_map["pixSizePhys"] = std::to_string(pixSizePhys);
  
  gerlumph::BaseProfile* base = NULL;
  base = gerlumph::FactoryProfile::getInstance()->createProfileFromHalfRadius(main_map);
  base->writeImageFITS("gaussian_rhalf.fits",1);

  double* stored = (double*) calloc(base->Nx*base->Ny,sizeof(double));
  for(long j=0;j<base->Nx*base->Ny;j++){
    stored[j] = base->data[j];
  }



  // Read in a json file containing the profile parameters
  Json::Value lc_in;
  std::ifstream fin;
  fin.open("testCAM-i_LC_intrinsic.json",std::ifstream::in);
  fin >> lc_in;
  fin.close();

  
  std::ostringstream step;
  int N = lc_in[0]["time"].size();
  for(int i=0;i<N;i++){
    double mag = lc_in[0]["signal"][i].asDouble();
    double total_flux = pow(10.0,-0.4*(mag-ZP));
    double factor = total_flux/pow(pixSizePhys,2);
    
    for(long j=0;j<base->Nx*base->Ny;j++){
      base->data[j] = stored[j]*factor;
    }
    
    step.str("");
    step.clear();
    step << std::setw(4) << std::setfill('0') << i;
    std::string fname = step.str()+".fits";
    base->writeImageFITS(fname,1);
  }

  
  delete(base);
  
  return 0;
}
