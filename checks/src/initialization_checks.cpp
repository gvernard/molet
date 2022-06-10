#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include <algorithm>

#include "json/json.h"
#include "vkllib.hpp"
#include "instruments.hpp"


inline bool file_exists(const std::string& name){
  if( FILE *file=fopen(name.c_str(),"r") ){
    fclose(file);
    return true;
  } else {
    return false;
  }
}


int main(int argc,char* argv[]){
  std::string molet_home = argv[1];
  std::string case_path  = argv[2];
  std::string input_file = argv[3];

  std::string input_dir = case_path + "input_files/";

  
  Json::Value root;
  std::ifstream fin(input_file,std::ifstream::in);
  fin >> root;
  fin.close();


  bool point_source = false;
  if( root.isMember("point_source") ){
    point_source = true;
  }


  


  // ===================== Checks start here ============================
  bool check = false;

  
  // Check if given instrument match the ones in the Molet modules and get a list of their names.
  std::vector<std::string> instruments;
  for(int i=0;i<root["instruments"].size();i++){
    std::string name = root["instruments"][i]["name"].asString();    
    if( Instrument::checkInstrumentExists(name) ){
      instruments.push_back(name);
    } else {
      fprintf(stderr,"Instrument %s is unknown! Add it to the modules or change instrument.\n",name.c_str());
      check = true;
    }
  }

  
  // Check if the given Field of View is larger than the instruments' PSFs
  for(int i=0;i<instruments.size();i++){
    double fov_width  = root["instruments"][i]["field-of-view_xmax"].asDouble() - root["instruments"][i]["field-of-view_xmin"].asDouble();
    double fov_height = root["instruments"][i]["field-of-view_ymax"].asDouble() - root["instruments"][i]["field-of-view_ymin"].asDouble();

    double ZP = root["instruments"][i]["ZP"].asDouble();
    Instrument test(instruments[i],ZP,root["instruments"][i]["noise"]);
    double psf_width  = test.original_psf->xmax - test.original_psf->xmin;
    double psf_height = test.original_psf->ymax - test.original_psf->ymin;

    if( fov_width < psf_width || fov_height < psf_height ){
      fprintf(stderr,"The Field of View (%.4f by %.4f) given for instrument %s is smaller than its PSF (%.4f by %.4f)! Increase the given Field of View.\n",fov_width,fov_height,instruments[i].c_str(),psf_width,psf_height);
      check = true;
    }
  }

  
  // Check mass model
  for(int j=0;j<root["lenses"].size();j++){
    for(int k=0;k<root["lenses"][j]["mass_model"].size();j++){
      std::string mass_type = root["lenses"][j]["mass_model"][k]["type"].asString();
      if( mass_type == "pert" ){
	std::string filename = root["lenses"][j]["mass_model"][k]["pars"]["filepath"].asString();
	if( !file_exists(input_dir+filename) ){
	  fprintf(stderr,"Custom mass model '%s' for lens %d not found! It must be provided in %s.\n",filename.c_str(),j,input_dir.c_str());
	  check = true;
	}
      }
    }
  }

  
  // Check compact mass model of the lenses if a point source exists
  if( point_source ){
    for(int j=0;j<root["lenses"].size();j++){
      for(int k=0;k<root["lenses"][j]["compact_mass_model"].size();j++){
	std::string compact_type = root["lenses"][j]["compact_mass_model"][k]["type"].asString();
	if( compact_type == "custom" ){
	  std::string filename = root["lenses"][j]["compact_mass_model"][k]["pars"]["filepath"].asString();
	  if( !file_exists(input_dir+filename) ){
	    fprintf(stderr,"Custom compact mass model '%s' for lens %d not found! It must be provided in %s.\n",filename.c_str(),j,input_dir.c_str());
	    check = true;
	  }
	}
      }
    }
  }

  
  // Exit here if there is something wrong at the first stage of the checks
  if( check ){
    return 1;
  }



  

  // Loop over instruments and check light profiles
  for(int i=0;i<instruments.size();i++){
    std::string name = instruments[i];
    
    // Check if extended source exists
    for(int j=0;j<root["source"]["light_profile"][name].size();j++){
      std::string light_type = root["source"]["light_profile"][name][j]["type"].asString(); 
      if( light_type == "custom" ){
	std::string filename = root["source"]["light_profile"][name][j]["pars"]["filepath"].asString();
	if( !file_exists(input_dir+filename) ){
	  fprintf(stderr,"Custom light profile '%s' for the source in instrument %s not found! It must be provided in %s.\n",filename.c_str(),name.c_str(),input_dir.c_str());
	  check = true;
	}
      }
    }
      
    // Check if lens light exists
    for(int j=0;j<root["lenses"].size();j++){
      for(int k=0;k<root["lenses"][j]["light_profile"][name].size();j++){
	std::string light_type = root["lenses"][j]["light_profile"][name][k]["type"].asString();
	if( light_type == "custom" ){
	  std::string filename = root["lenses"][j]["light_profile"][name][k]["pars"]["filepath"].asString();
	  if( !file_exists(input_dir+filename) ){
	    fprintf(stderr,"Custom light profile '%s' for lens %d in instrument %s not found! It must be provided in %s.\n",filename.c_str(),j,name.c_str(),input_dir.c_str());
	    check = true;
	  }
	}
      }
    }
  }


  
  // If a point source exists, loop over instruments and check light curves
  if( point_source ){

    // Check light curves if a point source exists
    if( root["point_source"]["variability"]["intrinsic"]["type"] == "custom" ){
      for(int i=0;i<instruments.size();i++){
	if( !file_exists(input_dir+instruments[i]+"_LC_intrinsic.json") ){
	  fprintf(stderr,"Intrinsic light curves for instrument %s not found! They must be provided in %s.\n",instruments[i].c_str(),input_dir.c_str());
	  check = true;
	}
      }
    } else {
      // Check if mean_mags match the given instruments
      for(int i=0;i<instruments.size();i++){
	if( !root["point_source"]["variability"]["intrinsic"]["mean_mag"].isMember(instruments[i]) ){
	  fprintf(stderr,"Intrinsic mean magnitude not provided for instrument %s!\n",instruments[i].c_str());
	  check = true;
	}
      }
    }
      
    // if( unmicro == "custom" ){
    //   // Check if file exists
    // } else {
    //   // checks for generating unmicrolensed light curve on the fly
    // }
      
    if( root["point_source"]["variability"]["extrinsic"]["type"] == "custom" ){
      for(int i=0;i<instruments.size();i++){
	if( !file_exists(input_dir+instruments[i]+"_LC_extrinsic.json") ){
	  fprintf(stderr,"Extrinsic light curves for instrument %s not found! They must be provided in %s.\n",instruments[i].c_str(),input_dir.c_str());
	  check = true;
	}
      }
    }

  }
  // ====================================================================

  
  if( check ){
    return 1;
  }

  return 0;
}
