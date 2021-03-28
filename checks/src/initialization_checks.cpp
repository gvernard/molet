#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include <algorithm>

#include "json/json.h"


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


  // Get the list of existing instruments
  DIR* modules;
  struct dirent* ent;
  std::vector<std::string> molet_instruments;
  if( (modules=opendir((molet_home+"instrument_modules/modules/").c_str())) != NULL ){
    while( (ent=readdir(modules)) != NULL ){
      molet_instruments.push_back(ent->d_name);
    }
    closedir(modules);
  } else {
    /* could not open directory */
    perror("");
    return 1;
  }

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
    if( std::find(molet_instruments.begin(),molet_instruments.end(),name) != molet_instruments.end() ){
      instruments.push_back(name);
    } else {
      fprintf(stderr,"Instrument %s is unknown! Add it to the modules or change instrument.\n",name.c_str());
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