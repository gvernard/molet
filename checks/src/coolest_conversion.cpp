#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <algorithm>

#include "json/json.h"



Json::Value parseMass(Json::Value coolest_mass_model,std::string name){
  Json::Value profiles = Json::Value(Json::arrayValue);
  for(int i=0;i<coolest_mass_model.size();i++){
    Json::Value profile;
    
    // Check each mass model parameter and convert its name and value as necessary
    std::string mmodel_type = coolest_mass_model[i]["type"].asString();
    if( mmodel_type == "external_shear" ){
      profile["type"] = "external_shear";
      for( auto const& key : coolest_mass_model[i]["parameters"].getMemberNames()){
	if( key == "gamma_ext" ){
	  profile["pars"]["g"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "phi_ext" ){
	  profile["pars"]["phi"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "x0" ){
	  profile["pars"]["x0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "y0" ){
	  profile["pars"]["y0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }

    } else if( mmodel_type == "SIE" ){
      profile["type"] = "sie";
      for( auto const& key : coolest_mass_model[i]["parameters"].getMemberNames()){
	if( key == "center_x" ){
	  profile["pars"]["x0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "center_y" ){
	  profile["pars"]["y0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "phi" ){
	  profile["pars"]["pa"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "q" ){
	  profile["pars"]["q"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "theta_E" ){
	  double r_inter = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  double q = coolest_mass_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  profile["pars"]["theta_E"] = r_inter/sqrt(q);
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
      
    } else if( mmodel_type == "SPEMD" ){
      profile["type"] = "spemd";
      for( auto const& key : coolest_mass_model[i]["parameters"].getMemberNames()){
	if( key == "center_x" ){
	  profile["pars"]["x0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "center_y" ){
	  profile["pars"]["y0"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "gamma" ){
	  profile["pars"]["gam"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "phi" ){
	  profile["pars"]["pa"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "q" ){
	  profile["pars"]["q"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "theta_E" ){
	  double r_inter = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  double q = coolest_mass_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  profile["pars"]["theta_E"] = r_inter/sqrt(q);
	} else if( key == "s" ){
	  profile["pars"]["s"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
	    
    } else {
      fprintf(stderr,"Unknown mass model type '%s' for lensing entity '%s'!\n",mmodel_type.c_str(),name.c_str());
    }
    
    profiles.append( profile );
  }

  return profiles;
}


Json::Value parseLight(Json::Value coolest_light_model,std::string name){
  Json::Value profiles = Json::Value(Json::arrayValue);
  for(int i=0;i<coolest_light_model.size();i++){
    Json::Value profile;
    
    // Check each light model parameter and convert its name and value as necessary
    std::string lmodel_type = coolest_light_model[i]["type"].asString();
    if( lmodel_type == "Sersic" ){
      profile["type"] = "sersic";
      for( auto const& key : coolest_light_model[i]["parameters"].getMemberNames()){
	if( key == "I_eff" ){
	  profile["pars"]["i_eff"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();	  
	} else if( key == "R_eff" ){
	  double q = coolest_light_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  double r_inter = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  profile["pars"]["r_eff"] = r_inter/sqrt(q); 
	} else if( key == "center_x" ){
	  profile["pars"]["x0"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "center_y" ){
	  profile["pars"]["y0"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "n_sersic" ){
	  profile["pars"]["n"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "phi" ){
	  profile["pars"]["pa"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "q" ){
	  profile["pars"]["q"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
    } else {
      fprintf(stderr,"Unknown light model type '%s' for lensing entity '%s'!\n",lmodel_type.c_str(),name.c_str());
    }
    
    profiles.append( profile );
  }

  return profiles;
}



      
int main(int argc,char* argv[]){
  std::string case_path  = argv[1];
  std::string input_file = argv[2];
  Json::Value::Members jmembers;

  Json::Value root;
  std::ifstream fin(input_file,std::ifstream::in);
  fin >> root;
  fin.close();


  
  // Check if the input file is in the coolest standard, if not then return.
  if( root.isMember("standard") ){
    if( root["standard"] != "coolest" ){
      fprintf(stdout,"No conversion required.\n");
      return 0;
    }
  } else {
    fprintf(stdout,"No conversion required.\n");
    return 0;
  }



  // Get unique redshifts
  std::vector<double> redshifts;
  for(int i=0;i<root["lensing_entities"].size();i++){
    redshifts.push_back( root["lensing_entities"][i]["redshift"].asDouble() );
  }
  std::sort(redshifts.begin(),redshifts.end());
  auto last = std::unique(redshifts.begin(),redshifts.end());
  redshifts.erase(last,redshifts.end());


  
  // Sort lens entities into redshift planes
  Json::Value planes = Json::Value(Json::arrayValue);
  for(int r=0;r<redshifts.size();r++){
    Json::Value plane = Json::Value(Json::arrayValue);
    for(int i=0;i<root["lensing_entities"].size();i++){
      if( root["lensing_entities"][i]["redshift"].asDouble() == redshifts[r] ){
	plane.append( root["lensing_entities"][i] );
      }
    }
    planes.append( plane );
  }
  if( planes.size() > 2 ){
    fprintf(stderr,"More than one lens plane found. MOLET cannot handle this yet!\n");
    return 0;
  }
  

  
  // Cosmology
  Json::Value molet_cosmo;
  jmembers = root["cosmology"].getMemberNames();
  for(int i=0;i<jmembers.size();i++){
    if( jmembers[i] == "Om0" ){
      molet_cosmo["Wm0"] = root["cosmology"][jmembers[i]].asDouble();
    } else {
      molet_cosmo[jmembers[i]] = root["cosmology"][jmembers[i]];
    }
  }



  // Find instrument
  std::string myinstrument = "testCAM";
  

  
  // Loop over the lens planes and get the 'light' and 'mass' profiles of each 'lensing entity', while converting each profile to the MOLET conventions
  Json::Value molet_lenses = Json::Value(Json::arrayValue);
  for(int i=0;i<planes.size()-1;i++){
    Json::Value lens;
    lens["redshift"] = redshifts[i];

    Json::Value mass_model = Json::Value(Json::arrayValue);
    for(int j=0;j<planes[i].size();j++){
      Json::Value mmodel = parseMass(planes[i][j]["mass_model"],planes[i][j]["name"].asString());
      for(int k=0;k<mmodel.size();k++){
	mass_model.append( mmodel[k] );
      }
    }
    lens["mass_model"] = mass_model;


    Json::Value light_model;
    light_model[myinstrument] = Json::Value(Json::arrayValue);
    for(int j=0;j<planes[i].size();j++){
      Json::Value lmodel = parseLight(planes[i][j]["light_model"],planes[i][j]["name"].asString());
      for(int k=0;k<lmodel.size();k++){
	light_model[myinstrument].append( lmodel[k] );
      }
    }
    lens["light_profile"] = light_model;

    
    molet_lenses.append( lens );
  }
    

  
  // Get the 'light' profile of the source (last lens plane)
  Json::Value molet_source;
  Json::Value light_model;
  light_model[myinstrument] = Json::Value(Json::arrayValue);
  Json::Value last_plane = planes[planes.size()-1];
  for(int j=0;j<last_plane.size();j++){
    Json::Value lmodel = parseLight(last_plane[j]["light_model"],last_plane[j]["name"].asString());
    for(int k=0;k<lmodel.size();k++){
      light_model[myinstrument].append( lmodel[k] );
    }
  }
  molet_source["redshift"] = redshifts.back();
  molet_source["light_profile"] = light_model;


  
  // Point source




  
  // Create the file in the proper format acceptable by molet
  Json::Value molet_final;
  molet_final["cosmology"] = molet_cosmo;
  molet_final["lenses"] = molet_lenses;
  molet_final["source"] = molet_source;
  std::cout << molet_final << std::endl;


  
  return 0;
}



