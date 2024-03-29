#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <algorithm>

#include "json/json.h"

#include "instruments.hpp"

Json::Value parseMass(Json::Value coolest_mass_model,std::string name){
  Json::Value profiles = Json::Value(Json::arrayValue);
  for(int i=0;i<coolest_mass_model.size();i++){
    Json::Value profile;
    
    // Check each mass model parameter and convert its name and value as necessary
    std::string mmodel_type = coolest_mass_model[i]["type"].asString();
    if( mmodel_type == "ExternalShear" ){
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
      if( !profile["pars"].isMember("x0") ){
	profile["pars"]["x0"] = 0.0;
      }
      if( !profile["pars"].isMember("y0") ){
	profile["pars"]["y0"] = 0.0;
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
	  // double r_inter = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  // double q = coolest_mass_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  // profile["pars"]["theta_E"] = r_inter/sqrt(q);
	  profile["pars"]["theta_E"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
      
    } else if( mmodel_type == "PEMD" ){
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
	  // double r_inter = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  // double q = coolest_mass_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  // profile["pars"]["theta_E"] = r_inter/sqrt(q);
	  profile["pars"]["theta_E"] = coolest_mass_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
      profile["pars"]["s"] = 0.000001;

    } else if( mmodel_type == "PixelatedRegularGrid" ){
      profile["type"] = "pert";
      for( auto const& key : coolest_mass_model[i]["pixels"].getMemberNames()){
	if( key == "fits_file" ){
	  profile["pars"]["filepath"] = coolest_mass_model[i]["pixels"][key]["path"].asString();
	} else if( key == "field_of_view_x" ){
	  profile["pars"]["xmin"] = coolest_mass_model[i]["pixels"][key][0].asDouble();
	  profile["pars"]["xmax"] = coolest_mass_model[i]["pixels"][key][1].asDouble();
	} else if( key == "field_of_view_y" ){
	  profile["pars"]["ymin"] = coolest_mass_model[i]["pixels"][key][0].asDouble();
	  profile["pars"]["ymax"] = coolest_mass_model[i]["pixels"][key][1].asDouble();
	} else if( key == "num_pix_x" ){
	  profile["pars"]["Nx"] = coolest_mass_model[i]["pixels"][key].asInt();
	} else if( key == "num_pix_y" ){
	  profile["pars"]["Ny"] = coolest_mass_model[i]["pixels"][key].asInt();
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


Json::Value parseLight(Json::Value coolest_light_model,std::string name,std::string case_path){
  Json::Value profiles = Json::Value(Json::arrayValue);
  for(int i=0;i<coolest_light_model.size();i++){
    Json::Value profile;
    
    // Check each light model parameter and convert its name and value as necessary
    std::string lmodel_type = coolest_light_model[i]["type"].asString();
    if( lmodel_type == "Sersic" ){
      profile["type"] = "sersic";
      for( auto const& key : coolest_light_model[i]["parameters"].getMemberNames()){
	if( key == "I_eff" ){
	  double q = coolest_light_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  double coolest_i_eff = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	  profile["pars"]["i_eff"] = q*coolest_i_eff;
	} else if( key == "theta_eff" ){
	  double q = coolest_light_model[i]["parameters"]["q"]["point_estimate"]["value"].asDouble();
	  double coolest_r_eff = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble(); // intermediate axis
	  profile["pars"]["r_eff"] = coolest_r_eff*sqrt(q);
	} else if( key == "center_x" ){
	  profile["pars"]["x0"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "center_y" ){
	  profile["pars"]["y0"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "n" ){
	  profile["pars"]["n"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "phi" ){
	  profile["pars"]["pa"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else if( key == "q" ){
	  profile["pars"]["q"] = coolest_light_model[i]["parameters"][key]["point_estimate"]["value"].asDouble();
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
      profiles.append( profile );

    } else if( lmodel_type == "PixelatedRegularGrid" ){
      profile["type"] = "custom";   
      for( auto const& key : coolest_light_model[i]["parameters"]["pixels"].getMemberNames()){
	if( key == "fits_file" ){
	  profile["pars"]["filepath"] = coolest_light_model[i]["parameters"]["pixels"][key]["path"].asString();
	} else if( key == "field_of_view_x" ){
	  profile["pars"]["xmin"] = coolest_light_model[i]["parameters"]["pixels"][key][0].asDouble();
	  profile["pars"]["xmax"] = coolest_light_model[i]["parameters"]["pixels"][key][1].asDouble();
	} else if( key == "field_of_view_y" ){
	  profile["pars"]["ymin"] = coolest_light_model[i]["parameters"]["pixels"][key][0].asDouble();
	  profile["pars"]["ymax"] = coolest_light_model[i]["parameters"]["pixels"][key][1].asDouble();
	} else if( key == "num_pix_x" ){
	  profile["pars"]["Nx"] = coolest_light_model[i]["parameters"]["pixels"][key].asInt();
	} else if( key == "num_pix_y" ){
	  profile["pars"]["Ny"] = coolest_light_model[i]["parameters"]["pixels"][key].asInt();
	} else if( key == "id" ){
	  continue;
	} else {
	  fprintf(stderr,"Unknown parameter name '%s' for lensing entity '%s'!\n",key.c_str(),name.c_str());
	}
      }
      // Convert image pixels from COOLEST units [e/s] to MOLET units [e/(s arcsec^2].
      double xmax = profile["pars"]["xmax"].asDouble();
      double xmin = profile["pars"]["xmin"].asDouble();
      double ymax = profile["pars"]["ymax"].asDouble();
      double ymin = profile["pars"]["ymin"].asDouble();
      int Nx = profile["pars"]["Nx"].asInt();
      int Ny = profile["pars"]["Ny"].asInt();
      double pix_area = (xmax-xmin)*(ymax-ymin)/(Nx*Ny);
      std::string src_fname = profile["pars"]["filepath"].asString();
      int lastindex = src_fname.find_last_of("."); 
      std::string new_src_fname = src_fname.substr(0,lastindex) + "_molet_units.fits";
      
      // Read, convert, and write source
      vkl::RectGrid input_src(Nx,Ny,xmin,xmax,ymin,ymax,case_path+"/"+src_fname);
      for(int i=0;i<input_src.Nz;i++){
	input_src.z[i] /= pix_area;
      }
      vkl::FitsInterface::writeFits(Nx,Ny,input_src.z,case_path+"/"+new_src_fname);
      profile["pars"]["filepath"] = new_src_fname;
      
      profiles.append( profile );

    } else if( lmodel_type == "LensedPS" ){
      // Do nothing
      
    } else {
      fprintf(stderr,"HERE");
      fprintf(stderr,"Unknown light model type '%s' for lensing entity '%s'!\n",lmodel_type.c_str(),name.c_str());
    }
    
  }

  return profiles;
}



Json::Value getPointSource(Json::Value last_plane){
  Json::Value ps_profile;
  bool exit_flag = false;
  
  for(int j=0;j<last_plane.size();j++){
  
    if( last_plane[j].isMember("light_model") ){
      
      for(int i=0;i<last_plane[j]["light_model"].size();i++){    
	std::string lmodel_type = last_plane[j]["light_model"][i]["type"].asString();
	if( lmodel_type == "LensedPS" ){
	  std::string flag_contains = last_plane[j]["light_model"][i]["parameters"]["flag_contains"]["point_estimate"].asString();
	  bool flag_coupled = last_plane[j]["light_model"][i]["parameters"]["flat_coupled"]["point_estimate"].asBool();
	  if( flag_contains == "true" || (flag_contains == "both" && flag_coupled == true) ){
	    ps_profile["x0"] = last_plane[j]["light_model"][i]["parameters"]["x_true"]["point_estimate"]["value"].asDouble();
	    ps_profile["y0"] = last_plane[j]["light_model"][i]["parameters"]["y_true"]["point_estimate"]["value"].asDouble();
	    ps_profile["M_tot_unlensed"] = last_plane[j]["light_model"][i]["parameters"]["m_true"]["point_estimate"]["value"].asDouble();
	  }
	  exit_flag = true;
	  break;
	}
      }
      
    }

    if( exit_flag ){
      break;
    }
  }
  
  return ps_profile;
}

      
int main(int argc,char* argv[]){
  std::string molet_home = argv[1];
  std::string case_path  = argv[2];
  std::string input_file = argv[3];
  Json::Value::Members jmembers;

  Json::Value root;
  std::ifstream fin(input_file,std::ifstream::in);
  fin >> root;
  fin.close();


  
  // Check if the input file is in the coolest standard, if not then return.
  if( root.isMember("standard") ){
    if( root["standard"] != "COOLEST" ){
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



  // Create new instrument
  std::string instrument_name = root["instrument"]["name"].asString() + "-" + root["instrument"]["band"].asString();
  if( Instrument::checkInstrumentExists(instrument_name) ){
    fprintf(stdout,"Instrument module '%s' already exists and will be used.\n",instrument_name.c_str());
    fprintf(stdout,"If you still want to provide it just rename the existing one in %s.\n",(molet_home+"instrument_modules/modules/").c_str());
  } else {

    // Check if name contains forbitten characters, e.g. '/'

    fprintf(stdout,"Creating new instrument: '%s'\n",instrument_name.c_str());
    Json::Value new_instrument;
    new_instrument["name"] = root["instrument"]["name"];
    new_instrument["band"] = root["instrument"]["band"];
    new_instrument["readout"] = root["instrument"]["readout_noise"];
    new_instrument["resolution"] = root["instrument"]["pixel_size"];

    Json::Value psf;
    psf["pix_x"]  = root["instrument"]["psf"]["pixels"]["num_pix_x"];
    psf["pix_y"]  = root["instrument"]["psf"]["pixels"]["num_pix_y"];
    psf["width"]  = root["instrument"]["psf"]["pixels"]["field_of_view_x"][1].asDouble() - root["instrument"]["psf"]["pixels"]["field_of_view_x"][0].asDouble();
    psf["height"] = root["instrument"]["psf"]["pixels"]["field_of_view_y"][1].asDouble() - root["instrument"]["psf"]["pixels"]["field_of_view_y"][0].asDouble();
    new_instrument["psf"] = psf;

    if( root["instrument"].isMember("wavelength") && root["instrument"].isMember("throughput") ){
      new_instrument["wavelength"] = root["instrument"]["wavelength"];
      new_instrument["throughput"] = root["instrument"]["throughput"];
    } else {
      Json::Value wavelength = Json::Value(Json::arrayValue);
      Json::Value throughput = Json::Value(Json::arrayValue);
      wavelength.append(500);
      wavelength.append(600);
      throughput.append(1);
      throughput.append(1);
      new_instrument["wavelength"] = wavelength;
      new_instrument["throughput"] = throughput;
    }
    
    Instrument::createNewInstrument(new_instrument,case_path+"psf.fits");
  }


  // Observation options
  Json::Value molet_instruments = Json::Value(Json::arrayValue);
  Json::Value instrument;
  instrument["name"] = instrument_name;
  instrument["field-of-view_xmin"] = root["observation"]["pixels"]["field_of_view_x"][0];
  instrument["field-of-view_xmax"] = root["observation"]["pixels"]["field_of_view_x"][1];
  instrument["field-of-view_ymin"] = root["observation"]["pixels"]["field_of_view_y"][0];
  instrument["field-of-view_ymax"] = root["observation"]["pixels"]["field_of_view_y"][1];
  instrument["ZP"] = root["observation"]["mag_zero_point"];
  std::string noise_type = root["observation"]["noise"]["type"].asString();
  Json::Value noise;
  if( noise_type == "InstrumentalNoise" ){
    noise["type"] = "PoissonNoise";
    noise["texp"] = root["observation"]["exposure_time"];
    noise["Msb"]  = root["observation"]["mag_sky_brightness"]; 
  } else if( noise_type == "UniformGaussianNoise" ){
    noise["type"] = "UniformGaussian";
    noise["sigma"] = root["observation"]["noise"]["std_dev"];
  } else if( noise_type == "NoiseMap" ){
    noise["type"] = "NoiseMap";
  } else if( noise_type == "NoiseRealization" ){
    noise["type"] = "NoiseRealization";
    // something
  } else {
    fprintf(stderr,"Incompatible noise type '%s'!\n",noise_type.c_str());
  }
  instrument["noise"] = noise;
  molet_instruments.append( instrument );

  
  
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
    light_model[instrument_name] = Json::Value(Json::arrayValue);
    for(int j=0;j<planes[i].size();j++){
      Json::Value lmodel = parseLight(planes[i][j]["light_model"],planes[i][j]["name"].asString(),case_path);
      for(int k=0;k<lmodel.size();k++){
	light_model[instrument_name].append( lmodel[k] );
      }
    }
    lens["light_profile"] = light_model;

    molet_lenses.append( lens );
  }
    

  
  // Get the 'light' profile of the source (last lens plane)
  Json::Value molet_source;
  Json::Value light_model;
  light_model[instrument_name] = Json::Value(Json::arrayValue);
  Json::Value last_plane = planes[planes.size()-1];
  for(int j=0;j<last_plane.size();j++){
    Json::Value lmodel = parseLight(last_plane[j]["light_model"],last_plane[j]["name"].asString(),case_path);
    for(int k=0;k<lmodel.size();k++){
      light_model[instrument_name].append( lmodel[k] );
    }
  }
  molet_source["redshift"] = redshifts.back();
  molet_source["light_profile"] = light_model;


  
  // Point source
  Json::Value molet_ps = getPointSource(last_plane);
  if( molet_ps != Json::nullValue ){


   
    Json::Value meta = root["meta"];

    if( meta.isMember("triangle_size") ){
      molet_ps["triangle_size"] = meta["triangle_size"].asDouble();
    } else {
      molet_ps["triangle_size"] = 0.3;
    }

    // Check for variability



    // Add mass-to-light in every lens
    for(int i=0;i<molet_lenses.size();i++){
      for(int j=0;j<molet_lenses[i]["light_profile"][instrument_name].size();j++){
	molet_lenses[i]["light_profile"][instrument_name][j]["mass-to-light"] = meta["mass_to_light"];
      }
    }

  }
  //std::cout << molet_ps << std::endl;

  

  // Create the file in the proper format acceptable by molet and write output
  Json::Value molet_final;
  molet_final["cosmology"] = molet_cosmo;
  molet_final["lenses"] = molet_lenses;
  molet_final["source"] = molet_source;
  molet_final["instruments"] = molet_instruments;
  if( molet_ps != Json::nullValue ){
    molet_final["point_source"] = molet_ps;
  }
  //std::cout << molet_final << std::endl;
  //std::cout << molet_instruments << std::endl;

  std::ofstream molet_input(case_path+"/molet_input.json");
  Json::StreamWriterBuilder wbuilder;
  wbuilder.settings_["precision"] = 6;
  wbuilder.settings_["indentation"] = "    ";
  std::unique_ptr<Json::StreamWriter> writer(wbuilder.newStreamWriter());
  writer->write(molet_final,&molet_input);
  molet_input.close();

  
  
  return 0;
}
