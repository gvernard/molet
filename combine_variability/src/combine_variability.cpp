




int main(int argc,char* argv[]){

  /*
    Requires:
    - main input
    - multiple_images.json
    - tobs.json
    - intrinsic light curves
    - extrinsic light curves
  */


  //=============== BEGIN:PARSE INPUT ===========================================
  std::ifstream fin;

  std::string root_file = argv[1];
  std::string in_path   = argv[2];
  std::string out_path  = argv[3];

  
  // Read the main projection parameters
  Json::Value root;
  fin.open(root_file,std::ifstream::in);                                         // INPUT FILE 1: main MOLET input file
  fin >> root;
  fin.close();
  double zs = root["source"]["redshift"].asDouble();
  std::string ex_type = root["point_source"]["variability"]["extrinsic"]["type"].asString();
  std::string in_type = root["point_source"]["variability"]["intrinsic"]["type"].asString();
  

  // Read the multiple images' parameters from JSON
  Json::Value images;
  fin.open(out_path+"output/multiple_images.json",std::ifstream::in);            // INPUT FILE 2: properties of the multiple images
  fin >> images;
  fin.close();

  
  // Get maximum image time delay
  double td_max = 0.0;
  for(int q=0;q<images.size();q++){
    double td = images[q]["dt"].asDouble();
    if( td > td_max ){
      td_max = td;
    }
  }

  
  // Get tobs_min and tobs_max (taking into account all instruments)
  Json::Value tobs_json;
  fin.open(out_path+"output/tobs.json",std::ifstream::in);                       // INPUT FILE 3: time vector properties
  fin >> tobs_json;
  fin.close();
  double tobs_max = tobs_json["tobs_max"].asDouble();
  double tobs_min = tobs_json["tobs_min"].asDouble();


  // Create the common 'continuous' time vector: daily cadence
  // The extrinsic light curves have to be calculated in exactly the same range as tcont.
  int Ndays = (int) (ceil(tobs_max) - floor(tobs_min));
  std::vector<double> tcont(Ndays);
  for(int t=0;t<Ndays;t++){
    tcont[t] = tobs_min + t;
  }
  //================= END:PARSE INPUT ===========================================  

  


  
  // But further work is done when it comes to time varying images in three nested loops:
  // - one over all the intrinsic variability light curves,
  // - one over all the extrinsic variability light curves,
  // - and one over the observed time.



      
  // Loop over the instruments
  // ===================================================================================================================
  // ===================================================================================================================
  for(int b=0;b<root["instruments"].size();b++){
    const Json::Value instrument = root["instruments"][b];
    std::string instrument_name = root["instruments"][b]["name"].asString();
    Instrument mycam(instrument_name,root["instruments"][b]["ZP"].asDouble(),root["instruments"][b]["noise"],true);  
    std::vector<double> tobs;
    for(int t=0;t<instrument["time"].size();t++){
      tobs.push_back(instrument["time"][t].asDouble());
    }

    
    //======================================== Intrinsic variability ========================================================
    std::vector<LightCurve*> LC_intrinsic;
    // Read intrinsic light curve(s) from JSON and apply conversions: from rest frame to observer's frame, from mag to intensity and scale if needed.
    Json::Value intrinsic_lc_json = readLightCurvesJson("intrinsic",in_type,instrument_name,in_path,out_path);       // INPUT FILE 4 (per instrument): intrinsic light curves
    double scale_factor = 1.0;
    if( root["point_source"]["variability"]["intrinsic"].isMember("scale_factor") ){
      scale_factor = root["point_source"]["variability"]["intrinsic"]["scale_factor"].asDouble();
    }
    LC_intrinsic = conversions(intrinsic_lc_json,zs,scale_factor,mycam.ZP);
    //=======================================================================================================================

      
      
    //======================================== Unmicrolensed flux ===========================================================
    // Check for unmicrolensed variability, read light curve(s) from JSON and apply conversions: from rest frame to observer's frame, from mag to intensity and scale if needed.
    bool unmicro = false;
    double unmicro_ratio;
    std::vector<LightCurve*> LC_unmicro;
    if( root["point_source"]["variability"].isMember("unmicro") ){
      unmicro = true;
      
      /*
      // This is to read in a custom unmicrolensed light curve
      Json::Value unmicro_lc_json = readLightCurvesJson("unmicro",root["point_source"]["variability"]["unmicro"]["type"].asString(),instrument_name,in_path,out_path);
      double scale_factor = 1.0;
      if( root["point_source"]["variability"]["unmicro"].isMember("scale_factor") ){
      scale_factor = root["point_source"]["variability"]["unmicro"]["scale_factor"].asDouble();
      } 
      LC_unmicro = conversions(unmicro_lc_json,zs,scale_factor,mycam.ZP);
      */
      
      LC_unmicro.resize(LC_intrinsic.size());
      unmicro_ratio = root["point_source"]["variability"]["unmicro"][instrument_name]["flux_ratio"].asDouble();
      BaseLagKernel* lag = nullptr;
      
      // Set the lag function kernel
      std::string lag_type = root["point_source"]["variability"]["unmicro"][instrument_name]["type"].asString();
      if( lag_type == "top-hat" ){
	double radius = root["point_source"]["variability"]["unmicro"][instrument_name]["pars"]["radius"].asDouble();
	lag = new TopHatKernel(radius);
      } else if( lag_type == "delta" ){
	double t_peak = root["point_source"]["variability"]["unmicro"][instrument_name]["pars"]["t_peak"].asDouble();
	lag = new DeltaKernel(t_peak);
      } else {
	fprintf(stderr,"Unknown Lag Kernel model: %s\n",lag_type.c_str());
	return 1;
      }
      
      // Get the unmicrolensed light curve that corresponds to each intrinsic curve
      for(int lc_in=0;lc_in<LC_intrinsic.size();lc_in++){
	LC_unmicro[lc_in] = new LightCurve(LC_intrinsic[lc_in]->time);
	
	// Evaluate the lag kernel on the same time vector as the given intrinsic light curve 
	int N = LC_intrinsic[lc_in]->time.size();
	std::vector<double> kernel(N);
	lag->getKernel(LC_intrinsic[lc_in]->time,kernel);
	
	// Perform the convolution
	std::vector<double> unmicro_signal(N);
	for(int n=0;n<N;n++){
	  double sum = 0.0;
	  for(int m=0;m<N;m++){
	    int index_in = n - m; // index to the intrinsic light curve
	    if( index_in < 0 ){
	      index_in += N; // wrap the convolution to the end of the intrinsic light curve
	    }
	    sum += kernel[m]*LC_intrinsic[lc_in]->signal[index_in];
	  }
	  LC_unmicro[lc_in]->signal[n] = sum;
	}
      }
      
      outputLightCurvesJson(LC_unmicro,mycam.ZP,out_path+"output/"+instrument_name+"_LC_unmicro.json");
    }
    //=======================================================================================================================

      

    //=============================================== Microlensing ==========================================================
    // Read extrinsic light curve(s) from JSON
    Json::Value extrinsic_lc = readLightCurvesJson("extrinsic",ex_type,instrument_name,in_path,out_path);            // INPUT FILE 5 (per instrument): extrinsic light curves
    
    // Split the indices of the multiple images to those that have microlensing light curves and those that not.
    std::vector<int> images_micro;
    std::vector<int> images_no_micro;
    for(int q=0;q<extrinsic_lc.size();q++){
      if( extrinsic_lc[q].size() > 0 ){
	images_micro.push_back(q);
      } else {
	images_no_micro.push_back(q);
      }
    }      
    int N_ex = extrinsic_lc[images_micro[0]].size();
    
    // No need to dilate the time because the effective velocity has length units in the source plane and time in the oberver's frame.
    // Also, the starting time is tobs_min for quasars and tin_min for supernovae, so no need to set it either.
    std::vector< std::vector<LightCurve*> > LC_extrinsic(images.size());
    for(int i=0;i<images_micro.size();i++){
      int q = images_micro[i];
      LC_extrinsic[q].resize(N_ex);
      for(int j=0;j<N_ex;j++){
	LC_extrinsic[q][j] = new LightCurve(extrinsic_lc[q][j]);
      }
    }
    //=======================================================================================================================      



    // Loop over intrinsic light curves
    //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
    for(int lc_in=0;lc_in<LC_intrinsic.size();lc_in++){
      // Loop over extrinsic light curves
      //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
      for(int lc_ex=0;lc_ex<N_ex;lc_ex++){ // use the first multiple image to get the number of extrinsic light curves
	
	// Output directory
	//char buffer[15];
	char buffer[30];
	sprintf(buffer,"mock_%04d_%04d",lc_in,lc_ex);
	std::string mock = buffer;
	
	// *********************** Product: Observed continuous and sampled light curves ***********************
	// cont_LC and samp_LC contain the light curves for ALL the images (i.e. with or without microlensing).
	// The first is based on the continuous time, tcont, and the second on the observed (sampled) time, tobs.
	// Units are: t in [days], and signal in [flux] NOT mag.
	std::vector<LightCurve*> cont_LC(images.size());
	std::vector<LightCurve*> samp_LC(images.size());
	for(int q=0;q<images.size();q++){
	  cont_LC[q] = new LightCurve(tcont);
	  samp_LC[q] = new LightCurve(tobs);
	}
	
	// // Images with microlensing: Combine extrinsic and intrinsic light curves with the given time delay at the same starting time
	// for(int i=0;i<images_micro.size();i++){
	//   int q = images_micro[i];
	//   double td = td_max - images[q]["dt"].asDouble();
	//   double macro_mag = fabs(images[q]["mag"].asDouble());
	//   combineSupernovaInExSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],cont_LC[q]);
	//   combineSupernovaInExSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],samp_LC[q]);
	// }
	
	// // Images without microlensing: Just shift the intrinsic light curve by the time delay
	// for(int i=0;i<images_no_micro.size();i++){
	//   int q = images_no_micro[i];
	//   double td = td_max - images[q]["dt"].asDouble();
	//   double macro_mag = fabs(images[q]["mag"].asDouble());
	//   justSupernovaInSignal(td,macro_mag,tcont,LC_intrinsic[lc_in],cont_LC[q]);
	//   justSupernovaInSignal(td,macro_mag,tobs,LC_intrinsic[lc_in],samp_LC[q]);
	// }	    
	
	
	if( ex_type == "custom" || ex_type == "moving_fixed_source" || ex_type == "moving_fixed_source_custom" ){
	  
	  // Calculate the combined light curve for each image with microlensing
	  for(int i=0;i<images_micro.size();i++){
	    int q = images_micro[i];
	    double td = td_max - images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    if( unmicro ){
	      combineInExUnSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],LC_unmicro[lc_in],cont_LC[q],unmicro_ratio);
	      combineInExUnSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],LC_unmicro[lc_in],samp_LC[q],unmicro_ratio);
	    } else {
	      combineInExSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],cont_LC[q]);
	      combineInExSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],samp_LC[q]);
	    }
	  }
	  
	  // Calculate the combined light curve for each image that doesn't have any microlensing
	  for(int i=0;i<images_no_micro.size();i++){
	    int q = images_no_micro[i];
	    double td = td_max - images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    if( unmicro ){
	      combineInUnSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_unmicro[lc_in],cont_LC[q],unmicro_ratio);
	      combineInUnSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_unmicro[lc_in],samp_LC[q],unmicro_ratio);
	    } else {
	      justOneSignal(td,macro_mag,tcont,LC_intrinsic[lc_in],cont_LC[q]);
	      justOneSignal(td,macro_mag,tobs,LC_intrinsic[lc_in],samp_LC[q]);
	    }
	  }
	  
	} else if( ex_type == "expanding_source" ){
	  
	  // Calculate the combined light curve for each image with microlensing
	  for(int i=0;i<images_micro.size();i++){
	    int q = images_micro[i];
	    double td = images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    combineSupernovaInExSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],cont_LC[q]);
	    combineSupernovaInExSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],samp_LC[q]);
	  }
	  
	  // Calculate the combined light curve for each image that doesn't have any microlensing
	  for(int i=0;i<images_no_micro.size();i++){
	    int q = images_no_micro[i];
	    double td = images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    justSupernovaInSignal(td,macro_mag,tcont,LC_intrinsic[lc_in],cont_LC[q]);
	    justSupernovaInSignal(td,macro_mag,tobs,LC_intrinsic[lc_in],samp_LC[q]);
	  }
	  
	} else if( ex_type == "moving_variable_source" ){
	  
	  // Calculate the combined light curve for each image with microlensing
	  for(int i=0;i<images_micro.size();i++){
	    int q = images_micro[i];
	    double td = td_max - images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    combineInExSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],cont_LC[q]);
	    combineInExSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],samp_LC[q]);
	  }
	  
	  // Calculate the combined light curve for each image that doesn't have any microlensing
	  for(int i=0;i<images_no_micro.size();i++){
	    int q = images_no_micro[i];
	    double td = td_max - images[q]["dt"].asDouble();
	    double macro_mag = fabs(images[q]["mag"].asDouble());
	    justOneSignal(td,macro_mag,tcont,LC_intrinsic[lc_in],cont_LC[q]);
	    justOneSignal(td,macro_mag,tobs,LC_intrinsic[lc_in],samp_LC[q]);
	  }
	  
	} else {
	  
	  // This is not really needed because it should have been taken care of in the checks.
	  fprintf(stderr,"Unknown variability model: %s\n",ex_type.c_str());
	  return 1;
	  
	}
	
	
	// Write json light curves and clean up
	outputLightCurvesJson(cont_LC,mycam.ZP,out_path+mock+"/"+instrument_name+"_LC_continuous.json");
	outputLightCurvesJson(samp_LC,mycam.ZP,out_path+mock+"/"+instrument_name+"_LC_sampled.json");
	for(int q=0;q<images.size();q++){
	  delete(cont_LC[q]);
	  delete(samp_LC[q]);
	}
	// *********************** End of product **************************************************************

	
      }
      // Loop over extrinsic light curves ends here
      //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
    }
    // Loop over intrinsic light curves ends here
    //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
    


    
    for(int lc_in=0;lc_in<LC_intrinsic.size();lc_in++){
      delete(LC_intrinsic[lc_in]);
    }
    
    if( unmicro ){
      for(int lc_in=0;lc_in<LC_intrinsic.size();lc_in++){
	delete(LC_unmicro[lc_in]);
      }
    }
    
    for(int i=0;i<images_micro.size();i++){
      int q = images_micro[i];
      for(int lc_ex=0;lc_ex<N_ex;lc_ex++){
	delete(LC_extrinsic[q][lc_ex]);
      }
    }


  }
  // Loop over the instruments ends here
  // ===================================================================================================================
  // ===================================================================================================================
  
  

  return 0;
}
