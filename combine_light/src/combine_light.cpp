#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>

#include "json/json.h"

#include "vkllib.hpp"

#include "auxiliary_functions.hpp"
#include "instruments.hpp"
#include "noise.hpp"

int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);                                   // INPUT 1: main MOLET input file
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];
  double zs = root["source"]["redshift"].asDouble();


  


  //=============== Initialization for time varying light ONLY ==================
  // Calculate the common continuous observed time and perform checks with the intrinsic and unmicro lensed light curves
  Json::Value images;
  double td_max = 0.0;
  std::vector<double> tcont;
  if( root.isMember("point_source") ){
    // Read the multiple images' parameters from JSON
    fin.open(out_path+"output/multiple_images.json",std::ifstream::in);  // INPUT FILE 2 (only for point source): properties of the multiple images
    fin >> images;
    fin.close();
    
    // Get maximum image time delay
    for(int q=0;q<images.size();q++){
      double td = images[q]["dt"].asDouble();
      if( td > td_max ){
	td_max = td;
      }
    }
  
    // Read tobs_min and tobs_max
    Json::Value tobs_json;
    fin.open(out_path+"output/tobs.json",std::ifstream::in);             // INPUT FILE 3 (only for point source): time vector properties
    fin >> tobs_json;
    fin.close();
    double tobs_max = tobs_json["tobs_max"].asDouble();
    double tobs_min = tobs_json["tobs_min"].asDouble();
    
    // Create a 'continuous' time vector: daily cadence
    // The extrinsic light curves have to be calculated in exactly the same range as tcont.
    int Ndays = (int) (ceil(tobs_max) - floor(tobs_min));
    tcont.resize(Ndays);
    for(int t=0;t<Ndays;t++){
      tcont[t] = tobs_min + t;
    }
  }
  //=============================================================================
  






  
  
  // Loop over the instruments
  // ===================================================================================================================
  // ===================================================================================================================
  for(int b=0;b<root["instruments"].size();b++){
    const Json::Value instrument = root["instruments"][b];
    std::string instrument_name = root["instruments"][b]["name"].asString();
    Instrument mycam(instrument_name,root["instruments"][b]["ZP"].asDouble(),root["instruments"][b]["noise"]);



    
    //=============== START: INITIALIZE FINAL IMAGE GRIDS ====================
    // Set output image plane in super-resolution
    double xmin = root["instruments"][b]["field-of-view_xmin"].asDouble();
    double xmax = root["instruments"][b]["field-of-view_xmax"].asDouble();
    double ymin = root["instruments"][b]["field-of-view_ymin"].asDouble();
    double ymax = root["instruments"][b]["field-of-view_ymax"].asDouble();
    int res_x = static_cast<int>(ceil((xmax-xmin)/mycam.resolution));
    int res_y = static_cast<int>(ceil((ymax-ymin)/mycam.resolution));
    int super_factor = 10;
    if( root["output_options"].isMember("super_factor") ){
      super_factor = root["output_options"]["super_factor"].asInt();
    }
    int super_res_x = super_factor*res_x;
    int super_res_y = super_factor*res_y;
    vkl::RectGrid supersim(super_res_x,super_res_y,xmin,xmax,ymin,ymax);    

    // Get the psf in super-resolution, crop it, and create convolution kernel
    mycam.preparePSF(&supersim,0.999);
    //vkl::FitsInterface::writeFits(mycam.scaled_psf->Nx,mycam.scaled_psf->Ny,mycam.scaled_psf->z,out_path + "output/supersampled_psf.fits");
    //vkl::FitsInterface::writeFits(mycam.cropped_psf->Nx,mycam.cropped_psf->Ny,mycam.cropped_psf->z,out_path + "output/cropped_psf.fits");
    //vkl::FitsInterface::writeFits(supersim.Nx,supersim.Ny,mycam.kernel,out_path + "output/kernel_psf.fits");

    vkl::RectGrid super_extended = vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/"+instrument_name+"_lensed_image_super.fits");        // INPUT FILE 4: super-resolved lensed source
    vkl::RectGrid super_lens_light = vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/"+instrument_name+"_lens_light_super.fits");        // INPUT FILE 5: super-resolved lens light

    // Create static image (lens light and extended lensed source)
    vkl::RectGrid obs_static = createObsStatic(&mycam,&super_extended,&super_lens_light,res_x,res_y,out_path); // This is in units of flux!

    // Assign noise grid
    mycam.noise->setGrid(&obs_static);
    
    // All the static light components have been created.
    // The resulting observed image (not super-resolved) is "obs_static".
    // If there is no time dimension required, the code just outputs the obs_static image and stops.
    // But further work is done when it comes to time varying images in three nested loops:
    // - one over all the intrinsic variability light curves,
    // - one over all the extrinsic variability light curves,
    // - and one over the observed time.
    //===============   END: INITIALIZE FINAL IMAGE GRIDS ====================
    



    //=============== CREATE THE TIME VARYING LIGHT ====================
    if( root.isMember("point_source") ){

      std::vector<double> tobs;
      for(int t=0;t<instrument["time"].size();t++){
	tobs.push_back(instrument["time"][t].asDouble());
      }
      std::string ex_type = root["point_source"]["variability"]["extrinsic"]["type"].asString();
      
      

      //======================================== Intrinsic variability ========================================================
      std::vector<LightCurve*> LC_intrinsic;
      // Read intrinsic light curve(s) from JSON and apply conversions: from rest frame to observer's frame, from mag to intensity and scale if needed.
      Json::Value intrinsic_lc_json = readLightCurvesJson("intrinsic",root["point_source"]["variability"]["intrinsic"]["type"].asString(),instrument_name,in_path,out_path);
      double scale_factor = 1.0;
      if( root["point_source"]["variability"]["intrinsic"].isMember("scale_factor") ){
	scale_factor = root["point_source"]["variability"]["intrinsic"]["scale_factor"].asDouble();
      }
      LC_intrinsic = conversions(intrinsic_lc_json,zs,scale_factor,mycam.ZP);
      //=======================================================================================================================

      
      
      //======================================== Unmicrolensed flux ===========================================================
      // Check for unmicrolensed variability, read light curve(s) from JSON and apply conversions: from rest frame to observer's frame, from mag to intensity and scale if needed.
      bool unmicro = false;
      std::vector<LightCurve*> LC_unmicro;
      if( root["point_source"]["variability"].isMember("unmicro") ){
	unmicro = true;
	Json::Value unmicro_lc_json = readLightCurvesJson("unmicro",root["point_source"]["variability"]["unmicro"]["type"].asString(),instrument_name,in_path,out_path);
	double scale_factor = 1.0;
	if( root["point_source"]["variability"]["unmicro"].isMember("scale_factor") ){
	  scale_factor = root["point_source"]["variability"]["unmicro"]["scale_factor"].asDouble();
	} 
	LC_unmicro = conversions(unmicro_lc_json,zs,scale_factor,mycam.ZP);
      }
      //=======================================================================================================================            

      

      //=============================================== Microlensing ==========================================================      
      // Read extrinsic light curve(s) from JSON
      Json::Value extrinsic_lc = readLightCurvesJson("extrinsic",root["point_source"]["variability"]["extrinsic"]["type"].asString(),instrument_name,in_path,out_path);

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



      //==================================== Configure the PSF for the point source ===========================================
      // Perturb the PSF at each image location
      std::vector<Instrument*> instrument_list(images.size());
      for(int q=0;q<images.size();q++){
	//TransformPSF* dum = new TransformPSF(images[q]["x"].asDouble(),images[q]["y"].asDouble(),0.0,false,false);
	//transPSF[q] = dum;
	instrument_list[q] = &mycam; // just copy the same instrument pointer per image
      }
      // Set the PSF related offsets for each image
      std::vector<offsetPSF> PSFoffsets(images.size());
      for(int q=0;q<images.size();q++){
	PSFoffsets[q] = instrument_list[q]->offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&supersim);
      }
      FILE* fh = fopen((out_path+"output/psf_locations.dat").c_str(),"w");
      for(int q=0;q<images.size();q++){
	PSFoffsets[q].printFrame(fh,supersim.Nx,supersim.Ny,supersim.xmin,supersim.xmax,supersim.ymin,supersim.ymax);
      }
      fclose(fh);
      // Calculate the appropriate PSF sums
      std::vector<double> psf_partial_sums(images.size());
      for(int q=0;q<images.size();q++){
	psf_partial_sums[q] = instrument_list[q]->sumPSF(&PSFoffsets[q]);
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
		combineInExUnSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],LC_unmicro[lc_in],cont_LC[q]);
		combineInExUnSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_extrinsic[q][lc_ex],LC_unmicro[lc_in],samp_LC[q]);
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
		combineInUnSignals(td,macro_mag,tcont,LC_intrinsic[lc_in],LC_unmicro[lc_in],cont_LC[q]);
		combineInUnSignals(td,macro_mag,tobs,LC_intrinsic[lc_in],LC_unmicro[lc_in],samp_LC[q]);
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
	  outputLightCurvesJson(cont_LC,out_path+mock+"/"+instrument_name+"_LC_continuous.json");
	  outputLightCurvesJson(samp_LC,out_path+mock+"/"+instrument_name+"_LC_sampled.json");
	  for(int q=0;q<images.size();q++){
	    delete(cont_LC[q]);
	    // we don't delete samp_LC yet because it is needed in case of outputing cutouts
	  }
	  // *********************** End of product **************************************************************
	    
	    
	  // *********************** Product: Observed sampled cut-outs (images) *****************************
	  if( root["point_source"]["output_cutouts"].asBool() ){
	    writeAllCutouts(tobs,images,samp_LC,&super_extended,&super_lens_light,&mycam,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y,mock,instrument_name,out_path);
	  }
	  // *********************** End of product **************************************************	    
	    
	  for(int q=0;q<images.size();q++){
	    delete(samp_LC[q]);
	  }
	}
	// Loop over extrinsic light curves ends here
	//0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
      }
      // Loop over intrinsic light curves ends here
      //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
	
      

      // *********************** Product: Just one cutout with only the macromagnification *****************************
      // For each image we use the same brightness value.
      // This value defines a flux from the equation F = 10^(-0.4*(M-ZP)). For M=ZP-7.5 (-5) I get a flux of 1000 (100).
      double total_image_flux = 0.0;
      double area = supersim.step_x*supersim.step_y;
      double single_image_flux = 10000.0; // In units of integrated (total) flux per image
      std::vector<double> image_signal(images.size());
      for(int q=0;q<images.size();q++){
	double macro_mag = fabs(images[q]["mag"].asDouble());
	image_signal[q] = single_image_flux*macro_mag;
	total_image_flux += image_signal[q]*area; // The area of a PS image is one pixel
      }
      double total_image_flux_mag = -2.5*log10(total_image_flux) + mycam.ZP;
      std::cout << "Total PS flux from all images: " << total_image_flux << " " << total_image_flux_mag << std::endl;
      vkl::RectGrid obs_ps_light = addPSLight(&obs_static,&supersim,image_signal,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y,mycam.ZP); // This is in flux units!
      convertFromFlux(&obs_ps_light,mycam.ZP);

      // Adding time-dependent noise here
      mycam.noise->initializeFromData(&obs_ps_light);
      mycam.noise->calculateNoise();
      //mycam.noise->addNoise(&obs_ps_light);
      mycam.noise->outputNoiseProperties(out_path + "output/",instrument_name+"_ps_macro");

      // Output the observed static image with a PS with macromagnification only
      writeCutout(&obs_ps_light,out_path+"output/OBS_"+instrument_name+"_ps_macro.fits");
      //vkl::FitsInterface::writeFits(obs_static.Nx,obs_static.Ny,obs_static.z,fname);
      // *********************** End of product ************************************************************************


      
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
    //================= END:CREATE THE TIME VARYING LIGHT ====================

    


    //=============== START: CREATE STATIC LIGHT ====================


    double total_flux,total_flux_mag;
    obs_static.integrate(total_flux,total_flux_mag,mycam.ZP);
    std::cout << "Static: " << total_flux << " " << total_flux_mag << std::endl;

      
    // A single image without any point source
    convertFromFlux(&obs_static,mycam.ZP);
    
    // Adding noise here and output realization
    mycam.noise->initializeFromData(&obs_static);
    mycam.noise->calculateNoise();
    //mycam.noise->addNoise(&obs_static); // <------ BE CAREFUL: Noise is added directly to the obs_static pixels.
    mycam.noise->outputNoiseProperties(out_path + "output/",instrument_name);
    
    // Output the observed static image
    writeCutout(&obs_static,out_path + "output/OBS_" + instrument_name + ".fits");
    //===============   END: CREATE STATIC LIGHT ====================



    
  }
  // Loop over the instruments ends here
  // ===================================================================================================================
  // ===================================================================================================================

  
  return 0;
}
