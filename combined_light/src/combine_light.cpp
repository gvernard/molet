#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "vkllib.hpp"

#include "auxiliary_functions.hpp"
#include "mask_functions.hpp"
#include "instruments.hpp"
#include "noise.hpp"

int main(int argc,char* argv[]){

  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[2];
  std::string out_path = argv[3];

  std::string cut_out_scale;
  if( root.isMember("output_options") ){
    cut_out_scale = root["output_options"]["cut_outs"]["scale"].asString();
  } else {
    cut_out_scale = "mag";
  }

  
  // Loop over the instruments
  // ===================================================================================================================
  // ===================================================================================================================
  for(int b=0;b<root["instruments"].size();b++){
    const Json::Value instrument = root["instruments"][b];
    std::string instrument_name = root["instruments"][b]["name"].asString();
    Instrument mycam(instrument_name,root["instruments"][b]["noise"]);
    
    // Set output image plane in super-resolution
    double xmin = root["instruments"][b]["field-of-view_xmin"].asDouble();
    double xmax = root["instruments"][b]["field-of-view_xmax"].asDouble();
    double ymin = root["instruments"][b]["field-of-view_ymin"].asDouble();
    double ymax = root["instruments"][b]["field-of-view_ymax"].asDouble();
    int res_x = static_cast<int>(ceil((xmax-xmin)/mycam.resolution));
    int res_y = static_cast<int>(ceil((ymax-ymin)/mycam.resolution));
    int super_res_x = 10*res_x;
    int super_res_y = 10*res_y;
    RectGrid mysim(super_res_x,super_res_y,xmin,xmax,ymin,ymax);
    

    
    // Get the psf in super-resolution, crop it, and create convolution kernel
    mycam.interpolatePSF(&mysim);
    mycam.cropPSF(0.99);
    mycam.createKernel(mysim.Nx,mysim.Ny);
    
    
    // Create the fixed extended lensed light
    RectGrid* extended = new RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/lensed_image_super.fits");
    mycam.convolve(extended);
    //extended->writeImage(output+"psf_lensed_image_super.fits");
    
    // Create the fixed lens galaxy light
    RectGrid* lens_light = new RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax,out_path+"output/lens_light_super.fits");
    mycam.convolve(lens_light);
    //lens_light->writeImage(output+"psf_lens_light_super.fits");
    

    // Combined light of the observed base image (binned from 'super' to observed resolution)
    RectGrid* base = new RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax); 
    for(int i=0;i<base->Nz;i++){
      base->z[i] = lens_light->z[i] + extended->z[i];
    }
    delete(extended);
    delete(lens_light);
    RectGrid* obs_base = base->embeddedNewGrid(res_x,res_y,"integrate");
    delete(base);



    // All the static light components have been created.
    // The resulting observed image (not super-resolved) is "obs_base".
    // If there is no time dimension required, the code just outputs the obs_base image and stops.
    // But further work is done when it comes to time varying images in three nested loops:
    // - one over all the intrinsic variability light curves,
    // - one over all the extrinsic variability light curves,
    // - and one over the observed time.

    

    if( !root.isMember("point_source") ){
      //=============== CREATE A SINGLE STATIC IMAGE ====================

      // Adding noise here
      mycam.noise->addNoise(obs_base);
      
      // Convert to magnitudes
      if( cut_out_scale == "mag" ){
	for(int i=0;i<obs_base->Nz;i++){
	  obs_base->z[i] = -2.5*log10(obs_base->z[i]);
	}
      }

      // Output the observed base image
      FitsInterface::writeFits(obs_base->Nx,obs_base->Ny,obs_base->z,out_path + "output/OBS_" + instrument_name + ".fits");
      delete(obs_base);
      
    } else {
      //=============== CREATE THE TIME VARYING LIGHT ====================
      
      // Read the multiple images' parameters from JSON
      Json::Value images;
      fin.open(out_path+"output/multiple_images.json",std::ifstream::in);
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

      // Get observed time vector
      std::vector<double> tobs;
      for(int t=0;t<instrument["time"].size();t++){
	tobs.push_back(instrument["time"][t].asDouble());
      }
      double tobs_t0   = tobs[0];
      double tobs_tmax = tobs.back();
      double tobs_Dt   = tobs_tmax - tobs_t0;

      // Create a 'continuous' time vector: daily cadence
      int Ndays = (int) ceil(tobs_Dt);
      std::vector<double> tcont(Ndays);
      for(int t=0;t<Ndays;t++){
	tcont[t] = tobs_t0 + t;
      }

      // Quick check on time delay and observing time compatibility
      if( td_max > tobs_tmax ){
	printf("Observing period (%f days) shorter than the maximum time delay (%f days).\n",tobs_tmax,td_max);
	printf("Increase the observing time period!!!\n");
	return 1;
      }

      // Read intrinsic light curve(s) from JSON
      Json::Value intrinsic_lc;
      if( root["point_source"]["variability"]["intrinsic"]["type"].asString() == "custom" ){
	fin.open(in_path+"/input_files/"+instrument_name+"_LC_intrinsic.json",std::ifstream::in);
      } else {
	fin.open(out_path+"output/"+instrument_name+"_LC_intrinsic.json",std::ifstream::in);
      }
      fin >> intrinsic_lc;
      fin.close();
      int N_in = intrinsic_lc.size();

      // Process the intrinsic light curves
      std::vector<LightCurve*> LC_intrinsic(N_in);
      for(int lc_in=0;lc_in<N_in;lc_in++){
	LC_intrinsic[lc_in] = new LightCurve(intrinsic_lc[lc_in]);

	// Convert from magnitudes to intensities
	for(int i=0;i<LC_intrinsic[lc_in]->signal.size();i++){
	  LC_intrinsic[lc_in]->signal[i] = pow(10.0,-0.4*LC_intrinsic[lc_in]->signal[i]);
	}
	
	// Check time limitations
	double tmax_intrinsic = LC_intrinsic[lc_in]->time.back();
	if( (td_max+tobs_tmax) > tmax_intrinsic ){
	  printf("Intrinsic light curve %i duration (%f days) is shorter than the maximum time delay plus the observing period (%f + %f days).\n",lc_in,tmax_intrinsic,td_max,tobs_tmax);
	  int i=0;
	  while( tobs[i] < (tmax_intrinsic-td_max) ){
	    i++;
	  }
	  tobs.resize(i);
	  printf("Observing period is truncated to %f days!!!\n",tobs_tmax);
	}
      }

      // Check for unmicrolensed variability and read unmicrolensed light curves from JSON
      bool unmicro = false;
      if( root["point_source"]["variability"].isMember("unmicro") ){
	unmicro = true;
      }

      std::vector<LightCurve*> LC_unmicro;
      if( unmicro ){
	Json::Value unmicro_lc;
	fin.open(in_path+"/input_files/"+instrument_name+"_LC_unmicro.json",std::ifstream::in);
	fin >> unmicro_lc;
	fin.close();
	int N_un = unmicro_lc.size();
	
	if( N_un != N_in ){
	  fprintf(stderr,"Number of intrinsic (%d) and unmicrolensed (%d) light curves should be the same!\n",N_in,N_un);
	  return -1;
	}
	
	// Process the unmicrolensed light curves
	LC_unmicro.resize(N_in);
	for(int lc_in=0;lc_in<N_in;lc_in++){
	  LC_unmicro[lc_in] = new LightCurve(unmicro_lc[lc_in]);
	  
	  // Convert from magnitudes to intensities
	  for(int i=0;i<LC_unmicro[lc_in]->signal.size();i++){
	    LC_unmicro[lc_in]->signal[i] = pow(10.0,-0.4*LC_unmicro[lc_in]->signal[i]);
	  }
	  
	  // Check time limitations
	  double tmax_unmicro = LC_unmicro[lc_in]->time.back();
	  if( (td_max+tobs_tmax) > tmax_unmicro ){
	    printf("Unmicrolensed light curve %i duration (%f days) is shorter than the maximum time delay plus the observing period (%f + %f days).\n",lc_in,tmax_unmicro,td_max,tobs_tmax);
	    int i=0;
	    while( tobs[i] < (tmax_unmicro-td_max) ){
	      i++;
	    }
	    tobs.resize(i);
	    printf("Observing period is truncated to %f days!!!\n",tobs_tmax);
	  }
	}
      }
      
      // Read extrinsic light curve(s) from JSON
      Json::Value extrinsic_lc;
      if( root["point_source"]["variability"]["extrinsic"]["type"].asString() == "custom" ){
	fin.open(in_path+"/input_files/"+instrument_name+"_LC_extrinsic.json",std::ifstream::in);
      } else {
	fin.open(out_path+"output/"+instrument_name+"_LC_extrinsic.json",std::ifstream::in);
      }
      fin >> extrinsic_lc;
      fin.close();
      int N_ex;
      for(int q=0;q<extrinsic_lc.size();q++){
	if( extrinsic_lc[q].size() > 0 ){
	  N_ex = extrinsic_lc[q].size();
	  break;
	}
      }

      // Process the extrinsic light curves
      std::vector< std::vector<LightCurve*> > LC_extrinsic(images.size());
      for(int q=0;q<images.size();q++){
	LC_extrinsic[q].resize(N_ex);
      }
      for(int q=0;q<images.size();q++){
	for(int lc_ex=0;lc_ex<N_ex;lc_ex++){
	  if( extrinsic_lc[q].size() > 0 ){
	    LC_extrinsic[q][lc_ex] = new LightCurve(extrinsic_lc[q][lc_ex]);
	    // Ex. light curves' time begins at 0, so I need to add t0 (of both tobs and tcont) so that the time vectors match
	    for(int t=0;t<LC_extrinsic[q][lc_ex]->time.size();t++){
	      LC_extrinsic[q][lc_ex]->time[t] += tobs_t0;
	    }
	  } else {
	    LC_extrinsic[q][lc_ex] = new LightCurve();
	  }
	}
      }

      
      // Configure the PSF for the point source
      // Perturb the PSF at each image location
      std::vector<Instrument*> Instrument_list(images.size());
      for(int q=0;q<images.size();q++){
	//TransformPSF* dum = new TransformPSF(images[q]["x"].asDouble(),images[q]["y"].asDouble(),0.0,false,false);
	//transPSF[q] = dum;
	Instrument_list[q] = &mycam; // just copy the same instrument pointer per image
      }
      // Set the PSF related offsets for each image
      std::vector<offsetPSF> PSFoffsets(images.size());
      for(int q=0;q<images.size();q++){
	PSFoffsets[q] = Instrument_list[q]->offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&mysim);
      }
      FILE* fh = fopen((out_path+"output/psf_locations.dat").c_str(),"w");
      for(int q=0;q<images.size();q++){
	PSFoffsets[q].printFrame(fh,mysim.Nx,mysim.Ny,mysim.width,mysim.height);
      }
      fclose(fh);
      // Calculate the appropriate PSF sums
      std::vector<double> psf_partial_sum(images.size());
      for(int q=0;q<images.size();q++){
	double sum = 0.0;
	for(int i=0;i<PSFoffsets[q].ni;i++){
	  for(int j=0;j<PSFoffsets[q].nj;j++){
	    int index_psf = i*Instrument_list[q]->cropped_psf->Nx + j;
	    sum += Instrument_list[q]->cropped_psf->z[index_psf];
	  }
	}
	psf_partial_sum[q] = sum;
      }
      
      
      // Loop over intrinsic light curves
      //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
      for(int lc_in=0;lc_in<N_in;lc_in++){
	// Loop over extrinsic light curves
	//0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
	for(int lc_ex=0;lc_ex<N_ex;lc_ex++){ // use the first multiple image to get the number of extrinsic light curves

	  // Output directory
	  char buffer[15];
	  sprintf(buffer,"mock_%04d_%04d",lc_in,lc_ex);
	  std::string mock = buffer;
	  //std::cout << mock << std::endl;

	  // *********************** Product: Observed continuous light curves ***********************
	  std::vector<LightCurve*> cont_LC(images.size());
	  for(int q=0;q<images.size();q++){
	    cont_LC[q] = new LightCurve(tcont);
	  }

	  // Calculate the combined light curve for each image
	  for(int q=0;q<images.size();q++){
	    double macro_mag = abs(images[q]["mag"].asDouble());
	    LightCurve* cont_LC_intrinsic = new LightCurve(tcont);
	    LC_intrinsic[lc_in]->interpolate(cont_LC_intrinsic,td_max - images[q]["dt"].asDouble());

	    if( unmicro ){
	      // === Combining three signals: intrinsic, intrinsic unmicrolensed, and extrinsic
	      LightCurve* cont_LC_unmicro = new LightCurve(tcont);
	      LC_unmicro[lc_in]->interpolate(cont_LC_unmicro,td_max - images[q]["dt"].asDouble());
	      
	      if( LC_extrinsic[q][lc_ex]->time.size() > 0 ){ // Check if multiple image does not have a corresponding extrinsic light curve (i.e. a maximum image without a magnification map)
		LC_extrinsic[q][lc_ex]->interpolate(cont_LC[q],0.0);
		for(int t=0;t<tcont.size();t++){
		  cont_LC[q]->signal[t] = macro_mag*(cont_LC[q]->signal[t]*cont_LC_intrinsic->signal[t] + cont_LC_unmicro->signal[t]);
		}
	      } else {
		for(int t=0;t<tcont.size();t++){
		  cont_LC[q]->signal[t] = macro_mag*(cont_LC_intrinsic->signal[t] + cont_LC_unmicro->signal[t]); // this line includes only the intrinsic signal and excludes microlensing
		}
	      }

	      delete(cont_LC_unmicro);
	    } else {
	      // === Combining two signals: intrinsic and extrinsic
	      if( LC_extrinsic[q][lc_ex]->time.size() > 0 ){ // Check if multiple image does not have a corresponding extrinsic light curve (i.e. a maximum image without a magnification map)
		LC_extrinsic[q][lc_ex]->interpolate(cont_LC[q],0.0);
		for(int t=0;t<tcont.size();t++){
		  cont_LC[q]->signal[t] = cont_LC[q]->signal[t] * macro_mag * cont_LC_intrinsic->signal[t];
		}
	      } else {
		for(int t=0;t<tcont.size();t++){
		  cont_LC[q]->signal[t] = macro_mag * cont_LC_intrinsic->signal[t]; // this line includes only the intrinsic signal and excludes microlensing
		}
	      }
	    }

	    delete(cont_LC_intrinsic);
	  }

	  // Write json light curves
	  outputLightCurvesJson(cont_LC,out_path+mock+"/"+instrument_name+"_LC_continuous.json");

	  // Clean up
	  for(int q=0;q<images.size();q++){
	    delete(cont_LC[q]);
	  }
	  // *********************** End of product **************************************************



	  // *********************** Product: Observed sampled light curves **************************
	  std::vector<LightCurve*> samp_LC(images.size());
	  for(int q=0;q<images.size();q++){
	    samp_LC[q] = new LightCurve(tobs);
	  }
	  
	  // Calculate the combined light curve for each image
	  for(int q=0;q<images.size();q++){
	    double macro_mag = abs(images[q]["mag"].asDouble());
	    LightCurve* samp_LC_intrinsic = new LightCurve(tobs);
	    LC_intrinsic[lc_in]->interpolate(samp_LC_intrinsic,td_max - images[q]["dt"].asDouble());
	    
	    if( unmicro ){
	      // === Combining three signals: intrinsic, intrinsic unmicrolensed, and extrinsic
	      LightCurve* samp_LC_unmicro = new LightCurve(tobs);
	      LC_unmicro[lc_in]->interpolate(samp_LC_unmicro,td_max - images[q]["dt"].asDouble());

	      if( LC_extrinsic[q][lc_ex]->time.size() > 0 ){ // Check if multiple image does not have a corresponding extrinsic light curve (i.e. a maximum image without a magnification map)
		LC_extrinsic[q][lc_ex]->interpolate(samp_LC[q],0.0);
		for(int t=0;t<tobs.size();t++){
		  samp_LC[q]->signal[t] = macro_mag*(samp_LC[q]->signal[t]*samp_LC_intrinsic->signal[t] + samp_LC_unmicro->signal[t]);
		}
	      } else {
		for(int t=0;t<tobs.size();t++){
		  samp_LC[q]->signal[t] = macro_mag*(samp_LC_intrinsic->signal[t] + samp_LC_unmicro->signal[t]); // this line includes only the intrinsic signal and excludes microlensing
		}
	      }

	      delete(samp_LC_unmicro);
	    } else {
	      // === Combining two signals: intrinsic and extrinsic
	      if( LC_extrinsic[q][lc_ex]->time.size() > 0 ){ // Check if multiple image does not have a corresponding extrinsic light curve (i.e. a maximum image without a magnification map)
		LC_extrinsic[q][lc_ex]->interpolate(samp_LC[q],0.0);
		for(int t=0;t<tobs.size();t++){
		  samp_LC[q]->signal[t] = samp_LC[q]->signal[t] * macro_mag * samp_LC_intrinsic->signal[t];
		}
	      } else {
		for(int t=0;t<tobs.size();t++){
		  samp_LC[q]->signal[t] = macro_mag * samp_LC_intrinsic->signal[t]; // this line includes only the intrinsic signal and excludes microlensing
		}
	      }    
	    }

	    delete(samp_LC_intrinsic);
	  }
	  
	  // Write json light curves
	  outputLightCurvesJson(samp_LC,out_path+mock+"/"+instrument_name+"_LC_sampled.json");
	  // *********************** End of product **************************************************


	  
	  
	  // *********************** Product: Observed sampled cut-outs (images) *****************************
	  if( root["point_source"]["output_cutouts"].asBool() ){
	    for(int t=0;t<tobs.size();t++){

	      // Loop over the truncated PSF (through PSF_offsets) for each image, and add their light to the pp_light image that contains all the point source light.
	      RectGrid pp_light(super_res_x,super_res_y,xmin,xmax,ymin,ymax); // this has to be in intensity units in order to be able to add the different light components
	      for(int q=0;q<images.size();q++){
		for(int i=0;i<PSFoffsets[q].ni;i++){
		  for(int j=0;j<PSFoffsets[q].nj;j++){
		    int index_img = PSFoffsets[q].offset_image + i*pp_light.Nx + j;
		    int index_psf = PSFoffsets[q].offset_cropped + i*Instrument_list[q]->cropped_psf->Nx + j;
		    //pp_light.z[index_img] += 1.0;
		    pp_light.z[index_img] += samp_LC[q]->signal[t]*Instrument_list[q]->cropped_psf->z[index_psf]/psf_partial_sum[q];
		  }
		}
	      }

	      /*
	      char buffer[4];
	      sprintf(buffer,"%03d",t);
	      std::string timestep = buffer;
	      FitsInterface::writeFits(pp_light.Nx,pp_light.Ny,pp_light.z,out_path+mock+"/OBS_"+instrument_name+"_"+timestep+".fits");
	      */
	      
	      // Check the expected brightness of the multiple images vs the image

		//double sum = 0.0;
		//for(int i=0;i<pp_light.Nm;i++){
		//sum += pp_light.img[i];
		//}
		//double fac = inf_dx*inf_dy;
		//double true_sum = 0.0;
		//for(int q=0;q<images.size();q++){
		//true_sum += img_signal[q][t];
		//}
		//printf("True: %15.10f (%15.10f)  Numerical: %15.10f (%15.10f)\n",true_sum,-2.5*log10(true_sum),sum,-2.5*log10(sum));

	      // Output Point Source image only
	      //char baf[4];
	      //sprintf(baf,"%03d",t);
	      //std::string timesteep = baf;
	      //pp_light.writeImage(out_path+mock+"/PS_"+instrument_name+"_"+timesteep+".fits");


	      // Bin image from 'super' to observed resolution
	      RectGrid* obs_img = pp_light.embeddedNewGrid(res_x,res_y,"additive");
	      
	      // Adding time-dependent noise here
	      mycam.noise->addNoise(obs_img);
	      
	      // Finalize output (e.g convert to magnitudes) and write
	      if( cut_out_scale == "mag" ){
		for(int i=0;i<obs_img->Nz;i++){
		  obs_img->z[i] = -2.5*log10(obs_img->z[i]);
		}
	      }
	      char buffer[4];
	      sprintf(buffer,"%03d",t);
	      std::string timestep = buffer;
	      FitsInterface::writeFits(obs_img->Nx,obs_img->Ny,obs_img->z,out_path+mock+"/OBS_"+instrument_name+"_"+timestep+".fits");
	      delete(obs_img);

	    }
	  }
	  // *********************** End of product **************************************************	    

	  // Do some cleanup
	  for(int q=0;q<images.size();q++){
	    delete(samp_LC[q]);
	  }

	  //std::cout << "done" << std::endl;
	}
	// Loop over extrinsic light curves ends here
	//0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=
      }
      // Loop over intrinsic light curves ends here
      //0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=

      
      for(int lc_in=0;lc_in<N_in;lc_in++){
	delete(LC_intrinsic[lc_in]);
      }

      if( unmicro ){
	for(int lc_in=0;lc_in<N_in;lc_in++){
	  delete(LC_unmicro[lc_in]);
	}
      }

      for(int q=0;q<images.size();q++){
	for(int lc_ex=0;lc_ex<N_ex;lc_ex++){
	  delete(LC_extrinsic[q][lc_ex]);
	}
      }

    }
    //================= END:CREATE THE TIME VARYING LIGHT ====================
    
    
    
  }
  // Loop over the instruments ends here
  // ===================================================================================================================
  // ===================================================================================================================

  
  return 0;
}
