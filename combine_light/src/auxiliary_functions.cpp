#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <dirent.h>
#include <algorithm>

#include "auxiliary_functions.hpp"

#include "vkllib.hpp"
#include "instruments.hpp"
#include "noise.hpp"

// START:LIGHTCURVE ========================================================================================
LightCurve::LightCurve(const Json::Value lc){
  for(int t=0;t<lc["time"].size();t++){
    this->time.push_back( lc["time"][t].asDouble() );
    this->signal.push_back( lc["signal"][t].asDouble() );
    this->dsignal.push_back( lc["dsignal"][t].asDouble() );
  }
}

LightCurve::LightCurve(std::vector<double> t){
  this->time = t;
  this->signal.resize(t.size());
  this->dsignal.resize(t.size());
}

LightCurve::LightCurve(std::vector<double> t,std::vector<double> s,std::vector<double> ds){
  this->time   = t;
  this->signal = s;
  this->signal = ds;
}

// I think the following function is not used somewhere anymore
void LightCurve::interpolate(std::vector<double> obs_time,double delay,double* interpolated){
  // Assumption 1: obs_time and interpolated have the same length
  // Assumption 2: the light curve is at least as long as obs_time (in actual time values)
  int i = 0;
  int j  = 1;
  while( i < obs_time.size() ){
    double t = obs_time[i] + delay;
    if( this->time[j] > t ){
      //interpolate
      interpolated[i] = this->signal[j-1]+(t-this->time[j-1])*(this->signal[j]-this->signal[j-1])/(this->time[j]-this->time[j-1]);
      i++;
    } else {
      j++;
    }
  }
}

void LightCurve::interpolate(LightCurve* interpolated,double delay){
  // Assumption 1: obs_time and interpolated have the same length
  // Assumption 2: the light curve is at least as long as obs_time (in actual time values)
  int i = 0;
  int j = 1;
  while( i < interpolated->time.size() ){
    double t = interpolated->time[i] + delay;
    if( this->time[j] > t ){
      //interpolate
      interpolated->signal[i]  = this->signal[j-1]+(t-this->time[j-1])*(this->signal[j]-this->signal[j-1])/(this->time[j]-this->time[j-1]);
      interpolated->dsignal[i] = this->dsignal[j-1]+(t-this->time[j-1])*(this->dsignal[j]-this->dsignal[j-1])/(this->time[j]-this->time[j-1]);
      i++;
    } else {
      j++;
    }
  }
}

Json::Value LightCurve::jsonOut(){
  Json::Value json_lc,json_time,json_signal,json_dsignal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(this->signal[i]);
    json_dsignal.append(this->dsignal[i]);
  }
  json_lc["time"]    = json_time;
  json_lc["signal"]  = json_signal;
  json_lc["dsignal"] = json_dsignal;
  return json_lc;
}

Json::Value LightCurve::jsonOutMag(double ZP){
  Json::Value json_lc,json_time,json_signal,json_dsignal;
  for(int i=0;i<this->time.size();i++){
    json_time.append(this->time[i]);
    json_signal.append(-2.5*log10(this->signal[i]) + ZP);
    // This is a symmetric error derived from the asymmetric logarithmic errors as: dm = [(m+)+(m-)]/2
    json_dsignal.append( -2.5*log10(sqrt(1.0-pow(this->dsignal[i]/this->signal[i],2))) + ZP );
  }
  json_lc["time"]    = json_time;
  json_lc["signal"]  = json_signal;
  json_lc["dsignal"] = json_dsignal;
  return json_lc;
}
// END:LIGHTCURVE ========================================================================================






// START:LIGHTCURVE MANIPULATION FUNCTIONS ========================================================================================
Json::Value readLightCurvesJson(std::string lc_type,std::string type,std::string instrument_name,std::string in_path,std::string out_path){
  std::ifstream fin;
  Json::Value lcs_json;
  if( type == "custom" ){
    fin.open(in_path+"/input_files/"+instrument_name+"_LC_"+lc_type+".json",std::ifstream::in);
  } else {
    fin.open(out_path+"output/"+instrument_name+"_LC_"+lc_type+".json",std::ifstream::in);
  }
  fin >> lcs_json;
  fin.close();
  return lcs_json;
}

std::vector<LightCurve*> conversions(Json::Value lcs_json,double zs,double scale,double ZP){
  int N = lcs_json.size();
  //  double fac = 1.0 + zs;
  std::vector<LightCurve*> lcs(N);
  for(int i=0;i<N;i++){
    lcs[i] = new LightCurve(lcs_json[i]);    
    for(int j=0;j<lcs[i]->signal.size();j++){
      lcs[i]->signal[j] = scale*pow(10.0,-0.4*(lcs[i]->signal[j]-ZP)); // Convert from magnitudes to intensities and scale by a factor if necessary
      //lcs[i]->signal[j] = pow(10.0,-0.4*(lcs[i]->signal[j]+scale)); // Add constant to magnitude (Dmag) and convert from magnitudes to intensities
      //lcs[i]->time[j] *= fac; // Convert time to the observer's frame
    }
  }	
  return lcs;
}


std::vector<LightCurve*> conversions_SN(Json::Value lcs_json,double zs,double scale,double ZP,double tobs_min,double td_max){
  int N = lcs_json.size();
  double fac = 1.0 + zs;
  std::vector<LightCurve*> lcs(N);
  for(int i=0;i<N;i++){

    // Check starting time of intrinsic light curve
    if( lcs_json[i]["time"][0] > (tobs_min - td_max) ){
      // Extend the light curves backwards in time to reach at least "tobs_min - td_max" - the start of the observations minus the maximum time delay.

      // Extend by two entries and make 10 mag dimmer
      double t00 = tobs_min - td_max - 1;
      double t01 = lcs_json[i]["time"][0].asDouble() - 0.1;
      double sig = lcs_json[i]["signal"][0].asDouble() + 10;

      Json::Value dum; 
      dum["time"] = Json::Value(Json::arrayValue);
      dum["signal"] = Json::Value(Json::arrayValue);
      dum["dsignal"] = Json::Value(Json::arrayValue);
      dum["time"].append(t00);
      dum["time"].append(t01);
      dum["signal"].append(sig);
      dum["signal"].append(sig);
      dum["dsignal"].append(0);
      dum["dsignal"].append(0);
      for(int j=0;j<lcs_json[i]["time"].size();j++){
	dum["time"].append(lcs_json[i]["time"][j]);
	dum["signal"].append(lcs_json[i]["signal"][j]);
	dum["dsignal"].append(lcs_json[i]["dsignal"][j]);
      }
      lcs_json[i] = dum;
    }

    lcs[i] = new LightCurve(lcs_json[i]);    
    for(int j=0;j<lcs[i]->signal.size();j++){
      lcs[i]->signal[j] = scale*pow(10.0,-0.4*(lcs[i]->signal[j]-ZP)); // Convert from magnitudes to intensities and scale by a factor if necessary
      //lcs[i]->signal[j] = pow(10.0,-0.4*(lcs[i]->signal[j]+scale)); // Add constant to magnitude (Dmag) and convert from magnitudes to intensities
      //lcs[i]->time[j] *= fac; // Convert time to the observer's frame
    }
  }	
  return lcs;
}

void justSupernovaInSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* target){
  // expand intrinsic light curve
  std::vector<double> tmp_time(4+LC_intrinsic->time.size());
  LightCurve* tmp_in = new LightCurve(tmp_time);
  double floor_val_in = LC_intrinsic->signal[0]/10.0; // light curves are in flux, not in magnitudes
  double time_shift = td; // the shift is measured from time_max backwards
  tailorLightCurve(LC_intrinsic,tmp_in,time_shift,time[0],time.back(),floor_val_in);
 
  LightCurve* base = new LightCurve(time);
  tmp_in->interpolate(base,0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*base->signal[t]; // this line includes only the intrinsic signal and excludes microlensing
    target->dsignal[t] = 0.0;
  }
  delete(base);
  delete(tmp_in);
}

void combineSupernovaInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target){
  double time_min = time[0]; 
  double time_max = time.back(); 
  double time_shift,floor_val;
  
  std::vector<double> tmp_time_in(4+LC_intrinsic->time.size());
  LightCurve* tmp_in = new LightCurve(tmp_time_in);
  floor_val = LC_intrinsic->signal[0]/10.0; // light curves are in flux, not in magnitudes
  time_shift = td; // the shift is measured from time_max backwards
  tailorLightCurve(LC_intrinsic,tmp_in,time_shift,time_min,time_max,floor_val);

  std::vector<double> tmp_time_ex(4+LC_extrinsic->time.size());
  LightCurve* tmp_ex = new LightCurve(tmp_time_ex);
  floor_val = 1.0; // no microlensing magnification 
  time_shift = LC_intrinsic->time[0] + td; // shifting the microlensing light curves to absolute time in order to match the intrinsic ones.
  tailorLightCurve(LC_extrinsic,tmp_ex,time_shift,time_min,time_max,floor_val);

  LightCurve* base = new LightCurve(time);
  tmp_in->interpolate(base,0);
  tmp_ex->interpolate(target,0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = base->signal[t]*macro_mag*target->signal[t];
    target->dsignal[t] = base->dsignal[t]*macro_mag*target->signal[t];
  }
  delete(base);
  delete(tmp_in);
  delete(tmp_ex);  
}

void tailorLightCurve(LightCurve* lc_input,LightCurve* lc_output,double td,double tstart,double tend,double floor_val){
  int Nout = lc_output->time.size();
  int Nin = lc_input->time.size();
  // First set the time
  lc_output->time[0] = tstart - 1;
  lc_output->time[1] = (lc_input->time[0] + td) - 0.1; // just 0.1 days earlier should be small enough
  lc_output->time[Nout-2] = (lc_input->time[Nin-1] + td) + 0.1;
  lc_output->time[Nout-1] = tend + 1;
  for(int i=0;i<Nin;i++){
    lc_output->time[2+i] = lc_input->time[i] + td;
  }
  // Then set the signal
  lc_output->signal[0] = floor_val;
  lc_output->signal[1] = floor_val;
  lc_output->signal[Nout-2] = floor_val;
  lc_output->signal[Nout-1] = floor_val;
  lc_output->dsignal[0] = 0.0;
  lc_output->dsignal[1] = 0.0;
  lc_output->dsignal[Nout-2] = 0.0;
  lc_output->dsignal[Nout-1] = 0.0;
  for(int i=0;i<Nin;i++){
    lc_output->signal[2+i]  = lc_input->signal[i];
    lc_output->dsignal[2+i] = lc_input->dsignal[i];
  }
}

void combineInExSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* target){
  // === Combining two signals: intrinsic and extrinsic
  LightCurve* base = new LightCurve(time); // base is the same as target, they are both initialized with tcont/time
  LC_intrinsic->interpolate(base,td);
  LC_extrinsic->interpolate(target,0.0);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = target->signal[t]*macro_mag*base->signal[t];
    target->dsignal[t] = target->dsignal[t]*macro_mag*base->signal[t]; // no error from macro-mag, intrinsic, and unmicrolensed flux, only from the microlensing mag
  }
  delete(base);
}

void combineInExUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_extrinsic,LightCurve* LC_unmicro,LightCurve* target,double unmicro_ratio){
  // === Combining three signals: intrinsic, intrinsic unmicrolensed, and extrinsic
  LightCurve* base = new LightCurve(time);
  LightCurve* base_unmicro = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  LC_unmicro->interpolate(base_unmicro,td);
  LC_extrinsic->interpolate(target,0.0);
  double factor = 1.0 - unmicro_ratio;
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*(factor*target->signal[t]*base->signal[t] + unmicro_ratio*base_unmicro->signal[t]);
    target->dsignal[t] = macro_mag*target->dsignal[t]*base->signal[t]; // no error from macro-mag, intrinsic, and unmicrolensed flux, only from the microlensing mag    
  }
  delete(base);
  delete(base_unmicro);
}

void combineInUnSignals(double td,double macro_mag,std::vector<double> time,LightCurve* LC_intrinsic,LightCurve* LC_unmicro,LightCurve* target,double unmicro_ratio){
  // === Combining two signals: intrinsic and unmicrolensed
  LightCurve* base = new LightCurve(time);
  LightCurve* base_unmicro = new LightCurve(time);
  LC_intrinsic->interpolate(base,td);
  LC_unmicro->interpolate(base_unmicro,0.0);
  double factor = 1.0 - unmicro_ratio;
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*(factor*base->signal[t] + unmicro_ratio*base_unmicro->signal[t]);
    target->dsignal[t] = 0.0;
  }
  delete(base);
  delete(base_unmicro);
}

void justOneSignal(double td,double macro_mag,std::vector<double> time,LightCurve* LC,LightCurve* target){
  // === Just one signal
  LightCurve* base = new LightCurve(time);
  LC->interpolate(base,td);
  for(int t=0;t<time.size();t++){
    target->signal[t]  = macro_mag*base->signal[t]; // this line includes only the intrinsic signal and excludes microlensing
    target->dsignal[t] = 0.0;
  }
  delete(base);
}

void outputLightCurvesJson(std::vector<LightCurve*> lcs,double ZP,std::string filename){
  Json::Value lcs_json;
  for(int q=0;q<lcs.size();q++){
    lcs_json.append( lcs[q]->jsonOutMag(ZP) );
    //lcs_json.append( lcs[q]->jsonOut() );
  }
  std::ofstream lcs_file(filename);
  lcs_file << lcs_json;
  lcs_file.close();									      
}
// END:LIGHTCURVE MANIPULATION FUNCTIONS ========================================================================================





// START:IMAGE MANIPULATION FUNCTIONS ========================================================================================
vkl::RectGrid createObsStatic(Instrument* mycam,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,int res_x,int res_y,double& F_conv_extended,double& F_conv_lens,bool convolve_lens){  
  // supersim is just carrying the resolution and extent of the image.
  int super_res_x = super_extended->Nx;
  int super_res_y = super_extended->Ny;
  double xmin = super_extended->xmin;
  double xmax = super_extended->xmax;
  double ymin = super_extended->ymin;
  double ymax = super_extended->ymax;
  
  // Get the fixed extended lensed light. This is in units of flux!
  vkl::RectGrid convolved_extended = *super_extended;
  mycam->convolve(super_extended,&convolved_extended);
  //writeCutout(&convolved_extended,out_path+"output/psf_lensed_image_super.fits");
  
  // Get the fixed lens galaxy light. This is in units of flux!
  vkl::RectGrid convolved_lens_light = *super_lens_light;
  if( convolve_lens ){
    mycam->convolve(super_lens_light,&convolved_lens_light);
    //writeCutout(&convolved_lens_light,out_path+"output/psf_lens_light_super.fits");  
  }
  
  // Calculate convolved fluxes
  F_conv_lens = convolved_lens_light.integrate();
  F_conv_extended = convolved_extended.integrate();
  
  // Combined light of the observed base image (integrating flux density from 'super' to observed resolution)
  vkl::RectGrid* base = new vkl::RectGrid(super_res_x,super_res_y,xmin,xmax,ymin,ymax); 
  for(int i=0;i<base->Nz;i++){
    base->z[i] = convolved_lens_light.z[i] + convolved_extended.z[i];
  }
  vkl::RectGrid obs_static = base->embeddedNewGrid(res_x,res_y,"average"); // This is in units of flux!
  delete(base);

  return obs_static;
}

vkl::RectGrid createObsPS(vkl::RectGrid* supersim,std::vector<double> image_signal,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y){
  // Loop over the truncated PSF (through PSF_offsets) for each image, and add their light to the ps_light image that contains all the point source light.
  // All vectors must have the same size, equal to the number of multiple images.

  // supersim is just carrying the resolution and extent of the image.
  vkl::RectGrid super_ps_light(supersim->Nx,supersim->Ny,supersim->xmin,supersim->xmax,supersim->ymin,supersim->ymax);
  
  for(int q=0;q<PSFoffsets.size();q++){
    for(int i=0;i<PSFoffsets[q].ni;i++){
      for(int j=0;j<PSFoffsets[q].nj;j++){
	int index_img = PSFoffsets[q].offset_image + i*super_ps_light.Nx + j;
	int index_psf = PSFoffsets[q].offset_cropped + i*instrument_list[q]->cropped_psf->Nx + j;
	//super_ps_light.z[index_img] += 1.0;
	super_ps_light.z[index_img] += image_signal[q]*instrument_list[q]->cropped_psf->z[index_psf]/psf_partial_sums[q];
      }
    }
  }

  // Bin image from 'super' to observed resolution
  vkl::RectGrid obs_ps_light = super_ps_light.embeddedNewGrid(res_x,res_y,"average"); // This is (and should be) still in flux units!!!
  
  return obs_ps_light;  
}

void writeCutout(vkl::RectGrid* obs,std::string fname){
  std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
  std::vector<std::string> values{std::to_string(obs->xmin),std::to_string(obs->xmax),std::to_string(obs->ymin),std::to_string(obs->ymax)};
  std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
  vkl::FitsInterface::writeFits(obs->Nx,obs->Ny,obs->z,keys,values,descriptions,fname);
}


void writeAllCutouts(std::vector<double> tobs,Json::Value images,std::vector<LightCurve*> samp_LC,vkl::RectGrid* super_extended,vkl::RectGrid* super_lens_light,Instrument* mycam,std::vector<offsetPSF>& PSFoffsets,std::vector<Instrument*>& instrument_list,std::vector<double>& psf_partial_sums,int res_x,int res_y,std::string mock,bool convolve_lens,std::string out_path){

  // Create list of PSF file names +++++++++++++++++++	    
  //std::vector<std::string> psf_fnames = getFileNames(tobs,argv[4]); // user-provided function to get the file names of the PSFs as a function of timestep t
  //++++++++++++++++++++++++++++++++++++++++++++++++++	    
  
  for(int t=0;t<tobs.size();t++){
    
    // Time-dependent PSF goes here ++++++++++++++++++++
    // //std::cout << psf_fnames[t] << std::endl;
    // mycam.replacePSF(psf_fnames[t]);
    // mycam.preparePSF(&supersim,0.999);
    // for(int q=0;q<instrument_list.size();q++){
    // 	PSFoffsets[q] = instrument_list[q]->offsetPSFtoPosition(images[q]["x"].asDouble(),images[q]["y"].asDouble(),&supersim);
    // }
    // for(int q=0;q<instrument_list.size();q++){
    // 	psf_partial_sums[q] = instrument_list[q]->sumPSF(&PSFoffsets[q]);
    // }
    // vkl::RectGrid obs_base_t = createObsBase(&mycam,&supersim,res_x,res_y,out_path);
    // ptr_obs_base = &obs_base_t;
    //++++++++++++++++++++++++++++++++++++++++++++++++++

    // The following line is only for a non-time-varying PSF
    double F_conv_extended,F_conv_lens;
    vkl::RectGrid obs_static = createObsStatic(mycam,super_extended,super_lens_light,res_x,res_y,F_conv_extended,F_conv_lens,convolve_lens); // This is in flux units!
    
    // construct a vector with the point source brightness in each image at the given time step
    std::vector<double> image_signal(images.size());
    for(int q=0;q<images.size();q++){
      image_signal[q] = samp_LC[q]->signal[t];
    }
    vkl::RectGrid obs_ps_light = createObsPS(super_extended,image_signal,PSFoffsets,instrument_list,psf_partial_sums,res_x,res_y);  // This is in flux units!

    // Calculate total flux of the PS light only
    //double F_conv_ps,F_conv_ps_mag;
    //obs_ps_light.integrate(F_conv_ps,F_conv_ps_mag,mycam->ZP);
    
    // The two grids, obs_ps_light and obs_static, have the same resolution, so I can simply add their fluxes
    for(int i=0;i<obs_ps_light.Nz;i++){
      obs_ps_light.z[i] += obs_static.z[i]; // <------ BE CAREFUL: 'obs_ps_light' now also contains the static light
    }
    
    // Adding a simple time-dependent noise here, no need to save the noise realizations
    convertGridFromFlux(&obs_ps_light,mycam->ZP);
    mycam->noise->initializeFromData(&obs_ps_light);
    mycam->noise->calculateNoise();
    mycam->noise->addNoise(&obs_ps_light);
    
    char buffer[100];
    sprintf(buffer,"%s%s/OBS_%s_%03d.fits",out_path.c_str(),mock.c_str(),mycam->name.c_str(),t);
    std::string fname(buffer);
    writeCutout(&obs_ps_light,fname);
  }
}

void convertGridFromFlux(vkl::RectGrid* obs,double ZP){
  for(int i=0;i<obs->Nz;i++){
    obs->z[i] = convertFromFlux(obs->z[i],ZP); // Convert to electrons/(s arcsec^2)
  }
}

double convertFromFlux(double flux,double ZP){
  if( flux == 0.0 ){
    return 0.0; // Flux is zero (cannot convert to electrons/(s arcsec^2))
  } else {
    return -2.5*log10(flux) + ZP; // Convert to electrons/(s arcsec^2)
  }
}
// END:IMAGE MANIPULATION FUNCTIONS ==========================================================================================




// START:TRANSFORM PSF =====================================================================================
TransformPSF::TransformPSF(double a,double b,double c,bool d,bool e){
  this->x0 = a;
  this->y0 = b;
  this->rot = c * 0.01745329251;// converting to rad
  this->flip_x = d;
  this->flip_y = e;

  this->cosrot = cos(this->rot);
  this->sinrot = sin(this->rot);
}

void TransformPSF::applyTransform(double xin,double yin,double& xout,double& yout){
  xout =  (xin-this->x0)*this->cosrot + (yin-this->y0)*this->sinrot;
  yout = -(xin-this->x0)*this->sinrot + (yin-this->y0)*this->cosrot;
  if( this->flip_x ){
    xout *= -1.0;
  }
  if( this->flip_y ){
    yout *= -1.0;
  }
}

double TransformPSF::interpolateValue(double x,double y,PSF* mypsf){
  return 1.0;
}
// END:TRANSFORM PSF =====================================================================================





// START:TIME LAG KERNELS =====================================================================================
void DeltaKernel::getKernel(std::vector<double> time,std::vector<double>& kernel){
  int index = 0;
  for(int i=0;i<time.size();i++){
    kernel[i] = 0.0;
    if( (time[i] - time[0]) < this->t_peak ){ // Need to take the difference with the first element because t_peak is an interval Delta t.
      index = i;
    }
  }
  kernel[index] = 1.0;
  std::cout << "Delta lag kernel: Delta location is: " << (time[index] - time[0]) << " (given was: " << this->t_peak << ") days" << std::endl;
}

TopHatKernel::TopHatKernel(double radius){
  this->t_limit = radius/this->speed_of_light;  // radius in [10^14 cm], speed of light in [10^14 cm / day]
}
void TopHatKernel::getKernel(std::vector<double> time,std::vector<double>& kernel){
  int counter = 0;
  for(int i=0;i<time.size();i++){
    if( (time[i] - time[0]) < this->t_limit ){ // Need to take the difference with the first element because t_peak is an interval Delta t.
      kernel[i] = 1.0;
      counter++;
    } else {
      kernel[i] = 0.0;
    }
  }

  for(int i=0;i<time.size();i++){
    kernel[i] /= counter;
  }
}
// END:TIME LAG KERNELS =====================================================================================





// user provided function for getting the file names of a time dependent PSF
std::vector<std::string> getFileNames(std::vector<double> tobs,std::string path){
  DIR* dir;
  struct dirent* diread;
  std::vector<std::string> files;
  if( (dir=opendir(path.c_str())) != nullptr ){
    while( (diread=readdir(dir)) != nullptr ){
      files.push_back(diread->d_name);
    }
    closedir(dir);
  } else {
    // directory not found
  }
  std::sort(files.begin(),files.end(),std::greater<std::string>()); // alphanumerically sorted in descending order
  files.pop_back(); // removing '.' and '..'
  files.pop_back();
  
  // for(auto file:files){
  //   std::cout << file << std::endl;
  // }
  
  std::vector<std::string> names(tobs.size());
  for(int t=0;t<tobs.size();t++){
    int index = t % files.size();
    names[t] = path+files[index];
    //names[t] = "/home/george/Desktop/My Papers/unpublished/2021_molet/extract_with_cosmouline/psfs/3_ECAM.2019-03-01T04:54:57.000_psf_abcdef.fits"; // reference image
  }
  return names;
}

