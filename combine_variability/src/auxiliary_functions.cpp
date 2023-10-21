#include "auxiliary_functions.hpp"



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



