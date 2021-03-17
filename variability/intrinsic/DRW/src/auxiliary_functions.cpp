#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "json/json.h"

#include "auxiliary_functions.hpp"

std::vector< std::vector<double> > getDRWLightCurve(std::vector<double> time,double zs,double Mi,double lrest,double mean_mag,int N_in,int seed){
  // Coeeficients to get the SF and the characteistic time scale from eq. 7 of MacLeod et al. 2010 (from their Tab. 1).
  std::vector<double> sf_coeffs{-0.51,-0.479,0.131,0.18,0};
  std::vector<double> tau_coeffs{2.4,0.17,0.03,0.21,0};
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uni(-1.0,1.0);
  std::vector< std::vector<double> > all_signals(N_in);

  // Here we use only the mean of the balck hole mass distribution given in eq. 11 of MacLeod et al. 2010.
  double logMBH = 2.0 - 0.27*Mi; // in M_solar
  double log_sf  = sf_coeffs[0]  + sf_coeffs[1]*log10(lrest/400.0)  + sf_coeffs[2]*(Mi+23)  + sf_coeffs[3]*(logMBH-9.0)  + sf_coeffs[4]*log10(1.0+zs);
  double log_tau = tau_coeffs[0] + tau_coeffs[1]*log10(lrest/400.0) + tau_coeffs[2]*(Mi+23) + tau_coeffs[3]*(logMBH-9.0) + tau_coeffs[4]*log10(1.0+zs);
  double sf  = pow(10.0,log_sf);
  double tau = pow(10.0,log_tau);

  //  sf  = 0.1;
  //  tau = 20;

  std::cout << "sf: " << sf << std::endl;
  std::cout << "tau: " << tau << std::endl;

  
  // Calculate the light curve
  double s2,u,v,z,value;
  for(int j=0;j<N_in;j++){
    std::vector<double> signal(time.size());
    for(int t=0;t<time.size();t++){
      
      double normal_mean = exp(-time[t]/tau)*mean_mag + mean_mag*(1.0-exp(-time[t]/tau));
      double normal_var  = 0.5*pow(sf,2)*(1.0-exp(-2*time[t]/tau));
      // Draw from Box-Muller transform
      s2 = 1.1;
      while(s2 > 1.0){
	u = uni(rng);
	v = uni(rng);
	s2 = u*u + v*v;
      }
      z = v*sqrt(-2.0*log(s2)/s2);	
      value = sqrt(normal_var)*z + normal_mean;
      signal[t] = value;
    }
    all_signals[j] = (signal);
  }

  return all_signals;
}


void writeLightCurves(std::vector<double> time,std::vector< std::vector<double> > all_lc_signals,std::string outfile){  
  Json::Value lcs;
  for(int j=0;j<all_lc_signals.size();j++){
    Json::Value lc;
    Json::Value j_time;
    Json::Value j_signal;
    for(int t=0;t<time.size();t++){
      j_time.append(time[t]);
      j_signal.append(all_lc_signals[j][t]);
    }
    lc["time"] = j_time;
    lc["signal"] = j_signal;
    lcs.append(lc);
  }

  std::ofstream file_lc(outfile);
  file_lc << lcs;
  file_lc.close();
}
