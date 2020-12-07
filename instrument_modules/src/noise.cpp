#include <algorithm>
#include <stdlib.h>

#include "vkllib.hpp"

#include "noise.hpp"

// START: NoNoise ====================================
NoNoise::NoNoise(){}
void NoNoise::addNoise(RectGrid* mydata){}
// END: NoNoise ======================================

// START: UniformGaussian ============================
UniformGaussian::UniformGaussian(double sn){
  this->sn = sn;
}
void UniformGaussian::addNoise(RectGrid* mydata){
  this->seed += 2; // increment seed at each call

  double maxdata = *std::max_element(mydata->z,mydata->z+mydata->Nz);
  double sigma = maxdata/this->sn;
  srand48(seed);

  double min_noise = sigma; // just a starting value
  double z1,z2,u1,u2,noise;
  //Applying the Box-Muller transformation
  for(int i=0;i<mydata->Nz;i++){
    u1 = drand48();
    u2 = drand48();
    z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);
    //    z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);    

    noise = z1*sigma;
    mydata->z[i] += noise;
    if( noise < min_noise ){
      min_noise = noise;
    }
  }

  // renormalize by adding the minimum (negative) noise value
  for(int i=0;i<mydata->Nz;i++){
    mydata->z[i] += abs(min_noise);
  }
}
// END: UniformGaussian ==============================
