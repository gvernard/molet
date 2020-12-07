#ifndef NOISE_HPP
#define NOISE_HPP

#define _USE_MATH_DEFINES

#include <string>

#include "json/json.h"

class RectGrid;

class BaseNoise {
public:
  int seed = 123;
  BaseNoise(){};
  ~BaseNoise(){};
  virtual void addNoise(RectGrid* mydata) = 0;
};

class NoNoise: public BaseNoise {
public:
  NoNoise();
  void addNoise(RectGrid* mydata);
};

class UniformGaussian: public BaseNoise {
public:
  const double two_pi = 2.0*M_PI;
  double sn; // signal to noise ratio
  UniformGaussian(double sn);
  void addNoise(RectGrid* mydata);
};

class FactoryNoiseModel{//This is a singleton class.
public:
  FactoryNoiseModel(FactoryNoiseModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryNoiseModel const&) = delete;

  static FactoryNoiseModel* getInstance(){
    static FactoryNoiseModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseNoise* createNoiseModel(Json::Value noise_pars){
    std::string type = noise_pars["type"].asString();
    if( type == "NoNoise" ){
      return new NoNoise();
    } else if( type == "UniformGaussian" ){
      double sn = noise_pars["sn"].asDouble();
      return new UniformGaussian(sn);
    } else {
      return new NoNoise();
    }
  }

private:
  FactoryNoiseModel(){};
};


#endif /* NOISE_HPP */
