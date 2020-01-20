#ifndef COVARIANCE_KERNELS_HPP
#define COVARIANCE_KERNELS_HPP

#include <string>
#include <vector>
#include <iostream>

class Nlpar;

class BaseCovKernel {
public:
  std::string type;
  double cmax;
  
  BaseCovKernel(){};
  BaseCovKernel(const BaseCovKernel& other){
    this->type = other.type;
    this->cmax = other.cmax;
  };
  ~BaseCovKernel(){};
  
  virtual BaseCovKernel* clone() = 0;
  virtual double getCovariance(double r) = 0;
  virtual double getCovarianceSelf() = 0;
  virtual void setParameters(std::vector<Nlpar*> pars) = 0;
  virtual void printParameters() = 0;
};


class GaussKernel: public BaseCovKernel {
public:
  double sdev;

  GaussKernel(std::vector<Nlpar*> pars);
  GaussKernel(const GaussKernel& other);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
  void printParameters();
  GaussKernel* clone(){
    return new GaussKernel(*this);
  };
private:
 double fac;
};

class ModGaussKernel: public BaseCovKernel {
public:
  double sdev;

  ModGaussKernel(std::vector<Nlpar*> pars);
  ModGaussKernel(const ModGaussKernel& other);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
  void printParameters();
  ModGaussKernel* clone(){
    return new ModGaussKernel(*this);
  };
};

class ExpGaussKernel: public BaseCovKernel {
public:
  double expo;
  double sdev;

  ExpGaussKernel(std::vector<Nlpar*> pars);
  ExpGaussKernel(const ExpGaussKernel& other);
  double getCovariance(double r);
  double getCovarianceSelf();
  void setParameters(std::vector<Nlpar*> pars);
  void printParameters();
  ExpGaussKernel* clone(){
    return new ExpGaussKernel(*this);
  };
private:
  double fac;
};

class FactoryCovKernel {//This is a singleton class.
public:
  FactoryCovKernel(FactoryCovKernel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryCovKernel const&) = delete;

  static FactoryCovKernel* getInstance(){
    static FactoryCovKernel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseCovKernel* createCovKernel(const std::string kernel_type,std::vector<Nlpar*> pars){
    if( kernel_type == "modgauss" ){
      return new ModGaussKernel(pars);
    } else if( kernel_type == "gauss" ){
      return new GaussKernel(pars);
    } else if( kernel_type == "expgauss" ){
      return new ExpGaussKernel(pars);
    } else {
      return NULL;
    }
  }

private:
  FactoryCovKernel(){};
};


#endif /* COVARIANCE_KERNELS_HPP */
