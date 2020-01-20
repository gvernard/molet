#ifndef NON_LINEAR_PARS_HPP
#define NON_LINEAR_PARS_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>

#include "json/json.h"

// In this file: BasePrior, Uni, Gauss, Exp, and FactoryPrior
//               Nlpar

class BasePrior;

class Nlpar {
public:
  std::string nam;
  int fix;
  int per;
  double val;
  double err;
  double min;
  double max;
  double ran;
  BasePrior* pri = NULL;


  Nlpar(){};
  Nlpar(std::string a,int b,int c,double d,double e,double f,double g);
  ~Nlpar(){
    free(pri);
  };

  std::string getName();
  double getValue();
  int getActive();
  void setNewPrior(BasePrior* prior);
  void printNlpar();

  static std::vector<Nlpar*> nlparsFromJsonVector(const Json::Value myjson);
  static std::vector<std::string> getVectorNames(std::vector<Nlpar*> pars);
  static std::vector<double> getVectorValues(std::vector<Nlpar*> pars);
  static std::vector<int> getVectorActive(std::vector<Nlpar*> pars);
  static double getValueByName(std::string nam,std::vector<Nlpar*> pars);
  static Nlpar* getParByName(std::string nam,std::vector<Nlpar*> pars);
  static bool getSampleReg(std::vector<Nlpar*> pars);
  static std::string removeSuffix(std::vector<Nlpar*> pars);
  static void addSuffix(std::vector<Nlpar*> pars,std::string suffix);

  static std::map<std::string,double> getSigmaIntervals(const std::vector<double>& x,const std::vector<double>& p,int sigma_interval);
  static bool compareFunc(const std::pair<double,double> &a,const std::pair<double,double> &b);
  static double calculateInterval(const std::vector<double>& x,const std::vector<double>& p,double limit);
};


class BasePrior {
public:
  std::string type;
  Nlpar* mother;

  BasePrior(){};
  ~BasePrior(){};

  virtual double prior(double x) = 0;
  virtual double fromUnitCube(double u) = 0;
  virtual void printPars() = 0;
  virtual std::map<std::string,double> getPars() = 0;
};


class Uni: public BasePrior {
public:
  Uni(Nlpar* p);
  ~Uni(){};
  
  double prior(double x);
  double fromUnitCube(double u);
  void printPars();
  std::map<std::string,double> getPars();

private:
  double p0; // the constant probability
};


class Gauss: public BasePrior {
public:
  Gauss(Nlpar* p,double a,double b);
  ~Gauss(){};
  
  double prior(double x);
  double fromUnitCube(double u);
  void printPars();
  std::map<std::string,double> getPars();

private:
  double mean;
  double sdev;
  double two_sdev2; // 2 times the sdev squared
  double fac;       // 1/(s*sqrt(2*pi)) the normalization factor of a normal distribution
  double den;       // the denominator of a truncated normal distribution (difference between two cumulative distributions)

  double F(double x);
};


class Exp: public BasePrior {
public:
  Exp(Nlpar* p,double a);
  ~Exp(){};

  double prior(double x);
  double fromUnitCube(double u);
  void printPars();
  std::map<std::string,double> getPars();

private:
  double beta;
  double fac;
};


class Log: public BasePrior {
public:
  Log(Nlpar* p);
  ~Log(){};

  double prior(double x);
  double fromUnitCube(double u);
  void printPars();
  std::map<std::string,double> getPars();

private:
  double fac;
};

class PowerLaw: public BasePrior {
public:
  PowerLaw(Nlpar* p,int beta);
  ~PowerLaw(){};

  double prior(double x);
  double fromUnitCube(double u);
  void printPars();
  std::map<std::string,double> getPars();

private:
  int beta;
  double fac;
};

//This is a singleton class.
class FactoryPrior {
public:
  FactoryPrior(FactoryPrior const&) = delete;//Stop the compiler generating methods of copying the object.
  void operator=(FactoryPrior const&) = delete;

  static FactoryPrior* getInstance(){
    static FactoryPrior dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BasePrior* createPrior(Nlpar* mother,std::map<std::string,std::string> prior){
    if( prior["type"] == "uni" ){ 
      return new Uni(mother);
    } else if( prior["type"] == "gauss" ){
      return new Gauss(mother,stof(prior["mean"]),stof(prior["sdev"]));
    } else if( prior["type"] == "exp" ){
      return new Exp(mother,stof(prior["beta"]));
    } else if( prior["type"] == "log10" ){
      return new Log(mother);
    } else if( prior["type"] == "plaw" ){
      return new PowerLaw(mother,stoi(prior["beta"]));
    } else {
      return NULL;
    }
  }
  
private:
  FactoryPrior(){};
};

#endif /* NON_LINEAR_PARS_HPP */
