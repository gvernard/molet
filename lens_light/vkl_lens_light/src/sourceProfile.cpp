#include "sourceProfile.hpp"

#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>

#include <CCfits/CCfits>

//Abstract class: BaseProfile
//===============================================================================================================
void BaseProfile::profile(int Sj,int Si,double* sx,double* sy,double* s){
  for(int i=0;i<Si;i++){
    for(int j=0;j<Sj;j++){
      s[i*Sj+j] += this->value(sx[i*Sj+j],sy[i*Sj+j]);
    }
  }
}

void BaseProfile::writeProfile(std::string filename,double half_range){
  // produces a square image centered at (0,0) on the source plane

  // create grid of source brightness profile
  std::valarray<double> array(output_res*output_res);
  double dpix = 2.0*half_range/output_res;
  for(int i=0;i<output_res;i++){
    for(int j=0;j<output_res;j++){
      double x = j*dpix - half_range;
      double y = i*dpix - half_range;
      //array[(output_res-1-i)*output_res+j] = this->value(x,y);
      array[i*output_res+j] = this->value(x,y);
    }
  }

  // Total brightness
  double sum = 0.0;
  for(long i=0;i<output_res*output_res;i++){
    sum += array[i]*pow(2*half_range/output_res,2);
  }
  //  std::cout << sum << std::endl;
  
  //Write FITS:
  long naxis    = 2;
  long naxes[2] = {(long) output_res,(long) output_res};
  long Ntot = (long) output_res*output_res;

  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  pFits.reset( new CCfits::FITS("!"+filename,DOUBLE_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) output_res);
  std::string newName("NEW-EXTENSION");
  CCfits::ExtHDU* imageExt = pFits->addImage(newName,DOUBLE_IMG,extAx);

  long fpixel(1);
  imageExt->write(fpixel,(long) Ntot,array);
  pFits->pHDU().addKey("WIDTH",2.0*half_range,"width of the image"); 
  pFits->pHDU().addKey("HEIGHT",2.0*half_range,"height of the image"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
}





//Derived class from BaseAnalyticFunction: Sersic
//===============================================================================================================
Sersic::Sersic(std::map<std::string,double> pars){
  this->type          = "sersic";
  this->pars["n"]     = pars["n"];
  this->pars["r_eff"] = pars["r_eff"];
  this->pars["i_eff"] = pars["i_eff"];
  this->pars["q"]     = pars["q"];
  this->pars["x0"]    = pars["x0"];
  this->pars["y0"]    = pars["y0"];
  this->pars["pa"]    = pars["pa"];
}

double Sersic::function_value(double x,double y){
  double bn = 1.9992*this->pars["n"] - 0.3271;//From Capaccioli 1989
  double u,v,r,fac2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  r = sqrt(this->pars["q"]*this->pars["q"]*u*u + v*v);
  fac2 = pow(r/this->pars["r_eff"],1./this->pars["n"]);
  return this->pars["i_eff"]*exp(-bn*fac2 - 1);
}

std::vector<double> Sersic::extent(){
  double dx = 3*this->pars["_reff"]*cos(this->pars["pa"]*this->fac);
  double xmin = this->pars["x0"] - dx;
  double xmax = this->pars["x0"] + dx;
  double dy = 3*this->pars["r_eff"]*sin(this->pars["pa"]*this->fac);
  double ymin = this->pars["y0"] - dy;
  double ymax = this->pars["y0"] + dy;
  std::vector<double> ranges = {xmin,xmax,ymin,ymax};
  return ranges;
}

//Derived class from BaseAnalyticFunction: Gauss
//===============================================================================================================
proGauss::proGauss(std::map<std::string,double> pars){
  this->type          = "gauss";
  this->pars["r_eff"] = pars["r_eff"];
  this->pars["i_eff"] = pars["i_eff"];
  this->pars["q"]     = pars["q"];
  this->pars["x0"]    = pars["x0"];
  this->pars["y0"]    = pars["y0"];
  this->pars["pa"]    = pars["pa"];
}

double proGauss::function_value(double x,double y){
  double u,v,r2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);  
  double sdev   = 2*this->pars["r_eff"]*this->pars["r_eff"];
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  //    u =   x*cosphi + y*sinphi;
  //    v = - x*sinphi + y*cosphi;
  r2 = (this->pars["q"]*this->pars["q"]*u*u + v*v)/sdev;
  //    return (this->ieff*exp(-r2)/(sqrt(sdev*3.14159)));
  return this->pars["i_eff"]*exp(-r2);
}

std::vector<double> proGauss::extent(){
  double dimg = 3.0*this->pars["r_eff"];
  double xmin = this->pars["x0"] - dimg;
  double xmax = this->pars["x0"] + dimg;
  double ymin = this->pars["y0"] - dimg;
  double ymax = this->pars["y0"] + dimg;
  std::vector<double> ranges = {xmin,xmax,ymin,ymax};
  return ranges;
}







//Derived class from BaseProfile: Analytic
//===============================================================================================================
Analytic::Analytic(std::vector<std::string> names,std::vector<std::map<std::string,double> > par_maps){
  this->type = "analytic";
  this->output_res = 500;
  for(int i=0;i<names.size();i++){
    BaseAnalyticFunction* function = FactoryAnalyticFunction::getInstance()->createAnalyticFunction(names[i],par_maps[i]);
    this->components.push_back( function );
  }
}

double Analytic::value(double x,double y){
  double value = 0.0;
  for(int i=0;i<this->components.size();i++){
    value += this->components[i]->function_value(x,y);
  }
  return value;
}

void Analytic::outputProfile(std::string filename){
  double half_range = 0.0;
  for(int i=0;i<this->components.size();i++){
    std::vector<double> ranges = this->components[i]->extent();
    for(int j=0;j<ranges.size();j++){
      if( fabs(ranges[j]) > half_range ){
	half_range = fabs(ranges[j]);
      }
    }
  }
  this->writeProfile(filename,half_range);
}







//Derived class from BaseProfile: fromFITS
//===============================================================================================================
fromFITS::fromFITS(std::string filename,int Ni,int Nj,double height,double width,double x0,double y0){
  this->type = "fromfits";
  this->Ni = Ni;
  this->Nj = Nj;
  this->height = height;
  this->width  = width;
  this->x0 = x0;
  this->y0 = y0;
  this->mySource = new ImagePlane(filename,Ni,Nj,height,width);
  scaleProfile();
  // ImagePlane sets the coordinate origin in the center of the image, so I can re-position it here.
  for(int i=0;i<this->mySource->Nm;i++){
    this->mySource->x[i] += x0;
    this->mySource->y[i] += y0;
  }
  this->mySource->xmin += x0;
  this->mySource->xmax += x0;
  this->mySource->ymin += y0;
  this->mySource->ymax += y0;
}

double fromFITS::value(double x,double y){
  if( this->mySource->xmin < x && x < this->mySource->xmax && this->mySource->ymin < y && y < this->mySource->ymax ){
    // Source and Image grids MUST be the same, therefore I just need to match the right pixels (no interpolation)
    int i = (int) floor((this->mySource->ymax - y)*this->mySource->Ni/this->mySource->height); // y-axis is reflected
    int j = (int) floor((x - this->mySource->xmin)*this->mySource->Nj/this->mySource->width);
    //int j = (int) floor((this->mySource->xmin - x)*this->mySource->Nj/this->mySource->width);
    return this->mySource->img[i*this->mySource->Nj + j];
  } else {
    return 0;
  }
}

void fromFITS::outputProfile(std::string filename){
  double xmin = this->x0 - this->width/2.0;
  double xmax = this->x0 + this->width/2.0;
  double ymin = this->y0 - this->height/2.0;
  double ymax = this->y0 + this->height/2.0;
  std::vector<double> ranges = {xmin,xmax,ymin,ymax};
  double half_range = 0.0;
  for(int i=0;i<ranges.size();i++){
    if( fabs(ranges[i]) > half_range ){
      half_range = fabs(ranges[i]);
    }
  }
  this->output_res = 2.0*half_range*this->Nj/this->width;
  this->writeProfile(filename,half_range);
}

void fromFITS::scaleProfile(){
  double max = 0.0;
  for(int i=0;i<this->mySource->Nm;i++){
    if( this->mySource->img[i] > max ){
      max = this->mySource->img[i];
    }
  }
  for(int i=0;i<this->mySource->Nm;i++){
    this->mySource->img[i] /= max;    
  }
}
