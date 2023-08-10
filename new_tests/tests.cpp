#define CATCH_CONFIG_MAIN
#include <stdio.h>
#include <stdlib.h>

#include <CCfits/CCfits>
#include <valarray>
#include "json/json.h"

#include "catch2/catch.hpp"
//#include <catch2/matchers/catch_matchers_floating_point.hpp>



std::vector<double> readFits(int& N,std::string filepath){
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  
  int Ny = image.axis(0);
  int Nx = image.axis(1);
  N = Nx*Ny;
  std::vector<double> img(N);
  
  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int i=0;i<N;i++){
    img[i] = tmp[i];
  }

  return img;
}


double maximum_difference(std::vector<double> img1,std::vector<double> img2){
  double diff,max_diff = 0.0;
  for(int i=0;i<img1.size();i++){
    diff = fabs(img1[i]-img2[i]);
    if( diff > max_diff ){
      max_diff = diff;
    }
  }
  return max_diff;
}


void compare_image(std::string test1,std::string test2,double tol,std::string image_name){
  int N1,N2;

  std::vector<double> img1 = readFits(N1,test1+"output/"+image_name);
  std::vector<double> img2 = readFits(N2,test2+"output/"+image_name);
  bool dimensions = ( N1 == N2 );
  REQUIRE( N1 == N2 );

  if( dimensions ){
    REQUIRE_THAT( maximum_difference(img1,img2),Catch::Matchers::WithinAbs(0.0,tol) );
  }
}


TEST_CASE("Basic test", "[general]")
{
  std::string test1 = "general/test_A/";
  std::string test2 = "quasar/test_A/";
  double tol = 0.0001;
  

  SECTION("Comparing static images"){
    compare_image(test1,test2,tol,"OBS_testCAM-i.fits");
  }
  SECTION("Comparing static images"){
    compare_image(test1,test2,tol,"OBS_testCAM-i_ps_macro.fits");
  }


  
}
