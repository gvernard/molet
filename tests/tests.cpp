#define CATCH_CONFIG_MAIN
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <random>

#include <CCfits/CCfits>
#include "json/json.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_test_case_info.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>


//=================================== UTILITY FUNCTIONS ===================================
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


double maximum_image_difference(std::vector<double> img1,std::vector<double> img2){
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
    REQUIRE_THAT( maximum_image_difference(img1,img2),Catch::Matchers::WithinAbs(0.0,tol) );
  }
}


Json::Value readJson(std::string filepath){
  std::ifstream fin;
  Json::Value json;
  fin.open(filepath,std::ifstream::in);
  fin >> json;
  fin.close();
  return json;
}

void compare_LC(std::string test1,std::string test2,double tol,std::string LC_name){
  Json::Value lc1 = readJson(test1+LC_name);
  Json::Value lc2 = readJson(test2+LC_name);

  REQUIRE( lc1.size() == lc2.size() );

  for(int q=0;q<lc1.size();q++){
    REQUIRE( lc1[q]["time"].size() == lc2[q]["time"].size() );
    REQUIRE( lc1[q]["signal"].size() == lc2[q]["signal"].size() );
    REQUIRE( lc1[q]["dsignal"].size() == lc2[q]["dsignal"].size() );
  }
  
  double diff_t,diff_time    = 0.0;
  double diff_s,diff_signal  = 0.0;
  double diff_d,diff_dsignal = 0.0;
  for(int q=0;q<lc1.size();q++){
    for(int i=0;i<lc1[q]["time"].size();i++){
      diff_t = fabs(lc1[q]["time"][i].asDouble()-lc2[q]["time"][i].asDouble());
      diff_s = fabs(lc1[q]["signal"][i].asDouble()-lc2[q]["signal"][i].asDouble());
      diff_d = fabs(lc1[q]["dsignal"][i].asDouble()-lc2[q]["dsignal"][i].asDouble());
      if( diff_t > diff_time ){
	diff_time = diff_t;
      }
      if( diff_s > diff_signal ){
	diff_signal = diff_s;
      }
      if( diff_d > diff_dsignal ){
	diff_dsignal = diff_d;
      }
    }
  }
  REQUIRE_THAT( diff_time,Catch::Matchers::WithinAbs(0.0,tol) );
  REQUIRE_THAT( diff_signal,Catch::Matchers::WithinAbs(0.0,tol) );
  REQUIRE_THAT( diff_dsignal,Catch::Matchers::WithinAbs(0.0,tol) );
}



void check_fluxes(std::string example,std::string instrument,double tol,bool check_PS){
  Json::Value json;
  std::ifstream fin(example+"/output/fluxes.json",std::ifstream::in);
  fin >> json;
  fin.close();
  double sum;
  
  double final_static = json["final_static"][instrument]["flux"].asDouble();
  double lens = json["lens_flux"][instrument]["total_convolved"]["flux"].asDouble();
  double extended = json["lensed_source_flux"][instrument]["total_convolved"]["flux"].asDouble();
  REQUIRE_THAT( (lens+extended),Catch::Matchers::WithinAbs(final_static,tol) );

  if( check_PS ){
    double final_with_ps = json["final_with_ps"][instrument]["flux"].asDouble();
    double lensed_ps = json["lensed_ps_flux"][instrument]["flux"].asDouble();
    REQUIRE_THAT( (lens+extended+lensed_ps),Catch::Matchers::WithinAbs(final_with_ps,tol) );
  }  
}


void run_molet(std::string example,std::string input="/molet_input.json"){
  std::string molet = "../molet_driver";
  char line_before[107],line_after[107];
  const char *padding="######################################################";

  int len = example.length();
  int padLen = 55 - len;
  sprintf(line_before,"%.*s %s %*.*s",50,padding,example.c_str(),padLen,padLen,padding);
  sprintf(line_after,"%.*s%.*s%.*s",50,padding,len+2,padding,padLen,padding);

  std::cout << line_before << std::endl;
  int exit = system((molet+" "+example+"/"+input).c_str());
  std::cout << line_after << std::endl;

  REQUIRE( exit == 0 );
}

//=================================== UTILITY FUNCTIONS ===================================




//=================================== EVENT LISTENERS ===================================
class testCaseListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void testCaseStarting(Catch::TestCaseInfo const& info) override {
      std::cout << "Test case: " << info.name << std::endl;
      
      if( info.name == "quasar" ){
	run_molet("general/test_A");
      } else if( info.name == "static" ){
	run_molet("general/test_D");
      }
    }
};
CATCH_REGISTER_LISTENER(testCaseListener)


class testSectionListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void sectionStarting(Catch::SectionInfo const& info) override {
      std::cout << "Section: " << info.name << std::endl;

      if( info.name == "QA" ){
	// Quasar first
	run_molet("quasar/test_A");
      } else if( info.name == "QB" ){
	run_molet("quasar/test_B");
      } else if( info.name == "QC" ){
	run_molet("quasar/test_C");
      } else if( info.name == "QD" ){
	run_molet("quasar/test_D");
	// Then static
      } else if( info.name == "SA" ){
	run_molet("static/test_A");
      } else if( info.name == "SB" ){
	run_molet("static/test_B");
      } else if( info.name == "SC" ){
	run_molet("static/test_C","coolest.json");
      } else if( info.name == "SD" ){
	run_molet("static/test_D");
      } else if( info.name == "SE" ){
	run_molet("static/test_E");
      } else if( info.name == "SF" ){
	run_molet("static/test_F");
      } else {
	// Do nothing
      }
    }
};
CATCH_REGISTER_LISTENER(testSectionListener)
//=================================== EVENT LISTENERS ===================================




double img_tol  = 0.0001;
double flux_tol = 0.01;
double lc_tol   = 0.01;
double lc_tol2  = 0.1;



TEST_CASE("static")
{
  std::string base_test = "general/test_D/";
  std::string test;

  
  SECTION("SA","Same as general/test_D but with a SPEMD mass model instead of an SIE."){
    std::cout << std::endl;
    test = "static/test_A/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_noiseless.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,false);
  }

  
  SECTION("SB","Same as general/test_D (and static/test_A) but with the Gaussian sources replaced by equivalent Sersic ones."){
    std::cout << std::endl;
    test = "static/test_B/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_noiseless.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,false);
  }

  
  SECTION("SC","Same as general/test_D (and static/test_A, static/test_B) but with input given in the COOLEST standard."){
    std::cout << std::endl;
    test = "static/test_C/MOLET_RUN/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_noiseless.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,false);
  }

  
  SECTION("SD","Is the same as general/test_D but the Gaussian sources are replaced with a pixelated source of the same total flux."){
    std::cout << std::endl;
    test = "static/test_D/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,0.65,"OBS_"+mycam+"_noiseless.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,false);
  }

  
  SECTION("SE","Is the same as general/test_D but the lens potential includes perturbations."){
    std::cout << std::endl;
    test = "static/test_E/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,0.5,"OBS_"+mycam+"_noiseless.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,false);
  }


  SECTION("SF","Is the same as general/test_D but there is a smaller satellite galaxy to the lens."){
    std::cout << std::endl;
    test = "static/test_F/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,2*flux_tol,false);
  }  
}


  

TEST_CASE("quasar")
{
  std::string base_test = "general/test_A/";
  std::string test;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distr(0,9);
  int digit1,digit2;
  

  SECTION("QA","Using fixed extrinsic light curves calculated in general/test_A."){
    std::cout << std::endl;
    test = "quasar/test_A/";
    std::string mycam = "testCAM-i";
    
    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+".fits");
    std::cout << "Comparing PS-macro images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_ps_macro.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,true);

    std::cout << "Comparing intrinsic light curves" << std::endl;
    compare_LC(base_test,test,lc_tol,"input_files/"+mycam+"_LC_intrinsic.json");
    digit1 = distr(gen);
    digit2 = distr(gen);
    std::string mock = "mock_0000_00"+std::to_string(digit1)+std::to_string(digit2);
    std::cout << "Comparing continuous light curves for " << mock  << std::endl;
    compare_LC(base_test,test,lc_tol,mock+"/"+mycam+"_LC_continuous.json");
    std::cout << "Comparing sampled light curves for " << mock << std::endl;
    compare_LC(base_test,test,lc_tol,mock+"/"+mycam+"_LC_sampled.json");
  }
 

  SECTION("QB","Same as test_A, but using a compact mass profile instead of a mass-to-light ratio."){
    std::cout << std::endl;
    test = "quasar/test_B/";
    std::string mycam = "testCAM-i";

    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+".fits");
    std::cout << "Comparing PS-macro images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_ps_macro.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,true);

    std::cout << "Comparing intrinsic light curves" << std::endl;
    compare_LC(base_test,test,lc_tol,"input_files/"+mycam+"_LC_intrinsic.json");
  }


  SECTION("QC","Using the output microlensing light curves of test_A as input and adding an unmicrolensed component."){
    std::cout << std::endl;
    test = "quasar/test_C/";
    std::string mycam = "testCAM-i";

    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+".fits");
    std::cout << "Comparing PS-macro images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_ps_macro.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,true);
    
    std::cout << "Comparing intrinsic light curves" << std::endl;
    compare_LC(base_test,test,lc_tol,"input_files/"+mycam+"_LC_intrinsic.json");
  }


  SECTION("QD","Same as test_A, but a custom profile is used for the accretion disc - the same as the analytic one of test_A but cast into a pixelated grid."){
    std::cout << std::endl;
    test = "quasar/test_D/";
    std::string mycam = "testCAM-i";

    std::cout << "Comparing static images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+".fits");
    std::cout << "Comparing PS-macro images" << std::endl;
    compare_image(base_test,test,img_tol,"OBS_"+mycam+"_ps_macro.fits");
    std::cout << "Checking fluxes" << std::endl;
    check_fluxes(test,mycam,flux_tol,true);

    std::cout << "Comparing intrinsic light curves" << std::endl;
    compare_LC(base_test,test,lc_tol,"input_files/"+mycam+"_LC_intrinsic.json");
    digit1 = distr(gen);
    digit2 = distr(gen);
    std::string mock = "mock_0000_00"+std::to_string(digit1)+std::to_string(digit2);
    std::cout << "Comparing continuous light curves for " << mock  << std::endl;
    compare_LC(base_test,test,lc_tol2,mock+"/"+mycam+"_LC_continuous.json");
    std::cout << "Comparing sampled light curves for " << mock << std::endl;
    compare_LC(base_test,test,lc_tol2,mock+"/"+mycam+"_LC_sampled.json");
  }
}

