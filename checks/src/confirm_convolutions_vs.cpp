#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "json/json.h"

#include "gerlumph.hpp"

int main(int argc,char* argv[]){
  
  /*
    Requires:
    - gerlumph_maps.json
    - multiple_images.json
    - tobs.json
  */


    //=============== BEGIN:INITIALIZE =======================
  std::ifstream fin;
  
  // Read the main parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();
  
  std::string out_path = argv[2];
  std::string output = out_path + "output/";
    
  // Read matching gerlumph map parameters
  Json::Value maps;
  fin.open(output+"gerlumph_maps.json",std::ifstream::in);
  fin >> maps;
  fin.close();

  // Read the multiple images' parameters from JSON to get the maximum time delay
  Json::Value multiple_images;
  fin.open(output+"multiple_images.json",std::ifstream::in);
  fin >> multiple_images;
  fin.close();
  // Get maximum image time delay
  double td_max = 0.0;
  for(int q=0;q<multiple_images.size();q++){
    double td = multiple_images[q]["dt"].asDouble();
    if( td > td_max ){
      td_max = td;
    }
  }
  
  // Read tobs_min and tobs_max
  Json::Value tobs_json;
  fin.open(output+"tobs.json",std::ifstream::in);
  fin >> tobs_json;
  fin.close();
  double tobs_max = tobs_json["tobs_max"].asDouble();
  double tobs_min = tobs_json["tobs_min"].asDouble();
  double duration = tobs_max - tobs_min;
  
  // Number of filters
  int Nfilters = root["instruments"].size();
  //================= END:INITIALIZE =======================




  //=============== BEGIN:REPORT NUMBER OF CONVOLUTIONS AND ASK FOR CONFIRMATION ===============
  int* Nsnaps = (int*) calloc(Nfilters,sizeof(int));
  for(int m=0;m<maps.size();m++){

    // Every multiple image has a different number of snapshots falling within the observed time frame because of the different time delays 
    if( maps[m]["id"].asString() != "none" ){
      for(int k=0;k<Nfilters;k++){
	std::string instrument_name = root["instruments"][k]["name"].asString();
	
	// Get time duration, start, and end, and set light curve number of sample points
	const Json::Value jtime = root["point_source"]["variability"]["extrinsic"][instrument_name]["time"];
	int Nsnap_tot = jtime.size();
	double Dtime = tobs_max-tobs_min;
	double tsrcmin = td_max - multiple_images[m]["dt"].asDouble();
	double tsrcmax = tsrcmin + Dtime;
	int jj_start;
	for(int jj=1;jj<Nsnap_tot;jj++){
	  if( jtime[jj].asDouble() > tsrcmin ){
	    jj_start = jj - 1;
	    break;
	  }
	}
	int jj_end;
	for(int jj=jj_start;jj<Nsnap_tot;jj++){
	  if( jtime[jj].asDouble() > tsrcmax ){
	    jj_end = jj;
	    break;
	  }
	}
	Nsnaps[k] += jj_end - jj_start + 1;
      }
    }
    
  }



    
  printf(">>>>>> Please pay ATTENTION to the following: <<<<<<\n");
  printf("   The number of convolutions per filter is:\n");
  int convs_tot = 0;
  for(int k=0;k<Nfilters;k++){
    std::string instrument_name = root["instruments"][k]["name"].asString();
    printf(" %s: %d \n",instrument_name.c_str(),Nsnaps[k]);
    convs_tot += Nsnaps[k];
  }
  printf("   The total number of convolutions will be %d.\n",convs_tot);
  double conv_fac = convs_tot/3600.0;
  printf("   It will approximately take %.2f hours (%d mins) on the CPU and %.2f hours (%d mins) on the GPU.\n",conv_fac*13.0,(int)std::ceil(conv_fac*13.0*60),conv_fac*3.0,(int)std::ceil(conv_fac*3.0*60)); // convolution on a CPU lasts 13s and on the GPU 3s (approximate numbers)
  printf("   (if you want to compille gerlumphpp for GPUs see here: https://github.com/gvernard/gerlumphpp)\n");
  printf("   If that is too long, consider shortening the observing time for each instrument, or increasing the fractional_increase parameter, or decreasing the size_cutoff parameter.\n");
  char choice;
  bool run = true;
  do{
    std::cout << "-> Would you like to proceed with the calculations? (answer 'y' or 'n')" << std::endl;
    std::cin >> choice;
    choice = std::tolower(choice);//Put your letter to its lower case
  } while( choice != 'n' && choice != 'y' );
  if( choice =='n' ){
    return 1;
  }
  //================= END:REPORT NUMBER OF CONVOLUTIONS AND ASK FOR CONFIRMATION ===============
  

  

  return 0;
}
