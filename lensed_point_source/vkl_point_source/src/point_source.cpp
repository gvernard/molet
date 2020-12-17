#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "polygons.hpp"
#include "pointImage.hpp"

#include "vkllib.hpp"
#include "instruments.hpp"



int main(int argc,char* argv[]){

  /*
    Requires:
    - angular_diameter_distances.json
  */
  
  //=============== BEGIN:PARSE INPUT =======================
  std::ifstream fin;
  Json::Value::Members jmembers;

  // Read the main projection parameters
  Json::Value root;
  fin.open(argv[1],std::ifstream::in);
  fin >> root;
  fin.close();

  std::string in_path = argv[1];
  std::string input   = in_path+"input_files/";

  std::string out_path = argv[2];
  std::string output = out_path+"output/";

  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(output+"angular_diameter_distances.json",std::ifstream::in);
  fin >> cosmo;
  fin.close();

  // Initialize image plane
  double xmin  = root["instruments"][0]["field-of-view_xmin"].asDouble();
  double xmax  = root["instruments"][0]["field-of-view_xmax"].asDouble();
  double ymin  = root["instruments"][0]["field-of-view_ymin"].asDouble();
  double ymax  = root["instruments"][0]["field-of-view_ymax"].asDouble();
  double res    = Instrument::getResolution(root["instruments"][0]["name"].asString());
  //================= END:PARSE INPUT =======================





  //=============== BEGIN:CREATE THE LENSES ====================
  const Json::Value jlens = root["lenses"][0];

  // Initialize mass model collection
  CollectionMassModels mass_collection = JsonParsers::parse_mass_model(jlens["mass_model"],input);

  // Scale dpsi mass models if necessary
  for(int k=0;k<jlens["mass_model"].size();k++){
    if( jlens["mass_model"][k]["pars"].isMember("scale_factor") ){
      double scale_factor = jlens["mass_model"][k]["pars"]["scale_factor"].asDouble();
      Pert* pert = static_cast<Pert*> (mass_collection.models[k]);
      for(int m=0;m<pert->Sm;m++){
	pert->z[m] *= scale_factor;
      }
      pert->updateDerivatives();
    }
  }
  //================= END:CREATE THE LENSES ====================






  //=============== BEGIN:FIND NUMBER OF IMAGES AND LOCATION =======================
  point point_source = {root["point_source"]["x0"].asDouble(),root["point_source"]["y0"].asDouble()};

  // Create and deflect image plane
  std::vector<RectGrid*> planes;
  RectGrid* img = new RectGrid(10,10,xmin,xmax,ymin,ymax);
  planes.push_back(img);
  std::vector<double> xc;
  std::vector<double> yc;
  std::vector<double> rc;
  double final_scale = res/100.0;
  bool condition = true;

  while( condition ){
    std::vector<double> xc_tmp;
    std::vector<double> yc_tmp;
    std::vector<double> rc_tmp;

    for(int p=0;p<planes.size();p++){
      // Deflect image plane
      double* tmp_defl_x = (double*) malloc(planes[p]->Nz*sizeof(double));
      double* tmp_defl_y = (double*) malloc(planes[p]->Nz*sizeof(double));
      
      for(int i=0;i<planes[p]->Ny;i++){
	for(int j=0;j<planes[p]->Nx;j++){
	  mass_collection.all_defl(planes[p]->center_x[j],planes[p]->center_y[i],tmp_defl_x[i*planes[p]->Nx+j],tmp_defl_y[i*planes[p]->Nx+j]);
	}
      }
      
      // Create triangle indices based on the image plane pixel indices
      std::vector<itriangle> triangles = imagePlaneToTriangleIndices(planes[p]);
      
      // Find which deflected triangles contain the point source
      std::vector<int> match;
      for(int k=0;k<triangles.size();k++){
	int indA = triangles[k].ya*planes[p]->Nx+triangles[k].xa;
	int indB = triangles[k].yb*planes[p]->Nx+triangles[k].xb;
	int indC = triangles[k].yc*planes[p]->Nx+triangles[k].xc;
	point p1 = {tmp_defl_x[indA],tmp_defl_y[indA]};
	point p2 = {tmp_defl_x[indB],tmp_defl_y[indB]};
	point p3 = {tmp_defl_x[indC],tmp_defl_y[indC]};
	
	if( pointInTriangle(point_source,p1,p2,p3) ){
	  match.push_back(k);
	}
      }
      free(tmp_defl_x);
      free(tmp_defl_y);
      
      // Get the center and radius of the circumcircle of each image triangle
      double xdum,ydum,rdum;
      for(int i=0;i<match.size();i++){
	point A = {planes[p]->center_x[triangles[match[i]].xa],planes[p]->center_y[triangles[match[i]].ya]};
	point B = {planes[p]->center_x[triangles[match[i]].xb],planes[p]->center_y[triangles[match[i]].yb]};
	point C = {planes[p]->center_x[triangles[match[i]].xc],planes[p]->center_y[triangles[match[i]].yc]};
	circumcircle(A,B,C,xdum,ydum,rdum);
	xc_tmp.push_back(xdum);
	yc_tmp.push_back(ydum);
	rc_tmp.push_back(rdum);
      }
      
      // Maybe Write rectangular of the image plane
    }

    // Evaluate stopping criterion: all image planes must be smaller than some fraction of a pixel
    // Otherwise create new image planes and redefine planes vector
    int counter = 0;
    for(int i=0;i<rc_tmp.size();i++){
      if( rc_tmp[i] > final_scale ){
	counter++;
      }
    }
    if( counter > 0 ){
      // Set new smalle image planes around xc,yc
      for(int p=0;p<planes.size();p++){
	delete(planes[p]);
      }
      planes.resize(xc_tmp.size());
      for(int p=0;p<planes.size();p++){
	RectGrid* img = new RectGrid(10,10,xc_tmp[p]-rc_tmp[p],xc_tmp[p]+rc_tmp[p],yc_tmp[p]-rc_tmp[p],yc_tmp[p]+rc_tmp[p]);
	planes[p] = img;
      }
      //std::cout << "Zooming in... (planes " << planes.size() << ")" <<  std::endl;
    } else {
      // Exit loop and write multiple image positions
      copy(xc_tmp.begin(),xc_tmp.end(),back_inserter(xc)); 
      copy(yc_tmp.begin(),yc_tmp.end(),back_inserter(yc)); 
      copy(rc_tmp.begin(),rc_tmp.end(),back_inserter(rc)); 
      for(int p=0;p<planes.size();p++){
	delete(planes[p]);
      }
      condition = false;
    }   

  }
  
  // Filter images by location
  std::vector<double> xc_final;
  std::vector<double> yc_final;
  std::vector<double> rc_final;
  for(int i=0;i<xc.size()-1;i++){
    bool check = true;
    for(int j=i+1;j<xc.size();j++){
      double d = hypot(xc[i]-xc[j],yc[i]-yc[j]);
      if( d < final_scale ){
	check = false;
      }
    }
    if( check ){
      xc_final.push_back(xc[i]);
      yc_final.push_back(yc[i]);
      rc_final.push_back(rc[i]);
    }
  }
  xc_final.push_back(xc.back());
  yc_final.push_back(yc.back());
  rc_final.push_back(rc.back());
  
  std::vector<pointImage*> multipleImages(xc_final.size());
  for(int i=0;i<xc_final.size();i++){
    pointImage* img = new pointImage(xc_final[i],yc_final[i],2*rc_final[i],2*rc_final[i],0,0,0,-999,0,0);
    multipleImages[i] = img;
  }
  //================= END:FIND NUMBER OF IMAGES AND LOCATION =======================




  
  
  //=============== BEGIN:CORRESPONDING KAPPA, GAMMA, AND TIME DELAY =======================
  for(int i=0;i<multipleImages.size();i++){
    double x = multipleImages[i]->x;
    double y = multipleImages[i]->y;
    multipleImages[i]->k = mass_collection.all_kappa(x,y);
    double gamma_mag,gamma_phi;
    mass_collection.all_gamma(x,y,gamma_mag,gamma_phi);
    multipleImages[i]->g    = gamma_mag;
    multipleImages[i]->phig = gamma_phi/0.01745329251 - 90.0; // in degrees east-of-north;
    multipleImages[i]->mag  = 1.0/mass_collection.detJacobian(x,y);
  }
    
  // Calculate time delays
  std::vector<double> delays(multipleImages.size());
  for(int i=0;i<multipleImages.size();i++){
    double x = multipleImages[i]->x;
    double y = multipleImages[i]->y;
    double psi_tot = mass_collection.all_psi(x,y);
    double time = 0.5*(pow(point_source.x-x,2) + pow(point_source.y-y,2)) - psi_tot;
    //double time = 0.5*(pow(point_source.x-x,2) + pow(point_source.y-y,2));
    delays[i] = time;
  }

  double d_min = delays[0];
  for(int i=1;i<delays.size();i++){
    if( delays[i] < d_min ){
      d_min = delays[i];
    }
  }

  double factor = 0.0281*(1.0+jlens["redshift"].asDouble())*cosmo[0]["Dl"].asDouble()*cosmo[0]["Ds"].asDouble()/(cosmo[0]["Dls"].asDouble()); // in days
  //*** this factor has to be mutliplied by rad^2, i.e. converted from arcsec^2 that are the units of the potential and the other time delay term.
  for(int i=0;i<multipleImages.size();i++){
    multipleImages[i]->dt = (delays[i] - d_min)*factor;
  }
  //================= END:CORRESPONDING KAPPA, GAMMA, AND TIME DELAY =======================





  //=============== BEGIN:OUTPUT =======================
  // Multiple images
  Json::Value json_images;
  for(int i=0;i<multipleImages.size();i++){
    Json::Value image;
    image["x"]    = multipleImages[i]->x;
    image["y"]    = multipleImages[i]->y;
    image["dx"]   = multipleImages[i]->dx;
    image["dy"]   = multipleImages[i]->dy;
    image["k"]    = multipleImages[i]->k;
    image["g"]    = multipleImages[i]->g;
    image["phig"] = multipleImages[i]->phig;
    image["s"]    = multipleImages[i]->s;
    image["mag"]  = multipleImages[i]->mag;
    image["dt"]   = multipleImages[i]->dt;
    json_images.append(image);
  }
  std::ofstream file_images(output+"multiple_images.json");
  file_images << json_images;
  file_images.close();

  for(int i=0;i<multipleImages.size();i++){
    delete(multipleImages[i]);
  }
  //================= END:OUTPUT =======================


  
  
  return 0;
}
