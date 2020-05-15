#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "json/json.h"

#include "contour-tracing.hpp"
#include "polygons.hpp"
#include "simplify_caustics.hpp"
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

  std::string out_path = argv[2];
  std::string output = out_path+"output/";
  
  // Read the cosmological parameters
  Json::Value cosmo;
  fin.open(output+"angular_diameter_distances.json",std::ifstream::in);
  fin >> cosmo;
  fin.close();

  
  // Initialize image plane
  double width  = root["instruments"][0]["field-of-view_x"].asDouble();
  double height = root["instruments"][0]["field-of-view_y"].asDouble();
  double res    = Instrument::getResolution(root["instruments"][0]["name"].asString());
  //================= END:PARSE INPUT =======================


  


  //=============== BEGIN:CREATE THE LENSES ====================
  const Json::Value jlens = root["lenses"][0];

  // Initialize mass model physical parameters
  jmembers = jlens["external_shear"].getMemberNames();
  
  std::vector<Nlpar*> ext_pars;
  for(int i=0;i<jmembers.size();i++){
    double value = jlens["external_shear"][jmembers[i]].asDouble();
    ext_pars.push_back( new Nlpar(jmembers[i],0,0,value,0,0,0) );
  }
  CollectionMassModels* mycollection = new CollectionMassModels(ext_pars);
  for(int i=0;i<ext_pars.size();i++){ delete(ext_pars[i]); }
  ext_pars.clear();

  // Initialize main mass model
  mycollection->models.resize(jlens["mass_model"].size());
  for(int k=0;k<jlens["mass_model"].size();k++){
    std::string mmodel = jlens["mass_model"][k]["type"].asString();

    if( mmodel == "custom" ){

      jmembers = jlens["mass_model"][k]["pars"].getMemberNames();
      std::map<std::string,std::string> pars;
      for(int i=0;i<jmembers.size();i++){
	pars[jmembers[i]] = jlens["mass_model"][k]["pars"][jmembers[i]].asString();
      }
      mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(mmodel,pars);

    } else if ( mmodel == "eagle" ){
      
    } else {
      
      jmembers = jlens["mass_model"][k]["pars"].getMemberNames();
      std::vector<Nlpar*> pars;
      for(int i=0;i<jmembers.size();i++){
	pars.push_back( new Nlpar(jmembers[i],0,0,jlens["mass_model"][k]["pars"][jmembers[i]].asDouble(),0,0,0) ); // only nam and val have meaning in this call
      }
      mycollection->models[k] = FactoryMassModel::getInstance()->createMassModel(mmodel,pars,cosmo[0]["Dls"].asDouble(),cosmo[0]["Ds"].asDouble());
    }
  }

  //  for(int i=0;i<mycollection->models.size();i++){
  //    mycollection->models[i]->printMassPars();
  //  }
  //  mycollection->printPhysPars();
  //================= END:CREATE THE LENSES ====================


  
  

  //=============== BEGIN:GET CRITICAL LINES AND CAUSTICS =======================
  int super_res_x = 10*( static_cast<int>(ceil(width/res)) );
  int super_res_y = 10*( static_cast<int>(ceil(height/res)) );
  ImagePlane detA(super_res_x,super_res_y,width,height);

  mycollection->detJacobian(&detA,&detA);
  for(int i=0;i<detA.Nm;i++){
    if( detA.img[i] > 0 ){
      detA.img[i] = 0;
    } else {
      detA.img[i] = 1;
    }
  }

  // Get detA contours on the image plane
  std::vector<Contour*> contours;
  mooreNeighborTracing(&detA,contours);

  /*
  // Remove every second point from the contours to make them smoother
  for(int i=0;i<contours.size();i++){
    bool toggle = false;
    std::vector<double> used;
    std::vector<double> unused;
    used.resize(0);
    unused.resize(0);
    std::partition_copy(contours[i]->x.begin(),contours[i]->x.end(),std::back_inserter(used),std::back_inserter(unused),[&toggle](int) { return toggle = !toggle; });
    contours[i]->x.swap(used);
    used.resize(0);
    unused.resize(0);
    std::partition_copy(contours[i]->y.begin(),contours[i]->y.end(),std::back_inserter(used),std::back_inserter(unused),[&toggle](int) { return toggle = !toggle; });
    contours[i]->y.swap(used);
  }
  */
  
  // Add first point as last and close the polygon
  for(int i=0;i<contours.size();i++){
    contours[i]->x.push_back( contours[i]->x[0] );
    contours[i]->y.push_back( contours[i]->y[0] );
  }

  // Create caustic contours, but don't fill them yet
  std::vector<Contour*> caustics(contours.size());
  for(int i=0;i<contours.size();i++){
    Contour* mycontour = new Contour();
    mycontour->x.resize(contours[i]->x.size());
    mycontour->y.resize(contours[i]->y.size());
    caustics[i] = mycontour;
  } 

  // Deflect the contours to create the caustics
  double xdefl,ydefl;
  for(int i=0;i<contours.size();i++){
    for(int j=0;j<contours[i]->x.size();j++){
      mycollection->all_defl(contours[i]->x[j],contours[i]->y[j],xdefl,ydefl);
      caustics[i]->x[j] = xdefl;
      caustics[i]->y[j] = ydefl;
    }
  }

  /*
  // Simplify the caustics (purely for visual purposes)
  std::vector<Contour*> simplified(contours.size());
  for(int i=0;i<simplified.size();i++){
    Contour* mycontour = new Contour();
    simplified[i] = mycontour;
  }
  simplifyPolygon(caustics,simplified);

  Json::Value json_simple;
  for(int i=0;i<simplified.size();i++){
    Json::Value simp_x;
    Json::Value simp_y;
    for(int j=0;j<simplified[i]->x.size();j++){
      simp_x.append(simplified[i]->x[j]);
      simp_y.append(simplified[i]->y[j]);
    }
    Json::Value simple;
    simple["x"] = simp_x;
    simple["y"] = simp_y;
    json_simple.append(simple);
  }
  for(int i=0;i<simplified.size();i++){
    delete(simplified[i]);
  }
  */
  
  // Create the json output object
  Json::Value json_caustics;
  Json::Value json_criticals;
  for(int i=0;i<contours.size();i++){
    Json::Value crit_x;
    Json::Value crit_y;
    Json::Value cau_x;
    Json::Value cau_y;
    for(int j=0;j<contours[i]->x.size();j++){
      crit_x.append(contours[i]->x[j]);
      crit_y.append(contours[i]->y[j]);
      cau_x.append(caustics[i]->x[j]);
      cau_y.append(caustics[i]->y[j]);
    }

    Json::Value caustic;
    caustic["x"] = cau_x;
    caustic["y"] = cau_y;
    json_caustics.append(caustic);
    Json::Value critical;
    critical["x"] = crit_x;
    critical["y"] = crit_y;
    json_criticals.append(critical);
  }
  
  for(int i=0;i<contours.size();i++){
    delete(contours[i]);
    delete(caustics[i]);
  }
  //================= END:GET CRITICAL LINES AND CAUSTICS =======================
  






  //=============== BEGIN:FIND NUMBER OF IMAGES AND LOCATION =======================
  point point_source = {root["point_source"]["x0"].asDouble(),root["point_source"]["y0"].asDouble()};

  // Create and deflect image plane
  std::vector<ImagePlane*> planes;
  ImagePlane* img = new ImagePlane(10,10,width,height);
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
      for(int i=0;i<planes[p]->Nm;i++){
	mycollection->all_defl(planes[p]->x[i],planes[p]->y[i],planes[p]->defl_x[i],planes[p]->defl_y[i]);
      }
      
      // Create triangle indices based on the image plane pixel indices
      std::vector<itriangle> triangles = imagePlaneToTriangleIndices(planes[p]);
      
      // Find which deflected triangles contain the point source
      std::vector<int> match;
      for(int i=0;i<triangles.size();i++){
	point p1 = {planes[p]->defl_x[triangles[i].ia],planes[p]->defl_y[triangles[i].ia]};
	point p2 = {planes[p]->defl_x[triangles[i].ib],planes[p]->defl_y[triangles[i].ib]};
	point p3 = {planes[p]->defl_x[triangles[i].ic],planes[p]->defl_y[triangles[i].ic]};
	
	if( pointInTriangle(point_source,p1,p2,p3) ){
	  match.push_back(i);
	}
      }
      
      // Get the center and radius of the circumcircle of each image triangle
      double xdum,ydum,rdum;
      for(int i=0;i<match.size();i++){
	point A = {planes[p]->x[triangles[match[i]].ia],planes[p]->y[triangles[match[i]].ia]};
	point B = {planes[p]->x[triangles[match[i]].ib],planes[p]->y[triangles[match[i]].ib]};
	point C = {planes[p]->x[triangles[match[i]].ic],planes[p]->y[triangles[match[i]].ic]};
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
	ImagePlane* img = new ImagePlane(10,10,2*rc_tmp[p],2*rc_tmp[p],xc_tmp[p],yc_tmp[p]);
	planes[p] = img;
      }
      //      std::cout << "Zooming in... (planes " << planes.size() << ")" <<  std::endl;
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
    multipleImages[i]->k = mycollection->all_kappa(x,y);
    double gamma_mag,gamma_phi;
    mycollection->all_gamma(x,y,gamma_mag,gamma_phi);
    multipleImages[i]->g    = gamma_mag;
    multipleImages[i]->phig = gamma_phi/0.01745329251 - 90.0; // in degrees east-of-north;
    multipleImages[i]->mag  = 1.0/mycollection->detJacobian(x,y);
  }
    
  // Calculate time delays
  std::vector<double> delays(multipleImages.size());
  for(int i=0;i<multipleImages.size();i++){
    double x = multipleImages[i]->x;
    double y = multipleImages[i]->y;
    double psi_tot = mycollection->all_psi(x,y);
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
  // Image plane magnification (0:positive, 1:negative)
  detA.writeImage(output + "detA.fits");

  // Caustics
  std::ofstream file_caustics(output+"caustics.json");
  file_caustics << json_caustics;
  file_caustics.close();
  //std::ofstream file_caustics(output+"caustics.json");
  //file_caustics << json_simple;
  //file_caustics.close();

  // Critical curves
  std::ofstream file_criticals(output+"criticals.json");
  file_criticals << json_criticals;
  file_criticals.close();

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
