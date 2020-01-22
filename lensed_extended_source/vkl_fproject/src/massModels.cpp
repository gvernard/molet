#include "massModels.hpp"

#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include "imagePlane.hpp"
#include "sourcePlane.hpp"
#include "tableDefinition.hpp"


//Abstract class: BaseMassModel
//===============================================================================================================
void BaseMassModel::setMassPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    this->mpars[nlpars[i]->nam] = nlpars[i]->val;
  }
}

void BaseMassModel::printMassPars(){
  typedef std::map<std::string,double>::iterator some_type;
  for(some_type iterator=this->mpars.begin();iterator!=this->mpars.end();iterator++){
    std::cout << iterator->first << " " << iterator->second << std::endl;
    //    printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
  }
}


//Derived class from BaseMassModel: Sie (Singular Isothermal Ellipsoid)
//===============================================================================================================
Sie::Sie(std::vector<Nlpar*> nlpars){
  this->n = 5;
  this->type = "sie";
  setMassPars(nlpars);
}

void Sie::defl(double xin,double yin,double& xout,double& yout){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  if( fabs(x_t) < 0.0001 && fabs(y_t) < 0.0001 ){
    if( std::signbit(x_t) ){
      x_t = -0.0001;
    } else {
      x_t =  0.0001;
    }
    if( std::signbit(y_t) ){
      y_t = -0.0001;
    } else {
      y_t =  0.0001;
    }
  }

  double fac   = 1.0-q*q;
  double omega = q*q*x_t*x_t + y_t*y_t; // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double fac2  = sqrt(fac/omega);

  double ax_t = (b/sqrt(fac))*atan(x_t*fac2);
  double ay_t = (b/sqrt(fac))*atanh(y_t*fac2);
  
  //rotate back according to position angle, no need to translate (this is equivalent to rotating by -pa using the same equations as above)
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

//Derived class from BaseMassModel: Spemd (Softened Power-law Elliptical Mass Density)
//===============================================================================================================
Spemd::Spemd(std::vector<Nlpar*> nlpars){
  this->n = 7;
  this->type = "spemd";
  setMassPars(nlpars);
}

void Spemd::defl(double xin,double yin,double& xout,double& yout){
  double q  = this->mpars["q"];
  double e  = this->mpars["e"];
  double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e)/(q*2.0);
  double pa = this->mpars["pa"] * 0.01745329251;//in rad
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double s2 = this->mpars["s"] * this->mpars["s"];
  double defl[2] = {0.0,0.0};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  fastelldefl_(&x_t,&y_t,&b,&e,&q,&s2,defl);

  double ax_t = defl[0];
  double ay_t = defl[1];

  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

//Derived class from BaseMassModel: Pert (perturbations on a grid)
//===============================================================================================================
Pert::Pert(int a,int b,double c,double d,std::string reg){
  this->type = "pert";
  this->dpsi = new FixedSource(a,b,c,d,reg);

  this->dpsi_dx = (double*) calloc(this->dpsi->Sm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Sm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Sm;
  this->Bdev.Tj = this->dpsi->Sm;

  this->di = this->dpsi->height/(this->dpsi->Si);
  this->dj = this->dpsi->width/(this->dpsi->Sj);

  //  createBdev();
}

Pert::Pert(int a,int b,ImagePlane* image,std::string reg){
  this->type = "pert";
  this->dpsi = new FixedSource(a,b,image->xmin,image->xmax,image->ymin,image->ymax,reg);
  this->dpsi->inMask(image);

  this->dpsi_dx = (double*) calloc(this->dpsi->Sm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Sm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Sm;
  this->Bdev.Tj = this->dpsi->Sm;

  this->di = this->dpsi->height/(this->dpsi->Si);
  this->dj = this->dpsi->width/(this->dpsi->Sj);

  //  createBdev();
}

Pert::Pert(std::string filepath,int a,int b,double c,double d,std::string reg){
  this->type = "pert";
  this->dpsi = new FixedSource(a,b,c,d,reg);
  
  this->dpsi_dx = (double*) calloc(this->dpsi->Sm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Sm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Sm;
  this->Bdev.Tj = this->dpsi->Sm;

  this->di = this->dpsi->height/(this->dpsi->Si);
  this->dj = this->dpsi->width/(this->dpsi->Sj);

  ImagePlane* image = new ImagePlane(filepath,a,b,c,d);
  replaceDpsi(image->img);
  delete(image);

  updatePert();
  //  createBdev();
}

Pert::Pert(FixedSource* new_dpsi){
  this->type = "pert";
  this->dpsi = new FixedSource(*new_dpsi);

  this->dpsi_dx = (double*) calloc(this->dpsi->Sm,sizeof(double));
  this->dpsi_dy = (double*) calloc(this->dpsi->Sm,sizeof(double));

  this->Bdev.Ti = 2*this->dpsi->Sm;
  this->Bdev.Tj = this->dpsi->Sm;

  this->di = this->dpsi->height/(this->dpsi->Si);
  this->dj = this->dpsi->width/(this->dpsi->Sj);

  updatePert();
  //  createBdev();
}

void Pert::replaceDpsi(double* new_dpsi){
  for(int i=0;i<this->dpsi->Sm;i++){
    this->dpsi->src[i] = new_dpsi[i];
  }
}

void Pert::addDpsi(double* corrections){
  for(int i=0;i<this->dpsi->Sm;i++){
    this->dpsi->src[i] += corrections[i];
  }
}

void Pert::updatePert(){
  int Si = this->dpsi->Si;
  int Sj = this->dpsi->Sj;

  // Calculate derivatives:
  // first row
  this->dpsi_dx[0] = (this->dpsi->src[1]-this->dpsi->src[0])/dj;
  this->dpsi_dy[0] = (this->dpsi->src[Sj]-this->dpsi->src[0])/di;
  for(int j=1;j<Sj-1;j++){
    this->dpsi_dx[j] = (this->dpsi->src[j+1]-this->dpsi->src[j-1])/(2*dj);
    this->dpsi_dy[j] = (this->dpsi->src[Sj+j]-this->dpsi->src[j])/di;
  }
  this->dpsi_dx[Sj-1] = (this->dpsi->src[Sj-1]-this->dpsi->src[Sj-2])/dj;
  this->dpsi_dy[Sj-1] = (this->dpsi->src[Sj+Sj-1]-this->dpsi->src[Sj-1])/di;

  // in-between rows
  for(int i=1;i<Si-1;i++){
    this->dpsi_dx[i*Sj+0] = (this->dpsi->src[i*Sj+1]-this->dpsi->src[i*Sj+0])/dj;
    this->dpsi_dy[i*Sj+0] = (this->dpsi->src[(i+1)*Sj+0]-this->dpsi->src[(i-1)*Sj+0])/(2*di);
    for(int j=1;j<Sj-1;j++){
      this->dpsi_dx[i*Sj+j] = (this->dpsi->src[i*Sj+j+1]-this->dpsi->src[i*Sj+j-1])/(2*dj);
      this->dpsi_dy[i*Sj+j] = (this->dpsi->src[(i+1)*Sj+j]-this->dpsi->src[(i-1)*Sj+j])/(2*di);
    }
    this->dpsi_dx[i*Sj+Sj-1] = (this->dpsi->src[i*Sj+Sj-1]-this->dpsi->src[i*Sj+Sj-2])/dj;
    this->dpsi_dy[i*Sj+Sj-1] = (this->dpsi->src[(i+1)*Sj+Sj-1]-this->dpsi->src[(i-1)*Sj+Sj-1])/(2*di);
  }

  // last row
  this->dpsi_dx[(Si-1)*Sj+0] = (this->dpsi->src[(Si-1)*Sj+1]-this->dpsi->src[(Si-1)*Sj+0])/dj;
  this->dpsi_dy[(Si-1)*Sj+0] = (this->dpsi->src[(Si-1)*Sj+0]-this->dpsi->src[(Si-2)*Sj+0])/di;
  for(int j=1;j<Sj-1;j++){
    this->dpsi_dx[(Si-1)*Sj+j] = (this->dpsi->src[(Si-1)*Sj+j+1]-this->dpsi->src[(Si-1)*Sj+j-1])/(2*dj);
    this->dpsi_dy[(Si-1)*Sj+j] = (this->dpsi->src[(Si-1)*Sj+j]-this->dpsi->src[(Si-2)*Sj+j])/di;
  }
  this->dpsi_dx[(Si-1)*Sj+Sj-1] = (this->dpsi->src[(Si-1)*Sj+Sj-1]-this->dpsi->src[(Si-1)*Sj+Sj-2])/dj;
  this->dpsi_dy[(Si-1)*Sj+Sj-1] = (this->dpsi->src[(Si-1)*Sj+Sj-1]-this->dpsi->src[(Si-2)*Sj+Sj-1])/di;

}

void Pert::createBdev(){
  // This function creates the table of finite difference coefficients for the x and y derivatives of dpsi.
  // The result is a Ndpsi x Ndpsi sparse matrix.
  // The coefficients are the same as in function Pert::updatePert
  int Sj = this->dpsi->Sj;
  int Si = this->dpsi->Si;

  std::vector<mytriplet> tmp;//need to make sure that the Ddpsi triplet vector is a new one

  // Calculate derivatives:
  // first row
  tmp.push_back({  0,  0,         -1.0/dj });
  tmp.push_back({  0,  1,          1.0/dj });
  tmp.push_back({  1,  0,         -1.0/di });
  tmp.push_back({  1, Sj,          1.0/di });
  for(int j=1;j<Sj-1;j++){
    tmp.push_back({   2*j,  j-1, -1.0/(2.0*dj) });
    tmp.push_back({   2*j,  j+1,  1.0/(2.0*dj) });
    tmp.push_back({ 2*j+1,    j,       -1.0/di });
    tmp.push_back({ 2*j+1, j+Sj,        1.0/di });
  }
  tmp.push_back({   2*(Sj-1),    Sj-2,      -1.0/dj });
  tmp.push_back({   2*(Sj-1),    Sj-1,       1.0/dj });
  tmp.push_back({ 2*(Sj-1)+1,    Sj-1,      -1.0/di });
  tmp.push_back({ 2*(Sj-1)+1, Sj+Sj-1,       1.0/di });

  // in-between rows
  for(int i=1;i<Si-1;i++){
    tmp.push_back({   2*(i*Sj),     i*Sj,       -1.0/dj });
    tmp.push_back({   2*(i*Sj),   i*Sj+1,        1.0/dj });
    tmp.push_back({ 2*(i*Sj)+1, (i-1)*Sj, -1.0/(2.0*di) });
    tmp.push_back({ 2*(i*Sj)+1, (i+1)*Sj,  1.0/(2.0*di) });
    for(int j=1;j<Sj-1;j++){
      tmp.push_back({   2*(i*Sj+j),   i*Sj+j-1, -1.0/(2.0*dj) });
      tmp.push_back({   2*(i*Sj+j),   i*Sj+j+1,  1.0/(2.0*dj) });
      tmp.push_back({ 2*(i*Sj+j)+1, (i-1)*Sj+j, -1.0/(2.0*di) });
      tmp.push_back({ 2*(i*Sj+j)+1, (i+1)*Sj+j,  1.0/(2.0*di) });
    }
    tmp.push_back({   2*(i*Sj+Sj-1),     i*Sj+Sj-2,       -1.0/dj });
    tmp.push_back({   2*(i*Sj+Sj-1),     i*Sj+Sj-1,        1.0/dj });
    tmp.push_back({ 2*(i*Sj+Sj-1)+1, (i-1)*Sj+Sj-1, -1.0/(2.0*di) });
    tmp.push_back({ 2*(i*Sj+Sj-1)+1, (i+1)*Sj+Sj-1,  1.0/(2.0*di) });
  }

  // last row
  tmp.push_back({   2*(Si-1)*Sj,   (Si-1)*Sj,       -1.0/dj });
  tmp.push_back({   2*(Si-1)*Sj, (Si-1)*Sj+1,        1.0/dj });
  tmp.push_back({ 2*(Si-1)*Sj+1,   (Si-2)*Sj,       -1.0/di });
  tmp.push_back({ 2*(Si-1)*Sj+1,   (Si-1)*Sj,        1.0/di });
  for(int j=1;j<Sj-1;j++){
    tmp.push_back({   2*((Si-1)*Sj+j), (Si-1)*Sj+j-1, -1.0/(2.0*dj) });
    tmp.push_back({   2*((Si-1)*Sj+j), (Si-1)*Sj+j+1,  1.0/(2.0*dj) });
    tmp.push_back({ 2*((Si-1)*Sj+j)+1,   (Si-2)*Sj+j,       -1.0/di });
    tmp.push_back({ 2*((Si-1)*Sj+j)+1,   (Si-1)*Sj+j,        1.0/di });
  }
  tmp.push_back({   2*((Si-1)*Sj+Sj-1), (Si-1)*Sj+Sj-2,   -1.0/dj });
  tmp.push_back({   2*((Si-1)*Sj+Sj-1), (Si-1)*Sj+Sj-1,    1.0/dj });
  tmp.push_back({ 2*((Si-1)*Sj+Sj-1)+1, (Si-2)*Sj+Sj-1,   -1.0/di });
  tmp.push_back({ 2*((Si-1)*Sj+Sj-1)+1, (Si-1)*Sj+Sj-1,    1.0/di });

  this->Bdev.tri.swap(tmp);
}

void Pert::createAint(ImagePlane* data){
  this->Aint.Ti = 2*data->Nm;
  this->Aint.Tj = 2*this->dpsi->Sm;

  std::vector<mytriplet> tmp;

  int i,j;
  double xa,ya,xb,yb,w00,w10,w01,w11,f00,f10,f01,f11;
  double den = this->di*this->dj;

  int Sj = this->dpsi->Sj;

  for(int k=0;k<data->Nm;k++){
    int j = floor( (data->x[k]-this->dpsi->xmin)/this->dj );
    int i = floor( (data->y[k]-this->dpsi->ymin)/this->di );  

    if( j == this->dpsi->Sj-1 ){
      j = j-2;
    }
    if( i == this->dpsi->Si-1 ){
      i = i-2;
    }

    ya  = data->y[k] - this->dpsi->y[i+1];
    yb  = this->dpsi->y[i] - data->y[k];
    xa  = data->x[k] - this->dpsi->x[j];
    xb  = this->dpsi->x[j+1] - data->x[k];

    w00 = xb*ya/den;
    w10 = xb*yb/den;
    w01 = xa*ya/den;
    w11 = xa*yb/den;

    tmp.push_back({ 2*k,       2*(i*Sj+j),  w00 });
    tmp.push_back({ 2*k,   2*((i+1)*Sj+j),  w10 });
    tmp.push_back({ 2*k,     2*(i*Sj+j+1),  w01 });
    tmp.push_back({ 2*k, 2*((i+1)*Sj+j+1),  w11 });

    tmp.push_back({ 2*k+1,       2*(i*Sj+j)+1,  w00 });
    tmp.push_back({ 2*k+1,   2*((i+1)*Sj+j)+1,  w10 });
    tmp.push_back({ 2*k+1,     2*(i*Sj+j+1)+1,  w01 });
    tmp.push_back({ 2*k+1, 2*((i+1)*Sj+j+1)+1,  w11 });
  }

  this->Aint.tri.swap(tmp);
}


void Pert::createCrosses(ImagePlane* image){
  int i0,j0;
  double xa,ya,xb,yb,w00,w10,w01,w11,f00,f10,f01,f11;
  double den = this->di*this->dj;
  double four_i[4]; // to hold the indices of the four vertices of a dpsi pixel
  double four_j[4]; // same as above
  double four_w[4]; // to hold the interpolation weights between the ray and each dpsi pixel vertex (in the image plane)
  int rel_ind[3];   // to hold the finite difference coefficients for 
  double coeffs[3];
  int cross_index;

  for(int h=0;h<image->Nm;h++){
    // Step 1: find the i-j indices of the top-left point in the dpsi grid, in the dpsi pixel where the ray lies (in the image plane)
    j0 = floor( (image->x[h] - this->dpsi->xmin)/this->dj );
    i0 = floor( (this->dpsi->ymax - image->y[h])/this->di ); // y-axis is reflected (indices i start from the top)
    if( j0 >= this->dpsi->Sj-1 ){
      j0 = this->dpsi->Sj - 2;
    }
    if( i0 >= this->dpsi->Si-1 ){
      i0 = this->dpsi->Si - 2;
    }
    four_i[0] = i0;
    four_i[1] = i0;
    four_i[2] = i0+1;
    four_i[3] = i0+1;
    four_j[0] = j0;
    four_j[1] = j0+1;
    four_j[2] = j0;
    four_j[3] = j0+1;

    // Step 2: get the interpolation weights for the ray in the dpsi pixel
    ya = image->y[h]                           - this->dpsi->y[(i0+1)*this->dpsi->Sj];
    yb = this->dpsi->y[i0*this->dpsi->Sj]      - image->y[h];
    xa = image->x[h]                           - this->dpsi->x[i0*this->dpsi->Sj+j0];
    xb = this->dpsi->x[i0*this->dpsi->Sj+j0+1] - image->x[h];
    four_w[0] = xb*ya/den;
    four_w[1] = xa*ya/den;
    four_w[2] = xb*yb/den;
    four_w[3] = xa*yb/den;

    // Step 3: loop over the dpsi pixel vertices
    //         cross_coeffs_x and cross_coeff_y have a particular structure, different in each case,
    //         which plays a role in how they are combined to create the final entries in the D_s(s_p)*D_psi table
    Cross* a_cross = new Cross(four_i[0],four_j[0]);
    for(int k=0;k<4;k++){
      // X derivatives
      this->derivativeDirection(four_j[k],this->dpsi->Sj-1,this->dj,rel_ind,coeffs);
      for(int q=0;q<3;q++){
	cross_index = (1 + (four_i[k]-i0)*4 + (four_j[k]-j0)) + rel_ind[q];
	a_cross->coeff_x[cross_index] += (coeffs[q]*four_w[k]);
      }

      // Y derivatives
      this->derivativeDirection(four_i[k],this->dpsi->Si-1,this->di,rel_ind,coeffs);
      for(int q=0;q<3;q++){
	cross_index = (1 + (four_j[k]-j0)*4 + (four_i[k]-i0)) + rel_ind[q];
	a_cross->coeff_y[cross_index] += (coeffs[q]*four_w[k]);
      }
    }
    image->crosses[h] = a_cross;
  }

}


/*
void Pert::createInterpolationWeights(ImagePlane* image){
  // This function is almost identical to the one that creates the interpolation weights for a source class.
  // The are two differences: the non-deflected x,y of the image plane are used, and the weights go in the dpsi_cells structure of the image class.
  double xp,yp;
  double dx     = this->x[1] - this->x[0];
  double dy     = this->y[this->Sj] - this->y[0];
  double norm   = 1.0/(dx*dy);
  int ic        = 0;
  int jc        = 0;
  double g1=0,g2=0,g3=0,g4=0;

  for(int i=0;i<image->Nm;i++){
    xp = image->x[i] - this->dpsi->xmin;
    yp = image->y[i] - this->dpsi->ymin;
    
    if( this->dpsi->pointInPolygon(xp,yp) ){
      //Indices corresponding to the bottom left pixel
      jc = (int) floor(xp/dx);
      ic = (int) floor(yp/dy);
      
      //Now interpolate between neighbouring pixels and add entries to L matrix.
      //The interpolation function could return an array of column indices in the row of the L matrix, and the corresponding weights.
      //The following is for bi-linear interpolation:
      g1 = (ic+1)*dy - yp;
      g2 = (jc+1)*dx - xp;
      g3 = yp - ic*dy;
      g4 = xp - jc*dx;
      
      delete(image->dpsi_cells[i]);
      InterpolationCell* cell = new InterpolationCell(4);
      cell->ind[0] = (this->Si-2-ic)*this->Sj + jc;
      cell->ind[1] = (this->Si-2-ic)*this->Sj + jc+1;
      cell->ind[2] = (this->Si-2-ic+1)*this->Sj + jc;
      cell->ind[3] = (this->Si-2-ic+1)*this->Sj + jc+1;
      cell->wei[0] = g1*g2*norm;
      cell->wei[1] = g1*g4*norm;
      cell->wei[2] = g3*g2*norm;
      cell->wei[3] = g3*g4*norm;
      image->dpsi_cells[i] = cell;
    } else {
      delete(image->dpsi_cells[i]);
      InterpolationCell* cell = new InterpolationCell(1);
      cell->ind[0] = 0;
      cell->wei[0] = 0.0;
      image->dpsi_cells[i] = cell;
    }

  }
}
*/

//private
void Pert::derivativeDirection(int q,int qmax,double den,int* rel_ind,double* coeff){
  if( q == 0 ){
    // forward
    rel_ind[0] = 0;
    rel_ind[1] = 1;
    rel_ind[2] = 2;
    coeff[0] = (-3.0/2.0)/den;
    coeff[1] = (2.0)/den;
    coeff[2] = (-1.0/2.0)/den;
  } else if( q == qmax ){
    // backward
    rel_ind[0] = 0;
    rel_ind[1] = -1;
    rel_ind[2] = -2;
    coeff[0] = (3.0/2.0)/den;
    coeff[1] = (-2.0)/den;
    coeff[2] = (1.0/2.0)/den;
  } else {
    // central
    rel_ind[0] = -1;
    rel_ind[1] = 0;
    rel_ind[2] = 1;
    coeff[0] = (-1.0/2.0)/den;
    coeff[1] = 0.0;
    coeff[2] = (1.0/2.0)/den;
  }
}

void Pert::defl(double xin,double yin,double& xout,double& yout){
  
  double xp = xin - this->dpsi->xmin;
  double yp = this->dpsi->ymax - yin;

  //Indices corresponding to the top left pixel
  int jc = (int) floor(xp/this->dj);
  int ic = (int) floor(yp/this->di);

  double g1 = yp - ic*this->di;
  double g2 = (ic+1)*this->di - yp;
  double g3 = xp - jc*this->dj;
  double g4 = (jc+1)*this->dj - xp;
  double norm = 1.0/(this->di*this->dj);

  //changing to weights for the top left pixel
  double w00 = g2*g4;
  double w01 = g2*g3;
  double w10 = g1*g4;
  double w11 = g1*g3;

  double f00,f01,f10,f11;

  f00 = this->dpsi_dx[ic*this->dpsi->Sj + jc];
  f01 = this->dpsi_dx[ic*this->dpsi->Sj + jc+1];
  f10 = this->dpsi_dx[(ic+1)*this->dpsi->Sj + jc];
  f11 = this->dpsi_dx[(ic+1)*this->dpsi->Sj + jc+1];
  double ax = (f00*w00 + f10*w10 + f01*w01 + f11*w11)*norm;

  f00 = this->dpsi_dy[ic*this->dpsi->Sj + jc];
  f01 = this->dpsi_dy[ic*this->dpsi->Sj + jc+1];
  f10 = this->dpsi_dy[(ic+1)*this->dpsi->Sj + jc];
  f11 = this->dpsi_dy[(ic+1)*this->dpsi->Sj + jc+1];
  double ay = (f00*w00 + f10*w10 + f01*w01 + f11*w11)*norm;

  xout = ax;
  yout = ay;
}

/*
void Pert::tableDefl(int Nm,double* xdefl,double* ydefl){
  Eigen::SparseMatrix<double> A(2*Nm,2*this->dpsi->Sm);
  A.reserve(Eigen::VectorXi::Constant(2*Nm,2));//overestimating the mask matrix number of entries per row
  for(int i=0;i<this->Aint.tri.size();i++){  A.insert(this->Aint.tri[i].i,this->Aint.tri[i].j) = this->Aint.tri[i].v;  }

  Eigen::SparseMatrix<double> B(2*this->dpsi->Sm,this->dpsi->Sm);
  B.reserve(Eigen::VectorXi::Constant(2*Nm,2));//overestimating the mask matrix number of entries per row
  for(int i=0;i<this->Bdev.tri.size();i++){  B.insert(this->Bdev.tri[i].i,this->Bdev.tri[i].j) = this->Bdev.tri[i].v;  }

  Eigen::Map<Eigen::VectorXd> dpsi(this->dpsi->src,this->dpsi->Sm);

  Eigen::VectorXd defl(2*Nm);

  defl = A*B*dpsi;

  for(int i=0;i<Nm;i++){
    xdefl[i] = defl[2*i];
    ydefl[i] = defl[2*i+1];
  }

  A.resize(0,0);
  B.resize(0,0);
}
*/

void Pert::getConvergence(ImagePlane* kappa){
  double ddx,ddy;
  int Ni = kappa->Ni;
  int Nj = kappa->Nj;

  double dx2 = pow( fabs(kappa->x[1]  - kappa->x[0]), 2);
  double dy2 = pow( fabs(kappa->y[Nj] - kappa->y[0]), 2);


  // Calculate second derivatives:

  /*
  // first row
  ddx = 1.0*this->dpsi->src[0] - 2.0*this->dpsi->src[1] + 1.0*this->dpsi->src[2];
  ddy = 1.0*this->dpsi->src[0] - 2.0*this->dpsi->src[Nj] + 1.0*this->dpsi->src[2*Nj];
  kappa->img[0] = 0.5*(ddx + ddy);
  for(int j=1;j<Nj-1;j++){
    ddx = 1.0*this->dpsi->src[j-1] - 2.0*this->dpsi->src[j] + 1.0*this->dpsi->src[j+1];
    ddy = 1.0*this->dpsi->src[j] - 2.0*this->dpsi->src[Nj+j] + 1.0*this->dpsi->src[2*Nj+j];
    kappa->img[j] = 0.5*(ddx + ddy);
  }
  ddx = 1.0*this->dpsi->src[Nj-3] - 2.0*this->dpsi->src[Nj-2] + 1.0*this->dpsi->src[Nj-1];
  ddy = 1.0*this->dpsi->src[Nj-1] - 2.0*this->dpsi->src[2*Nj-1] + 1.0*this->dpsi->src[3*Nj-1];
  kappa->img[Nj-1] = 0.5*(ddx + ddy);
  */
  for(int i=0;i<2;i++){
    for(int j=0;j<Nj;j++){
      kappa->img[i*Nj+j] = 0.0;
    }
  }

  // in-between rows
  for(int i=2;i<Ni-2;i++){
    /*
    ddx = 1.0*this->dpsi->src[i*Nj] - 2.0*this->dpsi->src[i*Nj+1] + 1.0*this->dpsi->src[i*Nj+2];
    ddy = 1.0*this->dpsi->src[(i-1)*Nj] - 2.0*this->dpsi->src[i*Nj] + 1.0*this->dpsi->src[(i+1)*Nj];
    kappa->img[i*Nj] = 0.5*(ddx + ddy);
    */
    kappa->img[i*Nj] = 0.0;
    kappa->img[i*Nj+1] = 0.0;
    for(int j=2;j<Nj-2;j++){
      ddx = (1.0*this->dpsi->src[i*Nj+j-1]   - 2.0*this->dpsi->src[i*Nj+j] + 1.0*this->dpsi->src[i*Nj+j+1]   )/dx2;
      ddy = (1.0*this->dpsi->src[(i-1)*Nj+j] - 2.0*this->dpsi->src[i*Nj+j] + 1.0*this->dpsi->src[(i+1)*Nj+j] )/dy2;
      kappa->img[i*Nj+j] = 0.5*(ddx + ddy);
    }
    /*
    ddx = 1.0*this->dpsi->src[i*Nj+Nj-3] - 2.0*this->dpsi->src[i*Nj+Nj-2] + 1.0*this->dpsi->src[i*Nj+Nj-1];
    ddy = 1.0*this->dpsi->src[(i-1)*Nj+Nj-1] - 2.0*this->dpsi->src[i*Nj+Nj-1] + 1.0*this->dpsi->src[(i+1)*Nj+Nj-1];
    kappa->img[i*Nj+Nj-1] = 0.5*(ddx + ddy);
    */
    kappa->img[i*Nj+Nj-2] = 0.0;
    kappa->img[i*Nj+Nj-1] = 0.0;
  }

  /*
  // last row
  ddx = 1.0*this->dpsi->src[(Ni-1)*Nj] - 2.0*this->dpsi->src[(Ni-1)*Nj+1] + 1.0*this->dpsi->src[(Ni-1)*Nj+2];
  ddy = 1.0*this->dpsi->src[(Ni-1)*Nj] - 2.0*this->dpsi->src[(Ni-2)*Nj] + 1.0*this->dpsi->src[(Ni-3)*Nj];
  kappa->img[(Ni-1)*Nj] = 0.5*(ddx + ddy);
  for(int j=1;j<Nj-1;j++){
    ddx = 1.0*this->dpsi->src[(Ni-1)*Nj+j-1] - 2.0*this->dpsi->src[(Ni-1)*Nj+j] + 1.0*this->dpsi->src[(Ni-1)*Nj+j+1];
    ddy = 1.0*this->dpsi->src[(Ni-1)*Nj+j] - 2.0*this->dpsi->src[(Ni-2)*Nj+j] + 1.0*this->dpsi->src[(Ni-3)*Nj+j];
    kappa->img[(Ni-1)*Nj+j] = 0.5*(ddx + ddy);
  }
  ddx = 1.0*this->dpsi->src[(Ni-1)*Nj+Nj-3] - 2.0*this->dpsi->src[(Ni-1)*Nj+Nj-2] + 1.0*this->dpsi->src[(Ni-1)*Nj+Nj-1];
  ddy = 1.0*this->dpsi->src[(Ni-1)*Nj+Nj-1] - 2.0*this->dpsi->src[(Ni-2)*Nj+Nj-1] + 1.0*this->dpsi->src[(Ni-3)*Nj+Nj-1];
  kappa->img[(Ni-1)*Nj+Nj-1] = 0.5*(ddx + ddy);
  */
  for(int i=Ni-2;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      kappa->img[i*Nj+j] = 0.0;
    }
  }


}





//Class: CollectionMassModels
//===============================================================================================================
CollectionMassModels::CollectionMassModels(){
  this->mpars["g1"] = 0.0;
  this->mpars["g2"] = 0.0;
}
CollectionMassModels::CollectionMassModels(std::vector<Nlpar*> nlpars){
  this->setPhysicalPars(nlpars);
};
CollectionMassModels::~CollectionMassModels(){
  for(int i=0;i<this->models.size();i++){
    delete(this->models[i]);
  }
  mpars.clear();
};

void CollectionMassModels::setPhysicalPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    this->mpars[nlpars[i]->nam] = nlpars[i]->val;
  }
  this->mpars["phi"] *= 0.01745329251;
  this->mpars["g1"]  = this->mpars["g"]*cos(2*this->mpars["phi"]);
  this->mpars["g2"]  = this->mpars["g"]*sin(2*this->mpars["phi"]);
}

void CollectionMassModels::all_defl(double xin,double yin,double& xout,double& yout){
  double ax   = 0.0;
  double ay   = 0.0;
  double dumx = 0.0;
  double dumy = 0.0;
  for(int i=0;i<this->models.size();i++){
    this->models[i]->defl(xin,yin,dumx,dumy);
    ax += dumx;
    ay += dumy;
  }
  xout = (1.0-this->mpars["g1"])*xin - this->mpars["g2"]*yin - ax;
  yout = (1.0+this->mpars["g1"])*yin - this->mpars["g2"]*xin - ay;
}

void CollectionMassModels::all_defl(ImagePlane* image){
  double dumx = 0.0;
  double dumy = 0.0;
  double xin,yin;
  
  for(int j=0;j<image->Nm;j++){
    xin = image->x[j];
    yin = image->y[j];
    double ax   = 0.0;
    double ay   = 0.0;
    for(int i=0;i<this->models.size();i++){
      this->models[i]->defl(xin,yin,dumx,dumy);
      ax += dumx;
      ay += dumy;
    }
    image->defl_x[j] = (1.0-this->mpars["g1"])*xin - this->mpars["g2"]*yin - ax;
    image->defl_y[j] = (1.0+this->mpars["g1"])*yin - this->mpars["g2"]*xin - ay;
  }
}
