#include <vector>

#include "polygons.hpp"
#include "massModels.hpp"
#include "imagePlane.hpp"


std::vector<triangle> imagePlaneToTriangles(RectGrid* image){
  std::vector<triangle> triangles;

  for(int i=0;i<image->Ny-1;i++){
    for(int j=0;j<image->Nx-1;j++){
      triangle tri;
      tri.a.x = image->center_x[i*image->Nx+j];
      tri.a.y = image->center_y[i*image->Nx+j];
      tri.b.x = image->center_x[i*image->Nx+j+1];
      tri.b.y = image->center_y[i*image->Nx+j+1];
      tri.c.x = image->center_x[(i+1)*image->Nx+j+1];
      tri.c.y = image->center_y[(i+1)*image->Nx+j+1];
      triangles.push_back(tri);

      tri.a.x = image->center_x[i*image->Nx+j];
      tri.a.y = image->center_y[i*image->Nx+j];
      tri.b.x = image->center_x[(i+1)*image->Nx+j+1];
      tri.b.y = image->center_y[(i+1)*image->Nx+j+1];
      tri.c.x = image->center_x[(i+1)*image->Nx+j];
      tri.c.y = image->center_y[(i+1)*image->Nx+j];
      triangles.push_back(tri);
    }
  }
   
  return triangles;
}

std::vector<itriangle> imagePlaneToTriangleIndices(RectGrid* image){
  std::vector<itriangle> itriangles;

  for(int i=0;i<image->Ny-1;i++){
    for(int j=0;j<image->Nx-1;j++){
      itriangle tri;
      tri.xa = j;
      tri.ya = i;
      tri.xb = j;
      tri.yb = i+1;
      tri.xc = j+1;
      tri.yc = i;
      itriangles.push_back(tri);

      tri.xa = j+1;
      tri.ya = i+1;
      tri.xb = j;
      tri.yb = i+1;
      tri.xc = j+1;
      tri.yc = i;
      itriangles.push_back(tri);
    }
  }
   
  return itriangles;
}

void deflectTriangles(const std::vector<triangle>& triangles_in,std::vector<triangle>& triangles_out,CollectionMassModels* mycollection){
  for(int i=0;i<triangles_in.size();i++){
    mycollection->all_defl(triangles_in[i].a.x,triangles_in[i].a.y,triangles_out[i].a.x,triangles_out[i].a.y);
    mycollection->all_defl(triangles_in[i].b.x,triangles_in[i].b.y,triangles_out[i].b.x,triangles_out[i].b.y);
    mycollection->all_defl(triangles_in[i].c.x,triangles_in[i].c.y,triangles_out[i].c.x,triangles_out[i].c.y);
  }
}


bool pointInTriangle(point p0,point p1,point p2,point p3){
  double termA = p2.y - p3.y;
  double termB = p1.y - p3.y;
  double termC = p3.x - p2.x;
  double termD = p1.x - p3.x;
  double termX = p0.x - p3.x;
  double termY = p0.y - p3.y;

  double den = termA*termD + termC*termB;
  
  double L1 = (termA*termX + termC*termY)/den;
  double L2 = (-termB*termX + termD*termY)/den;
  double L3 = 1.0 - L1 - L2;

  if( L1 < 0 || L2 < 0 || L3 < 0 ){
    return false;
  } else {
    return true;
  }
}



double determinant3x3(std::vector<double> row1,std::vector<double> row2,std::vector<double> row3){
  double termA =  row1[0]*(row2[1]*row3[2] - row3[1]*row2[2]);
  double termB = -row1[1]*(row2[0]*row3[2] - row3[0]*row2[2]);
  double termC =  row1[2]*(row2[0]*row3[1] - row3[0]*row2[1]);
  return termA + termB + termC;
}

void circumcircle(point A,point B,point C,double& xc,double& yc,double& r){
  std::vector<double> row1;
  std::vector<double> row2;
  std::vector<double> row3;
  double A2 = A.x*A.x + A.y*A.y;
  double B2 = B.x*B.x + B.y*B.y;
  double C2 = C.x*C.x + C.y*C.y;

  row1 = {A2,A.y,1};
  row2 = {B2,B.y,1};
  row3 = {C2,C.y,1};
  double sx = 0.5*determinant3x3(row1,row2,row3);

  row1 = {A.x,A2,1};
  row2 = {B.x,B2,1};
  row3 = {C.x,C2,1};
  double sy = 0.5*determinant3x3(row1,row2,row3);

  row1 = {A.x,A.y,1};
  row2 = {B.x,B.y,1};
  row3 = {C.x,C.y,1};
  double aa = determinant3x3(row1,row2,row3);

  row1 = {A.x,A.y,A2};
  row2 = {B.x,B.y,B2};
  row3 = {C.x,C.y,C2};
  double bb = determinant3x3(row1,row2,row3);

  xc = sx/aa;
  yc = sy/aa;
  r  = sqrt( bb/aa + (sx*sx+sy*sy)/pow(aa,2) );  
}
















int pnpoly(int nvert,double* vertx,double* verty,double testx,double testy){
  int i, j, c = 0;
  for(i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }

  return c;
}

void findBarycenter(int len,double* x,double* y,double& xc,double& yc){
  double area = 0.;
  for(int i=0;i<len-1;i++){
    area += (x[i]*y[i+1]-x[i+1]*y[i]);
  }
  area += (x[len-1]*y[0]-x[0]*y[len-1]);
  area /= 2.;

  double dum;
  xc = 0.;
  yc = 0.;
  for(int i=0;i<len-1;i++){
    dum = (x[i]*y[i+1]-x[i+1]*y[i]);
    xc += (x[i]+x[i+1])*dum;
    yc += (y[i]+y[i+1])*dum;
  }
  dum = (x[len-1]*y[0]-x[0]*y[len-1]);
  xc += (x[len-1]+x[0])*dum;
  yc += (y[len-1]+y[0])*dum;

  xc /= 6*area;
  yc /= 6*area;
}




std::vector<triangle> srcTriangles(int Ni,int Nj,double w,double h,double* x,double* y,CollectionMassModels* mycollection){

  triangle dum;
  std::vector<triangle> triangles;

  int i0 = floor(Ni/2);
  int j0 = floor(Nj/2);
  double di = w/Ni;
  double dj = h/Nj;


  double xx = 0;
  double yy = 0;
  double* defl_x = (double*) calloc(Ni*Nj,sizeof(double));
  double* defl_y = (double*) calloc(Ni*Nj,sizeof(double));

  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){

      xx =  (i-i0)*di;
      yy = -(j-j0)*dj;//reflect y-axis

      mycollection->all_defl(xx,yy,defl_x[j*Ni+i],defl_y[j*Ni+i]);
    }
  }

  for(int j=0;j<Nj-1;j++){
    for(int i=0;i<Ni-1;i++){
      dum.a.x = defl_x[j*Ni+i];
      dum.a.y = defl_y[j*Ni+i];
      dum.b.x = defl_x[j*Ni+i+1];
      dum.b.y = defl_y[j*Ni+i+1];
      dum.c.x = defl_x[(j+1)*Ni+i+1];
      dum.c.y = defl_y[(j+1)*Ni+i+1];
      triangles.push_back(dum);

      dum.b.x = defl_x[(j+1)*Ni+i+1];
      dum.b.y = defl_y[(j+1)*Ni+i+1];
      dum.c.x = defl_x[(j+1)*Ni+i];
      dum.c.y = defl_y[(j+1)*Ni+i];
      triangles.push_back(dum);
    }
  }

  return triangles;
}




std::vector<triangle> createTriangles(int Ni,int Nj,double w,double h){
  triangle dum;
  std::vector<triangle> triangles;

  int i0 = floor(Ni/2);
  int j0 = floor(Nj/2);
  double di = w/Ni;
  double dj = h/Nj;

  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){
      dum.a.x =  (i-i0)*di;
      dum.a.y = -(j-j0)*dj;//reflect y-axis
      dum.b.x =  (i+1-i0)*di;
      dum.b.y = -(j-j0)*dj;//reflect y-axis
      dum.c.x =  (i+1-i0)*di;
      dum.c.y = -(j+1-j0)*dj;//reflect y-axis
      triangles.push_back(dum);

      dum.b.x =  (i+1-i0)*di;
      dum.b.y = -(j+1-j0)*dj;//reflect y-axis
      dum.c.x =  (i-i0)*di;
      dum.c.y = -(j+1-j0)*dj;//reflect y-axis
      triangles.push_back(dum);
    }
  }

  return triangles;
}


std::vector<triangle> selectTriangles(std::vector<triangle> triangles,double rlim){
  std::vector<triangle> new_triangles;

  double xc,yc,r;
  for(int i=0;i<triangles.size();i++){
    xc = (triangles[i].a.x + triangles[i].b.x + triangles[i].c.x)/3.;
    yc = (triangles[i].a.y + triangles[i].b.y + triangles[i].c.y)/3.;

    r = sqrt(xc*xc + yc*yc);
    if( r > rlim ){
      new_triangles.push_back(triangles[i]);
    }
  }

  return new_triangles;
}


