#include "imagePlane.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <CCfits/CCfits>

#include "tableDefinition.hpp"


//ImagePlane class implementation
//============================================================================================
ImagePlane::ImagePlane(const std::string filepath,int i,int j,double h,double w){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  //Check if sizes agree
  std::valarray<float> contents(image.axis(0)*image.axis(1));
  this->readFits(filepath,contents);

  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));

  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = height/Ni;
  double dj = width/Nj;

  // xmax,xmin are the x coordinates of the leftmost and rightmost pixels.
  // REMEMBER: xmax-xmin != width (Nj*dj).
  // Similarly for y.
  xmin = -j0*dj;
  xmax = (this->Nj-j0)*dj;
  ymin = -i0*di;
  ymax = (this->Ni-i0)*di;

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      x[i*Nj+j]   =  (j-j0)*dj;
      y[i*Nj+j]   = -(i-i0)*di;//reflect y-axis
      img[i*Nj+j] = contents[i*Nj+j];
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(const std::string filepath,int i,int j,double h,double w,double x0,double y0){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  //Check if sizes agree
  std::valarray<float> contents(image.axis(0)*image.axis(1));
  this->readFits(filepath,contents);

  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));

  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = height/Ni;
  double dj = width/Nj;

  // xmax,xmin are the x coordinates of the leftmost and rightmost pixels.
  // REMEMBER: xmax-xmin != width (Nj*dj).
  // Similarly for y.
  xmin = -j0*dj + x0;
  xmax = (this->Nj-1-j0)*dj + x0;
  ymin = -i0*di + y0;
  ymax = (this->Ni-1-i0)*di + y0;

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      x[i*Nj+j]   = xmin + j*dj;
      y[i*Nj+j]   = ymax - i*di;//reflect y-axis
      img[i*Nj+j] = contents[i*Nj+j];
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(int i,int j,double h,double w){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));

  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = height/Ni;
  double dj = width/Nj;

  // xmax,xmin are the x coordinates of the leftmost and rightmost pixels.
  // REMEMBER: xmax-xmin != width (Nj*dj).
  // Similarly for y.
  xmin = -j0*dj;
  xmax = (this->Nj-1-j0)*dj;
  ymin = -i0*di;
  ymax = (this->Ni-1-i0)*di;

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      x[i*Nj+j]   =  (j-j0)*dj;
      y[i*Nj+j]   = -(i-i0)*di;//reflect y-axis
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(int i,int j,double h,double w,double x0,double y0){
  Ni     = i;
  Nj     = j;
  Nm     = Ni*Nj;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));

  int i0    = floor(Ni/2);
  int j0    = floor(Nj/2);
  double di = height/Ni;
  double dj = width/Nj;

  // xmax,xmin are the x coordinates of the leftmost and rightmost pixels.
  // REMEMBER: xmax-xmin != width (Nj*dj).
  // Similarly for y.
  xmin = -j0*dj + x0;
  xmax = (this->Nj-1-j0)*dj + x0;
  ymin = -i0*di + y0;
  ymax = (this->Ni-1-i0)*di + y0;

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      x[i*Nj+j]   = xmin + j*dj;
      y[i*Nj+j]   = ymax - i*di;//reflect y-axis
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(int i,double w,double h){
  Ni     = 0;
  Nj     = 0;
  Nm     = i;
  width  = w;
  height = h;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
}

ImagePlane::ImagePlane(const ImagePlane& image){
  Ni = image.Ni;
  Nj = image.Nj;
  Nm = image.Nm;
  width  = image.width;
  height = image.height;
  xmin = image.xmin;
  xmax = image.xmax;
  ymin = image.ymin;
  ymax = image.ymax;
  img     = (double*) calloc(Nm,sizeof(double));
  x       = (double*) calloc(Nm,sizeof(double));
  y       = (double*) calloc(Nm,sizeof(double));
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  for(int i=0;i<Nm;i++){
    img[i] = image.img[i];
    x[i] = image.x[i];
    y[i] = image.y[i];
    cells[i] = NULL;
    crosses[i] = NULL;
    dpsi_cells[i] = NULL;
  }
}

ImagePlane::~ImagePlane(){
  free(img);
  free(x);
  free(y);
  free(defl_x);
  free(defl_y);
  free(active);
  for(int i=0;i<this->Nm;i++){
    delete(cells[i]);
    delete(crosses[i]);
    delete(dpsi_cells[i]);
  }
  free(cells);
  free(crosses);
  free(dpsi_cells);
}


void ImagePlane::writeBin(const std::string filename){
  std::ofstream out(filename,std::ios::out|std::ios::binary);
  for(int i=0;i<this->Nm;i++){
    out.write((const char*) (&this->img[i]),sizeof(double));
  }
  out.close();
}

void ImagePlane::writeImage(const std::string filename){  
  long naxis    = 2;
  long naxes[2] = {(long) this->Ni,(long) this->Nj};
  long Ntot = (long) this->Nm;
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::unique_ptr<CCfits::FITS> pFits( new CCfits::FITS("!"+filename,FLOAT_IMG,naxis,naxes) );

  std::vector<long> extAx(2,(long) this->Ni);
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  //Need Ni and Nj as index counters to flip image
  std::valarray<float> array(Ntot);
  long count = 0;
  for(int j=0;j<this->Nj;j++){
    for(int i=0;i<this->Ni;i++){
      array[(this->Nj-1-j)*this->Ni+i] = (float) (this->img[count]);
      count++;
    }
  }
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,array);
  //  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().addKey("WIDTH",this->width,"width of the image in arcsec");
  pFits->pHDU().addKey("HEIGHT",this->height,"height of the image in arcsec");
  pFits->pHDU().write(fpixel,Ntot,array); 
  //  std::cout << pFits->pHDU() << std::endl;
}

void ImagePlane::readFits(const std::string filepath,std::valarray<float>& contents){
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  
  int Nj = image.axis(0);
  int Ni = image.axis(1);

  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){
      contents[j*Ni+i] = tmp[(Nj-j-1)*Ni+i];
    }
  }
}

void ImagePlane::readS(const std::string filepath){
  if( filepath == "0" ){
    
    for(int i=0;i<this->Nm;i++){
      this->S.tri.push_back({i,i,1.});
    }
    this->Nmask = this->Nm;

  } else {

    std::valarray<float> contents(this->Nm);
    this->readFits(filepath,contents);

    this->Nmask = 0.0;
    for(int i=0;i<this->Nm;i++){
      if( contents[i] == 1 ){
	this->S.tri.push_back({i,i,1.0});
	this->Nmask++;
      } else {
	this->S.tri.push_back({i,i,0.0});
      }
    }

  }
}

void ImagePlane::maskData(std::map<int,int> lookup,ImagePlane* masked){
  typedef std::map<int,int>::iterator it;
  for(it myiterator=lookup.begin();myiterator!=lookup.end();myiterator++){
    masked->img[ myiterator->second ] = this->img[ myiterator->first ];
    masked->x[ myiterator->second ]   = this->x[ myiterator->first ];
    masked->y[ myiterator->second ]   = this->y[ myiterator->first ];
  }
}

/*
void ImagePlane::setMaskedC(mytable* Cout,mytable* S,mytable* C){
  //this function does the same as multiplying algebraically tables S and C

  int* Sfull = (int*) calloc(this->Nm,sizeof(int));
  for(int k=0;k<this->S->tri.size();k++){
    Sfull[ this->S->tri[k].i ] = 1;
  }

  int i0,j0;
  double val;
  for(int k=0;k<C->tri.size();k++){
    i0  = C->tri[k].i;
    j0  = C->tri[k].j;
    val = C->tri[k].v;

    if( Sfull[i0] == 1 && Sfull[j0] == 1 ){
      Cout->tri.push_back({this->lookup[i0],this->lookup[j0],val});
    }
    
  }

  free(Sfull);
}
*/

void ImagePlane::readC(const std::string flag,const std::string filepath){
  double value;
  this->noise_flag = flag;

  if( flag == "uniform" ){

    //There is only one element in the file, repeated along the diagonal
    std::ifstream infile(filepath);
    infile >> value;
    for(int i=0;i<this->Nm;i++){
      this->C.tri.push_back({i,i,value});
    }
    infile.close();

  } else if( flag == "map" ){

    //There are exactly Ni x Nj (diagonal) elements in the file and in table C
    std::string extension = filepath.substr(filepath.find(".")+1,std::string::npos);

    if( extension == "dat" ){

      std::ifstream infile(filepath);
      int i;
      while( true ){
	infile >> i >> i >> value;
	if( infile.eof() ) break;
	this->C.tri.push_back({i,i,value});
      }
      infile.close();

    } else if( extension == "fits" ){
      
      std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
      CCfits::PHDU& image = pInfile->pHDU();
      image.readAllKeys();
      std::valarray<float> contents(image.axis(0)*image.axis(1));
      this->readFits(filepath,contents);
      
      for(int i=0;i<this->Ni*this->Nj;i++){
	this->C.tri.push_back({i,i,contents[i]});
	//	  std::cout << 1./pow(contents[i*this->Nj+j],2) << std::endl;
      }

    }


  } else if( flag == "correlated" ){

    //There are at least Ni x Nj (diagonal plus off-diagonal) elements in the file and in table C
    std::ifstream infile(filepath);
    int i,j;
    while( true ){
      infile >> i >> j >> value;
      if( infile.eof() ) break;
      this->C.tri.push_back({i,j,value});
    }
    infile.close();

  } else {

  }

}

void ImagePlane::readB(const std::string filepath,int i,int j,int ci,int cj){
  if( filepath == "0" ){

    //Create a diagonal matrix of the given dimensions.
    for(int i=0;i<this->B.Ti;i++){
      this->B.tri.push_back({i,i,1.});
    }

  } else {

    int Pi     = i;
    int Pj     = j;
    int Ncropx = ci;
    int Ncropy = cj;

    if( Pi%2 == 0 ){
      if( ci%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( ci%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }
    if( Pj%2 == 0 ){
      if( cj%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( cj%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }



    // Read PSF from file
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
    CCfits::PHDU& image = pInfile->pHDU();
    image.readAllKeys();
    //Check if sizes agree
    std::valarray<float> contents(image.axis(0)*image.axis(1));
    this->readFits(filepath,contents);

    if( Pj != image.axis(0) ){
      std::cout << std::endl << std::endl << "PROBLEM!!! PSF DIMENSIONS DON'T MATCH" << std::endl << std::endl;
    }
    if( Pi != image.axis(1) ){
      std::cout << std::endl << std::endl << "PROBLEM!!! PSF DIMENSIONS DON'T MATCH" << std::endl << std::endl;
    }



    // Report on the peak location of the PSF
    double vmax = contents[0];
    double imax = 0;
    for(int i=1;i<Pj*Pi;i++){
      if( contents[i] > vmax ){
	vmax = contents[i];
	imax = i;
      }
    }

    if( Pi%2 == 0 && Pj%2 == 0 ){
      int Pi2 = Pi/2;
      int Pj2 = Pj/2;
      int indices[] = {(Pi2-1)*Pj+Pj2-1,(Pi2-1)*Pj+Pj2,Pi2*Pj+Pj2-1,Pi2*Pj+Pj2};
      int loc[]     = {0,0,0,0};
      bool flag     =  true;
      for(int i=0;i<4;i++){
	if( indices[i] == imax ){
	  loc[i] = 1;
	  flag = false;
	}
      }
      if( flag ){
	//	std::cout << "Particular form of the PSF: it is even-even but the peak lies outside the four central pixels" << std::endl;
      } else {
	//	std::cout << std::endl << "PSF is even-even with the peak located at: " << std::endl;
	//	printf("%2d%2d\n%2d%2d\n",loc[0],loc[1],loc[2],loc[3]);
      }
    } else if( Pi%2 != 0 && Pj%2 != 0 ){
      int Pi2 = (Pi-1)/2;
      int Pj2 = (Pj-1)/2;
      int centre = Pi2*Pj + Pj2;
      if( centre != imax ){
	//	std::cout << "Particular form of the PSF: it is odd-odd but the peak is not the central pixel" << std::endl;	
      }
    } else {
      //      std::cout << "Particular form of the PSF: it it neither even-even nor odd-odd" << std::endl;
    }




    // Set cropped PSF and normalize
    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset_y = (Pi - Ncropy)/2;
    int offset_x = (Pj - Ncropx)/2;
    int offset_tot = (offset_y)*Pj + offset_x;
    double sum = 0.0;
    for(int i=0;i<Ncropy;i++){
      for (int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = contents[offset_tot+i*Pj+j];
	sum += blur[i*Ncropx+j];
      }
    }
    double fac = 1.0/sum;
    for(int i=0;i<Ncropx*Ncropy;i++){
      blur[i] *= fac;
    }    
    



    // Set quad-kernel and odd/even functions
    int Nquadx,Nquady;
    void (ImagePlane::*xfunc)(int,int,int,int,int&,int&,int&);
    void (ImagePlane::*yfunc)(int,int,int,int,int&,int&,int&);
    if( Ncropx%2 == 0){
      Nquadx = Ncropx/2;
      xfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquadx = ceil(Ncropx/2.0);
      xfunc = &ImagePlane::setCroppedLimitsOdd;
    }
    if( Ncropy%2 == 0){
      Nquady = Ncropy/2;
      yfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquady = ceil(Ncropy/2.0);
      yfunc = &ImagePlane::setCroppedLimitsOdd;
    }

    // Create the blurring matrix
    int Nleft,Nright,Ntop,Nbottom,crop_offsetx,crop_offsety,crop_offset;
    int ic,jc;
    double val;
    for(int i=0;i<this->Ni;i++){
      for(int j=0;j<this->Nj;j++){

	(this->*(xfunc))(j,Ncropx,this->Nj,Nquadx,Nleft,Nright,crop_offsetx);
	(this->*(yfunc))(i,Ncropy,this->Ni,Nquady,Ntop,Nbottom,crop_offsety);
	crop_offset = crop_offsety*Ncropx + crop_offsetx;

	for(int ii=i-Ntop;ii<i+Nbottom;ii++){
	  ic = ii - i + Ntop;
	  for(int jj=j-Nleft;jj<j+Nright;jj++){
	    jc = jj - j + Nleft;
  
	    val = blur[crop_offset + ic*Ncropx + jc];

	    this->B.tri.push_back({i*this->Nj+j,     ii*this->Nj+jj,     val });
	  }
	}

      }
    }

    free(blur);

  }
}
  
void ImagePlane::setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < Nquad ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad;
    Npost  = Nquad;
    offset = 0;
  }
}

void ImagePlane::setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < (Nquad - 1) ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - 1 - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad - 1;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad-1;
    Npost  = Nquad;
    offset = 0;
  }
}


//void ImagePlane::printCross(int k,mytable Ds){
void ImagePlane::printCross(int k){
  double coeffs[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double dsx = 1.0;//Ds.tri[2*k].v;
  double dsy = 1.0;//Ds.tri[2*k+1].v;

  coeffs[0]  = this->crosses[k]->coeff_y[0]*dsy;
  coeffs[1]  = this->crosses[k]->coeff_y[4]*dsy;
  coeffs[2]  = this->crosses[k]->coeff_x[0]*dsx;
  coeffs[3]  = this->crosses[k]->coeff_y[1]*dsy + this->crosses[k]->coeff_x[1]*dsx;
  coeffs[4]  = this->crosses[k]->coeff_y[5]*dsy + this->crosses[k]->coeff_x[2]*dsx;
  coeffs[5]  = this->crosses[k]->coeff_x[3]*dsx;
  coeffs[6]  = this->crosses[k]->coeff_x[4]*dsx;
  coeffs[7]  = this->crosses[k]->coeff_y[2]*dsy + this->crosses[k]->coeff_x[5]*dsx;
  coeffs[8]  = this->crosses[k]->coeff_y[6]*dsy + this->crosses[k]->coeff_x[6]*dsx;
  coeffs[9]  = this->crosses[k]->coeff_x[7]*dsx;
  coeffs[10] = this->crosses[k]->coeff_y[3]*dsy;
  coeffs[11] = this->crosses[k]->coeff_y[7]*dsy;

  printf("\n");
  printf("%10s %10.3f %10.3f %10s\n"," ",coeffs[0],coeffs[1]," ");
  printf("%10.3f %10.3f %10.3f %10.3f\n",coeffs[2],coeffs[3],coeffs[4],coeffs[5]);
  printf("%10.3f %10.3f %10.3f %10.3f\n",coeffs[6],coeffs[7],coeffs[8],coeffs[9]);
  printf("%10s %10.3f %10.3f %10s\n"," ",coeffs[10],coeffs[11]," ");
}

void ImagePlane::lowerResRebinAdditive(ImagePlane* newImage){
  double inf_dx = this->width/this->Nj;
  double inf_dy = this->height/this->Ni;
  double new_dx = newImage->width/newImage->Nj;
  double new_dy = newImage->height/newImage->Ni;
  for(int i=0;i<this->Ni;i++){
    int ii = (int) floor(i*inf_dy/new_dy);
    for(int j=0;j<this->Nj;j++){
      int jj = (int) floor(j*inf_dx/new_dx);
      newImage->img[ii*newImage->Nj + jj] += this->img[i*this->Nj + j];
    }
  }
}

void ImagePlane::lowerResRebinIntegrate(ImagePlane* newImage){
  // Integrating over equally sized elements is equivalent to getting the average
  int* counts = (int*) calloc(newImage->Nm,sizeof(int));
  double inf_dx = this->width/this->Nj;
  double inf_dy = this->height/this->Ni;
  double new_dx = newImage->width/newImage->Nj;
  double new_dy = newImage->height/newImage->Ni;
  for(int i=0;i<this->Ni;i++){
    int ii = (int) floor(i*inf_dy/new_dy);
    for(int j=0;j<this->Nj;j++){
      int jj = (int) floor(j*inf_dx/new_dx);
      newImage->img[ii*newImage->Nj + jj] += this->img[i*this->Nj + j];
      counts[ii*newImage->Nj + jj] += 1;
    }
  }
  for(int i=0;i<newImage->Nm;i++){
    newImage->img[i] = newImage->img[i]/counts[i];
  }
  free(counts);
}
