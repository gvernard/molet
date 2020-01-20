#ifndef FITS_HEADER_HPP
#define FITS_HEADER_HPP

#include <string>
#include <CCfits/CCfits>

void addFitsHeader(const std::string filepath,std::map<std::string,std::string> splane){
  std::string filename = filepath;


  // Read FITS file
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  std::valarray<float> contents(image.axis(0)*image.axis(1));
  image.read(contents);




  // Write FITS file
  long naxis    = 2;
  long naxes[2] = {(long) image.axis(0),(long) image.axis(1)};
  long Ntot = (long) image.axis(0)*image.axis(1);
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::unique_ptr<CCfits::FITS> pFits( new CCfits::FITS("!"+filename,FLOAT_IMG,naxis,naxes) );

  std::vector<long> extAx(2,(long) image.axis(0));
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,contents);
  pFits->pHDU().addKey("SIZE",std::stof(splane["size"]),"width of the image"); 
  pFits->pHDU().addKey("WIDTH",std::stof(splane["size"]),"width of the image"); 
  pFits->pHDU().addKey("HEIGHT",std::stof(splane["size"]),"height of the image"); 
  pFits->pHDU().write(fpixel,Ntot,contents); 
}


#endif /* FITS_HEADER_HPP */
