//=========================================================================================================
/**
 * This algorithm is called Moore Neighbor Tracing
 * An explanation of the algorithm can be found here:
 * http://www.thebigblob.com/moore-neighbor-tracing-algorithm-in-c/
 * http://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/moore.html
 *
 * @author Erik Smistad <smistad@idi.ntnu.no>
 */
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

#include "polygons.hpp"
#include "contour-tracing.hpp"
#include "vkllib.hpp"

void mooreNeighborTracing(RectGrid* image,std::vector<Contour*>& contours){
  const double EMPTY = 0;
  const double FILLED = 1;

  bool inside = false;
  int pos = 0;
  int Nx  = image->Nx;
  int Ny = image->Ny;
  
  // Need to start by padding the image by 1 pixel
  RectGrid paddedImage(Nx+2,Ny+2,image->xmin-image->step_x,image->xmax+image->step_x,image->ymin-image->step_y,image->ymax+image->step_y);
  padImage(image,&paddedImage,EMPTY);
  
  // Allocate a new image as a 1D array
  RectGrid borderImage = paddedImage;

  
  // Set entire image to EMPTY
  for(int y=0;y<(Ny+2);y++){
    for(int x=0;x<(Nx+2);x++){
      borderImage.z[x + y*(Nx+2)] = EMPTY;
    }
  }

  for(int y=0;y<(Ny+2);y++){
    for(int x=0;x<(Nx+2);x++){
      pos = x + y*(Nx+2);
      
      // Scan for FILLED pixel
      if( borderImage.z[pos] == FILLED && !inside ){	      // Entering an already discovered border
	inside = true;
      } else if( paddedImage.z[pos] == FILLED && inside ){  // Already discovered border point
	continue;
      } else if( paddedImage.z[pos] == EMPTY && inside ){   // Leaving a border
	inside = false;
      } else if( paddedImage.z[pos] == FILLED && !inside ){ // Undiscovered border point
	borderImage.z[pos] = FILLED; // Mark the start pixel
	int checkLocationNr = 1;	// The neighbor number of the location we want to check for a new border point
	int checkPosition;		// The corresponding absolute array address of checkLocationNr
	int newCheckLocationNr; 	// Variable that holds the neighborhood position we want to check if we find a new border at checkLocationNr
	int startPos = pos;		// Set start position
	int counter = 0; 		// Counter is used for the jacobi stop criterion
	int counter2 = 0; 		// Counter2 is used to determine if the point we have discovered is one single point
	
	// Defines the neighborhood offset position from current position and the neighborhood
	// position we want to check next if we find a new border at checkLocationNr
	int neighborhood[8][2] = {
	  {-1,7},
	  {-3-Nx,7},
	  {-Nx-2,1},
	  {-1-Nx,1},
	  {1,3},
	  {3+Nx,3},
	  {Nx+2,5},
	  {1+Nx,5}
	};
	// Trace around the neighborhood
	Contour* mycontour = new Contour();
	while( true ){
	  checkPosition = pos + neighborhood[checkLocationNr-1][0];
	  newCheckLocationNr = neighborhood[checkLocationNr-1][1];
	  
	  if( paddedImage.z[checkPosition] == FILLED ){ // Next border point found
	    if( checkPosition == startPos ){
	      counter ++;
	      
	      // Stopping criterion (jacob)
	      if( newCheckLocationNr == 1 || counter >= 3 ){
		// Close loop
		inside = true; // Since we are starting the search at where we first started we must set inside to true
		break;
	      }
	    }
	      
	    checkLocationNr = newCheckLocationNr;      // Update which neighborhood position we should check next
	    pos = checkPosition;
	    counter2 = 0; 			       // Reset the counter that keeps track of how many neighbors we have visited
	    borderImage.z[checkPosition] = FILLED;   // Set the border pixel
	    mycontour->x.push_back( borderImage.center_x[checkPosition] );
	    mycontour->y.push_back( borderImage.center_y[checkPosition] );
	  } else {
	    // Rotate clockwise in the neighborhood
	    checkLocationNr = 1 + (checkLocationNr % 8);
	    if( counter2 > 8 ){
	      // If counter2 is above 8 we have traced around the neighborhood and
	      // therefore the border is a single black pixel and we can exit
	      counter2 = 0; 
	      break; 
	    } else { 
	      counter2++;
	    } 
	  }
	}
	contours.push_back( mycontour );
      }
    }
  }

}

/**
 * Pads an image represented by a 1D pixel array with 1 pixel with a color
 * specified by paddingColor
 */
void padImage(RectGrid* image,RectGrid* paddedImage,double paddingColor){
  for(int x=0;x<image->Nx+2;x++){
    for(int y=0;y<image->Ny+2;y++){
      if( x == 0 || y == 0 || x == image->Nx+1 || y == image->Ny+1 ){
	paddedImage->z[x + y*(image->Nx+2)] = paddingColor;
      } else {
	paddedImage->z[x+y*(image->Nx+2)] = image->z[x-1 + (y-1)*image->Nx];
      }
    }
  }
}
//=========================================================================================================

