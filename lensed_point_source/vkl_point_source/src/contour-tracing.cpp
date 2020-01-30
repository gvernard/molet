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
#include "imagePlane.hpp"

void mooreNeighborTracing(ImagePlane* image,std::vector<Contour*>& contours){
  const double EMPTY = 0;
  const double FILLED = 1;

  bool inside = false;
  int pos = 0;
  int width  = image->Nj;
  int height = image->Ni;
  double phys_width = image->width;
  double phys_height = image->height;
  
  // Need to start by padding the image by 1 pixel
  ImagePlane paddedImage(width+2,height+2,phys_width+2.0*phys_width/width,phys_height+2.0*phys_height/height);
  padImage(image,&paddedImage,EMPTY);
  
  // Allocate a new image as a 1D array
  ImagePlane borderImage = paddedImage;

  
  // Set entire image to EMPTY
  for(int y=0;y<(height+2);y++){
    for(int x=0;x<(width+2);x++){
      borderImage.img[x + y*(width+2)] = EMPTY;
    }
  }

  for(int y=0;y<(height+2);y++){
    for(int x=0;x<(width+2);x++){
      pos = x + y*(width+2);
      
      // Scan for FILLED pixel
      if( borderImage.img[pos] == FILLED && !inside ){	      // Entering an already discovered border
	inside = true;
      } else if( paddedImage.img[pos] == FILLED && inside ){  // Already discovered border point
	continue;
      } else if( paddedImage.img[pos] == EMPTY && inside ){   // Leaving a border
	inside = false;
      } else if( paddedImage.img[pos] == FILLED && !inside ){ // Undiscovered border point
	borderImage.img[pos] = FILLED; // Mark the start pixel
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
	  {-3-width,7},
	  {-width-2,1},
	  {-1-width,1},
	  {1,3},
	  {3+width,3},
	  {width+2,5},
	  {1+width,5}
	};
	// Trace around the neighborhood
	Contour* mycontour = new Contour();
	while( true ){
	  checkPosition = pos + neighborhood[checkLocationNr-1][0];
	  newCheckLocationNr = neighborhood[checkLocationNr-1][1];
	  
	  if( paddedImage.img[checkPosition] == FILLED ){ // Next border point found
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
	    borderImage.img[checkPosition] = FILLED;   // Set the border pixel
	    mycontour->x.push_back( borderImage.x[checkPosition] );
	    mycontour->y.push_back( borderImage.y[checkPosition] );
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
void padImage(ImagePlane* image,ImagePlane* paddedImage,double paddingColor){
  for(int x=0;x<image->Nj+2;x++){
    for(int y=0;y<image->Ni+2;y++){
      if( x == 0 || y == 0 || x == image->Nj+1 || y == image->Ni+1 ){
	paddedImage->img[x + y*(image->Nj+2)] = paddingColor;
      } else {
	paddedImage->img[x+y*(image->Nj+2)] = image->img[x-1 + (y-1)*image->Nj];
      }
    }
  }
}
//=========================================================================================================



//=========================================================================================================
// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 
// a point is defined by its coordinates {int x, y;}

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
double isLeft(point P0,point P1,point P2){
  return ( (P1.x - P0.x) * (P2.y - P0.y) - (P2.x -  P0.x) * (P1.y - P0.y) );
}


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int windingNumber(point P,Contour* C){
  int wn = 0;    // the  winding number counter
  //  std::reverse(C->x.begin(),C->x.end());
  //  std::reverse(C->y.begin(),C->y.end());

  
  // loop through all edges of the polygon
  for(int i=0;i<C->x.size()-1;i++){   // edge from V[i] to  V[i+1]
    if( C->y[i] <= P.y ){             // start y <= P.y
      if( C->y[i+1] > P.y ){          // an upward crossing
	point C1 = {C->x[i],C->y[i]};
	point C2 = {C->x[i+1],C->y[i+1]};
	if( isLeft(C1,C2,P) > 0 ){  // P left of  edge
	  ++wn;                     // have  a valid up intersect
	}
      }
    } else {                        // start y > P.y (no test needed)
      if( C->y[i+1] <= P.y ){       // a downward crossing
	point C1 = {C->x[i],C->y[i]};
	point C2 = {C->x[i+1],C->y[i+1]};
	if( isLeft(C1,C2,P) < 0 ){  // P right of  edge
	  --wn;                     // have  a valid down intersect
	}
      }
    }
  }
  return wn;
}
//=========================================================================================================
