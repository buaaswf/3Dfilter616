#pragma once
//#include <wtypes.h>
//#include <afx.h>
#include"DetectMemLeak.h"
#include <stdlib.h>
#define BOOL int
#define TRUE 1
#define FALSE 0
class trilateralfilter
{
public:
	trilateralfilter(void);
	~trilateralfilter(void);
};

//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_

/*


*/


#define PIXTYPE unsigned char
#define M_EXP 2.7182818284590452353602874713527

//=====================
// Forward Declarations	
//=====================
class Raw;
class Raw3D; 

/**
	A 3D array of pixels, consisting of an array of PIXTYPE objects. Note
	that the Raw1D objects are not required to be the same size, though
	the sizer() function will create them as all the same size.
**/
class Raw  {
private:   			//-----------------DATA----------------- 
    	int xsize;		// # of pixels per scanline,
    	int ysize;		// # of scanlines in this Raw.
		int zsize;
    	PIXTYPE *y;		// 1D array of PIXTYPE that are accessed as a 2D array.
		
public:				//---------------init fcns-------------
    	Raw(int,int,int,PIXTYPE *);	
		Raw(void);// constructor for 'empty' Raws
    	~Raw(void);		// destructor; releases memory
   	void sizer(int ixsize, int iysize,int izsize);		// get mem for rectangle of pixels
   	void sizer(Raw* src);					// get same amt. of mem as 'src'
	int getXsize(void) {return xsize;};		// get # pixels per scanline
	int getYsize(void) {return ysize;};		// get # of scanlines.
	int getZsize(void) {return zsize;};		//get # of image numbers
    int wipecopy(Raw* src);			// copy, even with size mismatch change from bool swf 2013 4 16
   	void put(int ix, int iy,int iz, PIXTYPE val) {	// write 'val' at location ix,iy.iz.
	y[ix + xsize*iy+xsize*ysize*iz] = val; 
	};
	inline PIXTYPE get(int ix, int iy,int iz) {	// read the value at ix,iy.
		int index=ix + xsize*iy+xsize*ysize*iz;
		return(y[index]); 
	};
	PIXTYPE getXYZ(int ixyz){		// read value at 1D address ixyz
		return y[ixyz];
	};
	void putXYZ(int ixyz,PIXTYPE val){// write value at 1D address ixy
		y[ixyz] = val;
	};
				//---------------Trilateral Filter fcns-------------

	//Trilateral filter consisting of gradient filter, adaptive neighborhood
	//computation and detail filter
	void TrilateralFilter(Raw* srcImg, PIXTYPE sigmaC); 

	//Computes X , Y and Z gradients of the input image
   	void ComputeGradients(Raw* pX, Raw* pY,Raw *pZ); 

	//Bilaterally filters  gradients pX and pY 
   	void BilateralGradientFilter(Raw* pX, Raw* pY,Raw*pZ, Raw* pSmoothX, 
		Raw* pSmoothY,Raw * pSmoothZ,float sigmaC,float sigmaS, int filterSize); //sigmaC,sugmaS is change  from PIXTYPE to float

	//Builds the stack of min-max image gradients; returns the range variance
   	PIXTYPE buildMinMaxImageStack(Raw* pX, Raw* pY, Raw* pZ,Raw3D* pMinStack,
		Raw3D* pMaxStack , int levelMax, float beta);		// beta is changed from PIXTYPE to float

	//Finds the adaptive neighborhood size (stack level) 
	//from the min-max gradient stack
   	void findAdaptiveRegion(Raw3D* pMinStack, Raw3D* pMaxStack, PIXTYPE R, int levelMax); 
						
	//Filters the detail signal and computes the final output image	
   	void DetailBilateralFilter(Raw* srcImg, Raw* pSmoothX, Raw* pSmoothY, Raw* pSmoothZ, 
		Raw* fTheta, float sigmaCTheta, float sigmaRTheta); // sigmaCTheta sigmaRTheta is changed from PIXTYPE to float
    	
};



/** 
	A 3D array of pixels, consisting of a 1D array of Raw objects. 
	Use this class for images, with one Raw object for each 
	a 'field' of that image, such as R,G,B,A, etc. 
**/
class Raw3D {
public:
    	Raw *z;	// dynam. allocated space for a set of Raw objects.
    	int rawNum;	// # of Raw objects stored.
																		
public:							
    	Raw3D(void);// 'empty' Raw3D constructor.
		Raw3D(int rawNum,Raw *src);//swf add for read data 
    	~Raw3D(void);	// destructor.
    	void sizer(int ixsize, int iysize, int izsize,int irawNum); // reserve memory
    	void sizer(Raw3D* src);			// get same amt. of mem as 'src
	int getrawNum(void) {				// return # of Raw's we hold;
		return(rawNum); 
	};
	int getYsize() {					// # of Raw1D's in zval-th Raw;
		return(z[0].getYsize()); 
	};
	int getXsize(){						// # of pixels on yval,zval-th line
		return(z[0].getXsize()); 
	};
	int getZsize(){
	return (z[0].getZsize());}

    	PIXTYPE get(int ix, int iy, int iz,int rawNum) {
    		return(z[rawNum].get(ix,iy,iz));	// write 'val' at location ix,iy,iz. 
    	};
    	void put(int ix, int iy, int iz,int rawNum,  PIXTYPE val) { 
	    	z[rawNum].put(ix,iy,iz,val);		//write 'val' at location ix,iy,iz.
    	};
	void wipecopy(Raw3D& src);			// copy, resize as needed.
};
void rawarray(int xsize,int ysize,int zsize,PIXTYPE *y);