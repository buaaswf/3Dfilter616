#include "vol_math_trilateralfilter.h"
#include<stdlib.h>
#include<iostream>
using namespace std;
#ifndef DataType unsigned char 
#define DataType unsigned char 
class SIZE
{
public: 
	int xsize;
	int ysize;
	int zsize;
		SIZE()
		{		
			 xsize=0;
			ysize=0;
			 zsize=0;

		};
		~SIZE()
		{
			cout<<"size is deconstuctor"<<endl;
		};


};

void Grad(SIZE sz, DataType* pGray,int *pGradX, int *pGradY, int *pMag);


#endif