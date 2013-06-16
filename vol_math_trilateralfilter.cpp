#include "vol_math_trilateralfilter.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#define u_char unsigned char
using namespace std;

trilateralfilter::trilateralfilter(void)
{
}


trilateralfilter::~trilateralfilter(void)
{
}
static float lgtt=log10(2.0f);
Raw::Raw(void)
{
	xsize=0;
	ysize=0;
	zsize=0;
	y=NULL;
}

Raw::Raw(int xsize,int ysize,int zsize,PIXTYPE *y)
{
	//float i=0.0f;
	this->xsize=xsize;
	this->ysize=ysize;
	this->zsize=zsize;
	this->y=y;
}

Raw::~Raw(void)
{
	if(this->y!=NULL)
		delete [] this->y;
	y=NULL;
	//cout<<"Raw is deconstruct"<<endl;
}
/* Raw3D::Raw3D(void)
{
z=new Raw();
zsize=0;

}*/
Raw3D::Raw3D(int rawNum,Raw *src)
{
	this->rawNum=rawNum;
	this->z=src;
}
Raw3D::Raw3D(void)
{
	z=0;
}

Raw3D::~Raw3D(void){
	if(this->z!=NULL)
		delete [] this->z;
	z=NULL;
	//cout<<"Raw3D is deconstruct"<<endl;
}
//added by me
void rawarray(int xsize,int ysize,int zsize,PIXTYPE *yy)
{
	int i=0,j=0,k=0,kk=0;
	int thread=100;
	for(kk=0;kk<(zsize/thread);kk++)
	{
	//Raw *Raw[10000];//=new Raw [zsize];
	//Raw Raw[zsize];
	//PIXTYPE *p=new PIXTYPE[135161];
	try
	{

			//p+=zsize;
		PIXTYPE *p=new PIXTYPE[xsize*ysize*thread];
		//	PIXTYPE *p=&yy[0];
		for(i=0;i<xsize*ysize*thread;i++)
		{
			//printf("%d",*(yy+i));
			*(p+i)=*(yy+i+kk*(xsize*ysize*thread));
			//p=p+i;
	
		}	
			try
			{
				
				Raw *rawData=new Raw(xsize,ysize,thread,p);
				rawData->TrilateralFilter(rawData,1);
				delete rawData;
			}
			catch(std::bad_alloc)
			{
				delete [] p;
				p=NULL;
				i--;
			}
		


			//delete Raw;
			printf("no  %d\n" ,i+1);

		//for i=..


	}//try ...
	catch (const bad_alloc& e)
	{
		printf("p alloc failed");
	}


	//delete []p;
	}//for(kk...

	//return Raw[zsize];
};
//FUNCTION CALL THAT IMPLEMENTS TRILATERAL FILTER
//=====================================================================================================
void Raw::TrilateralFilter(Raw* pSrcImg, PIXTYPE sigmaC)
	//=====================================================================================================
{
	Raw destImg; 			
	Raw fTheta; 			//stores Adaptive neighborhood size for each pixel
	Raw3D minGradientStack, maxGradientStack; 
	Raw xGradient, yGradient,zGradient; 	//X and Y gradients of the input image
	Raw xSmoothGradient, ySmoothGradient,zSmoothGradient; 	//Bilaterally filtered X and Y gradients of the input image
	int levX, levY,levZ, levelMax, filterSize;  //level = log2(xsize) or log2(ysize)
	float sigmaR, sigmaCTheta, sigmaRTheta,beta;
	PIXTYPE  R;
	//domain variance for the two filters: sigmaC, sigmaCTheta
	//range variance of the two filters: sigmaR, sigmaRTheta
	//R -- threshold to compute adaptive region (see paper for details)

	//Default internal Parameters
	sigmaCTheta = sigmaC; //Variance of the Domain Filter, the only user set parameter
	beta = (float)0.15; //Always set between 0.1 and 0.2
	filterSize = (int) sigmaC; 

	//Compute the image stack height
	/*levX=(int) (log10(float (pSrcImg->getXsize()))/log10(2.0f));
	levY=(int) (log10(float (pSrcImg->getYsize()))/log10(2.0f));
	*/

	levX=(int) (log10(float (pSrcImg->getXsize()))/lgtt);
	levY=(int) (log10(float (pSrcImg->getYsize()))/lgtt);
	levZ=(int) (log10(float (pSrcImg->getZsize()))/lgtt);
	if(levX < levY)
		levelMax = levX+1;
	else
		levelMax = levY+1;

	//Allocate memory for the Min-Max Image Stack
	minGradientStack.sizer(pSrcImg->getXsize(),pSrcImg->getYsize(),pSrcImg->getZsize(),levelMax);
	maxGradientStack.sizer(pSrcImg->getXsize(),pSrcImg->getYsize(),pSrcImg->getZsize(),levelMax);

	//Allocate memory for the gradient vectors and output image
	xGradient.sizer(pSrcImg); 
	yGradient.sizer(pSrcImg);
	zGradient.sizer(pSrcImg);
	xSmoothGradient.sizer(pSrcImg);
	ySmoothGradient.sizer(pSrcImg); 
	zSmoothGradient.sizer(pSrcImg); 
	fTheta.sizer(pSrcImg);
	destImg.sizer(pSrcImg);

	/**
	Compute Gradients using Forward Difference (Step 1)
	X gradient is stored in xGradient
	Y gradient is stored in yGradient
	**/
	pSrcImg->ComputeGradients(&xGradient,&yGradient,&zGradient);

	/**
	Builds the Min-Max Image Stack consisting of Image Gradients (Step 2).
	Also computes sigmaR, the range variance, (see equation 11 in the paper for further details).
	min and max gradients are stored in two separate stacks.
	Range variance sigmaR is calculated from equation 11 in the paper.
	Large ssq improves noise reduction, but also reduces outlier
	rejection, and may blur weaker boundaries of slight intensity
	changes.
	**/
	sigmaR = pSrcImg->buildMinMaxImageStack(&xGradient,&yGradient,&zGradient,&minGradientStack,&maxGradientStack,levelMax,beta);

	//Set the remaining internal parameters required for trilateral filter
	R = sigmaR;
	sigmaRTheta =sigmaR;//from  sigmaRTheta= R = sigmaR
	/**
	Bilaterally filter the X and Y gradients of the input image (Step 3, equation 4 and 5)
	to produce xSmoothGradient and ySmoothGradient.
	**/
	pSrcImg->BilateralGradientFilter(&xGradient, &yGradient,&zGradient, &xSmoothGradient, &ySmoothGradient, &zSmoothGradient,sigmaC, sigmaR, filterSize);

	/**
	Find the adaptive neighborhood fTheta for each pixel location (Step 4). fTheta size is
	given by stack level. The min-max gradient stacks and range threshold "R" are used for this calculation
	(see equation 10 in the paper for details).
	**/
	fTheta.findAdaptiveRegion(&minGradientStack, &maxGradientStack, R, levelMax);

	/**
	Performs bilateral filter on the detail signal (Step 5).
	See Equation 6, 7, 8 and 9.
	Output is stored in destimg (end result of equation 8, Section 3.1)
	**/
	destImg.DetailBilateralFilter(pSrcImg, &xSmoothGradient, &ySmoothGradient,&zSmoothGradient,&fTheta, sigmaCTheta, sigmaRTheta);

	//Copying the result to the output image
	//wipecopy(&destImg);

	FILE *p=fopen("F:\\3D.raw","ab+");
	fwrite(destImg.y,sizeof(PIXTYPE),281*481*100,p);
	fclose(p);
	fflush(stdout);

	printf("write is ok");
	/*delete [] xGradient.y;
	delete [] yGradient.y;
	delete [] xSmoothGradient.y;
	delete [] ySmoothGradient.y;
	delete [] minGradientStack.z;
	delete [] maxGradientStack.z;
	delete [] fTheta.y;
	delete [] destImg.y;*/


	/*
	fwrite(pointer,sizeof(T),length,nfile);
	fclose(nfile);*/

	}
	




/**
Computes the X and Y gradients using forward difference
X gradient is stored in pX 
Y gradient is stored in pY
**/

void Raw::ComputeGradients(Raw* pX, Raw *pY,Raw* pZ)
{

	int i, j,k, imax, jmax,kmax, jN, iE,kU;
	PIXTYPE Cval, Eval, Nval,Uval, gE, gN,gU;

	imax = getXsize();		// get image size,
	jmax = getYsize();
	kmax = getZsize();
	for(k=0; k<kmax; k++)
	{
		kU = k+1;
		if(kU>kmax-1) kU=kmax-1;
	for(j=0; j<jmax; j++)			// for each scanline,
	{
		jN = j+1;					// addr. of North neighbor;
		if(jN>jmax-1) jN=jmax-1;
		for(i=0;i<imax;i++)			// and each pixel on a scanline,
		{
			iE = i+1;
			if(iE>imax-1) iE=imax-1;	// addr. of East neighbor
			Cval = get(i,j,k);
			Eval = get(iE,j,k);
			Nval = get(i,jN,k);
			Uval = get(i,j,kU);
			gE = (PIXTYPE) (Eval-Cval); //gradient computation with forward difference
			gN = (PIXTYPE) (Nval-Cval);
			gU = (PIXTYPE) (Uval-Cval);
			//	if(gE!=0)printf("gE");
			pX->put(i,j,k,gE);		//copy the gradient values to pX and pY
			pY->put(i,j,k,gN);
			pZ->put(i,j,k,gU);
		}// for(i...
	}// for(j...
	}//for(k...

}

/***
Building the Min-Max Stack of Image Gradients
Input:
pX and pY -- pointers to X and Y gradients of original image
levelMax -- height of the image stack
beta -- user specified parameter used to compute the range variance
Output:
pMinStack and pMaxStack -- pointers to min and max stack of image gradients
Return Value;
Range variance (sigmaR), equation 11.
***/

PIXTYPE Raw::buildMinMaxImageStack(Raw* pX, Raw* pY,Raw *pZ, Raw3D* pMinStack,
	Raw3D* pMaxStack , int levelMax, float beta)
{
	int imax, jmax,kmax, i, j,k, lev, m, n,l;
	PIXTYPE min, max, minGrad = 1000000.0, maxGrad = -1000000.0, tmp, tmpMin, tmpMax, rangeVariance;

	imax=getXsize();	//get image size
	jmax=getYsize();
	kmax=getZsize();


	for(lev=0; lev < levelMax; lev++) 
	{
		for(k=0;k<kmax;k++)
		{
			for(j=0; j < jmax; j++) 
			{
				for(i=0; i < imax; i++) 
				{
					if( lev == 0) { //stack level 0 is the magnitude of the gradients of original image
						tmp = (PIXTYPE) sqrt((float) pX->get(i,j,k)*pX->get(i,j,k) + pY->get(i,j,k)*pY->get(i,j,k) +pZ->get(i,j,k)*pZ->get(i,j,k));
						if(maxGrad < tmp)
							maxGrad = tmp;
						if(minGrad > tmp)
							minGrad = tmp;
						max = min = tmp;
						pMinStack->put(i,j,k,0,min);
						pMaxStack->put(i,j,k,0,max);
					}
					if( lev > 0) { //Gradients at each level of the image stack is computed from the level below 
						min = pMinStack->get(i,j,k,lev-1);
						max = pMaxStack->get(i,j,k,lev-1);

						for(m=-1; m <= 1; m++) 
						{
							for(n=-1; n <= 1; n++) 
							{	
								for(l=-1;l<=1;l++)
								{
									//Computes the maximum and minimum gradient value in the neighborhood
									if((i+m) >=0 && (i+m) < imax && (j+n) >=0 && (j+n) < jmax &&(k+l)>=0&&(k+l)<kmax)
									{
										tmpMin = (PIXTYPE) pMinStack->get(i+m,j+n,k+l,lev-1);
										tmpMax = (PIXTYPE) pMaxStack->get(i+m,j+n,k+l,lev-1);
										if(min > tmpMin)
											min = tmpMin;
										if(max < tmpMax)
											max = tmpMax;
									}
								}
							}
						}

						pMinStack->put(i,j,k,lev,min);
						pMaxStack->put(i,j,k,lev,max);
					} 

				}	// for(i...
			}	// for(j...
		} //for (k...
	}	// for(lev...

	//range variance is computed as a ratio of difference between max and min gradient
	rangeVariance = beta*(maxGrad - minGrad); 
	return rangeVariance;

}


/**
Bilaterally filters the X and Y gradients of the input image.
Input:
pX -- X gradient of the input image
pY -- Y gradient of the input image
sigmaC -- domain variance of the bilateral filter
sigmaR -- range variance of the bilateral filter
filterSize -- size of the filter kernel
Output:
pSmoothX -- bilaterally filtered X gradient
pSmoothY -- bilaterally filtered Y gradient

**/
void Raw::BilateralGradientFilter(Raw* pX, Raw* pY,Raw*pZ, Raw* pSmoothX, Raw* pSmoothY,Raw* pSmoothZ, 
	float sigmaC, float sigmaR, int filterSize)
{
	int i,j,k,m,n,l,imax,jmax,kmax,halfSize;
	PIXTYPE tmpX, tmpY, tmpZ,posDiff, gradDiff, domainWeight, rangeWeight,g1, g2;
	float  normFactor;
	imax = getXsize(); //get image size
	jmax = getYsize(); 
	kmax = getZsize();
	halfSize = (int) (filterSize-1)/2; //size of the filter kernel


	for(i=0; i<imax; i++) //X scanline origion i=0
	{		
		for(j=0; j<jmax; j++) //Y scaline origion i=0
		{ 
			for(k=0;k<kmax;k++)//Z scanline origion k=0
			{
				normFactor=0.0;
				tmpX=0.0;
				tmpY=0.0;
				tmpZ=0.0;
				for(m = -halfSize; m<=halfSize; m++)
				{
					for (n = -halfSize; n<=halfSize; n++)
					{

						for(l = - halfSize; l<=halfSize; l++ )
						{
							posDiff=(PIXTYPE) (m*m+n*n+l*l); 
							//Compute the weight for the domain filter (domainWeight). The domain filter
							//is a Gaussian low pass filter
							domainWeight = (PIXTYPE) pow(M_EXP, (double) (-posDiff/(2*sigmaC*sigmaC*sigmaC)));
							if( (i+m) >= 0 && (i+m) <imax && (j+n) >=0 &&(j+n) < jmax &&(k+l)>=0&&(k+l)<kmax) {
								g1 = (PIXTYPE) (pow( float (pX->get(i+m,j+n,k+l)),2.0f) + pow(float(pY->get(i+m,j+n,k+l)),2.0f) );
								g2 = (PIXTYPE) (pow( float(pX->get(i,j,k)),2.0f) + pow(float(pY->get(i,j,k)),2.0f)+pow(float(pZ->get(i,j,k)),2.0f) );
								//g3 = (PIXTYPE) (pow( float(pX->get(i,j)),2.0f) + pow(float(pY->get(i,j)),2.0f) );
								//Compute the gradient difference between a pixel and its neighborhood pixel 
								gradDiff = (PIXTYPE) (sqrt(double(g1)) - sqrt(double(g2)));
								//Compute the weight for the range filter (rangeWeight). The range filter
								//is a Gaussian filter defined by the difference in gradient magnitude.
								if(sigmaR==0){sigmaR=0.1;}
								rangeWeight = (PIXTYPE) pow(M_EXP, (double) (-(gradDiff*gradDiff)/(2*sigmaR*sigmaR*sigmaR)));	
								tmpX += pX->get(i+m,j+n,k+l)*domainWeight*rangeWeight;
								tmpY += pY->get(i+m,j+n,k+l)*domainWeight*rangeWeight;
								tmpZ += pZ->get(i+m,j+n,k+l)*domainWeight*rangeWeight;
								//Bilateral filter normalized by normFactor (eq. 5, Section 3.1) 
								normFactor += domainWeight*rangeWeight;
								//printf("nornfactor====>>%d",normFactor);
							}
						}
					}	
					if(normFactor==0)normFactor=0.1;
					tmpX = tmpX/normFactor;
					/*if(tmpX||tmpY)
					printf("not 0"); for test */
					tmpY = tmpY/normFactor;
					tmpZ = tmpZ/normFactor;
					pSmoothX->put(i,j,k,tmpX);	//copy smooth gradients to pSmoothX and pSmoothY

					/*if(pSmoothX->get(i,j)!=tmpX)
					printf("pSmoothX->get(i,j)!=tmpX");*/
					pSmoothY->put(i,j,k,tmpY);
					pSmoothZ->put(i,j,k,tmpZ);
				}  //for(k...
			}	// for(j...
		}	// for(i...
	}

}

/**
Finds adaptive neighborhood for every pixel. 
Input:
pMinStack -- pointer to the min gradient stack
pmaxStack -- pointer to the max gradient stack
R -- threshold value for computing fTheta
levelmax -- maximum level of the image stack
Output:
fTheta -- stack level that satisfies equation 10
**/

void Raw::findAdaptiveRegion(Raw3D* pMinStack, Raw3D* pMaxStack, PIXTYPE R, int levelMax)
{
	int imax, jmax,kmax,i,j,k,lev;

	imax=getXsize();	//get image size
	jmax=getYsize();
	kmax=getZsize();
	for(k=0; k<kmax; k++)
	{
	for(j=0; j<jmax; j++) {
		for(i=0; i<imax;i++) {
			for(lev=0; lev < levelMax; lev++) {
				//Compute the adaptive neighboirhood based on the similarity of
				//the neighborhood gradients, equation 10, Section 3.2.
				if ( pMaxStack->get(i,j,k,lev) > (pMaxStack->get(i,j,k,0) + R) ||
					pMinStack->get(i,j,k,lev) < (pMaxStack->get(i,j,k,0) - R) )
					break;
			}
			put(i,j,k,(PIXTYPE) (lev-1) );
		}	// for(i...
	}	// for(j...
	}//for(k...
}



/**
Filters the detail signal and computes the output (2nd filtering pass for trilateral filter).
Input:
srcImg -- input image
pSmoothX -- bilaterally filtered X gradient of srcImg
pSmoothY -- bilaterally filtered Y gradient of srcImg
fTheta -- Adaptive neighborhood for each pixel of srcImg
sigmaCTheta -- domain variance of the bilateral filter
sigmaRTheta -- range variance of the bilateral filter
Output:
Trilaterally filtered input.
**/


void Raw::DetailBilateralFilter(Raw* srcImg, Raw* pSmoothX, Raw* pSmoothY, Raw* pSmoothZ,
	Raw* fTheta, float sigmaCTheta, float sigmaRTheta)
{

	int i,j,k,m,n,l,imax,jmax,kmax,halfSize;
	int countvar=0;
	PIXTYPE tmp, diff, detail;
	float  domainWeight, rangeWeight, normFactor;
	PIXTYPE coeffA, coeffB, coeffC,coeffD; //Plane Equation is z = coeffA.x + coeffB.y + coeffC.z+coeffD
	//coeffA = dI/dx, coeffB = dI/dy, coeffC = I at center pixel of the filter kernel

	imax=getXsize();	//get image size
	jmax=getYsize();	
	kmax=getZsize();
	for(i=0; i<imax; i++) { //X scankline...
		for(j=0; j<jmax; j++) {	//Y scanline...
			for(k=0;k<kmax;k++){ //zscanline..
			normFactor=0.0;
			tmp=0.0;
			//filter window width is calculated from fTheta
			//halfsize is half of the filter window width
			halfSize=(int) fTheta->get(i,j,k); 
			halfSize = (int) (pow(2.0f,halfSize)/2);
			//halfSize=halfSize*halfSize;
			//halfSize=1.5;
			if(halfSize>2){halfSize=2;}//halfsize=5
			//Coefficients defining the centerplane (equation6, section 3.1) is calculated
			//from the smoothed image gradients
			coeffA=pSmoothX->get(i,j,k); 
			assert(coeffA==pSmoothX->get(i,j,k));
			coeffB=pSmoothY->get(i,j,k);
			coeffC=pSmoothZ->get(i,j,k);
			coeffD=srcImg->get(i,j,k);
			
			for(m = -halfSize; m<=halfSize; m++) {
				for (n = -halfSize; n<=halfSize; n++) {
					for(l= -halfSize; l<=halfSize; l++)
					diff = (PIXTYPE) (m*m+n*n+l*l);
					//Compute the weight for the domain filter (domainWeight). The domain filter
					//is a Gaussian lowpass filter
					domainWeight = (PIXTYPE) pow(M_EXP, (double) (-diff/(2*sigmaCTheta*sigmaCTheta)));		
					if( (i+m) >= 0 && (i+m) <imax && (j+n) >=0 && (j+n) < jmax && (k+l) >=0&& (k+l)< kmax) {
						//Compute the detail signal (detail) based on the difference between a 
						//neighborhood pixel and the centerplane passing through the center-pixel 
						//of the filter window. See equation 7, section 3.1 for details
						detail=(PIXTYPE) (srcImg->get(i+m,j+n,k+l) - coeffA*m - coeffB*n - coeffC*l-coeffD);	
						if(detail!=0)//printf("detail!=0%d",detail); 
						//printf("detail====>>%d ",detail);
						countvar++;
						//Compute the weight for the range filter (rangeWeight). The range filter
						//is a Gaussian filter defined by the detail signal.
						if(sigmaRTheta==0){sigmaRTheta=0.1;}// 1===>0.1
						rangeWeight = (PIXTYPE) pow(M_EXP, (double) (-(detail*detail)/(2*sigmaRTheta*sigmaRTheta)));	
						if(rangeWeight==0)rangeWeight=0.1;
						if(domainWeight==0)domainWeight=0.1;
						//tmp+=detail*domainWeight*rangeWeight;
						tmp+=detail*domainWeight*rangeWeight;
						
						//Detail Bilateral filter normalized by normFactor (eq. 9, Section 3.1)
						normFactor+= domainWeight*rangeWeight;
					}
				}
			}
			if(normFactor==0) normFactor=0.1;
			tmp = tmp/normFactor;
			if(tmp!=0)
			{
				//printf("tmp=%d\n",tmp);
				tmp+=coeffD;
				if(tmp==coeffD)
					printf("tmp==get(i,j)");
			

			}
			else tmp= coeffD;
			//if(tmp!=0){if(tmp==coeffC)printf("tmp==get(i,j)");}
			put(i,j,k,tmp);//copy to the output

		}	
		// for(j...
	}	// for(i...
	} //for(k..
//	printf("i=%d,j=%d\n",i,j);
	printf("countvar=%d\n",countvar);
}


//================================================================================
//Some helper functions for Raw and Raw3D class.
//================================================================================

/**
Memory allocator for a rectangular 2D object. Takes an 'empty' Raw
object and reserves a rectangular field of ixsize by iysize pixels.
**/
void Raw::sizer(int ixsize, int iysize,int izsize) {
	if(y!=NULL)
		delete [] this->y;
	y=NULL;
	this->y = new PIXTYPE[ixsize*iysize*izsize];	// & allocate memory.
	xsize = ixsize;				// set new image size,
	ysize = iysize;
	zsize=izsize;
}

/**
Memory allocator for a 3D object.  Allocates space for an object of
same size as 'src'.
**/
void Raw::sizer(Raw* src) {
	int ix, iy,iz;

	ix = src->getXsize();
	iy = src->getYsize();
	iz =src->getZsize();
	sizer(ix,iy,iz);
}



/**
Copy all the pixels from 'src'.  If size of 'src' is different from
us, 'wipe' out our current pixels, and change our size to match 'src'.
If a 'wipe' was necessary, return 'FALSE', else return 'TRUE'.
**/

BOOL Raw::wipecopy(Raw* src) {
	BOOL out=TRUE;
	int i,imax;	

	if(getYsize() != src->getYsize() || getXsize()!=src->getXsize()) { // resize to fit 'src'.
		sizer(src);
		out=FALSE;
	}
	imax = getXsize()*getYsize();
	for(i=0; i<imax; i++) {
		putXYZ(i,src->getXYZ(i));
	}

	return(out);
}



/**
Sets the size and allocates memory. Discards any existing contents if
the 3D object is not already 'empty'
**/
void Raw3D::sizer(int ixsize, int iysize, int izsize,int rawNum) {
	int ii;
	if(z!=NULL)
		delete[]this->z;
	z = new Raw[rawNum];			// make room for the 2D objects,
	this->rawNum = 0;   //rawNum

	for(ii=0; ii< rawNum; ii++) 
		z[ii].sizer(ixsize,iysize,izsize);	// use the Raw sizer.	
}

void Raw3D::sizer(Raw3D* src)
{
	z=src->z;
	rawNum=src->rawNum;

}
/**
Copy all the pixels from 'src', resizing as needed. Do not change
name of this object.
**/ 
void Raw3D::wipecopy(Raw3D& src) {
	int k,kmax;

	if(&src==NULL)return;
	if(src.rawNum==0) return;		// ignore empty inputs.
	if(src.getrawNum()!=rawNum || src.getZsize() != z[0].getZsize() || src.getYsize() != z[0].getYsize() ||
		src.getXsize() != getXsize()) {
			sizer(&src);
	}
	kmax = getrawNum();
	for(k=0; k< kmax; k++) {		// copy each field;
		z[k].wipecopy(&(src.z[k]));
	}
}


