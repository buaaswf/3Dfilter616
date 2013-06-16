#include <stdio.h>
#include <stdlib.h>
#include "vol_math_trilateralfilter.h"
#include "image.h"
#include <crtdbg.h> 


//using namespace std;


int main()
{
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF; 
	_CrtSetDbgFlag( tmpFlag ); 
	image *img=new image(281,481,2501);
	//image *img1=new image(281,481,2501);
	img->readImage(img->buf,"F:\\lab\\VTKproj\\mig.raw",img->getlength());
	//img1->readImage(img1->buf,"F:\\out3.vol",img1->getlength());
	//int k=0;
	//int count=0;
	//for(k=281*481*189;k<281*481*2501;k++)
	//{
	//	
	//	
	//	PIXTYPE *P=img->buf;
	//	PIXTYPE *P1=img1->buf;
	//	if(k=281*481*189)
	//	{P=P+k;P1=P1+k;}
	//	else
	//	{P=P+1;P1=P1+1;}
	//	if(*P!=127||*P1!=127)
	//	printf("*p%d%d\n",(*P),*P1);
	//	if(*P1!=*P)
	//	{
	//		count++;
	//	}
	//	//else 
	//	
	////P++;P1++;
	//}
	////float f=count/(281*481*2501);
	//printf("%d",count);
	//img->readImage(img->buf,"E:\\mig.vol",img->getlength());
	//Raw * Raw1=new Raw(281,481,img->buf2float(img->buf));
	//Raw *Raw;//=new Raw();
	//Raw=rawarray(281,481,2501,img->buf2float(img->buf));
	rawarray(281,481,2501,img->buf);
	delete img;
	//delete img1;
	//Raw[0]->TrilateralFilter(Raw[0],1.0);
	//Raw3D *raw3d=new Raw3D(2501,Raw);
	//raw3d->sizer(281,481,2501);
	//raw3d->wipecopy(raw3d);
	//Raw1->put(1,2,3.0f);
	//Raw2->put(1,1,4.0f);
	//Raw1->ComputeGradients(Raw1,Raw2);

	system("pause");

}