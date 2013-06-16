#include "image.h"
//#include <mbstring.h>
#include <stdio.h>
#include <math.h>
# include "vol_math_trilateralfilter.h"
#include<iostream>
using namespace std;
image::image(int length,int height,int width)
{
	this->height=height;
	this->length=length;
	this->width=width;
	int length1=this->getlength();
	this->buf=(u_char*)malloc(length1);
}

void image::readImage(u_char* buf,char const *file ,int size)
{
	//this->height=height;
	//this->length=length;
	//this->width=width;
	//long size=size;
	FILE *op;
		op=fopen(file,"rb");
	if(op==NULL)
	{
		printf("open fail");
	}
	fread(buf,sizeof(u_char),size,op);
	fclose(op);
	printf("read is ok");
}
float* image::buf2float(u_char *buf)
{
	u_char *p;
	p=buf;
	int i=0;
	long length=this->getlength();
	//p=p[length]
	//int *imageval=new int[length];
	float *imgf=new float[length];
	while(p)
	{
		
		//*(buf+i)=*p;
		imgf[i]= (float)p[i] ;
		imgf+=i;
		p++;i++;
	}
	return imgf;
}	
//float image::checkresult()
image::~image(void)
{
	//free(buf);
	cout<<"free buf"<<endl;
}
