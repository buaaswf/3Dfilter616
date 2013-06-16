#include "vol_math_Canny.h"
#include"math.h"

void Grad(SIZE sz, DataType *pGray,int *pGradX, int *pGradY, int *pMag)
{
	int  z,y,x;

	//x����ķ�����
	for (z=1;z<sz.zsize-1;z++)
	{
		for(y=1;y<sz.ysize-1;y++)
		{
			for(x=1;x<sz.xsize-1;x++)
			{
				pGradX[y*sz.xsize +x] = (int)( pGray[y*sz.xsize+x+1]-pGray[y*sz.xsize+ x-1]  );
			}
		}
	}
	//y��������

	for(x=1;x<sz.xsize-1;x++)
	{
		for(y=1;y<sz.ysize-1;y++)
		{
			for (z=1;z<sz.zsize;z++)
			{
				pGradY[y*sz.xsize +x] = (int)(pGray[(y+1)*sz.xsize +x] - pGray[(y-1)*sz.xsize +x]);
			}
		}

		//z������



		//���ݶ�

		//�м����
		double dSqt1;
		double dSqt2;

		for(y=0; y<sz.xsize; y++)
		{
			for(x=0; x<sz.ysize; x++)
			{
				//���׷������ݶ�
				dSqt1 = pGradX[y*sz.xsize + x]*pGradX[y*sz.xsize + x];
				dSqt2 = pGradY[y*sz.xsize + x]*pGradY[y*sz.xsize + x];
				pMag[y*sz.xsize+x] = (int)(sqrt(dSqt1+dSqt2)+0.5);
			}
		}
	}
}