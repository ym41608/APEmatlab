#include "mex.h"
#include <math.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{

	/* Retrieve the input data */
    double tx_min  = mxGetScalar(prhs[0]);
    double tx_max  = mxGetScalar(prhs[1]);
    double ty_min  = mxGetScalar(prhs[2]);
    double ty_max  = mxGetScalar(prhs[3]);
    double tz_min  = mxGetScalar(prhs[4]);
    double tz_max  = mxGetScalar(prhs[5]);
    double rx_min  = mxGetScalar(prhs[6]);
    double rx_max  = mxGetScalar(prhs[7]);
    double rz_min  = mxGetScalar(prhs[8]);
    double rz_max  = mxGetScalar(prhs[9]);

	double tx_step = mxGetScalar(prhs[10]);
    double ty_step = mxGetScalar(prhs[11]);
	double tz_step = mxGetScalar(prhs[12]);
	double rx_step = mxGetScalar(prhs[13]);
	double rz0_step = mxGetScalar(prhs[14]);
	double rz1_step = mxGetScalar(prhs[15]);
    double Sxf = mxGetScalar(prhs[16]);
    double img_w = mxGetScalar(prhs[17]);
    double Syf = mxGetScalar(prhs[18]);
    double img_h = mxGetScalar(prhs[19]);
	double marker_w = mxGetScalar(prhs[20]);
	double marker_h = mxGetScalar(prhs[21]);

    double tz_mid = sqrt(tz_max*tz_min);
    
    // Pre calculate size
	double tz = tz_min;
    int numTranslate = 0;
	double length = sqrt(marker_w*marker_w + marker_h*marker_h);
    while (tz <= tz_max)
	{
        double rz1 = rz_min;
		while (rz1 <= rz_max)  // -pi pi
		{
            double rx = rx_min;
			while (rx >= -rx_max) // -pi pi
			{
                double rz0 = rz_min;
				while (rz0 <= rz_max)
				{
                    double tx = tx_min;
					while (tx <= tx_max)
					{
                        double ty = ty_min;
						while (ty < ty_max)
						{
							// inside stuff
                            if (fabs(tx) < fabs(img_w * tz / Sxf - marker_h) &&
                                fabs(ty) < fabs(img_h * tz / Syf - marker_h))// &&
                            //    fabs(rx) + fabs(ry) < 2.3)
                                numTranslate++;
                            ty += ty_step * (tz+length*sin(rx));//sqrt((tz+marker_w*sin(rx))*(tz-marker_w*sin(rx)));//
						}
						tx += tx_step * (tz+length*sin(rx));//sqrt((tz+marker_w*sin(rx))*(tz-marker_w*sin(rx)));//
					}
					rz0 += rz0_step*(tz_mid);// * sqrt((tz+0.5*sin(rx))*(tz-0.5*sin(rx)));
                    if (rx == 0)
                        rz0 = rz_max+1;
				}
                double sinValuey = 1 / (1/(2+sin(rx)) + rx_step) - 2;
                if (sinValuey <= 1 && sinValuey >= -1)
                    rx = asin(sinValuey);
                else
                    rx = -rx_max - 1;
			}
            rz1 += rz1_step*(tz_mid);// * tz;
		}
        tz += tz*tz * tz_step / (1 - tz_step*tz);
	}

    plhs[0] = mxCreateDoubleMatrix(6, numTranslate, mxREAL );  // 6 for variable for pose variable
    double *configs = mxGetPr(plhs[0]);//
    
	// MAIN LOOP
	int gridInd = 0;
    tz = tz_min;
    while (tz <= tz_max)
	{
        double rz1 = rz_min;
		while (rz1 <= rz_max)  // -pi pi
		{
            double rx = rx_min;
			while (rx >= -rx_max) // -pi pi
			{
                double rz0 = rz_min;
				while (rz0 <= rz_max)
				{
                    double tx = tx_min;
					while (tx <= tx_max)
					{
                        double ty = ty_min;
						while (ty < ty_max)
						{
							// inside stuff
                            if (fabs(tx) < fabs(img_w * tz / Sxf - marker_h) &&
                                fabs(ty) < fabs(img_h * tz / Syf - marker_h)) {
                                numTranslate++;
								configs[gridInd++] = tx;
                                configs[gridInd++] = ty;
                                configs[gridInd++] = tz;
                                configs[gridInd++] = -rx;
                                configs[gridInd++] = rz0;
                                configs[gridInd++] = rz1;
							}
                            ty += ty_step * (tz+length*sin(rx));//sqrt((tz+marker_w*sin(rx))*(tz-marker_w*sin(rx)));//
						}
						tx += tx_step * (tz+length*sin(rx));//sqrt((tz+marker_w*sin(rx))*(tz-marker_w*sin(rx)));//
					}
					rz0 += rz0_step*(tz_mid);// * sqrt((tz+0.5*sin(rx))*(tz-0.5*sin(rx)));
                    if (rx == 0)
                        rz0 = rz_max+1;
				}
                double sinValuey = 1 / (1/(2+sin(rx)) + rx_step) - 2;
                if (sinValuey <= 1 && sinValuey >= -1)
                    rx = asin(sinValuey);
                else
                    rx = -rx_max - 1;
			}
            rz1 += rz1_step*(tz_mid);// * tz;
		}
        tz += tz*tz * tz_step / (1 - tz_step*tz);
	}
	
	
	
	
}
