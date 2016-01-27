#include "mex.h"
#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{

	

	// inputs
	double *poses;
	double Sxf, Syf, x_w, y_h, marker_w, marker_h;
	int hI, wI;

	// outputs
	double *trans;
	int *insiders;

	// num. of poses
	int numPoses;
	numPoses = mxGetN(prhs[0]);

	// mxArray for outputs
	plhs[0] = mxCreateDoubleMatrix(16, numPoses, mxREAL );    // 16 for 4*4 matrix
	plhs[1] = mxCreateNumericMatrix(1, numPoses, mxINT32_CLASS, mxREAL);
	trans = mxGetPr(plhs[0]);
	insiders = (int*)mxGetPr(plhs[1]);

	// get input data
	poses = mxGetPr(prhs[0]);
	hI = double(mxGetScalar(prhs[1]));
	wI = double(mxGetScalar(prhs[2]));
    Sxf = mxGetScalar(prhs[3]);
    x_w = mxGetScalar(prhs[4]);
    Syf = mxGetScalar(prhs[5]);
    y_h = mxGetScalar(prhs[6]);
	marker_w = mxGetScalar(prhs[7]);
	marker_h = mxGetScalar(prhs[8]);


	// MAIN LOOP
    int i = 0;
	#pragma omp parallel for private(i) num_threads(8)
	for (i = 0 ; i < numPoses ; i++)
	{
        double r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,rx,rz0,rz1;
	
		tx = poses[6*i];
		ty = poses[6*i+1];
		tz = poses[6*i+2];
		rx = poses[6*i+3];
		rz0 = poses[6*i+4];
		rz1 = poses[6*i+5];
		
		//  z coordinate is y cross x   so add minus
		r11 =  cos(rz0)*cos(rz1) - sin(rz0)*cos(rx)*sin(rz1);
		r12 = -cos(rz0)*sin(rz1) - sin(rz0)*cos(rx)*cos(rz1);
		r13 = -(sin(rz0)*sin(rx));
		r21 =  sin(rz0)*cos(rz1) + cos(rz0)*cos(rx)*sin(rz1);
		r22 = -sin(rz0)*sin(rz1) + cos(rz0)*cos(rx)*cos(rz1);
		r23 = -(-cos(rz0)*sin(rx));
		r31 =  sin(rx)*sin(rz1);
		r32 =  sin(rx)*cos(rz1);
		r33 = -(cos(rx));
		
		// final transfomration
		trans[16*i]   = Sxf*r11 + x_w*r31;
		trans[16*i+1] = Sxf*r12 + x_w*r32;
		trans[16*i+2] = Sxf*r13 + x_w*r33;
		trans[16*i+3] = Sxf*tx  + x_w*tz;
		trans[16*i+4] = Syf*r21 + (y_h-1)*r31;
		trans[16*i+5] = Syf*r22 + (y_h-1)*r32;
		trans[16*i+6] = Syf*r23 + (y_h-1)*r33;
		trans[16*i+7] = Syf*ty  + (y_h-1)*tz;
		trans[16*i+8] = r31;
		trans[16*i+9] = r32;
		trans[16*i+10] = r33;
		trans[16*i+11] = tz;
		trans[16*i+12] = 0;
		trans[16*i+13] = 0;
		trans[16*i+14] = 0;
		trans[16*i+15] = 1;

		// reject transformations make marker out of boundary
		double c1x = (trans[16*i]   *(-marker_w) + trans[16*i+1]*(-marker_h) + trans[16*i+3]) / 
                     (trans[16*i+8] *(-marker_w) + trans[16*i+9]*(-marker_h) + trans[16*i+11]);
		double c1y = (trans[16*i+4] *(-marker_w) + trans[16*i+5]*(-marker_h) + trans[16*i+7]) / 
                     (trans[16*i+8] *(-marker_w) + trans[16*i+9]*(-marker_h) + trans[16*i+11]);
		double c2x = (trans[16*i]   *(+marker_w) + trans[16*i+1]*(-marker_h) + trans[16*i+3]) / 
                     (trans[16*i+8] *(+marker_w) + trans[16*i+9]*(-marker_h) + trans[16*i+11]);
		double c2y = (trans[16*i+4] *(+marker_w) + trans[16*i+5]*(-marker_h) + trans[16*i+7]) /
                     (trans[16*i+8] *(+marker_w) + trans[16*i+9]*(-marker_h) + trans[16*i+11]);
		double c3x = (trans[16*i]   *(+marker_w) + trans[16*i+1]*(+marker_h) + trans[16*i+3]) /
                     (trans[16*i+8] *(+marker_w) + trans[16*i+9]*(+marker_h) + trans[16*i+11]);
		double c3y = (trans[16*i+4] *(+marker_w) + trans[16*i+5]*(+marker_h) + trans[16*i+7]) / 
                     (trans[16*i+8] *(+marker_w) + trans[16*i+9]*(+marker_h) + trans[16*i+11]);
		double c4x = (trans[16*i]   *(-marker_w) + trans[16*i+1]*(+marker_h) + trans[16*i+3]) / 
                     (trans[16*i+8] *(-marker_w) + trans[16*i+9]*(+marker_h) + trans[16*i+11]);
		double c4y = (trans[16*i+4] *(-marker_w) + trans[16*i+5]*(+marker_h) + trans[16*i+7]) /
                     (trans[16*i+8] *(-marker_w) + trans[16*i+9]*(+marker_h) + trans[16*i+11]);
        int k = int((c1x>=0)&&(c1x<wI)&&
					(c1y>=0)&&(c1y<hI)&&
                    (c2x>=0)&&(c2x<wI)&&
					(c2y>=0)&&(c2y<hI)&&
                    (c3x>=0)&&(c3x<wI)&&
					(c3y>=0)&&(c3y<hI)&&
                    (c4x>=0)&&(c4x<wI)&&
					(c4y>=0)&&(c4y<hI));
       
        double distance1 = (c1x-c3x)*(c1x-c3x) + (c1y-c3y)*(c1y-c3y);
        double distance2 = (c2x-c4x)*(c2x-c4x) + (c2y-c4y)*(c2y-c4y);
        if (distance1 > 16 && distance2 > 16)
			insiders[i] = k;
		else
			insiders[i] = 0;
	}		
}

