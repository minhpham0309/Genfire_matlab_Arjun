// -- splinterp3.cpp --
// Created by AJ Pryor on 2/2/17.
//
#include <vector>
#include <complex>
#include "mex.h"
#include "splinterp.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    if (nrhs != 4) mexErrMsgTxt("Incorrect number of arguments. Syntax is Vq = splinterp3(V,Xq,Yq,Zq)");
    if (mxIsComplex(prhs[0])){
        float const *Matrix;
        float const *Matrix_r;
        float const *Matrix_i;
        float const *x;
        float const *y;
        float const *z;
        float *result_r;
        float *result_i;

        const mwSize* dims   = mxGetDimensions(prhs[0]);
        const mwSize nrows   = dims[0];
        const mwSize ncols   = dims[1];
        const mwSize nlayers = dims[2];
        if ( nrows==1 | ncols==1 | nlayers==1 )mexErrMsgTxt("Input data is not 3D.");

        const size_t ncols_out = mxGetN(prhs[1]);
        const size_t nrows_out = mxGetM(prhs[1]);
        
        
        const mwSize ndims_out  = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims_out  = mxGetDimensions(prhs[1]);
        size_t npoints = 1;
        for (auto i = 0; i < ndims_out; ++i) npoints*=dims_out[i];
        plhs[0] = mxCreateNumericArray(ndims_out, dims_out, mxSINGLE_CLASS, mxCOMPLEX);
//
        Matrix_r = (float*)mxGetPr(prhs[0]);
        Matrix_i = (float*)mxGetPi(prhs[0]);
        y        = (float*)mxGetPr(prhs[1]);
        x        = (float*)mxGetPr(prhs[2]);
        z        = (float*)mxGetPr(prhs[3]);
        result_r = (float*)mxGetPr(plhs[0]);
        result_i = (float*)mxGetPi(plhs[0]);
      
        splinterp::parallel_interp3_cx(splinterp::interp3_F_cx<float>,Matrix_r, Matrix_i, nrows, ncols, nlayers, x, y, z, npoints, result_r, result_i, 1);

    } 
    else {

        float const *Matrix;
        float const *x;
        float const *y;
        float const *z;
        float *result;
        
        const mwSize ndims   = mxGetNumberOfDimensions(prhs[0]);
        const mwSize* dims   = mxGetDimensions(prhs[0]);
        const mwSize nrows   = dims[0];
        const mwSize ncols   = dims[1];
        const mwSize nlayers = dims[2];
        if ( nrows==1 | ncols==1 | nlayers==1 )mexErrMsgTxt("Input data is not 3D.");
        
        const size_t ncols_out = mxGetN(prhs[1]);
        const size_t nrows_out = mxGetM(prhs[1]);
        
        const mwSize ndims_out  = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims_out  = mxGetDimensions(prhs[1]);
        size_t npoints = 1;
        for (auto i = 0; i < ndims_out; ++i) npoints*=dims_out[i];
        plhs[0] = mxCreateNumericArray(ndims_out, dims_out, mxSINGLE_CLASS, mxREAL);
        
        Matrix   = (float*)mxGetPr(prhs[0]);
        y        = (float*)mxGetPr(prhs[1]);
        x        = (float*)mxGetPr(prhs[2]);
        z        = (float*)mxGetPr(prhs[3]);
        result   = (float*)mxGetPr(plhs[0]);
     
        splinterp::parallel_interp3(splinterp::interp3_F<float>,Matrix, nrows, ncols, nlayers, x, y, z, npoints, result, 1);

   }
}

