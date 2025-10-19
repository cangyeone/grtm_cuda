/*
 * mathematics.h
 *
 *  Created on: 2022/12/1
 *      Author: Yu Ziye
 */

#ifndef MATHEMATICS_H_
#define MATHEMATICS_H_
#include <complex.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "cucomplex.hpp"
#define EPSILON (1e-10)
#define PI (3.1415926535897932384626433)
#define PI2 (3.1415926535897932384626433*2)
#define ZPI (make_cuDoubleComplex(PI, 0.0))
#define ZPI2 (make_cuDoubleComplex(PI * 2, 0.0))
#define cuZI (make_cuDoubleComplex(0.0, 1.0))
#define cuZ1 (make_cuDoubleComplex(1.0, 0.0))
#define cuZ2 (make_cuDoubleComplex(2.0, 0.0))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define ABS(x) (((x)>0)?(x):(-(x)))
#define float32 double 
#define float64 double 
#define cprint(x) printf("(%.11lf,%.11lf)\n", real(x), imag(x)); 
#define cprint2(x) real(x), imag(x)
#define cprint3(a, x) printf("a(%.10lf,%.10lf)\n", real(x), imag(x)); 
//typedef float float32; 
//typedef double float64; 

typedef complex<double> cplx;


float32 besselj(float32, int); 
float32 d_besselj(float32, int); 
//int cplx_minv_4(cplx a[4][4], cplx c[4][4]); 
//int cplx_minvc_4(cplx a[4][4]); 
//int cplx_minvc_2(cplx a[2][2]); 
//void cplx_mm_4(cplx a[4][4], cplx b[4][4], cplx c[4][4]);
//void cplx_mm_2(cplx a[2][2], cplx b[2][2], cplx c[2][2]);
//void cplx_mmc_4(cplx a[4][4], cplx b[4][4]);
//void cplx_mmmc_4(cplx a[4][4], cplx b[4][4], cplx c[4][4]); 
#ifdef DEBUG
#define CUSOLVER_CHECK(err) (HandlecusolverError(err, __FILE__, __LINE__))
#define CUDA_CHECK(err) (HandleError(err, __FILE__, __LINE__))
#define CUBLAS_CHECK(err) (HandleBlasError(err, __FILE__, __LINE__))
#else
#define CUSOLVER_CHECK(err) (err)
#define CUDA_CHECK(err) (err)
#define CUBLAS_CHECK(err) (err)
#endif
#endif 