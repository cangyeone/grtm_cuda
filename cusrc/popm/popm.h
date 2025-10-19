/*
 *   popm.h
 *
 *  Created on: 2022-11-20
 *      Author: Yu Ziye 
 *  Used to define the basic parameter of wavenumber method
 */

#ifndef POPM_H_
#define POPM_H_
#include "../scilib/mathematics.h"
#include <stdlib.h>
#include <math.h>
//#include "scilib/mydebug.h"

struct layer_parameter
{
    float32 h; 
    cplx a;
    cplx b;
    float32 a0; 
    float32 b0; 
    float32 rho; 
    float32 qp; 
    float32 qs; 
    cplx u;
};
struct parameter
{
    float32 h;
    cplx a;
    cplx b;
    cplx u;
    cplx v;
    cplx x;
    cplx r;
};

struct mrtc{
    cplx tdd1[2][2]; 
    cplx tuu1[2][2]; 
    cplx rdu1[2][2]; 
    cplx rud1[2][2]; 
    cplx tdd0; 
    cplx tuu0; 
    cplx rdu0; 
    cplx rud0; 
}; 

struct grtc{
    cplx rud1[2][2]; 
    cplx rdu1[2][2]; 
    cplx tu1[2][2]; 
    cplx td1[2][2]; 
    cplx rud0; 
    cplx rdu0; 
    cplx tu0; 
    cplx td0; 
}; 

struct solution_source{
    cplx su0[2], su1[2], su2[2]; 
    cplx sd0[2], sd1[2], sd2[2]; 
    cplx shu1, shd1, shu2, shd2;  
    cplx gshu, gshd; 
    cplx gu0[2], gu1[2], gd0[2], gd1[2]; 
}; 

struct force{
    cplx f[4][4]; 
    cplx f0; 
    cplx f1; 
    cplx f2; 
    cplx f3; 
}; 

struct int_parameter{
    int n_layer; 
    int s_layer; 
    int r_layer; 
    float32 img_freq; 
    float32 zs; 
    float32 zr; 
    float32 hs; 
    float32 hr; 
    float32 dist; 
    float32 kmin; 
    float32 kmax; 
    cplx df;
    int nf1, nf2; 
    float32 Twin; 
    float32 dt; 
    float32 t0; 
    cplx dw; 
    float32 dk; 
    float32 k_factor; 
    double taper; 
    int n_ptam; 
    int src_type; 
}; 
struct popm_parameter{
    cplx yu0, yd0;
    cplx yu1[2][2], yd1[2][2]; 
    cplx a0[2]; 
    cplx a1[2]; 
    cplx a2[2]; 
    cplx b0; 
    cplx b1; 
    cplx b2; 
    cplx yu0t, yd0t;
    cplx yu1t[2][2], yd1t[2][2]; 
    cplx a0t[2]; 
    cplx a1t[2]; 
    cplx a2t[2]; 
    cplx b0t; 
    cplx b1t; 
    cplx b2t; 

}; 

struct engine_matrix{
    cplx E1[4][4];
    cplx E2[2][2]; 
}; 

struct int_k_w{
    double k;
    cplx w; 
    int nk;
    int nw;  
}; 

struct int_out{
    cplx sums[30]; 
    cplx w; 
}; 

typedef struct popm_parameter poppar; 
typedef struct parameter parameter;
typedef struct layer_parameter layer_parameter;
typedef struct mrtc mrtc; 
typedef struct grtc grtc; 
typedef struct solution_source ssrc; 
typedef struct force force; 
typedef struct int_parameter intpar; 
typedef struct engine_matrix engine_matrix; 
typedef struct int_k_w int_k_w; 
typedef struct int_out int_out; 


//void mat_PSV_E(parameter par, float32 k, cplx w, cplx E[4][4]);
//void mat_SH_E(parameter par, float32 k, cplx w, cplx E[2][2]); 
//void mat_PSV_RT(cplx Eu[4][4], cplx Ed[4][4]);
//grtc *grtm(float32, cplx, parameter *, int, int); 
//poppar ypars(float32 k, cplx w, parameter *par, int n_layer, int s_layer, float32 zs, int r_layer, float32 zr, grtc *grt); 
//ssrc moment_source(float32 k, cplx w, parameter spar, int s_layer, float32 zs, grtc *grt); 
//ssrc single_source(float32 k, cplx w, parameter spar, int s_layer, float32 zs, grtc *grt); 
//void integrate(cplx, float32, layer_parameter *par, intpar ipar, cplx sums[7]); 
extern "C" {
void green(char *, char *, char *); 
}
//float32 cal_kmax(cplx, layer_parameter *, intpar, int); 
#endif