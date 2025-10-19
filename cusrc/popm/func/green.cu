#include "../popm.h" 
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cufft.h>

    cublasOperation_t char_to_cublas_trans(char trans)
    {
        cublasOperation_t cuTrans;
        switch (trans)
        {
        case 'n':
        case 'N':
            cuTrans = CUBLAS_OP_N;
            break;
        case 't':
        case 'T':
            cuTrans = CUBLAS_OP_T;
            break;
        case 'c':
        case 'C':
            cuTrans = CUBLAS_OP_C;
            break;
        default:
            exit(-1);
        }
        return cuTrans;
    }
/*
void int_func(float32 k, cplx w, layer_parameter *laypar, intpar ipar, cplx intp[14]){
    grtc *grt; 
    poppar ppar; 
    ssrc src; 
    cplx a0[2], a1[2], a2[2], b1, b2; 
    float32 j0, j1, j2, j0d, j1d, j2d; 
    cplx func[7]; 
    parameter *par; 
    int i; 
    cplx wa, wb; 

    par = (parameter *)malloc(sizeof(parameter)*ipar.n_layer); 
    
    for(int lay=0;lay<ipar.n_layer;lay++){
        par[lay].a = laypar[lay].a;
        par[lay].b = laypar[lay].b; 
        par[lay].h = laypar[lay].h; 
        par[lay].u = laypar[lay].u; 
        wa = w/par[lay].a; 
        wb = w/par[lay].b; 
        par[lay].r = csqrt(k*k-wa*wa); 
        par[lay].v = csqrt(k*k-wb*wb); 
        if(creal(par[lay].r)<0){
            par[lay].r = -par[lay].r; 
        }
        if(creal(par[lay].v)<0){
            par[lay].v = -par[lay].v; 
        }
        par[lay].x = k*k + par[lay].v * par[lay].v; 

    }
    grt = grtm(k, w, par, ipar.n_layer, ipar.s_layer);
    ppar = ypars(k, w, par, ipar.n_layer, ipar.s_layer, ipar.zs, ipar.r_layer, ipar.zr, grt); 
    src = single_source(k, w, par[ipar.s_layer], ipar.s_layer, ipar.zs, grt); 
    if(((ipar.r_layer==ipar.s_layer)&&(ipar.zr<ipar.zs))||(ipar.r_layer<ipar.s_layer)){
        b1 = ppar.yu0 * src.gshu; 
        for(i=0;i<2;i++){
            a0[i] = ppar.yu1[i][0] * src.gu0[0] + ppar.yu1[i][1] * src.gu0[1]; 
            a1[i] = ppar.yu1[i][0] * src.gu1[0] + ppar.yu1[i][1] * src.gu1[1]; 
        }
        
    }
    else if(((ipar.r_layer==ipar.s_layer)&&(ipar.zr>ipar.zs))||(ipar.r_layer>ipar.s_layer)){
        b1 = ppar.yd0 * src.gshd; 
        for(i=0;i<2;i++){
            a0[i] = ppar.yd1[i][0] * src.gd0[0] + ppar.yd1[i][1] * src.gd0[1]; 
            a1[i] = ppar.yd1[i][0] * src.gd1[0] + ppar.yd1[i][1] * src.gd1[1]; 
        }        
    }
    j0 = besselj(ipar.dist*k, 0); 
    j1 = besselj(ipar.dist*k, 1); 
    j2 = besselj(ipar.dist*k, 2); 
    j0d = d_besselj(ipar.dist*k, 0); 
    j1d = d_besselj(ipar.dist*k, 1); 
    j2d = d_besselj(ipar.dist*k, 2); 

    func[0] = b1; 
    func[1] = a1[0] * k; 
    func[2] = a0[0] * k; 
    func[3] = a1[0]; 
    func[4] = b1 * k; 
    func[5] = a1[1] * k; 
    func[6] = a0[1] * k; 

    
    intp[0] = func[0] * j1 / ipar.dist + func[1] * j1d; 
    intp[1] = func[2] * j0d; 
    intp[2] = func[3] * j1 / ipar.dist + func[4] * j1d; 
    intp[3] = func[5] * j1; 
    intp[4] = func[6] * j0; 


    free(par); 
    free(grt); 

}
*/
float32 cal_kmax(cplx w, layer_parameter *par, intpar ipar, int nf){
    float32 pmin, dzmax, kmax; 
    float32 hr, hs;  
    hr = ipar.hr; 
    hs = ipar.hs; 
    int lay; 
    if(ipar.r_layer==(ipar.s_layer-1)){
        pmin = MAX(abs(w/par[ipar.s_layer].b), abs(w/par[ipar.s_layer-1].b)); 
        dzmax = MIN(ipar.zs, par[ipar.r_layer].h-ipar.zr); 
    }
    else if (ipar.r_layer<ipar.s_layer-1){
        pmin = MAX(abs(w/par[ipar.s_layer].b), abs(w/par[ipar.s_layer-1].b)); 
         
        dzmax = MIN(ipar.zs, hs-hr); 
        for(lay=ipar.r_layer; lay<=ipar.s_layer-2; lay++){
            if(pmin<abs(w/par[lay+1].b))pmin = abs(w/par[lay+1].b); 
            if(dzmax>par[lay].h)dzmax = par[lay].h; 
        }
    }
    else if (ipar.r_layer==ipar.s_layer)
    {
        pmin = abs(w/par[ipar.s_layer].b); 
        dzmax = ABS(hr-hs); 
        //printf("pmin:%lf,dzmax:%lf\n", pmin, dzmax); 
    }
    else if (ipar.r_layer==(ipar.s_layer+1))
    {
        pmin = MAX(abs(w/par[ipar.s_layer].b), abs(w/par[ipar.s_layer+1].b)); 
        dzmax = MIN(par[ipar.s_layer].h-ipar.zs, ipar.zr); 
    }
    else if (ipar.r_layer>ipar.s_layer+1)
    {
        pmin = MAX(abs(w/par[ipar.s_layer].b), abs(w/par[ipar.s_layer+1].b)); 

        dzmax = MIN(par[ipar.s_layer].h-ipar.zs, ABS(hr-hs)); 
        for(lay=ipar.s_layer+1; lay<ipar.r_layer-1; lay++){
            if(pmin<abs(w/par[lay].b))pmin = abs(w/par[lay].b); 
            if(dzmax>par[lay].h)dzmax = par[lay].h; 
        }
    }
    
    if (dzmax<0.2){
        if(nf<10){
            kmax = MIN(pmin+2.0, 20. * pmin); 
        }
        else{
            kmax = MIN(pmin+2.0, 2.0 * pmin); 
        }
    }
    else{
        kmax = sqrt((3.0/dzmax)*(3.0/dzmax)+pmin*pmin);
    }
    if((kmax-pmin)<0.5)kmax=pmin+0.5; 
    return kmax; 
}

__global__ void set_pars(
            int_k_w *ikw, 
            layer_parameter *laypar, 
            parameter *par){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int bid = blockIdx.x; 
    int lay = threadIdx.x;  
    cplx wa, wb; 
    cplx w; 
    cplx ZI={0.0, 1.0}; 
    double k; 
    k = ikw[bid].k; 
    w = ikw[bid].w; 
    wa = (1.0+log(w/PI2)/(PI*laypar[lay].qs)+ZI/(2.0*laypar[lay].qs)); 
    wb = (1.0+log(w/PI2)/(PI*laypar[lay].qp)+ZI/(2.0*laypar[lay].qp)); 
    par[idx].a = laypar[lay].a0 * wb; 
    par[idx].b = laypar[lay].b0 * wa; 
    par[idx].h = laypar[lay].h; 
    par[idx].u = laypar[lay].u;  

        wa = w/par[idx].a; 
        wb = w/par[idx].b; 
        par[idx].r = sqrt(k*k-wa*wa); 
        par[idx].v = sqrt(k*k-wb*wb); 
        if(real(par[idx].r)<0){
            par[idx].r = -par[idx].r; 
        }
        if(real(par[idx].v)<0){
            par[idx].v = -par[idx].v; 
        }
        par[idx].x = k*k + par[idx].v * par[idx].v; 
}
__device__ void mat_PSV_E(parameter par, float32 k, cplx w, cplx E[4][4]){
    //calculate engine matrix E of parameter 
    E[0][0] = par.a * k / w; 
    E[0][1] = par.b * par.v / w; 
    E[0][2] = E[0][0];
    E[0][3] = E[0][1]; 
    E[1][0] = par.a * par.r / w; 
    E[1][1] = par.b * k / w; 
    E[1][2] = -E[1][0]; 
    E[1][3] = -E[1][1]; 
    E[2][0] = -2.0 * par.a * par.u * k * par.r / w; 
    E[2][1] = -par.b * par.u * par.x / w; 
    E[2][2] = -E[2][0];
    E[2][3] = -E[2][1]; 
    E[3][0] = -par.a * par.u * par.x / w; 
    E[3][1] = -2.0 * par.b * par.u * k * par.v / w; 
    E[3][2] = E[3][0]; 
    E[3][3] = E[3][1];
}
__device__ void mat_SH_E(parameter par, float32 k, cplx w, cplx E[2][2]){
    //calculate engine matrix E of parameter 
    E[0][0] = 1; 
    E[0][1] = 1; 
    E[1][0] = -par.u * par.v; 
    E[1][1] = -E[1][0]; 
}
__global__ void get_emat(
            int_k_w *ikw, 
            parameter *par, 
            engine_matrix *emat){

    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    double k = ikw[blockIdx.x].k;
    cplx w = ikw[blockIdx.x].w;  
    mat_PSV_E(par[idx], k, w, emat[idx].E1); 
    mat_SH_E( par[idx], k, w, emat[idx].E2); 

}

__global__ void get_emat_e1e2(
            engine_matrix *emat, 
            cuDoubleComplex *e1, 
            cuDoubleComplex *e2){

    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int i, j; 
    cplx E1[4][4], E2[4][4]; 
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            E1[i][j]     =  emat[idx+1].E1[i][j];
            E1[i+2][j]   =  emat[idx+1].E1[i+2][j]; 
            E1[i][j+2]   = -emat[idx].E1[i][j+2]; 
            E1[i+2][j+2] = -emat[idx].E1[i+2][j+2]; 
            E2[i][j]     =  emat[idx].E1[i][j];
            E2[i+2][j]   =  emat[idx].E1[i+2][j]; 
            E2[i][j+2]   = -emat[idx+1].E1[i][j+2]; 
            E2[i+2][j+2] = -emat[idx+1].E1[i+2][j+2];   
        }
    } 

    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            e1[(idx)*16+j*4+i] = make_cuDoubleComplex(real(E1[i][j]), imag(E1[i][j])); 
            e2[(idx)*16+j*4+i] = make_cuDoubleComplex(real(E2[i][j]), imag(E2[i][j])); 
        }
    }
}

__global__ void set_e_zero( 
            int lay, 
            cuDoubleComplex *e1, 
            int n_layer){

    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    int i; 
    for(i=0;i<16;i++){
        e1[(base+lay)*16+i] = make_cuDoubleComplex(0.0, 0.0);
    }
}

__global__ void set_e_one( 
            int lay, 
            cuDoubleComplex *e1, 
            int n_layer){

    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    int i; 

    for(i=0;i<16;i++){
        e1[(base+lay)*16+i] = make_cuDoubleComplex(1.0, 0.0);
    }
}

__global__ void print_mat(
            int_k_w *ikw, 
            int lay, 
            cuDoubleComplex *e1, 
            cuDoubleComplex **e2, 
            int n_layer, 
            int sel, 
            int nnk, 
            int nnw){

    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    int i, j; 
    double x, y; 
    if(ikw[idx].nk==nnk&&ikw[idx].nw==nnw){
        if(sel==0){
            for(i=0;i<4;i++){
                for(j=0;j<4;j++){
                    y = cuCimag(e1[(idx)*16+i*4+j]); 
                    x = cuCreal(e1[(idx)*16+i*4+j]); 
                    printf("INV:%d,%d,%d,%lf,%lf\n", lay, i, j, x, y); 
                }
            }
        }
        else{
            for(i=0;i<4;i++){
                for(j=0;j<4;j++){
                    y = cuCimag(e2[idx][i*4+j]); 
                    x = cuCreal(e2[idx][i*4+j]); 
                    printf("INV:%d,%d,%d,%lf,%lf\n", lay, i, j, x, y); 
                }
            }
        }

    }

}


__global__ void setting_array(cuDoubleComplex *dA, cuDoubleComplex *dB, cuDoubleComplex *dC, 
    cuDoubleComplex **dA_B, cuDoubleComplex **dB_B, cuDoubleComplex **dC_B,
                   int m){
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    dA_B[idx] = &dA[idx*m];
    dB_B[idx] = &dB[idx*m];
    dC_B[idx] = &dC[idx*m];
}

__global__ void get_mrt_sh(
            int_k_w *ikw, 
            parameter *par, 
            mrtc *mrt){
        //last layer has similiar format. 
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    cplx wa, wb; 

    wa = par[idx].v   * par[idx].u; 
    wb = par[idx+1].v * par[idx+1].u; 
    mrt[idx].tdd0 = 2. * wa / (wa + wb); 
    mrt[idx].rdu0 = (wa-wb) / (wa + wb); 
    mrt[idx].rud0 = - mrt[idx].rdu0; 
    mrt[idx].tuu0 = 2. * wb / (wa + wb);  
}
__inline__ __device__ cplx get_cplx_from_cucplx(cuDoubleComplex a){
    cplx c={cuCreal(a), cuCimag(a)}; 
    return c; 
}

__global__ void get_mrt_psv(
            cuDoubleComplex **e1, 
            mrtc *mrt){
        //last layer has similiar format. 
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int i, j; 
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            //mrt[base+lay].tdd1[i][j] = get_cplx_from_cucplx(e1[(base)*16+(j  )*4+(i  )]); 
            //mrt[base+lay].tuu1[i][j] = get_cplx_from_cucplx(e1[(base)*16+(j+2)*4+(i+2)]); 
            //mrt[base+lay].rdu1[i][j] = get_cplx_from_cucplx(e1[(base)*16+(j  )*4+(i+2)]); 
            //mrt[base+lay].rud1[i][j] = get_cplx_from_cucplx(e1[(base)*16+(j+2)*4+(i  )]); 
            mrt[idx].tdd1[i][j] = get_cplx_from_cucplx(e1[idx][(j  )*4+(i  )]); 
            mrt[idx].tuu1[i][j] = get_cplx_from_cucplx(e1[idx][(j+2)*4+(i+2)]); 
            mrt[idx].rdu1[i][j] = get_cplx_from_cucplx(e1[idx][(j  )*4+(i+2)]); 
            mrt[idx].rud1[i][j] = get_cplx_from_cucplx(e1[idx][(j+2)*4+(i  )]); 
            //if(ikw[idx].nk==50&&ikw[idx].nw==223){
            //    printf("tdd1%d:", lay); 
            //    cprint(mrt[base+lay].tdd1[i][j]); 
            //    printf("tuu1%d:", lay); 
            //    cprint(mrt[base+lay].tuu1[i][j]); 
            //    printf("rdu1%d:", lay); 
            //    cprint(mrt[base+lay].rdu1[i][j]); 
            //    printf("rud1%d:", lay); 
            //    cprint(mrt[base+lay].rud1[i][j]); 
            //}
        }
    }
}
__device__ void mat_inv_2(cplx a[2][2]){
    cplx d, b[2][2]; 
    b[0][0] = a[0][0];
    b[0][1] = a[0][1];
    b[1][0] = a[1][0];
    b[1][1] = a[1][1];

    d = a[0][0] * a[1][1] - a[0][1] * a[1][0]; 
    a[0][0] = b[1][1]/d; 
    a[1][1] = b[0][0]/d; 
    a[0][1] = -b[0][1]/d; 
    a[1][0] = -b[1][0]/d; 
}

__device__ void cplx_mm_2(cplx a[2][2], cplx b[2][2], cplx c[2][2])
{
    int i,j,k;
    int n = 2; 
    cplx temp;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            c[i][j] = 0.0; 
        }
    }
    for(k=0;k<n;k++)
    {
        for(i=0;i<n;i++)
        {
            temp=a[i][k];
            for(j=0;j<n;j++)
            {
                c[i][j]+=temp*b[k][j];
            }
        }
    }
}
__global__ void get_grt(
    int_k_w *ikw, 
    parameter *par, 
    engine_matrix *emat, 
    mrtc *mrt, 
    grtc *grt, 
    int s_layer, 
    int n_layer
    ){
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    
    
    int i, j; 
    cplx E1[4][4]; 
    cplx a[2][2], b[2][2], c[2][2]; 
    int lay; 
    mrtc tmrt; 
    cplx exu[2], exd[2]; 
    cplx tem1, tem2, tem3; 

    int nnk = ikw[idx].nk; 
    int nnw = ikw[idx].nw; 
    cplx w = ikw[idx].w; 

    for(i = 0; i < 4; i++){
        for(j = 0; j < 4; j++){
            E1[i][j] = emat[base+0].E1[i][j];
        }
    }
    // 自由界面条件
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            a[i][j] = E1[i+2][j]; 
            b[i][j] = E1[i+2][j+2]; 
        }
    }
    mat_inv_2(a); 
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            grt[base+0].rud1[i][j] = -(a[i][0]*b[0][j] + a[i][1]*b[1][j]);
        }
    }

    grt[base+0].rud0 = 1.0; 

    //call grt coeff 
    for(lay=0;lay<s_layer;lay++){
        exu[0] = exp(-par[base+lay].r*par[base+lay].h); 
        exu[1] = exp(-par[base+lay].v*par[base+lay].h); 
        tmrt = mrt[base+lay]; 

        tem1 = grt[base+lay].rud0 * exu[1]; 
        tem2 = tmrt.rdu0 * exu[1]; 
        tem3 = tmrt.tdd0 * exu[1]; 
        
        //SH
        grt[base+lay+1].tu0 = tmrt.tuu0 / (1. - tem1 * tem2);
        grt[base+lay+1].rud0 = tmrt.rud0 + tem3 * tem1 * grt[base+lay+1].tu0; 
        

        // PSV
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                a[i][j] = -(tmrt.rdu1[i][0] * exu[0] * grt[base+lay].rud1[0][j] 
                           +tmrt.rdu1[i][1] * exu[1] * grt[base+lay].rud1[1][j]) * exu[j]; 
                b[i][j] =  (tmrt.tdd1[i][0] * exu[0] * grt[base+lay].rud1[0][j] 
                           +tmrt.tdd1[i][1] * exu[1] * grt[base+lay].rud1[1][j]) * exu[j];  
            }
            a[i][i] = 1.0 + a[i][i]; 
        }
        
        //if(ikw[idx].nk==955&&ikw[idx].nw==222){
        //    printf("tempa1%d\n", lay);   
        //    for(int ii=0;ii<4;ii++){
        //            cprint(tmrt.rdu1[(int)(ii/2)][ii%2]); 
        //    }
        //}


        mat_inv_2(a); 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[base+lay+1].tu1[i][j] = a[i][0] * tmrt.tuu1[0][j] + a[i][1] * tmrt.tuu1[1][j]; 
                
            }
        }      




        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[base+lay+1].rud1[i][j] = tmrt.rud1[i][j] + b[i][0] * grt[base+lay+1].tu1[0][j] + b[i][1] * grt[base+lay+1].tu1[1][j];
            }
        }    
    }
    // 半无限空间
    tmrt = mrt[base+n_layer-2]; 
    grt[base+n_layer-1].rdu0 = tmrt.rdu0; 
    grt[base+n_layer-1].td0  = tmrt.tdd0; 
    
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            grt[base+n_layer-1].td1[i][j] = tmrt.tdd1[i][j]; 
            grt[base+n_layer-1].rdu1[i][j] = tmrt.rdu1[i][j]; 
        }
    }        
    for(lay=n_layer-3;lay>=s_layer;lay--){
        exu[0] = exp(-par[base+lay+1].r*par[base+lay+1].h); 
        exu[1] = exp(-par[base+lay+1].v*par[base+lay+1].h); 
        tmrt = mrt[base+lay];

        tem1 = tmrt.rud0 * exu[1]; 
        tem2 = grt[base+lay+2].rdu0 * exu[1]; 
        tem3 = tmrt.tuu0 * exu[1]; 
        //if(ikw[idx].nk==955&&ikw[idx].nw==222){
        //     printf("LL=%d", lay); 
        //    cprint(tmrt.rdu0); 
        //    cprint(tmrt.tuu0); 
        //}
        
        //SH
        grt[base+lay+1].td0 = tmrt.tdd0 / (1.0 - tem1 * tem2); 
        grt[base+lay+1].rdu0 = tmrt.rdu0 + tem3 * tem2 * grt[base+lay+1].td0; 
        //if(ikw[idx].nk==955&&ikw[idx].nw==222){
        //    cprint(tem1); 
        //    cprint(tem2); 
        //    cprint(tem3); 
        //    printf("CC:%d,(%lf,%lf),(%lf,%lf)\n", lay, cprint2(grt[base+lay+1].rdu0), cprint2(grt[base+lay+1].td0)); 
        //}
        // PSV
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                a[i][j] = -(tmrt.rud1[i][0] * exu[0] * grt[base+lay+2].rdu1[0][j] 
                           +tmrt.rud1[i][1] * exu[1] * grt[base+lay+2].rdu1[1][j]) * exu[j]; 
                b[i][j] =  (tmrt.tuu1[i][0] * exu[0] * grt[base+lay+2].rdu1[0][j] 
                           +tmrt.tuu1[i][1] * exu[1] * grt[base+lay+2].rdu1[1][j]) * exu[j];  
            }
            a[i][i] = 1. + a[i][i]; 
        }


        //if(ikw[idx].nk==955&&ikw[idx].nw==222){
        //    printf("tempa2%d\n", lay);    
        //    for(int ii=0;ii<4;ii++){
        //            cprint(a[(int)(ii/2)][ii%2]); 
        //    }
        //}


        mat_inv_2(a); 
        
        
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[base+lay+1].td1[i][j] = a[i][0] * tmrt.tdd1[0][j] + a[i][1] * tmrt.tdd1[1][j]; 
            }
        }      
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[base+lay+1].rdu1[i][j] = tmrt.rdu1[i][j] + b[i][0] * grt[base+lay+1].td1[0][j] + b[i][1] * grt[base+lay+1].td1[1][j]; 
            }
        }    
        //if(nnk==116&&nnw==222){
        //    printf("GRT%d\n", lay);
        //    cprint(w)
        //    cprint(grt[base+lay+1].rdu1[0][0]); 
        //    cprint(grt[base+lay+1].rdu1[1][1]); 
        //}
    }
    //if(ikw[idx].nk==13&&ikw[idx].nw==0){
    //    printf("GRTUDU1,%d\n", base); 
    //    cprint(grt[base+3].rdu0);
    //    cprint(grt[base+2].rud0); 
    //}
    
    //if(ikw[idx].nk==955&&ikw[idx].nw==222){
    //    printf("CUDAgrt.tu1:,%d,%d,%lf,(%lf,%lf)\n", ikw[idx].nk, ikw[idx].nw, ikw[idx].k, cprint2(ikw[idx].w));  
    //    for(lay=0;lay<4;lay++){
    //        printf("lay:%d\n", lay); 
    //        for(int ii=0;ii<4;ii++){
    //                cprint(grt[base+lay].td1[(int)(ii/2)][ii%2]); 
    //                cprint(grt[base+lay].tu1[(int)(ii/2)][ii%2]); 
    //        }
    //    }
    //}
    

}

__global__ void get_ypars_a_single(
    int_k_w *ikw, 
    engine_matrix *emat,  
    parameter *par, 
    grtc *grt, 
    ssrc *src, 
    popm_parameter *ppar, 
    double zr, 
    double zs, 
    int r_layer, 
    int s_layer, 
    int n_layer 
){
    cplx sexd[2], sexu[2], rexd[2], rexu[2], rexdt[2], rexut[2];
    cplx tmp1[2][2], tmp2[2][2]; 
    cplx tmp1t[2][2], tmp2t[2][2]; 
    cplx iden[2][2]; 
    cplx E[4][4]; 
    cplx u0t; 
    cplx d0t; 
    int i, j, lay; 
    cplx yu_psv[2][2], yd_psv[2][2], y_sh, yu0, yd0;
    cplx yu_psvt[2][2], yd_psvt[2][2], y_sht, yu0t, yd0t=0.0;
    cplx u[2][2], d[2][2], u0, d0, yu1[2][2], yd1[2][2]; 
    cplx ut[2][2], dt[2][2], ut0, dt0, yu1t[2][2], yd1t[2][2]; 
    cplx ex[2]; 
    cplx a0t[2], a1t[2], a2t[2], b1t, b2t; 
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    int nnk = ikw[idx].nk; 
    int nnw = ikw[idx].nw; 
    cplx w = ikw[idx].w; 
    //get grt matrix, need to be freed 
    //grt = grtm(k, w, par, n_layer, s_layer); 
    //Identify matrix
    iden[0][0] = iden[1][1] = 1.0; 
    iden[1][0] = iden[0][1] = 0.0; 

    //source exp pars. 
    sexd[0] = exp(-par[base+s_layer].r * (par[base+s_layer].h - zs)); 
    sexd[1] = exp(-par[base+s_layer].v * (par[base+s_layer].h - zs)); 
    sexu[0] = exp(-par[base+s_layer].r * (zs)); 
    sexu[1] = exp(-par[base+s_layer].v * (zs)); 

    
    rexd[0] = exp(-par[base+r_layer].r * zr); 
    rexd[1] = exp(-par[base+r_layer].v * zr); 

    rexdt[0] = - par[base+r_layer].r * exp(-par[base+r_layer].r * zr); 
    rexdt[1] = - par[base+r_layer].v * exp(-par[base+r_layer].v * zr); 
    if ((n_layer-1) == r_layer){
        rexu[0] = 0.0; 
        rexu[1] = 0.0; 
        rexut[0] = 0.0; 
        rexut[1] = 0.0; 
    }
    else{
        rexu[0] = exp(-par[base+r_layer].r * (par[base+r_layer].h - zr)); 
        rexu[1] = exp(-par[base+r_layer].v * (par[base+r_layer].h - zr));      
        rexut[0] = exp(-par[base+r_layer].r * (par[base+r_layer].h - zr)) * par[r_layer].r; 
        rexut[1] = exp(-par[base+r_layer].v * (par[base+r_layer].h - zr)) * par[r_layer].v;        
    }

    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            tmp1[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j];
            tmp2[i][j] = sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 

        }
    }
    //A12-cd:ERER 
    cplx_mm_2(tmp1, tmp2, yu_psv); 
    cplx_mm_2(tmp2, tmp1, yd_psv); 


    //A12-cd:1-ERER
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            yu_psv[i][j] = iden[i][j] - yu_psv[i][j]; 
            yd_psv[i][j] = iden[i][j] - yd_psv[i][j]; 
        }
    }
    //A12-cd:(1-ERER)^-1
    mat_inv_2(yu_psv); 
    mat_inv_2(yd_psv); 


    y_sh = 1./(1.-sexd[1]*grt[base+s_layer+1].rdu0*sexd[1]*sexu[1]*grt[base+s_layer].rud0*sexu[1]);  

    for(i = 0; i < 4; i++){
        for(j = 0; j < 4; j++){
            E[i][j] = emat[base+r_layer].E1[i][j];
        }
    }
    if((s_layer==r_layer)&&(zr<=zs)){
        rexd[0] = exp(-par[base+r_layer].r * (zr)); 
        rexd[1] = exp(-par[base+r_layer].v * (zr)); 
        rexu[0] = exp(-par[base+r_layer].r * (zs-zr)); 
        rexu[1] = exp(-par[base+r_layer].v * (zs-zr));     
        rexdt[0] = - par[base+r_layer].r * exp(-par[base+r_layer].r * (zr)); 
        rexdt[1] = - par[base+r_layer].v * exp(-par[base+r_layer].v * (zr)); 
        rexut[0] = par[base+r_layer].r * exp(-par[base+r_layer].r * (zs-zr)); 
        rexut[1] = par[base+r_layer].v * exp(-par[base+r_layer].v * (zs-zr));           
        
        u0 = rexd[1] * grt[base+s_layer-1+1].rud0 * sexu[1] + rexu[1]; 
        u0t = rexdt[1] * grt[base+s_layer-1+1].rud0 * sexu[1] + rexut[1] * par[base+r_layer].v;

        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[base+s_layer-1+1].rud1[0][j] * sexu[j] 
                        + E[i][1] * rexd[1] * grt[base+s_layer-1+1].rud1[1][j] * sexu[j] 
                        + E[i][j+2] * rexu[j];
                ut[i][j] = E[i][0] * rexdt[0] * grt[base+s_layer-1+1].rud1[0][j] * sexu[j] 
                        + E[i][1] * rexdt[1] * grt[base+s_layer-1+1].rud1[1][j] * sexu[j] 
                        + E[i][j+2] * rexut[j];  
            }
        }

        yu0 = u0 * y_sh; 
        yu0t = u0t * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = u[i][0] * yu_psv[0][j] + u[i][1] * yu_psv[1][j]; 
                yu1t[i][j] = ut[i][0] * yu_psv[0][j] + ut[i][1] * yu_psv[1][j]; 
            }
        }


    }
    else if ((s_layer==r_layer)&&(zr>zs))
    {
        rexd[0] = exp(-par[base+s_layer].r * (zr-zs)); 
        rexd[1] = exp(-par[base+s_layer].v * (zr-zs)); 
        rexu[0] = exp(-par[base+s_layer].r * (par[base+s_layer].h-zr)); 
        rexu[1] = exp(-par[base+s_layer].v * (par[base+s_layer].h-zr)); 
        rexdt[0] = - par[base+r_layer].r * rexd[0]; 
        rexdt[1] = - par[base+r_layer].v * rexd[1]; 
        rexut[0] = par[base+r_layer].r   * rexu[0]; 
        rexut[1] = par[base+r_layer].v   * rexu[1];  
        d0 = rexd[1] + rexu[1] * grt[base+s_layer+1].rdu0 * sexd[1]; 
        d0t = rexdt[1] + rexut[1] * grt[base+s_layer+1].rdu0 * sexd[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[base+s_layer+1].rdu1[0][j] * sexd[j] 
                        + E[i][3] * rexu[1] * grt[base+s_layer+1].rdu1[1][j] * sexd[j]; 
                dt[i][j] = E[i][j] * rexdt[j]
                        + E[i][2] * rexut[0] * grt[base+s_layer+1].rdu1[0][j] * sexd[j] 
                        + E[i][3] * rexut[1] * grt[base+s_layer+1].rdu1[1][j] * sexd[j]; 
            }
        }
        yd0 = d0 * y_sh; 
        yd0t = d0t * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = d[i][0] * yd_psv[0][j] + d[i][1] * yd_psv[1][j]; 
                yd1t[i][j] = dt[i][0] * yd_psv[0][j] + dt[i][1] * yd_psv[1][j]; 
            }
        }         
    }
    else if(r_layer<s_layer){
        ex[0] = exp(-par[base+r_layer].r*par[base+r_layer].h); 
        ex[1] = exp(-par[base+r_layer].v*par[base+r_layer].h); 
        u0 = rexd[1] * grt[base+r_layer].rud0 * ex[1] + rexu[1];
        u0t = rexdt[1] * grt[base+r_layer].rud0 * ex[1] + rexut[1];  
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[base+r_layer-1+1].rud1[0][j] * ex[j] 
                        + E[i][1] * rexd[1] * grt[base+r_layer-1+1].rud1[1][j] * ex[j] 
                        + E[i][j+2] * rexu[j]; 
                ut[i][j] = E[i][0] * rexdt[0] * grt[base+r_layer-1+1].rud1[0][j] * ex[j] 
                        + E[i][1] * rexdt[1] * grt[base+r_layer-1+1].rud1[1][j] * ex[j] 
                        + E[i][j+2] * rexut[j]; 
            }
        }
        


        for(lay=s_layer-1;lay>=r_layer;lay--){
            if(lay==(s_layer-1)){
                yu0 = grt[base+lay+1].tu0 * sexu[1] * y_sh; 
                //if(ikw[idx].nk==13&&ikw[idx].nw==0){
                //    printf("CC:\n"); 
                //    cprint(grt[base+lay+1].tu0);
                //    cprint(sexu[1]); 
                //    cprint(y_sh); 
                //    cprint(yu0)
                //}
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yu1[i][j] = grt[base+lay+1].tu1[i][0] * sexu[0] * yu_psv[0][j]
                                  + grt[base+lay+1].tu1[i][1] * sexu[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[base+lay+1].r * par[base+lay+1].h); 
                ex[1] = exp(-par[base+lay+1].v * par[base+lay+1].h); 


                yu0 = grt[base+lay+1].tu0 * ex[1] * yu0; 
                //if(ikw[idx].nk==13&&ikw[idx].nw==0){
                //    printf("DD:\n"); 
                //    cprint(grt[base+lay+1].tu0);
                //    cprint(yu0); 
                //    cprint(ex[1])
                //}
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[base+lay+1].tu1[i][0] * ex[0] * yu1[0][j]
                                   + grt[base+lay+1].tu1[i][1] * ex[1] * yu1[1][j]; 
                    }
                }




                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yu1[i][j] = tmp1[i][j]; 
                    }
                }              

            }
        }
        
        yu0 = u0 * yu0; 
        yu0t = u0t * yu0;  
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = u[i][0] * yu1[0][j] + u[i][1] * yu1[1][j]; 
                tmp1t[i][j] = ut[i][0] * yu1[0][j] + ut[i][1] * yu1[1][j];
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = tmp1[i][j]; 
                yu1t[i][j] = tmp1t[i][j]; 
            }
        }                    


    }
    else if (r_layer>s_layer)
    {
        ex[0] = exp(-par[base+r_layer-1].r*par[base+r_layer-1].h); 
        ex[1] = exp(-par[base+r_layer-1].v*par[base+r_layer-1].h); 
        d0 = rexu[1] * grt[base+r_layer+1].rdu0 * ex[1] + rexd[1]; 
        d0t = rexut[1] * grt[base+r_layer+1].rdu0 * ex[1] + rexdt[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[base+r_layer+1].rdu1[0][j] * ex[j] 
                        + E[i][3] * rexu[1] * grt[base+r_layer+1].rdu1[1][j] * ex[j]; 
                dt[i][j] = E[i][j] * rexdt[j]
                        + E[i][2] * rexut[0] * grt[base+r_layer+1].rdu1[0][j] * ex[j] 
                        + E[i][3] * rexut[1] * grt[base+r_layer+1].rdu1[1][j] * ex[j]; 
            }
        }
        for(lay=s_layer+1;lay<=r_layer;lay++){
            if(lay==(s_layer+1)){
                yd0 = grt[base+lay+1].td0 * sexd[1] * y_sh; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yd1[i][j] = grt[base+lay+1].td1[i][0] * sexd[0] * yu_psv[0][j]
                                  + grt[base+lay+1].td1[i][1] * sexd[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[base+lay-2].r * par[base+lay-2].h); 
                ex[1] = exp(-par[base+lay-2].v * par[base+lay-2].h); 
                yd0 = grt[base+lay+1].td0 * ex[1] * yd0; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[base+lay+1].td1[i][0] * ex[0] * yd1[0][j]
                                   + grt[base+lay+1].td1[i][1] * ex[1] * yd1[1][j]; 
                    }
                }
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yd1[i][j] = tmp1[i][j]; 
                    }
                }                                
            }
        }
        yd0 = d0 * yd0; 
        yd0t = d0 * yd0; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = d[i][0] * yd1[0][j] + d[i][1] * yd1[1][j]; 
                tmp1t[i][j] = dt[i][0] * yd1[0][j] + dt[i][1] * yd1[1][j]; 
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = tmp1[i][j]; 
                yd1t[i][j] = tmp1t[i][j]; 
            }
        }      
    }


    if(((r_layer==s_layer)&&(zr<zs))||(r_layer<s_layer)){
        ppar[idx].b1 = yu0 * src[idx].gshu; 
        ppar[idx].b1t = yu0t * src[idx].gshu; 
        for(i=0;i<2;i++){
            ppar[idx].a0[i] = yu1[i][0] * src[idx].gu0[0] + yu1[i][1] * src[idx].gu0[1]; 
            ppar[idx].a1[i] = yu1[i][0] * src[idx].gu1[0] + yu1[i][1] * src[idx].gu1[1]; 
            ppar[idx].a0t[i] = yu1t[i][0] * src[idx].gu0[0] + yu1t[i][1] * src[idx].gu0[1]; 
            ppar[idx].a1t[i] = yu1t[i][0] * src[idx].gu1[0] + yu1t[i][1] * src[idx].gu1[1];             
        }
        
    }
    else if(((r_layer==s_layer)&&(zr>zs))||(r_layer>s_layer)){
        ppar[idx].b1 = yd0 * src[idx].gshd; 
        ppar[idx].b1t = yd0t * src[idx].gshd; 
        for(i=0;i<2;i++){
            ppar[idx].a0[i] = yd1[i][0] * src[idx].gd0[0] + yd1[i][1] * src[idx].gd0[1]; 
            ppar[idx].a1[i] = yd1[i][0] * src[idx].gd1[0] + yd1[i][1] * src[idx].gd1[1]; 
            ppar[idx].a0t[i] = yd1t[i][0] * src[idx].gd0[0] + yd1t[i][1] * src[idx].gd0[1]; 
            ppar[idx].a1t[i] = yd1t[i][0] * src[idx].gd1[0] + yd1t[i][1] * src[idx].gd1[1]; 
        }        
    }
}


__global__ void get_ypars_a_tensor(
    int_k_w *ikw, 
    engine_matrix *emat,  
    parameter *par, 
    grtc *grt, 
    ssrc *src, 
    popm_parameter *ppar, 
    double zr, 
    double zs, 
    int r_layer, 
    int s_layer, 
    int n_layer 
){
    cplx sexd[2], sexu[2], rexd[2], rexu[2], rexdt[2], rexut[2];
    cplx tmp1[2][2], tmp2[2][2]; 
    cplx tmp1t[2][2], tmp2t[2][2]; 
    cplx iden[2][2]; 
    cplx E[4][4]; 
    cplx u0t; 
    cplx d0t; 
    int i, j, lay; 
    cplx yu_psv[2][2], yd_psv[2][2], y_sh, yu0, yd0;
    cplx yu_psvt[2][2], yd_psvt[2][2], y_sht, yu0t, yd0t=0.0;
    cplx u[2][2], d[2][2], u0, d0, yu1[2][2], yd1[2][2]; 
    cplx ut[2][2], dt[2][2], ut0, dt0, yu1t[2][2], yd1t[2][2]; 
    cplx ex[2]; 
    cplx a0t[2], a1t[2], a2t[2], b1t, b2t; 
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    
    //get grt matrix, need to be freed 
    //grt = grtm(k, w, par, n_layer, s_layer); 
    //Identify matrix
    iden[0][0] = iden[1][1] = 1.0; 
    iden[1][0] = iden[0][1] = 0.0; 

    //source exp pars. 
    sexd[0] = exp(-par[base+s_layer].r * (par[base+s_layer].h - zs)); 
    sexd[1] = exp(-par[base+s_layer].v * (par[base+s_layer].h - zs)); 
    sexu[0] = exp(-par[base+s_layer].r * (zs)); 
    sexu[1] = exp(-par[base+s_layer].v * (zs)); 

    
    rexd[0] = exp(-par[base+r_layer].r * zr); 
    rexd[1] = exp(-par[base+r_layer].v * zr); 

    rexdt[0] = - par[base+r_layer].r * exp(-par[base+r_layer].r * zr); 
    rexdt[1] = - par[base+r_layer].v * exp(-par[base+r_layer].v * zr); 
    if ((n_layer-1) == r_layer){
        rexu[0] = 0.0; 
        rexu[1] = 0.0; 
        rexut[0] = 0.0; 
        rexut[1] = 0.0; 
    }
    else{
        rexu[0] = exp(-par[base+r_layer].r * (par[base+r_layer].h - zr)); 
        rexu[1] = exp(-par[base+r_layer].v * (par[base+r_layer].h - zr));      
        rexut[0] = exp(-par[base+r_layer].r * (par[base+r_layer].h - zr)) * par[r_layer].r; 
        rexut[1] = exp(-par[base+r_layer].v * (par[base+r_layer].h - zr)) * par[r_layer].v;        
    }

    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            tmp1[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j];
            tmp2[i][j] = sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 

        }
    }
    //A12-cd:ERER 
    cplx_mm_2(tmp1, tmp2, yu_psv); 
    cplx_mm_2(tmp2, tmp1, yd_psv); 


    //A12-cd:1-ERER
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            yu_psv[i][j] = iden[i][j] - yu_psv[i][j]; 
            yd_psv[i][j] = iden[i][j] - yd_psv[i][j]; 
        }
    }
    //A12-cd:(1-ERER)^-1
    mat_inv_2(yu_psv); 
    mat_inv_2(yd_psv); 


    y_sh = 1./(1.-sexd[1]*grt[base+s_layer+1].rdu0*sexd[1]*sexu[1]*grt[base+s_layer].rud0*sexu[1]);  

    for(i = 0; i < 4; i++){
        for(j = 0; j < 4; j++){
            E[i][j] = emat[base+r_layer].E1[i][j];
        }
    }
    if((s_layer==r_layer)&&(zr<=zs)){
        rexd[0] = exp(-par[base+r_layer].r * (zr)); 
        rexd[1] = exp(-par[base+r_layer].v * (zr)); 
        rexu[0] = exp(-par[base+r_layer].r * (zs-zr)); 
        rexu[1] = exp(-par[base+r_layer].v * (zs-zr));     
        rexdt[0] = - par[base+r_layer].r * exp(-par[base+r_layer].r * (zr)); 
        rexdt[1] = - par[base+r_layer].v * exp(-par[base+r_layer].v * (zr)); 
        rexut[0] = par[base+r_layer].r * exp(-par[base+r_layer].r * (zs-zr)); 
        rexut[1] = par[base+r_layer].v * exp(-par[base+r_layer].v * (zs-zr));           
        
        u0 = rexd[1] * grt[base+s_layer-1+1].rud0 * sexu[1] + rexu[1]; 
        u0t = rexdt[1] * grt[base+s_layer-1+1].rud0 * sexu[1] + rexut[1] * par[base+r_layer].v;
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[base+s_layer-1+1].rud1[0][j] * sexu[j] 
                        + E[i][1] * rexd[1] * grt[base+s_layer-1+1].rud1[1][j] * sexu[j] 
                        + E[i][j+2] * rexu[j];
                ut[i][j] = E[i][0] * rexdt[0] * grt[base+s_layer-1+1].rud1[0][j] * sexu[j] 
                        + E[i][1] * rexdt[1] * grt[base+s_layer-1+1].rud1[1][j] * sexu[j] 
                        + E[i][j+2] * rexut[j];  
            }
        }

        yu0 = u0 * y_sh; 
        yu0t = u0t * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = u[i][0] * yu_psv[0][j] + u[i][1] * yu_psv[1][j]; 
                yu1t[i][j] = ut[i][0] * yu_psv[0][j] + ut[i][1] * yu_psv[1][j]; 
            }
        }
    }
    else if ((s_layer==r_layer)&&(zr>zs))
    {
        rexd[0] = exp(-par[base+s_layer].r * (zr-zs)); 
        rexd[1] = exp(-par[base+s_layer].v * (zr-zs)); 
        rexu[0] = exp(-par[base+s_layer].r * (par[base+s_layer].h-zr)); 
        rexu[1] = exp(-par[base+s_layer].v * (par[base+s_layer].h-zr)); 
        rexdt[0] = - par[base+r_layer].r * rexd[0]; 
        rexdt[1] = - par[base+r_layer].v * rexd[1]; 
        rexut[0] = par[base+r_layer].r   * rexu[0]; 
        rexut[1] = par[base+r_layer].v   * rexu[1];  
        d0 = rexd[1] + rexu[1] * grt[base+s_layer+1].rdu0 * sexd[1]; 
        d0t = rexdt[1] + rexut[1] * grt[base+s_layer+1].rdu0 * sexd[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[base+s_layer+1].rdu1[0][j] * sexd[j] 
                        + E[i][3] * rexu[1] * grt[base+s_layer+1].rdu1[1][j] * sexd[j]; 
                dt[i][j] = E[i][j] * rexdt[j]
                        + E[i][2] * rexut[0] * grt[base+s_layer+1].rdu1[0][j] * sexd[j] 
                        + E[i][3] * rexut[1] * grt[base+s_layer+1].rdu1[1][j] * sexd[j]; 
            }
        }
        yd0 = d0 * y_sh; 
        yd0t = d0t * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = d[i][0] * yd_psv[0][j] + d[i][1] * yd_psv[1][j]; 
                yd1t[i][j] = dt[i][0] * yd_psv[0][j] + dt[i][1] * yd_psv[1][j]; 
            }
        }         
    }
    else if(r_layer<s_layer){
        ex[0] = exp(-par[base+r_layer].r*par[base+r_layer].h); 
        ex[1] = exp(-par[base+r_layer].v*par[base+r_layer].h); 
        u0 = rexd[1] * grt[base+r_layer].rud0 * ex[1] + rexu[1];
        u0t = rexdt[1] * grt[base+r_layer].rud0 * ex[1] + rexut[1];  
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[base+r_layer-1+1].rud1[0][j] * ex[j] 
                        + E[i][1] * rexd[1] * grt[base+r_layer-1+1].rud1[1][j] * ex[j] 
                        + E[i][j+2] * rexu[j]; 
                ut[i][j] = E[i][0] * rexdt[0] * grt[base+r_layer-1+1].rud1[0][j] * ex[j] 
                        + E[i][1] * rexdt[1] * grt[base+r_layer-1+1].rud1[1][j] * ex[j] 
                        + E[i][j+2] * rexut[j]; 
            }
        }
        


        for(lay=s_layer-1;lay>=r_layer;lay--){
            if(lay==(s_layer-1)){
                yu0 = grt[base+lay+1].tu0 * sexu[1] * y_sh; 
                //if(ikw[idx].nk==13&&ikw[idx].nw==0){
                //    printf("CC:\n"); 
                //    cprint(grt[base+lay+1].tu0);
                //    cprint(sexu[1]); 
                //    cprint(y_sh); 
                //    cprint(yu0)
                //}
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yu1[i][j] = grt[base+lay+1].tu1[i][0] * sexu[0] * yu_psv[0][j]
                                  + grt[base+lay+1].tu1[i][1] * sexu[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[base+lay+1].r * par[base+lay+1].h); 
                ex[1] = exp(-par[base+lay+1].v * par[base+lay+1].h); 


                yu0 = grt[base+lay+1].tu0 * ex[1] * yu0; 
                //if(ikw[idx].nk==13&&ikw[idx].nw==0){
                //    printf("DD:\n"); 
                //    cprint(grt[base+lay+1].tu0);
                //    cprint(yu0); 
                //    cprint(ex[1])
                //}
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[base+lay+1].tu1[i][0] * ex[0] * yu1[0][j]
                                   + grt[base+lay+1].tu1[i][1] * ex[1] * yu1[1][j]; 
                    }
                }




                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yu1[i][j] = tmp1[i][j]; 
                    }
                }              

            }
        }
        
        yu0 = u0 * yu0; 
        yu0t = u0t * yu0;  
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = u[i][0] * yu1[0][j] + u[i][1] * yu1[1][j]; 
                tmp1t[i][j] = ut[i][0] * yu1[0][j] + ut[i][1] * yu1[1][j];
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = tmp1[i][j]; 
                yu1t[i][j] = tmp1t[i][j]; 
            }
        }                    


    }
    else if (r_layer>s_layer)
    {
        ex[0] = exp(-par[base+r_layer-1].r*par[base+r_layer-1].h); 
        ex[1] = exp(-par[base+r_layer-1].v*par[base+r_layer-1].h); 
        d0 = rexu[1] * grt[base+r_layer+1].rdu0 * ex[1] + rexd[1]; 
        d0t = rexut[1] * grt[base+r_layer+1].rdu0 * ex[1] + rexdt[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[base+r_layer+1].rdu1[0][j] * ex[j] 
                        + E[i][3] * rexu[1] * grt[base+r_layer+1].rdu1[1][j] * ex[j]; 
                dt[i][j] = E[i][j] * rexdt[j]
                        + E[i][2] * rexut[0] * grt[base+r_layer+1].rdu1[0][j] * ex[j] 
                        + E[i][3] * rexut[1] * grt[base+r_layer+1].rdu1[1][j] * ex[j]; 
            }
        }
        for(lay=s_layer+1;lay<=r_layer;lay++){
            if(lay==(s_layer+1)){
                yd0 = grt[base+lay+1].td0 * sexd[1] * y_sh; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yd1[i][j] = grt[base+lay+1].td1[i][0] * sexd[0] * yu_psv[0][j]
                                  + grt[base+lay+1].td1[i][1] * sexd[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[base+lay-2].r * par[base+lay-2].h); 
                ex[1] = exp(-par[base+lay-2].v * par[base+lay-2].h); 
                yd0 = grt[base+lay+1].td0 * ex[1] * yd0; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[base+lay+1].td1[i][0] * ex[0] * yd1[0][j]
                                   + grt[base+lay+1].td1[i][1] * ex[1] * yd1[1][j]; 
                    }
                }
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yd1[i][j] = tmp1[i][j]; 
                    }
                }                                
            }
        }
        yd0 = d0 * yd0; 
        yd0t = d0 * yd0; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = d[i][0] * yd1[0][j] + d[i][1] * yd1[1][j]; 
                tmp1t[i][j] = dt[i][0] * yd1[0][j] + dt[i][1] * yd1[1][j]; 
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = tmp1[i][j]; 
                yd1t[i][j] = tmp1t[i][j]; 
            }
        }      
    }

    if(((r_layer==s_layer)&&(zr<zs))||(r_layer<s_layer)){
        ppar[idx].b1 = yu0 * src[idx].shu1; 
        ppar[idx].b2 = yu0 * src[idx].shu2; 
        ppar[idx].b1t = yu0t * src[idx].shu1; 
        ppar[idx].b2t = yu0t * src[idx].shu2; 
        for(i=0;i<2;i++){
            ppar[idx].a0[i] = yu1[i][0] * src[idx].su0[0] + yu1[i][1] * src[idx].su0[1]; 
            ppar[idx].a1[i] = yu1[i][0] * src[idx].su1[0] + yu1[i][1] * src[idx].su1[1]; 
            ppar[idx].a2[i] = yu1[i][0] * src[idx].su2[0] + yu1[i][1] * src[idx].su2[1];     
            ppar[idx].a0t[i] = yu1t[i][0] * src[idx].su0[0] + yu1t[i][1] * src[idx].su0[1]; 
            ppar[idx].a1t[i] = yu1t[i][0] * src[idx].su1[0] + yu1t[i][1] * src[idx].su1[1]; 
            ppar[idx].a2t[i] = yu1t[i][0] * src[idx].su2[0] + yu1t[i][1] * src[idx].su2[1];                                                          
        }
    }
    else if(((r_layer==s_layer)&&(zr>zs))||(r_layer>s_layer)){
        ppar[idx].b1 = yd0 * src[idx].shd1;
        ppar[idx].b2 = yd0 * src[idx].shd2;
        ppar[idx].b1t = yd0t * src[idx].shd1;
        ppar[idx].b2t = yd0t * src[idx].shd2;
        for(i=0;i<2;i++){
            ppar[idx].a0[i] = yd1[i][0] * src[idx].sd0[0] + yd1[i][1] * src[idx].sd0[1]; 
            ppar[idx].a1[i] = yd1[i][0] * src[idx].sd1[0] + yd1[i][1] * src[idx].sd1[1]; 
            ppar[idx].a2[i] = yd1[i][0] * src[idx].sd2[0] + yd1[i][1] * src[idx].sd2[1]; 
            ppar[idx].a0t[i] = yd1t[i][0] * src[idx].sd0[0] + yd1t[i][1] * src[idx].sd0[1]; 
            ppar[idx].a1t[i] = yd1t[i][0] * src[idx].sd1[0] + yd1t[i][1] * src[idx].sd1[1]; 
            ppar[idx].a2t[i] = yd1t[i][0] * src[idx].sd2[0] + yd1t[i][1] * src[idx].sd2[1]; 
        }        
    }
}



__global__ void get_src_single(
    int_k_w *ikw, 
    engine_matrix *emat,  
    parameter *par, 
    grtc *grt, 
    ssrc *src, 
    double zs,  
    int s_layer, 
    int n_layer
){
    int i, j; 
    cplx ssh; 
    cplx tmp1, f[2][2]; 
    cplx sexd[2], sexu[2]; 
    cplx iden[2][2]; 
    cplx wru0[2][2], wrd0[2][2], wru1[2][2], wrd1[2][2]; 
    parameter spar; 
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    double k = ikw[idx].k;
    cplx w = ikw[idx].w;  
    iden[0][0] = iden[1][1] = 1.0; 
    iden[1][0] = iden[0][1] = 0.0; 
    
    spar = par[base+s_layer]; 
    
    sexd[0] = exp(-spar.r * (spar.h - zs)); 
    sexd[1] = exp(-spar.v * (spar.h - zs)); 
    sexu[0] = exp(-spar.r * (zs)); 
    sexu[1] = exp(-spar.v * (zs)); 
    
    tmp1 = spar.b / (2.0 * spar.u * spar.a * w * spar.r * spar.v); 
    f[0][0] = spar.b * spar.v * spar.r * tmp1; 
    f[0][1] = - spar.a * k * spar.r * tmp1; 
    f[1][0] = k * spar.b * spar.v * tmp1; 
    f[1][1] = - spar.a * spar.r * spar.v * tmp1; 


    ssh = 1. / (2. * spar.u * spar.v); 
    
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            //ER-I
            wru0[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j] - iden[i][j]; 
            //I-ER
            wrd0[i][j] = iden[i][j] - sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 
            //ER+I
            wru1[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j] + iden[i][j]; 
            //I+ER
            wrd1[i][j] = iden[i][j] + sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 
        }
    }

    for(i=0;i<2;i++){
        src[idx].gu0[i] = wru0[i][0] * f[0][0] + wru0[i][1] * f[0][1]; 
        src[idx].gd0[i] = wrd0[i][0] * f[0][0] + wrd0[i][1] * f[0][1]; 
        src[idx].gu1[i] = wru1[i][0] * f[1][0] + wru1[i][1] * f[1][1]; 
        src[idx].gd1[i] = wrd1[i][0] * f[1][0] + wrd1[i][1] * f[1][1]; 
    }
    src[idx].gshu = (1. + sexd[1] * grt[base+s_layer+1].rdu0 * sexd[1]) * ssh;
    src[idx].gshd = (sexu[1] * grt[base+s_layer-1+1].rud0 * sexu[1] + 1.) * ssh;
}

__global__ void get_src_tensor(
    int_k_w *ikw, 
    engine_matrix *emat,  
    parameter *par, 
    grtc *grt, 
    ssrc *src, 
    double zs,  
    int s_layer, 
    int n_layer
){
    int i, j; 
    cplx ssh; 
    cplx s0[2], s1[2], s2[2];
    cplx tmp1;
    cplx wru0[2][2], wru1[2][2], wrd0[2][2], wrd1[2][2];  
    cplx sexd[2], sexu[2]; 
    cplx iden[2][2]; 
    parameter spar; 
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int base = idx * (n_layer); 
    double k = ikw[idx].k;
    cplx w = ikw[idx].w;  
    iden[0][0] = iden[1][1] = 1.0; 
    iden[1][0] = iden[0][1] = 0.0; 
    
    spar = par[base+s_layer]; 
    
    sexd[0] = exp(-spar.r * (spar.h - zs)); 
    sexd[1] = exp(-spar.v * (spar.h - zs)); 
    sexu[0] = exp(-spar.r * (zs)); 
    sexu[1] = exp(-spar.v * (zs)); 
    
    tmp1 = spar.b / (2.0 * spar.u * spar.a * w); 


    ssh = 1. / (2. * spar.u * spar.v); 
    s2[0] = tmp1 * k / spar.r * k * spar.b; 
    s2[1] = - k * tmp1 / spar.r * spar.a * spar.r; 
    s1[0] = tmp1 / spar.v * (-2.0) * k * spar.v * spar.b; 
    s1[1] = tmp1 / spar.v * (k*k + spar.v * spar.v) * spar.a; 
    s0[0] = tmp1 * spar.r * spar.b; 
    s0[1] = tmp1 * spar.a * (-k); 
    //if(ikw[idx].nk==12&&ikw[idx].nw==223){
    //    cprint(spar.r); 
    //    cprint(spar.b); 
    //    printf("ssst%lf\n", k); 
    //    cprint(ssh); 
    //    cprint(s2[0]); 
    //    cprint(s2[1]); 
    //    cprint(s1[0]); 
    //    cprint(s1[1]); 
    //    cprint(s0[0]); 
    //    cprint(s0[1]); 
    //}
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            //ER-I
            wru0[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j] - iden[i][j]; 
            //I-ER
            wrd0[i][j] = iden[i][j] - sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 
            //ER+I
            wru1[i][j] = sexd[i] * grt[base+s_layer+1].rdu1[i][j] * sexd[j] + iden[i][j]; 
            //I+ER
            wrd1[i][j] = iden[i][j] + sexu[i] * grt[base+s_layer-1+1].rud1[i][j] * sexu[j]; 
        }
    }

    for(i=0;i<2;i++){
        src[idx].su0[i] = wru1[i][0] * s0[0] + wru1[i][1] * s0[1]; 
        src[idx].sd0[i] = wrd1[i][0] * s0[0] + wrd1[i][1] * s0[1]; 
        src[idx].su1[i] = wru0[i][0] * s1[0] + wru0[i][1] * s1[1]; 
        src[idx].sd1[i] = wrd0[i][0] * s1[0] + wrd0[i][1] * s1[1];  
        src[idx].su2[i] = wru1[i][0] * s2[0] + wru1[i][1] * s2[1]; 
        src[idx].sd2[i] = wrd1[i][0] * s2[0] + wrd1[i][1] * s2[1];                
    }
    src[idx].shu1 = (1. + sexd[1] * grt[base+s_layer+1].rdu0 * sexd[1]) * spar.v * ssh;
    src[idx].shd1 = (sexu[1] * grt[base+s_layer-1+1].rud0 * sexu[1] - 1.) * spar.v * ssh;
    src[idx].shu2 = (1. + sexd[1] * grt[base+s_layer+1].rdu0 * sexd[1]) * k * ssh;
    src[idx].shd2 = (sexu[1] * grt[base+s_layer-1+1].rud0 * sexu[1] - 1.) * k * ssh;
}




__host__ __device__  float32 besselj0(float32 x)
{
    float32 ax,x1,x2,theta,fct;
    float32 besselj0;
    float32 a[7]={1.00000000,-2.24999970,1.26562080,-0.31638660,0.04444790,-0.00394440,0.00021000};
    float32 b[7]={0.79788456,-0.00000077,-0.00552740,-0.00009512,0.00137237,-0.00072805,0.00014476};
    float32 c[7]={-0.78539816,-0.04166397,-0.00003954, 0.00262573,-0.00054125,-0.00029333,0.00013558};
    ax=abs(x);
    if(ax<=3)
    {
        x2=(ax/3.0)*(ax/3.0);
        besselj0=a[0]+x2*(a[1]+x2*(a[2]+x2*(a[3]+x2*(a[4]+x2*(a[5]+x2*a[6])))));
    }
    else
    {
        x1=3.0/ax;
        fct=b[0]+x1*(b[1]+x1*(b[2]+x1*(b[3]+x1*(b[4]+x1*(b[5]+x1*b[6])))));
        theta=ax+c[0]+x1*(c[1]+x1*(c[2]+x1*(c[3]+x1*(c[4]+x1*(c[5]+x1*c[6])))));
        besselj0=fct*cos(theta)/sqrt(ax);
    }
    return besselj0;
}

__host__ __device__   float32 besselj1(float32 x)
{
    float32 ax,x1,x2,theta,fct;
    float32 besselj1;
    float32 sign;
    float32 a[7]={0.50000000,-0.56249985,0.21093573,-0.03954289,0.00443319,-0.00031761,0.00001109};
    float32 b[7]={0.79788456,0.00000156,0.01659667,0.00017105,-0.00249511,0.00113653,-0.00020033};
    float32 c[7]={-2.35619449,0.12499612,0.00005650,-0.00637879,0.00074348,0.00079824,-0.00029166};
    ax=abs(x);
    if(ax<=3.0)
    {
        x2=(ax/3.0)*(ax/3.0);
        besselj1=x*(a[0]+x2*(a[1]+x2*(a[2]+x2*(a[3]+x2*(a[4]+x2*(a[5]+x2*a[6]))))));
    }
    else
    {
        x1=3.0/ax;
        fct=b[0]+x1*(b[1]+x1*(b[2]+x1*(b[3]+x1*(b[4]+x1*(b[5]+x1*b[6])))));
        theta=ax+c[0]+x1*(c[1]+x1*(c[2]+x1*(c[3]+x1*(c[4]+x1*(c[5]+x1*c[6])))));
        if(x>=0)sign=1.0;
        else sign=-1.0;
        besselj1=sign*fct*cos(theta)/sqrt(ax);
    }
    return besselj1;
}

__host__ __device__ float32 besselj(float32 x, int n)
{
      const int iacc=40;
      float32 bigno,bigni;
      bigno=10000000000.0;
      bigni=0.0000000001;
      int j,jsum,m,nn;
      float32 besselj;
      float32 ax,bj,bjm,bjp,sum,tox;
      nn=n;
      if(0==n)besselj=besselj0(x);
      else if(1==n)besselj=besselj1(x);
      else
      {
          if(n<0)n=-n;
          ax=abs(x);
          if(0==ax)besselj=0.0;
          else if(ax>n)
          {
              tox=2.0/ax;
              bjm=besselj0(ax);
              bj=besselj1(ax);
              for(j=0;j<n-1;j++)
              {
                  bjp=((float32)(j+1))*tox*bj-bjm;
                  bjm=bj;
                  bj=bjp;
              }
              besselj=bj;
          }
          else
          {
              tox=2.0/ax;
              m=2*((n+(int)(sqrt((float32)(iacc*n))))/2);
              besselj=0.0;
              jsum=0;
              sum=0.0;
              bjp=0.0;
              bj=1.0;
              for(j=m;j>=1;j--)
              {
                  bjm=((float32)(j))*tox*bj-bjp;
                  bjp=bj;
                  bj=bjm;
                  if(fabs(bj)>bigno)
                  {
                      bj=bj*bigni;
                      bjp=bjp*bigni;
                      besselj=besselj*bigni;
                      sum=sum*bigni;
                  }
                  if(jsum!=0)sum=sum+bj;
                  jsum=1-jsum;
                  if(j==n)besselj=bjp;
              }
              sum=2.0*sum-bj;
              besselj=besselj/sum;
          }
          if(x<0.0&&(n%2)==1)besselj=-besselj;
      }
      if(nn<0&&(n%2)==1)besselj=-besselj;
      return besselj;
}
__device__ __host__  float32 d_besselj(float32 x,int n)
{
    float32 bsl;
    float32 bigni=0.0000000001;
    if(fabs(x)<bigni)
    {
        bsl=(besselj(bigni,n)-besselj(0,n))/bigni;
        return bsl;
    }
    bsl=((float32)n)*besselj(x,n)/x-besselj(x,n+1);
    return  bsl;
}

__device__ __host__  float32 dd_besselj(float32 x,int n)
{
    float32 bsl;
    float32 bigni=0.0000000001;
    if(fabs(x)<bigni)
    {
        bsl=(d_besselj(bigni,n)-d_besselj(0,n))/bigni;
        return bsl;
    }
    bsl=((float32)n)*besselj(x,n-1)/x-besselj(x,n) * (n + n * n + x * x)/x/x;
    return  bsl;
}



__global__ void get_int_val_single(
    int_k_w *ikw, 
    popm_parameter *ppar, 
    int_out *iout, 
    double t0, 
    double dist
){
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int i; 
    double k = ikw[idx].k; 
    double j0 = besselj(dist*k, 0); 
    double j1 = besselj(dist*k, 1); 
    double j2 = besselj(dist*k, 2); 
    double j0d = d_besselj(dist*k, 0); 
    double j1d = d_besselj(dist*k, 1); 
    double j2d = d_besselj(dist*k, 2); 
    double j0dd = dd_besselj(dist*k, 0); 
    double j1dd = dd_besselj(dist*k, 1); 
    double j2dd = dd_besselj(dist*k, 2); 
    double j3dd = dd_besselj(dist*k, 3); 
    cplx intp[15]; 
    double phi; 


    intp[0] = ppar[idx].b1 * j1 / dist + ppar[idx].a1[0] * k * j1d; 
    intp[1] = ppar[idx].a0[0] * k * j0d; 
    intp[2] = ppar[idx].a1[0] * j1 / dist + ppar[idx].b1 * k * j1d; 
    intp[3] = ppar[idx].a1[1] * k * j1; 
    intp[4] = ppar[idx].a0[1] * k * j0; 

    //du/dz
    intp[5] = ppar[idx].b1t * j1 / dist + ppar[idx].a1t[0] * k * j1d; 
    intp[6] = ppar[idx].a0t[0] * k * j0d; 
    intp[7] = ppar[idx].a1t[0] * j1 / dist + ppar[idx].b1t * k * j1d; 
    intp[8] = ppar[idx].a1t[1] * k * j1; 
    intp[9] = ppar[idx].a0t[1] * k * j0; 

    //du/dr
    intp[10] = -ppar[idx].b1 * j1 / dist / dist + ppar[idx].b1 * j1d * k  + ppar[idx].a1[0] * k * j1dd * k; 
    intp[11] =  ppar[idx].a0[0] * k * j0dd * k; 
    intp[12] = -ppar[idx].a1[0] * j1 / dist/dist + ppar[idx].a1[0] * j1d * k / dist + ppar[idx].b1 * k * j1dd * k; 
    intp[13] =  ppar[idx].a1[1] * k * j1d * k; 
    intp[14] =  ppar[idx].a0[1] * k * j0d * k; 



    for(i=0;i<15;i++){
        iout[idx].sums[i] = intp[i]; 
    }

    phi = real(ikw[idx].w) * t0; 
    for(i=0;i<5;i++){
        //iout[idx].sums[i] = iout[idx].sums[i] * cplx(cos(phi), sin(phi)); 
    }
    //if(ikw[idx].nk==13&&ikw[idx].nw==222){
    //    printf("TESTIMG2\n"); 
    //    cprint(iout[idx].sums[0]); 
    //    cprint(iout[idx].sums[1]);
    //    cprint(iout[idx].sums[2]); 
    //    cprint(iout[idx].sums[3]);
    //    cprint(iout[idx].sums[4]);
    //}
    
}




__global__ void get_int_val_tensor(
    int_k_w *ikw, 
    popm_parameter *ppar, 
    int_out *iout, 
    double t0, 
    double dist
){
    int idx=blockIdx.x*blockDim.x + threadIdx.x;
    int i; 
    double k = ikw[idx].k; 
    double j0 = besselj(dist*k, 0); 
    double j1 = besselj(dist*k, 1); 
    double j2 = besselj(dist*k, 2); 
    double j0d = d_besselj(dist*k, 0); 
    double j1d = d_besselj(dist*k, 1); 
    double j2d = d_besselj(dist*k, 2); 
    double j0dd = dd_besselj(dist*k, 0); 
    double j1dd = dd_besselj(dist*k, 1); 
    double j2dd = dd_besselj(dist*k, 2); 
    double j3dd = dd_besselj(dist*k, 3); 
    cplx intp[30]; 
    double phi; 


        //Ur
        intp[0] = ppar[idx].b2 * 2.0 * j2 / dist + ppar[idx].a2[0] * k * j2d; 
        intp[1] = ppar[idx].b1 * j1 / dist + ppar[idx].a1[0] * k * j1d; 
        intp[2] = ppar[idx].a2[1] * k * j0d; 
        intp[3] = ppar[idx].a0[1] * k * j0d; 
        //Uo
        intp[4] = ppar[idx].b2 * k * j2d + 2.0 * ppar[idx].a2[0] * j2 / dist; 
        intp[5] = ppar[idx].b2 * k * j1d + ppar[idx].a1[0] * j1 / dist; 
        //Uz
        intp[6] = ppar[idx].a2[1] * k * j2; 
        intp[7] = ppar[idx].a1[1] * k * j1; 
        intp[8] = ppar[idx].a0[1] * k * j0; 
        intp[9] = ppar[idx].a2[1] * k * j0; 


        //dUr/dz
        intp[10] = ppar[idx].b2t * 2.0 * j2 / dist + ppar[idx].a2t[0] * k * j2d; 
        intp[11] = ppar[idx].b1t * j1 / dist + ppar[idx].a1t[0] * k * j1d; 
        intp[12] = ppar[idx].a2t[1] * k * j0d; 
        intp[13] = ppar[idx].a0t[1] * k * j0d; 
        //dUo/dz
        intp[14] = ppar[idx].b2t * k * j2d + 2.0 * ppar[idx].a2t[0] * j2 / dist; 
        intp[15] = ppar[idx].b2t * k * j1d + ppar[idx].a1t[0] * j1 / dist; 
        //dUz/dz
        intp[16] = ppar[idx].a2t[1] * k * j2; 
        intp[17] = ppar[idx].a1t[1] * k * j1; 
        intp[18] = ppar[idx].a0t[1] * k * j0; 
        intp[19] = ppar[idx].a2t[1] * k * j0; 

        //dUr/dr
        intp[10] = - ppar[idx].b2 * 2.0 * j2 / dist / dist + ppar[idx].b2 * 2.0 * j2d * k / dist + ppar[idx].a2[0] * k * j2dd * k; 
        intp[11] = -ppar[idx].b1 * j1 /dist /dist + ppar[idx].b1 * j1d * k / dist + ppar[idx].a1[0] * k * j1dd * k; 
        intp[12] = ppar[idx].a2[1] * k * j0dd * k; 
        intp[13] = ppar[idx].a0[1] * k * j0dd * k; 
        //dUo/dr
        intp[14] = ppar[idx].b2 * k * j2dd * k - 2.0 * ppar[idx].a2[0] * j2 / dist / dist + 2.0 * ppar[idx].a2[0] * j2d * k / dist; 
        intp[15] = ppar[idx].b2 * k * j1dd * k - ppar[idx].a1[0] * j1 / dist / dist + ppar[idx].a1[0] * j1d * k / dist; 
        //dUz/dr
        intp[16] = ppar[idx].a2[1] * k * j2d * k; 
        intp[17] = ppar[idx].a1[1] * k * j1d * k; 
        intp[18] = ppar[idx].a0[1] * k * j0d * k; 
        intp[19] = ppar[idx].a2[1] * k * j0d * k; 
    for(i=0;i<30;i++){
        iout[idx].sums[i] = intp[i]; 
    }
    //if(ikw[idx].nk==366&&ikw[idx].nw==222){
    //    printf("TESTIMG%d\n", ikw[idx].nk); 
    //    cprint(iout[idx].sums[0]); 
    //    cprint(iout[idx].sums[1]);
    //    cprint(iout[idx].sums[2]); 
    //    cprint(iout[idx].sums[3]);
    //    cprint(iout[idx].sums[4]);
    //}
    phi = real(ikw[idx].w) * t0; 
    for(i=0;i<5;i++){
        //iout[idx].sums[i] = iout[idx].sums[i] * cplx(cos(phi), sin(phi)); 
    }
    //if(ikw[idx].nk==13&&ikw[idx].nw==222){
    //    printf("TESTIMG2\n"); 
    //    cprint(iout[idx].sums[0]); 
    //    cprint(iout[idx].sums[1]);
    //    cprint(iout[idx].sums[2]); 
    //    cprint(iout[idx].sums[3]);
    //    cprint(iout[idx].sums[4]);
    //}
    
}


void integrate(int_k_w *ikw, int nkw, int nw, layer_parameter *par, intpar ipar, int_out *out){
    
    cplx wt;         
    
    cuDoubleComplex alpha={1.0,0.0};
    cuDoubleComplex beta={0.0,0.0};
    mrtc *d_mrt, tmrt; 
    grtc *d_grt; 
    int_k_w *d_ikw; 
    cplx tem1, tem2, tem3; 
    parameter *d_par; 
    popm_parameter *d_ppar; 
    layer_parameter *d_lay_par; 
    engine_matrix *d_emat; 
    int_out *d_iout; 
    force *d_fc; 
    ssrc *d_src; 
    int_out *h_iout; 
    cuDoubleComplex *d_E1; 
    cuDoubleComplex *d_E2, *d_E3; 
    cuDoubleComplex **d_dE1, **d_dE2, **d_dE3;
    
    cudaError_t cerr; 
    struct timespec rt1, rt2, rt3;  
    size_t time1, time2; 
    float32 rtime; 

    char transa = 'n', transb = 'n';
    int n_layer = ipar.n_layer; 


    //定义blas句柄

    cublasHandle_t blasHandle;
    int *info; 

    clock_gettime(CLOCK_MONOTONIC,&rt1);
    time1 = clock(); 
    cerr = cudaGetLastError(); 
    printf("cerr0=%s\n", cudaGetErrorString(cerr)); 
    CUBLAS_CHECK(cublasCreate(&blasHandle));
    CUDA_CHECK(cudaMalloc((void **)(&info), nkw  * sizeof(int) * n_layer));
    //RT系数
    CUDA_CHECK(cudaMalloc((void **)&d_mrt, (n_layer+0) * nkw * sizeof(mrtc)));
    CUDA_CHECK(cudaMalloc((void **)&d_grt, (n_layer+0) * nkw * sizeof(grtc)));
    CUDA_CHECK(cudaMalloc((void **)&d_emat, (n_layer+0) * nkw * sizeof(grtc) + 1));
    CUDA_CHECK(cudaMalloc((void **)&d_E1, nkw * sizeof(cuDoubleComplex)*16 * n_layer));
    CUDA_CHECK(cudaMalloc((void **)&d_E2, nkw * sizeof(cuDoubleComplex)*16 * n_layer));
    CUDA_CHECK(cudaMalloc((void **)&d_E3, nkw * sizeof(cuDoubleComplex)*16 * n_layer));
    CUDA_CHECK(cudaMalloc((void **)&d_dE1, nkw * sizeof(cuDoubleComplex *) * n_layer));
    CUDA_CHECK(cudaMalloc((void **)&d_dE2, nkw * sizeof(cuDoubleComplex *) * n_layer));
    CUDA_CHECK(cudaMalloc((void **)&d_dE3, nkw * sizeof(cuDoubleComplex *) * n_layer));
    //积分所有频点
    CUDA_CHECK(cudaMalloc((void **)&d_ikw, nkw * sizeof(int_k_w)));
    //拷贝参数
    CUDA_CHECK(cudaMalloc((void **)&d_par,   n_layer   * nkw *  sizeof(parameter) + 1));
    CUDA_CHECK(cudaMalloc((void **)&d_lay_par,   n_layer    *   sizeof(layer_parameter)));
    CUDA_CHECK(cudaMalloc((void **)&d_ppar,   nkw    *   sizeof(popm_parameter)));
    CUDA_CHECK(cudaMalloc((void **)&d_src,   nkw    *   sizeof(ssrc)));
    CUDA_CHECK(cudaMalloc((void **)&d_fc,   n_layer * nkw    *   sizeof(force)));
    CUDA_CHECK(cudaMemcpy(d_lay_par, par, n_layer * sizeof(layer_parameter), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ikw, ikw, nkw * sizeof(int_k_w), cudaMemcpyHostToDevice));
    cerr = cudaGetLastError(); 
    clock_gettime(CLOCK_MONOTONIC,&rt2);
    time2 = clock(); 
    rtime = (rt2.tv_sec-rt1.tv_sec) * 1.0 + (rt2.tv_nsec - rt1.tv_nsec) * 1e-9; 
    printf("Data allocated time %lf\n", (float)(time2-time1)/((float)CLOCKS_PER_SEC)); 
    printf("Data Allocated:%s, time=%lf\n", cudaGetErrorString(cerr), rtime); 
    cudaDeviceSynchronize();
    int nthreads = 32;
    int ngrids=(nkw+nthreads-1)/nthreads;
    //ngrids = 1;
    //nthreads = 64; 
    cudaStream_t s; cudaStreamCreate(&s);
    setting_array<<<nkw*n_layer, 1,0,s>>>(d_E1,d_E2,d_E3,d_dE1,d_dE2,d_dE3,16);
    
    //cudaDeviceSynchronize();
    cerr = cudaGetLastError(); 
    printf("Setting array:%s\n", cudaGetErrorString(cerr)); 
    //设置传播矩阵参数
    set_pars<<<nkw, n_layer,0,s>>>(d_ikw, d_lay_par, d_par);    
    //cudaDeviceSynchronize(); 
    cerr = cudaGetLastError(); 
    printf("Setting pars:%s\n", cudaGetErrorString(cerr)); 
    //计算本征值矩阵
    get_emat<<<nkw, n_layer,0,s>>>(d_ikw, d_par, d_emat); 
    //cudaDeviceSynchronize(); 
    cerr = cudaGetLastError(); 
    printf("Setting emat:%s\n", cudaGetErrorString(cerr)); 
    //计算传播矩阵 
    get_mrt_sh<<<nkw, n_layer,0,s>>>(d_ikw, d_par, d_mrt); 
    for(int lay=0;lay<n_layer-1;lay++){
        get_emat_e1e2<<<nkw, n_layer,0,s>>>(d_emat, d_E1, d_E2); 
        cerr = cudaGetLastError(); 
        CUBLAS_CHECK(cublasZmatinvBatched(blasHandle, 4, d_dE1, 4, d_dE3, 4, info, nkw*(n_layer)));  
     
        cerr = cudaGetLastError(); 
        CUBLAS_CHECK(cublasZgemmBatched(blasHandle, char_to_cublas_trans(transa), char_to_cublas_trans(transb), 4, 4, 4, &alpha, d_dE3, 4, d_dE2, 4, &beta, d_dE1, 4, nkw*(n_layer)));
        get_mrt_psv<<<nkw, n_layer,0,s>>>(d_dE1, d_mrt); 
    }
    get_grt<<<ngrids, nthreads,0,s>>>(d_ikw, d_par, d_emat, d_mrt, d_grt, ipar.s_layer, n_layer); 
    //cudaDeviceSynchronize(); 
    cerr = cudaGetLastError(); 
    clock_gettime(CLOCK_MONOTONIC,&rt2);
    rtime = (rt2.tv_sec-rt1.tv_sec) * 1.0 + (rt2.tv_nsec - rt1.tv_nsec) * 1e-9; 
    time2 = clock();
    printf("Setting time %lf\n", (float)(time2-time1)/((float)CLOCKS_PER_SEC)); 
    printf("Setting grt:%s, time=%lf\n", cudaGetErrorString(cerr), rtime); 
    printf("source type %d\n", ipar.src_type); 
    if(ipar.src_type==0){
        get_src_single<<<ngrids,nthreads,0,s>>>(d_ikw, d_emat, d_par, d_grt, d_src, ipar.zs, ipar.s_layer, n_layer); 
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting src:%s\n", cudaGetErrorString(cerr)); 
        get_ypars_a_single<<<ngrids,nthreads,0,s>>>(d_ikw, d_emat, d_par, d_grt, d_src, d_ppar, ipar.zr, ipar.zs, ipar.r_layer, ipar.s_layer, n_layer);
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting ypars_a:%s\n", cudaGetErrorString(cerr)); 
        CUDA_CHECK(cudaMalloc((void **)&d_iout,  nkw    *   sizeof(int_out)));

        get_int_val_single<<<ngrids,nthreads,0,s>>>(d_ikw, d_ppar, d_iout, ipar.t0, ipar.dist); 
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting int_val:%s\n", cudaGetErrorString(cerr)); 
    }
    else{
        get_src_tensor<<<ngrids,nthreads,0,s>>>(d_ikw, d_emat, d_par, d_grt, d_src, ipar.zs, ipar.s_layer, n_layer); 
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting src:%s\n", cudaGetErrorString(cerr)); 
        get_ypars_a_tensor<<<ngrids,nthreads,0,s>>>(d_ikw, d_emat, d_par, d_grt, d_src, d_ppar, ipar.zr, ipar.zs, ipar.r_layer, ipar.s_layer, n_layer);
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting ypars_a:%s\n", cudaGetErrorString(cerr)); 
        CUDA_CHECK(cudaMalloc((void **)&d_iout,  nkw    *   sizeof(int_out)));

        get_int_val_tensor<<<ngrids,nthreads,0,s>>>(d_ikw, d_ppar, d_iout, ipar.t0, ipar.dist); 
        //cudaDeviceSynchronize(); 
        cerr = cudaGetLastError(); 
        printf("Setting int_val:%s\n", cudaGetErrorString(cerr)); 
    }
    cublasSetStream(blasHandle, s);
    clock_gettime(CLOCK_MONOTONIC,&rt2);
    rtime = (rt2.tv_sec-rt1.tv_sec) * 1.0 + (rt2.tv_nsec - rt1.tv_nsec) * 1e-9; 
    printf("Finished time:%s, time=%lf\n", cudaGetErrorString(cerr), rtime);
    time2 = clock(); 
    printf("Finished time %lf\n", (float)(time2-time1)/((float)CLOCKS_PER_SEC)); 
    h_iout = (int_out *)malloc(sizeof(int_out)*nkw); 
    
    CUDA_CHECK(cudaMemcpy(h_iout, d_iout, nkw * sizeof(int_out), cudaMemcpyDeviceToHost));
    //int_out *out; 
    int_k_w tikw; 
    int count=0; 
    
    for(int jj=0;jj<nw*2;jj++){
        for(int kk=0;kk<30;kk++){
            out[jj].sums[kk] = 0.0; 
        }
    }
    for(int jj=0;jj<nkw;jj++){
        tikw = ikw[jj]; 
        //if(tikw.nw==222){
        //    for(int i=0;i<5;i++){
        //        fprintf(tfp, "%.11lf,%.11lf,", real(h_iout[jj].sums[i]), imag(h_iout[jj].sums[i])); 
        //    }
        //    fprintf(tfp, "\n"); 
        //    fflush(tfp);
        //}
        if(tikw.nk==0&&tikw.k==0.0){
            for(int kk=0;kk<30;kk++){
                out[tikw.nw].sums[kk] = h_iout[jj].sums[kk] * ipar.dk / 2.0; 
            }
        }
        if(tikw.nk==0&&tikw.k!=0.0){
            for(int kk=0;kk<30;kk++){
                out[tikw.nw].sums[kk] = h_iout[jj].sums[kk] * ipar.dk; 
            }
        }
        if(tikw.nk!=0){
            for(int kk=0;kk<30;kk++){
                out[tikw.nw].sums[kk] = out[tikw.nw].sums[kk] + h_iout[jj].sums[kk] * ipar.dk; 
            }
        }
        out[tikw.nw].w = tikw.w; 
    }
   
    double st0; 
    double ifftcoef; 
    double phi; 
    int wc = (int)nw * (1-ipar.taper);
    for(int iw=0;iw<(int)(nw);iw++){
        //printf("nw=%d\n", n_half_sample-nw); 
        for(int j=0;j<30;j++){
            phi = real(out[iw].w)*ipar.t0; 
            out[iw].sums[j] = out[iw].sums[j] * cplx(cos(phi), sin(phi));  
            st0 = 1; 
            if(iw>wc){
                st0 = 0.5*(1.+cos((iw-wc+1)*PI/(nw-wc+1))); 
            }
            ifftcoef = real(ipar.df) * nw / PI; 
            out[iw].sums[j] = st0 * out[iw].sums[j] * ifftcoef; 
        }
    }



    cplx invs; 
    for(int iw=0;iw<(int)(nw);iw++){
        //printf("nw=%d\n", n_half_sample-nw); 
        for(int i=0;i<30;i++){
            if(nw-iw>=nw){
                invs = 0.0; 
            }
            else{
                invs = out[nw-iw].sums[i]; 
            }
            invs = conj(invs); 
            out[nw+iw].sums[i] = invs; 
        }
    }


    clock_gettime(CLOCK_MONOTONIC,&rt2);
    float32 r = (rt2.tv_sec-rt1.tv_sec) * 1.0 + (rt2.tv_nsec - rt1.tv_nsec) * 1e-9; 
    printf("CUDA实际执行时间:%lf\n", r); 
    CUDA_CHECK(cudaFree(d_fc));
    CUDA_CHECK(cudaFree(d_src));
    CUDA_CHECK(cudaFree(d_ppar));
    CUDA_CHECK(cudaFree(d_mrt)); 
    CUDA_CHECK(cudaFree(d_emat)); 
    CUDA_CHECK(cudaFree(d_grt)); 
    CUDA_CHECK(cudaFree(d_par)); 
    CUDA_CHECK(cudaFree(d_lay_par)); 
    CUDA_CHECK(cudaFree(d_E1)); 
    CUDA_CHECK(cudaFree(d_E2)); 
    CUDA_CHECK(cudaFree(info)); 
    //cudaDeviceReset(); 
}

struct share_args{
    cplx w; 
    float32 k; 
    float32 kmax; 
    int_parameter ipar;
    layer_parameter *par;
    cplx sums[12]; 
}; 
typedef struct share_args share_args; 

//void *thread_func(void *args){
//    share_args *temp = (share_args *)args;
//    printf("TESTW%lf\n", temp->kmax); 
//    cprint(temp->w); 
//    integrate(temp->w, temp->kmax, temp->par, temp->ipar, temp->sums);
//}

void green(char *iname, char *oname_w, char *oname_t){
    layer_parameter *par;
    cplx w; 
   
    int n_layer;
    int i;  
    time_t t1, t2; 
    intpar ipar; 
    cplx sums[10]; 
    float32 vmax, vmin; 
    FILE *ifp; 
    float32 hr, hs; 
    int log2_n_sample; 

    ifp = fopen(iname, "r");


    float32 th1, ta01, tb01, trho1, tqp1, tqs1; 
    float32 th2, ta02, tb02, trho2, tqp2, tqs2; 
    int tcnt=0; 
    fscanf(ifp, "%lf\n", &hs);
    fscanf(ifp, "%lf\n", &hr); 
    fscanf(ifp, "%d\n", &n_layer); 
    printf("Nlayer%d\n", n_layer); 
    //Add one more layer to store source 
    par = (layer_parameter *)malloc(sizeof(layer_parameter)*(n_layer+1)); 
    fscanf(ifp, "%lf,%lf,%lf,%lf,%lf,%lf\n", &th1, &ta01, &tb01, &trho1, &tqp1, &tqs1);
    par[tcnt].h=th1; 
    par[tcnt].a0=ta01; 
    par[tcnt].b0=tb01; 
    par[tcnt].rho=trho1; 
    par[tcnt].qp=tqp1; 
    par[tcnt].qs=tqs1; 
    tcnt += 1; 
    for(i=1;i<n_layer;i++){
        fscanf(ifp, "%lf,%lf,%lf,%lf,%lf,%lf\n", &th2, &ta02, &tb02, &trho2, &tqp2, &tqs2);
        printf("Layer:%lf,%lf\n", th1, th2);
        if((th1<hs)&&(th2>hs)){
            par[tcnt].h=hs; 
            par[tcnt].a0=ta01; 
            par[tcnt].b0=tb01; 
            par[tcnt].rho=trho1; 
            par[tcnt].qp=tqp1; 
            par[tcnt].qs=tqs1; 
            tcnt += 1;
        }
        par[tcnt].h=th2; 
        par[tcnt].a0=ta02; 
        par[tcnt].b0=tb02; 
        par[tcnt].rho=trho2; 
        par[tcnt].qp=tqp2; 
        par[tcnt].qs=tqs2; 
        tcnt += 1; 
        th1 = th2; ta01=ta02; tb01=tb02; trho1=trho2; tqp1=tqp2; tqs1=tqs2;  
    }
    n_layer=n_layer + 1; 
    ipar.n_layer = n_layer; 
    vmax = 0; 
    vmin = 100000; 
    for(i=0;i<n_layer;i++){
        if(par[i].a0>vmax)vmax = par[i].a0; 
        if(par[i].b0<vmin)vmin = par[i].b0;
    }

    printf("\nLayer Parameter:\n"); 
    for(i=0;i<n_layer;i++){
        printf("%d,%lf,%lf,%lf,%lf,%lf,%lf\n", i, par[i].h, par[i].a0, par[i].b0, par[i].rho, par[i].qp, par[i].qs);
    }

    for(i=0;i<n_layer-1;i++){
        par[i].h = par[i+1].h - par[i].h;
    }
    for(i=0;i<n_layer;i++){
        par[i].u = par[i].b0 * par[i].b0 * par[i].rho; 
    }  

    fscanf(ifp, "%lf\n", &ipar.dt); 
    fscanf(ifp, "%lf\n", &ipar.dist); 
    fscanf(ifp, "%d\n", &ipar.nf1); 
    fscanf(ifp, "%d\n", &ipar.nf2); 
    fscanf(ifp, "%lf\n", &ipar.Twin);
    fscanf(ifp, "%d\n", &log2_n_sample); 
    fscanf(ifp, "%lf\n", &ipar.t0); 
    fscanf(ifp, "%lf\n", &ipar.taper); 
    fscanf(ifp, "%d\n", &ipar.n_ptam); 
    fscanf(ifp, "%d\n", &ipar.src_type); 
    printf("t0:%lf\n", ipar.t0); 
    printf("log2:%ld,%lf\n", log2_n_sample,hr); 
    ipar.k_factor = 1.1; 

    int n_half_sample = (int)pow(2, log2_n_sample-1);
    int n_sample = n_half_sample * 2; 
    float32 t_win = n_sample * ipar.dt; 
    float32 fmax = 1.0 / (2 * ipar.dt); 
    float32 df = 1/(n_sample * ipar.dt); 
    float32 img_freq = 2.0 / t_win; 
    float32 vavg = (vmax + vmin) * 0.5; 
    float32 max_dist = 0; //ipar.dist; 
    ipar.hs = hs; 
    ipar.hr = hr; 
    ipar.kmin = 0.0; 
    max_dist = 0; 
    float32 L = (vmax * t_win + max_dist +  vmax / vavg * sqrt(max_dist * max_dist + (hs-hr)*(hs-hr)) + 100);
    //float32 L = (vmax * t_win) + 100; 
    L = L * 2; 
    ipar.dk = 2 * PI / L; 
    for(i=0;i<n_layer;i++){
        if(hs<=par[i].h){
            ipar.s_layer = i; 
            ipar.zs = hs; 
            break; 
        }
        else{
            hs -= par[i].h; 
        }
    }

    for(i=0;i<n_layer;i++){
        if(hs<=par[i].h){
            ipar.r_layer = i; 
            ipar.zr = hr; 
            break; 
        }
        else{
            hr -= par[i].h; 
        }
    }


    printf("\nBasic parameters:\n"); 
    printf("Parameters:%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", vmax, t_win, max_dist, vmin, vavg, ipar.zs, ipar.zr); 
    printf("Time Window:%lf\n", t_win); 
    printf("Sampling points:%d\n", n_sample); 
    printf("Time delta:%lf\n", ipar.dt); 
    printf("Max freq.:%lf\n", fmax); 
    printf("Max vel.:%lf\n", vmax);
    printf("Min vel.:%lf\n", vmin);  
    printf("L:%lf\n", L); 
    printf("df:%lf\n", df); 
    printf("dk:%lf\n", ipar.dk); 
    printf("hs:%lf\n", hs); 
    printf("hr:%lf\n", hr); 
    
    

    //printf("par=%lf,%lf\n", hs, ipar.dist);

    printf("s_layer%d\n", ipar.s_layer);
    printf("s_layer_h%lf\n", ipar.zs);
    printf("rlayer%d\n", ipar.r_layer);
    printf("rlayer_h%lf\n", ipar.zr);

    printf("START INT2 %d,%d,%lf,%lf,%lf,%lf\n", ipar.r_layer, ipar.s_layer, ipar.zs, ipar.zr, hs, hr);
    cplx wm; 
    cplx cs, cp; 
    cplx wt; 
    cplx ZI={0.0, 1.0}; 
    float32 kmax; 
    struct timespec rt1, rt2; 
    
    share_args sargs; 
    //#pragma omp parallel for
    //printf(", )
    clock_gettime(CLOCK_MONOTONIC,&rt1);

    t1 = clock(); 
    int_k_w *ikw;
    int nkw=0;  
    ipar.df = df; 
    ipar.img_freq = img_freq; 
    for(int nw=0;nw<n_half_sample;nw++){
        w = cplx(2 * PI * df * nw, - img_freq); 
        wm = cplx(2*PI, 0.0); 
        for(int i=0;i<n_layer;i++){
            cs = (1.0+log(w/PI2)/(PI*par[i].qs)+ZI/(2.0*par[i].qs)); 
            cp = (1.0+log(w/PI2)/(PI*par[i].qp)+ZI/(2.0*par[i].qp)); 
            par[i].a = (par[i].a0 * cp); 
            par[i].b = (par[i].b0 * cs);  
        }
 
        if(nw<3){
            wt = PI2 * df * 3.0 - img_freq * ZI; 
            kmax = cal_kmax(wt, par, ipar, nw); 
        }
        else{
            kmax = cal_kmax(w, par, ipar, nw); 
        }
        kmax = ipar.k_factor * kmax; 
        nkw = nkw + (int)(kmax/ipar.dk) + 2;  
        //if(nw==222)printf("NK:%d,%f\n", (int)(kmax/ipar.dk) + 2, kmax); 
    }
    ikw = (int_k_w *)malloc(sizeof(int_k_w)*nkw); 
    int cnt=0; 
    for(int nw=0;nw<n_half_sample;nw++){
        w = cplx(2 * PI * df * nw, - img_freq); 
        wm = cplx(2*PI, 0.0); 
        for(int i=0;i<n_layer;i++){
            cs = (1.0+log(w/PI2)/(PI*par[i].qs)+ZI/(2.0*par[i].qs)); 
            cp = (1.0+log(w/PI2)/(PI*par[i].qp)+ZI/(2.0*par[i].qp)); 
            par[i].a = (par[i].a0 * cp); 
            par[i].b = (par[i].b0 * cs);  
        }
 
        if(nw<3){
            wt = PI2 * df * 3.0 - img_freq * ZI; 
            kmax = cal_kmax(wt, par, ipar, nw); 
        }
        else{
            kmax = cal_kmax(w, par, ipar, nw); 
        }
        kmax = ipar.k_factor * kmax; 
        //if(nw==0)printf("%d,%lf,%lf\n", (int)(kmax/ipar.dk)+2, kmax, ipar.dk); 
        for(int ik=0; ik<(int)(kmax/ipar.dk) + 2;ik++){
            ikw[cnt].k = ipar.dk * ik; 
            ikw[cnt].nk = ik; 
            ikw[cnt].w = w;   
            ikw[cnt].nw = nw; 
            cnt ++; 
        }
    }
    int_out *out; 
    out = (int_out *)malloc(n_half_sample*2 * sizeof(int_out));
    integrate(ikw, nkw, n_half_sample, par, ipar, out); 
    t2 = clock(); 
    clock_gettime(CLOCK_MONOTONIC,&rt2);
    printf("C99测试CPU时间:%lf\n", (float)(t2-t1)/((float)CLOCKS_PER_SEC)); 
    float32 r = (rt2.tv_sec-rt1.tv_sec) * 1.0 + (rt2.tv_nsec - rt1.tv_nsec) * 1e-9; 
    printf("C99测试实际时间:%lf\n", r); 





    //cuDoubleComplex *h_signal; 
    cufftDoubleComplex *h_freq; 
    cufftDoubleComplex *d_freq; 
    cufftDoubleComplex *h_signal; 
    cufftDoubleComplex *d_signal; 
    cufftHandle plan;
    double *signal_out, **sig_out; 
    
    signal_out = (double *)malloc(sizeof(double)*n_half_sample*2*30);
    sig_out = (double **)malloc(sizeof(double *)*30); 
    for(int i=0;i<30;i++){
        sig_out[i] = &signal_out[i*n_half_sample*2]; 
    } 
    
    cufftPlan1d(&plan, n_half_sample*2, CUFFT_Z2Z, 1);
    h_freq = (cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*(n_half_sample*2)); 
    h_signal = (cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*(n_half_sample*2)); 
    CUDA_CHECK(cudaMalloc((void **)&d_freq,   (n_half_sample*2)    *   sizeof(cufftDoubleComplex)));
    CUDA_CHECK(cudaMalloc((void **)&d_signal,   (n_half_sample*2)    *   sizeof(cufftDoubleComplex)));
    for(int c=0;c<30;c++){
        for(int i=0;i<n_half_sample*2;i++){
            h_freq[i].x = real(out[i].sums[c]); 
            h_freq[i].y = imag(out[i].sums[c]); 
        }
        cudaMemcpy(d_freq, h_freq, n_half_sample * 2 * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice); 
        cufftExecZ2Z(plan, d_freq, d_signal, CUFFT_INVERSE); 
        cudaMemcpy(h_signal, d_signal, n_half_sample * 2 * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
        for(int jj=0;jj<n_half_sample*2;jj++){
            sig_out[c][jj] = h_signal[jj].x/((double)n_half_sample)/2.0; 
            sig_out[c][jj] = sig_out[c][jj] * exp(((double)jj)*ipar.dt*ipar.img_freq) * exp(ipar.t0*ipar.img_freq); 
        }
    }
    
    //printf("REAL:%lf,\n", h_signal[15]); 
    FILE *f = fopen(oname_t, "w"); 
    for(int jj=0;jj<n_half_sample*2;jj++){
        for(int c=0;c<30;c++){
            fprintf(f, "%.18e,", sig_out[c][jj]); 
        }
        fprintf(f, "\n"); 
    }
    fclose(f); 

    FILE *fp = fopen(oname_w, "w"); 
    for(int jj=0;jj<n_half_sample*2;jj++){
        for(int i=0;i<30;i++){
            fprintf(fp, "%.18e,%.18e,", real(out[jj].sums[i]), imag(out[jj].sums[i])); 
        }
        fprintf(fp, "\n"); 
        fflush(fp);
    }



    /*
    for(int nf=ipar.nf1;nf<ipar.nf2;nf++){
        for(i=0;i<7;i++){
            sums[i] = 0.0; 
        }        
        t3 = clock(); 
        integrate(nf, par, ipar, sums); 
        t4 = clock(); 
        if(nf%10==0)printf("NF=%d, T=%lfs\n", nf, (double)(t4-t3)/CLOCKS_PER_SEC); 
        for(i=0;i<7;i++){
            fprintf(fp, "%f,%f,", creal(sums[i]), cimag(sums[i])); 
        }
        fprintf(fp, "\n"); 
        fflush(fp); 
    }
    fclose(fp); 
    t2 = clock(); 
    printf("Time=%lfs\n", (double)(t2-t1)/CLOCKS_PER_SEC); 
    */
}

