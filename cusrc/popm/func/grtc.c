/*
 *   popm.h
 *
 *  Created on: 2022-12-03
 *      Author: Yu Ziye 
 *  Used to calculate the GRTM and modified RTM. 
 */

#include "../popm.h" 
#include <stdlib.h>
grtc *grtm(float32 k, cplx w, parameter *par, int n_layer, int s_layer){
    int lay; 
    cplx wa, wb; 
    cplx E1[4][4], E2[4][4], a[2][2], b[2][2]; 
    cplx exu[2]; 
    int i, j; 
    mrtc *mrt, tmrt; 
    grtc *grt; 
    cplx tem1, tem2, tem3; 
    mrt = (mrtc *)malloc(sizeof(mrtc)*(n_layer-1)); 
    grt = (grtc *)malloc(sizeof(grtc)*(n_layer)); 

    //last layer has similiar format. 
    for(lay=0;lay<n_layer-1;lay++){
        mat_PSV_E(par[lay], k, w, E1); 
        mat_PSV_E(par[lay+1], k, w, E2); 
        mat_PSV_RT(E1, E2); 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                mrt[lay].tdd1[i][j] = E2[i][j]; 
                mrt[lay].tuu1[i][j] = E2[i+2][j+2]; 
                mrt[lay].rdu1[i][j] = E2[i+2][j]; 
                mrt[lay].rud1[i][j] = E2[i][j+2]; 
            }
        }
        wa = par[lay].v * par[lay].u; 
        wb = par[lay+1].v * par[lay+1].u; 
        mrt[lay].tdd0 = 2. * wa / (wa + wb); 
        mrt[lay].rdu0 = (wa-wb) / (wa + wb); 
        mrt[lay].rud0 = - mrt[lay].rdu0; 
        mrt[lay].tuu0 = 2. * wb / (wa + wb);   
    }
    

    // free surface condition 
    mat_PSV_E(par[0], k, w, E1); 
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            a[i][j] = E1[i+2][j]; 
            b[i][j] = E1[i+2][j+2]; 
        }
    }
    //cplx_minvc_2(a); 
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            grt[0].rud1[i][j] = -(a[i][0]*b[0][j] + a[i][1]*b[1][j]);
        }
    }
    grt[0].rud0 = 1.0; 
    //call grt coeff 
    for(lay=0;lay<s_layer;lay++){
        exu[0] = exp(-par[lay].r*par[lay].h); 
        exu[1] = exp(-par[lay].v*par[lay].h); 
        tmrt = mrt[lay]; 

        tem1 = grt[lay].rud0 * exu[1]; 
        tem2 = tmrt.rdu0 * exu[1]; 
        tem3 = tmrt.tdd0 * exu[1]; 
        
        //SH
        grt[lay+1].tu0 = tmrt.tuu0 / (1.0 - tem1 * tem2);
        grt[lay+1].rud0 = tmrt.rud0 + tem3 * tem1 * grt[lay+1].tu0; 
         
        // PSV
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                a[i][j] = -(tmrt.rdu1[i][0] * exu[0] * grt[lay].rud1[0][j] 
                           +tmrt.rdu1[i][1] * exu[1] * grt[lay].rud1[1][j]) * exu[j]; 
                b[i][j] =  (tmrt.tdd1[i][0] * exu[0] * grt[lay].rud1[0][j] 
                           +tmrt.tdd1[i][1] * exu[1] * grt[lay].rud1[1][j]) * exu[j];  
            }
            a[i][i] = 1.0 + a[i][i]; 
        }
        //cplx_minvc_2(a); 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[lay+1].tu1[i][j] = a[i][0] * tmrt.tuu1[0][j] + a[i][1] * tmrt.tuu1[1][j]; 
                
            }
        }      
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[lay+1].rud1[i][j] = tmrt.rud1[i][j] + b[i][0] * grt[lay+1].tu1[0][j] + b[i][1] * grt[lay+1].tu1[1][j]; 
            }
        }    
    }
    // 半无限空间
    tmrt = mrt[n_layer-2]; 
    grt[n_layer-1].rdu0 = tmrt.rdu0; 
    grt[n_layer-1].td0  = tmrt.tdd0; 
    
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            grt[n_layer-1].td1[i][j] = tmrt.tdd1[i][j]; 
            grt[n_layer-1].rdu1[i][j] = tmrt.rdu1[i][j]; 
        }
    }        
    for(lay=n_layer-3;lay>=s_layer;lay--){
        exu[0] = exp(-par[lay+1].r*par[lay+1].h); 
        exu[1] = exp(-par[lay+1].v*par[lay+1].h); 
        tmrt = mrt[lay];

        tem1 = tmrt.rud0 * exu[1]; 
        tem2 = grt[lay+2].rdu0 * exu[1]; 
        tem3 = tmrt.tuu0 * exu[1]; 

        
        //SH
        grt[lay+1].td0 = tmrt.tdd0 / (1.0 - tem1 * tem2); 
        grt[lay+1].rdu0 = tmrt.rdu0 + tem3 * tem2 * grt[lay+1].td0; 
        // PSV
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                a[i][j] = -(tmrt.rud1[i][0] * exu[0] * grt[lay+2].rdu1[0][j] 
                           +tmrt.rud1[i][1] * exu[1] * grt[lay+2].rdu1[1][j]) * exu[j]; 
                b[i][j] =  (tmrt.tuu1[i][0] * exu[0] * grt[lay+2].rdu1[0][j] 
                           +tmrt.tuu1[i][1] * exu[1] * grt[lay+2].rdu1[1][j]) * exu[j];  
            }
            a[i][i] = 1. + a[i][i]; 
        }


        //cplx_minvc_2(a); 
        
        
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[lay+1].td1[i][j] = a[i][0] * tmrt.tdd1[0][j] + a[i][1] * tmrt.tdd1[1][j]; 
            }
        }      
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                grt[lay+1].rdu1[i][j] = tmrt.rdu1[i][j] + b[i][0] * grt[lay+1].td1[0][j] + b[i][1] * grt[lay+1].td1[1][j]; 
            }
        }    
    }
    free(mrt); 
    return grt; 
}


poppar ypars(float32 k, cplx w, parameter *par, int n_layer, int s_layer, float32 zs, int r_layer, float32 zr, grtc *grt){
    cplx sexd[2], sexu[2], rexd[2], rexu[2];
    cplx tmp1[2][2], tmp2[2][2]; 
    cplx iden[2][2]; 
    cplx E[4][4]; 
    int i, j, lay; 
    cplx yu_psv[2][2], yd_psv[2][2], y_sh, yu0, yd0;
    cplx u[2][2], d[2][2], u0, d0, yu1[2][2], yd1[2][2]; 
    cplx ex[2]; 
    poppar ppar; 
    //get grt matrix, need to be freed 
    //grt = grtm(k, w, par, n_layer, s_layer); 
    //Identify matrix
    iden[0][0] = iden[1][1] = 1.0;
    iden[1][0] = iden[0][1] = 0.0; 

    //source exp pars. 
    sexd[0] = exp(-par[s_layer].r * (par[s_layer].h - zs)); 
    sexd[1] = exp(-par[s_layer].v * (par[s_layer].h - zs)); 
    sexu[0] = exp(-par[s_layer].r * (zs)); 
    sexu[1] = exp(-par[s_layer].v * (zs)); 
    
    
    rexd[0] = exp(-par[r_layer].r * zr); 
    rexd[1] = exp(-par[r_layer].v * zr); 


    if ((n_layer-1) == r_layer){
        rexu[0] = 0.0; 
        rexu[1] = 0.0; 
    }
    else{
        rexu[0] = exp(-par[r_layer].r * (par[r_layer].h - zr)); 
        rexu[1] = exp(-par[r_layer].v * (par[r_layer].h - zr));         
    }

    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            tmp1[i][j] = sexd[i] * grt[s_layer+1].rdu1[i][j] * sexd[j];
            tmp2[i][j] = sexu[i] * grt[s_layer-1+1].rud1[i][j] * sexu[j]; 
        }
    }
    //A12-cd:ERER 
    //cplx_mm_2(tmp1, tmp2, yu_psv); 
    //cplx_mm_2(tmp2, tmp1, yd_psv); 


    //A12-cd:1-ERER
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            yu_psv[i][j] = iden[i][j] - yu_psv[i][j]; 
            yd_psv[i][j] = iden[i][j] - yd_psv[i][j]; 
        }
    }
    //A12-cd:(1-ERER)^-1
    //cplx_minvc_2(yu_psv); 
    //cplx_minvc_2(yd_psv); 





    y_sh = 1./(1.-sexd[1]*grt[s_layer+1].rdu0*sexd[1]*sexu[1]*grt[s_layer].rud0*sexu[1]);  

    mat_PSV_E(par[r_layer], k, w, E); 
    if((s_layer==r_layer)&&(zr<=zs)){
        rexd[0] = exp(-par[r_layer].r * (zr)); 
        rexd[1] = exp(-par[r_layer].v * (zr)); 
        rexu[0] = exp(-par[r_layer].r * (zs-zr)); 
        rexu[1] = exp(-par[r_layer].v * (zs-zr));     
        u0 = rexd[1] * grt[s_layer-1+1].rud0 * sexu[1] + rexu[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[s_layer-1+1].rud1[0][j] * sexu[j] 
                        + E[i][1] * rexd[1] * grt[s_layer-1+1].rud1[1][j] * sexu[j] 
                        + E[i][j+2] * rexu[j]; 
            }
        }

        yu0 = u0 * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = u[i][0] * yu_psv[0][j] + u[i][1] * yu_psv[1][j]; 
            }
        }
    }
    else if ((s_layer==r_layer)&&(zr>zs))
    {
        rexd[0] = exp(-par[s_layer].r * (zr-zs)); 
        rexd[1] = exp(-par[s_layer].v * (zr-zs)); 
        rexu[0] = exp(-par[s_layer].r * (par[s_layer].h-zr)); 
        rexu[1] = exp(-par[s_layer].v * (par[s_layer].h-zr)); 
        d0 = rexd[1] + rexu[1] * grt[s_layer+1].rdu0 * sexd[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[s_layer+1].rdu1[0][j] * sexd[j] 
                        + E[i][3] * rexu[1] * grt[s_layer+1].rdu1[1][j] * sexd[j]; 
            }
        }
        yd0 = d0 * y_sh; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = d[i][0] * yd_psv[0][j] + d[i][1] * yd_psv[1][j]; 
            }
        }         
    }
    else if(r_layer<s_layer){
        ex[0] = exp(-par[r_layer].r*par[r_layer].h); 
        ex[1] = exp(-par[r_layer].v*par[r_layer].h); 
        u0 = rexd[1] * grt[r_layer].rud0 * ex[1] + rexu[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                u[i][j] = E[i][0] * rexd[0] * grt[r_layer-1+1].rud1[0][j] * ex[j] 
                        + E[i][1] * rexd[1] * grt[r_layer-1+1].rud1[1][j] * ex[j] 
                        + E[i][j+2] * rexu[j]; 
            }
        }
        
        for(lay=s_layer-1;lay>=r_layer;lay--){
            if(lay==(s_layer-1)){
                yu0 = grt[lay+1].tu0 * sexu[1] * y_sh; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yu1[i][j] = grt[lay+1].tu1[i][0] * sexu[0] * yu_psv[0][j]
                                  + grt[lay+1].tu1[i][1] * sexu[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[lay+1].r * par[lay+1].h); 
                ex[1] = exp(-par[lay+1].v * par[lay+1].h); 
                yu0 = grt[lay+1].tu0 * ex[1] * yu0; 

                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[lay+1].tu1[i][0] * ex[0] * yu1[0][j]
                                   + grt[lay+1].tu1[i][1] * ex[1] * yu1[1][j]; 
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
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = u[i][0] * yu1[0][j] + u[i][1] * yu1[1][j]; 
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yu1[i][j] = tmp1[i][j]; 
            }
        }                    
    }
    else if (r_layer>s_layer)
    {
        ex[0] = exp(-par[r_layer-1].r*par[r_layer-1].h); 
        ex[1] = exp(-par[r_layer-1].v*par[r_layer-1].h); 
        d0 = rexu[1] * grt[r_layer+1].rdu0 * ex[1] + rexd[1]; 
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                d[i][j] = E[i][j] * rexd[j]
                        + E[i][2] * rexu[0] * grt[r_layer+1].rdu1[0][j] * ex[j] 
                        + E[i][3] * rexu[1] * grt[r_layer+1].rdu1[1][j] * ex[j]; 
            }
        }
        for(lay=s_layer+1;lay<=r_layer;lay++){
            if(lay==(s_layer+1)){
                yd0 = grt[lay+1].td0 * sexd[1] * y_sh; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        yd1[i][j] = grt[lay+1].td1[i][0] * sexd[0] * yu_psv[0][j]
                                  + grt[lay+1].td1[i][1] * sexd[1] * yu_psv[1][j]; 
                    }
                }                
            }
            else{
                ex[0] = exp(-par[lay-2].r * par[lay-2].h); 
                ex[1] = exp(-par[lay-2].v * par[lay-2].h); 
                yd0 = grt[lay+1].td0 * ex[1] * yd0; 
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        tmp1[i][j] = grt[lay+1].td1[i][0] * ex[0] * yd1[0][j]
                                   + grt[lay+1].td1[i][1] * ex[1] * yd1[1][j]; 
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
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                tmp1[i][j] = d[i][0] * yd1[0][j] + d[i][1] * yd1[1][j]; 
            }
        }
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                yd1[i][j] = tmp1[i][j]; 
            }
        }      
    }
    ppar.yu0 = yu0; 
    ppar.yd0 = yd0; 
    for(i=0;i<2;i++){
        for(j=0;j<2;j++){
            ppar.yd1[i][j] = yd1[i][j]; 
            ppar.yu1[i][j] = yu1[i][j]; 
        }
    }      
    return ppar;  
}
