//#include "../popm.h" 
//#include "../../scilib/mathematics.h"
/*
void integrate(cplx w, float32 kmax, layer_parameter *par, intpar ipar, cplx sums[7]){
    
    cplx wt;         
    float32 k; 
    int nk; 
    cplx intp[14]; 
    int i, j;  
    int count; 
    k = ipar.kmin; 
    */
    /*
    //int_func(k, w, par, ipar, intp); 
    for(i=0;i<7;i++){
        sums[i] = 0.0; 
    }
    if(ipar.kmin==0.0){
        for(i=0;i<7;i++){
            sums[i] = intp[i] * ipar.dk / 2.0; 
        }
    }
    else{
        for(i=0;i<7;i++){
            sums[i] = intp[i] * ipar.dk; 
        }
    }
    nk = (int)(kmax/ipar.dk) + 1; 
    for(i=1;i<=nk;i++){
        k = ipar.kmin + i * ipar.dk; 
        //int_func(k, w, par, ipar, intp);
        for(j=0;j<7;j++){
            sums[j] += intp[j] * ipar.dk; 
        }
    }
    */
    /*
    count = 0; 
    float32 rp[10][3]={0}, ip[10][3]={0}, kk[10][3]={0};
    float32 peaks[10][10][2]={0}, sss[10][2]={0}; 
    size_t pcount[10][2]={0}; 
    int flag[10][2]={0}; 
    int truesum=0; 
    while(1){
        k = k + ipar.dk; 
        int_func(k, w, par, ipar, intp);
        break; 
        for(j=0;j<7;j++){
            sums[j] += intp[j] * ipar.dk; 
            peak_trough_averaging_method_base(
                count, 
                sums[j], 
                rp[j], 
                ip[j], 
                kk[j], 
                sss[j], 
                flag[j]
            ); 
            if(flag[j][0]==1 && pcount[j][0]<10){
                peaks[j][pcount[j][0]][0] = sss[j][0]; 
                pcount[j][0] += 1; 
            }
            if(flag[j][1]==1 && pcount[j][1]<10){
                peaks[j][pcount[j][1]][1] = sss[j][1]; 
                pcount[j][1] += 1; 
            }
            //if(j==5)printf("%lf,%lf,%lf|%d\n", real(rp[j][0]), real(rp[j][1]), real(rp[j][2]), ABS(rp[j][0]-rp[j][1])<EPSILON); 
        }
        truesum = 0; 

        if (truesum>=14*10)break; 
        count += 1; 
    }
    for(i=9; i>-1;i--){
        for(j=0;j<i;j++){
            for(int k=0;k<7;k++){
                for(int m=0;m<2;m++){
                    peaks[k][j][m] = (peaks[k][j][m]+peaks[k][j+1][m])/2.; 
                }
            }            
        }
    }
    for(j=0;j<7;j++){
        //sums[j] = peaks[j][0][0] + peaks[j][0][1] * 1.I;
    }    
}
*/