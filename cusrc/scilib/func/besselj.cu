/*
 *  Created on: 30/03/2015
 *      Author: Cangye
 */
#include <math.h>
#include "../mathematics.h"
/*
float32 besselj0(float32 x)
{
    float32 ax,x1,x2,theta,fct;
    float32 besselj0;
    float32 a[7]={1.00000000,-2.24999970,1.26562080,-0.31638660,0.04444790,-0.00394440,0.00021000};
    float32 b[7]={0.79788456,-0.00000077,-0.00552740,-0.00009512,0.00137237,-0.00072805,0.00014476};
    float32 c[7]={-0.78539816,-0.04166397,-0.00003954, 0.00262573,-0.00054125,-0.00029333,0.00013558};
    ax=fabs(x);
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

float32 besselj1(float32 x)
{
    float32 ax,x1,x2,theta,fct;
    float32 besselj1;
    float32 sign;
    float32 a[7]={0.50000000,-0.56249985,0.21093573,-0.03954289,0.00443319,-0.00031761,0.00001109};
    float32 b[7]={0.79788456,0.00000156,0.01659667,0.00017105,-0.00249511,0.00113653,-0.00020033};
    float32 c[7]={-2.35619449,0.12499612,0.00005650,-0.00637879,0.00074348,0.00079824,-0.00029166};
    ax=fabs(x);
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

float32 besselj(float32 x, int n)
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
          ax=fabs(x);
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
float32 d_besselj(float32 x,int n)
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
*/