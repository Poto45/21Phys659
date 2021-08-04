#include "tremaine_roz.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

tremaine_roz::tremaine_roz()
{
}

tremaine_roz::~tremaine_roz()
{
}

void tremaine_roz::setStructure(double a, double b, double lz, double eps)
{
 _a=a;
 _b=b;
 _lz=lz;
 _eps=eps;
}

void tremaine_roz::setStructure(double a, double b, double lx, double lz, double eps)
{
 _a=a;
 _b=b;
 _lx=lx;
 _lz=lz;
 _eps=eps;
}

void tremaine_roz::setNModes(int n)
{
 _n=n;
}

void tremaine_roz::setNEigenX(int nefx)
{
 _nefx=nefx;
}

void tremaine_roz::calcEigenFreqsInf()
{
 double f1,f2,cf,x1,x2,cx,approx;
 int m;
// solve eqn.: x*tan(x)=c
 for(int i=0; i<_n; ++i)
   {
    x1=i*PI;
    x2=x1+PI/2.;
    while(fabs(x2-x1)>SMALL)
      {
       cx=(x2+x1)/2.;
       f1=calcFunc1(x1);
       cf=calcFunc1(cx);
//       cout<<"x1: "<<x1<<" cx: "<<cx<<" x2: "<<x2<<" m: "<<m<<endl;
//       cout<<"f1: "<<f1<<" cf: "<<cf<<endl;
       if(f1*cf<0)
         {
          x2=cx;
         }
       else
         {
          x1=cx;
         }
      }
    _freqs[i]=cx*CSPEED/sqrt(_eps-1)/(_b-_a)/2./PI*CONV2GHZ;
    approx=(i+0.5)*PI*CSPEED/sqrt(_eps-1)/(_b-_a)/2./PI*CONV2GHZ;
//    _freqs[i]=approx;
    _lambdas[i]=CSPEED/_freqs[i]*CONV2GHZ;
    cout<<"i: "<<i<<" f (GHz): "<<_freqs[i]<<" approx: "<<approx<<" l (mm): "<<_lambdas[i]<<endl;
//    cout<<"i: "<<i<<" f (GHz): "<<_freqs[i]<<" l (mm): "<<_lambdas[i]<<endl;
   }
// cout<<_freqs[0]<<" "<<_freqs[1]<<endl;
}

void tremaine_roz::calcAmplitudesInf()
{
 double A;
 for(int i=0; i<_n; ++i)
   {
    A=sqrt(_eps-1)*(_b-_a)*2*PI*_freqs[i]/CONV2GHZ/CSPEED;
    A=1.0/sin(A);
    _amplit[i]=4*PI/(_a+_eps*A*A*(_b-_a));
    cout<<"i: "<<i<<" amplit: "<<_amplit[i]<<endl;
   }
}

void tremaine_roz::calcWakefieldsDelta(double zmin, double zmax, double z0)
{
 double z,dz,wake;
 dz=(zmax-zmin)/NDIV;
 for(int i=0; i<NDIV; ++i)
   {
    z=zmin+i*dz;
    wake=0;
    if(z<=z0)
      {
       for(int j=0; j<_n; ++j)
         {
          wake=wake+cos(2*PI*_freqs[j]/CONV2GHZ/CSPEED*(z-z0))*_amplit[j];
         }
      }
    _wakeZ[i]=wake;
    _z[i]=z;
   }
}

void tremaine_roz::calcWakefieldsConv(double zmin, double zmax, double sigz, double z0)
{
// zmax should not be >Lz
// zmin should not be <0 or >Lz
 if(zmax>_lz || zmin<0 || zmin>_lz)
   {
    cout<<"ERROR: choose 0 < zmin,max < Lz "<<endl;
    exit(0);
   }
 double z,dz,wake,zp,dzp,zpmin,zmaxint;
 int ni;
 dz=(zmax-zmin)/NDIV;
 zmaxint=z0+NSIGZ*sigz;
 if(zmaxint>_lz)
   zmaxint=_lz;
 for(int i=0; i<NDIV; ++i)
   {
    z=zmin+i*dz;
    wake=0;
    for(int k=0; k<_n; ++k)
      {
       ni=int((zmaxint-z)/dz)+1;
       for(int j=0; j<ni; ++j)
         {
          zp=z+j*dz;
          wake=wake+calcFunc2(zp,sigz,z0)*cos(2*PI*_freqs[k]/CONV2GHZ/CSPEED*
               (zp-z))*_amplit[k]*dz;
         }
      }
    _wakeZ[i]=wake;
    _z[i]=z;
   }
 _nzOut=NDIV;
}

void tremaine_roz::overlapWakefields(double zmin, double zmax, double sigz, 
                                     double z0, double delz, int nb)
{
 double z,dz,wake,zp,dzp,zpmin,zmaxint;
 int ni;
 dz=(zmax-zmin)/NDIV;
 zmaxint=z0+NSIGZ*sigz;
 if(zmaxint>_lz)
   zmaxint=_lz;
 for(int i=0; i<NDIV; ++i)
   {
    z=zmin+i*dz;
    wake=0;
    for(int k=0; k<_n; ++k)
      {
       ni=int((zmaxint-z)/dz)+1;
       for(int j=0; j<ni; ++j)
         {
          zp=z+j*dz;
          wake=wake+calcFunc3(zp,sigz,z0,delz,nb)*cos(2*PI*_freqs[k]/
               CONV2GHZ/CSPEED*(zp-z))*_amplit[k]*dz;
         }
      }
    _wakeZ[i]=wake;
    _z[i]=z;
   }
}

void tremaine_roz::normWakefieldsToGeV()
{
 for(int i=0; i<NDIV; ++i)
   {
    _wakeZ[i]=_wakeZ[i]/1000; // convert to GV/m from MV/m
    _wakeEX[i]=_wakeEX[i]/1000; // convert to GV/m from MV/m
    _wakeEY[i]=_wakeEY[i]/1000; // convert to GV/m from MV/m
   }
   
}

void tremaine_roz::convZtoMicrons()
{
 for(int i=0; i<NDIV; ++i)
   {
    _z[i]=_z[i]*1000; 
   }
}

void tremaine_roz::outputWakefields1D(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_nzOut; ++i)
    {
     fout<<_z[i]<<" "<<_wakeZ[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputWakefields1Dall(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_nzOut; ++i)
    {
     fout<<_z[i]<<" "<<_wakeEX[i]<<" "<<_wakeEY[i]<<" "<<_wakeZ[i]<<" "
           <<_wakeBX[i]<<" "<<_wakeBY[i]<<" "<<_wakeBZ[i]<<
           " "<<_FX[i]<<" "<<_FY[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputRhox(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_nx; ++i)
    {
     fout<<_x[i]<<" "<<_rhox[i]<<" "<<_recRhox[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputRhoy(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_ny; ++i)
    {
     fout<<_y[i]<<" "<<_rhoy[i]<<" "<<_recRhoy[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputRhoz(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_nz; ++i)
    {
     fout<<_zb[i]<<" "<<_rhoz[i]<<" "<<_recRhoz[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputRhoxy(char fname[80])
{
 fstream fout;
 fout.open(fname,ios::out);
 fout<<_nx<<" "<<_ny<<endl;
  for(int i=0; i<_nx; ++i)
    {
     for(int j=0; j<_ny; ++j)
       {
        fout<<_x[i]<<" "<<_y[j]<<" "<<_rhoxy[i][j]<<" "<<_recRhoXY[i][j]<<endl;
       }
    }
 fout.close();
}

void tremaine_roz::setQInf(double q)
{
 for(int i=0; i<NDIV; ++i)
   {
    _wakeZ[i]=_wakeZ[i]*9.0*q; // 9*10^9 -> 9.0 due to nC unit. results are in MV/m
// "-" sign is due to the fact that wakes decelerate the charge
   }
}

void tremaine_roz::setQ3D(double q)
{
 _q3d=q*1.0e-9;
}


void tremaine_roz::calcEigenFreqs3DSym()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,xas,_kx,fac;
 int m,found,parity,check,ng;
 ng=int(_n/_nefx);
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: tan(x)=c1*x/(c3*x*x-c2)
 c1=cosh(_kx*_a)/sinh(_kx*_a)/_kx/(_b-_a);
 c2=1.-1.0/_eps; // T-R
// c2=sinh(kx*(_b-_a))/cosh(kx*(_b-_a));
 c3=1./_eps/_kx/_kx/(_b-_a)/(_b-_a);
// cout<<"c's: "<<c1<<" "<<c2<<" "<<c3<<endl;
 xas=sqrt(c2/c3);
 found=0;
// cout<<"Xasym: "<<sqrt(c2/c3)<<endl;
 for(int i=0; found<ng; ++i)
   {
    x1=i*PI/2;
    x2=(i+1)*PI/2;
    if(i-2*int(i/2)==0)
      {
       parity=1; // i is even
      }
    else
      {
       parity=-1; // i is odd
      }
    if(xas<x1)
      {
       if(parity==-1)
         x2=x1;
      }
    else if(xas>x1 && xas<x2)
      {
       if(parity==1)
         {
          x1=xas+EPS;
         }
       else
         {
          x2=xas-EPS;
         }
      }
    else
      {
       if(parity==1)
         {
          x2=x1;
         }
      }
//       cout<<"x1: "<<x1<<" x2: "<<x2<<endl;
       while(fabs(x2-x1)>SMALL)
         {
          cx=(x2+x1)/2.;
          f1=calcFunc5(x1,c1,c2,c3);
          cf=calcFunc5(cx,c1,c2,c3);
//       cout<<"x1: "<<x1<<" cx: "<<cx<<" x2: "<<x2<<" m: "<<m<<endl;
//       cout<<"f1: "<<f1<<" cf: "<<cf<<endl;
          if(f1*cf<0)
            {
             x2=cx;
            }
          else
            {
             x1=cx;
            }
         }
//       cout<<"i: "<<i*PI/2<<" "<<(i+1)*PI/2<<" x: "<<cx<<" kx: "<<_kx<<endl;
//       cout<<"x1: "<<x1<<" x2: "<<x2<<endl;
       if(x1!=x2)
        {
         _freqsSym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasSym[l*ng+found]=CSPEED/_freqsSym[l*ng+found]*CONV2GHZ;
//         cout<<"found: "<<found<<" f (GHz): "<<_freqsSym[found]<<" l (mm): "<<_lambdasSym[found]<<endl;
//         cout<<"x: "<<cx<<" i: "<<i<<endl;
//    cout<<"i: "<<i<<" f (GHz): "<<_freqsSym[found]<<" l (mm): "<<_lambdasSym[found]<<endl;
//       cout<<"f: "<<_freqsSym[l*ng+found]<<" x: "<<cx<<endl;
//       cout<<"TR (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<endl;
       ++found;
        }
   } // end "i"
  } // end "l";
// cout<<_freqsSym[0]<<" "<<_freqsSym[1]<<endl;
}

void tremaine_roz::calcEigenFreqs3DAsym()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,xas,_kx,fac;
 int m,found,parity,check,ng;
 ng=_n/_nefx;
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: tan(x)=c1*x/(c3*x*x-c2)
 c1=sinh(_kx*_a)/cosh(_kx*_a)/_kx/(_b-_a);
 c2=1.-1.0/_eps/fac;
// c2=sinh(kx*(_b-_a))/cosh(kx*(_b-_a));
 c3=1./_eps/_kx/_kx/(_b-_a)/(_b-_a);
// cout<<"c's: "<<c1<<" "<<c2<<" "<<c3<<endl;
 xas=sqrt(c2/c3);
 found=0;
// cout<<"Xasym: "<<sqrt(c2/c3)<<endl;
 for(int i=0; found<ng; ++i)
   {
    x1=i*PI/2;
    x2=(i+1)*PI/2;
    if(i-2*int(i/2)==0)
      {
       parity=1; // i is even
      }
    else
      {
       parity=-1; // i is odd
      }
    if(xas<x1)
      {
       if(parity==-1)
         x2=x1;
      }
    else if(xas>x1 && xas<x2)
      {
       if(parity==1)
         {
          x1=xas+EPS;
         }
       else
         {
          x2=xas-EPS;
         }
      }
    else
      {
       if(parity==1)
         {
          x2=x1;
         }
      }
//       cout<<"x1: "<<x1<<" x2: "<<x2<<endl;
       while(fabs(x2-x1)>SMALL)
         {
          cx=(x2+x1)/2.;
          f1=calcFunc5(x1,c1,c2,c3);
          cf=calcFunc5(cx,c1,c2,c3);
//       cout<<"x1: "<<x1<<" cx: "<<cx<<" x2: "<<x2<<" m: "<<m<<endl;
//       cout<<"f1: "<<f1<<" cf: "<<cf<<endl;
          if(f1*cf<0)
            {
             x2=cx;
            }
          else
            {
             x1=cx;
            }
         }
       if(x1!=x2)
        {
         _freqsAsym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                           _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasAsym[l*ng+found]=CSPEED/_freqsAsym[l*ng+found]*CONV2GHZ;
//         cout<<"found: "<<found<<" f (GHz): "<<_freqsAsym[found]<<" l (mm): "<<_lambdasAsym[found]<<endl;
//         cout<<"x: "<<cx<<" i: "<<i<<endl;
//    cout<<"i: "<<i<<" f (GHz): "<<_freqsAsym[found]<<" l (mm): "<<_lambdasAsym[found]<<endl;
       ++found;
        }
   }
  } // end "l"
}

void tremaine_roz::calcEigenFreqs3DSymDM()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac,dx,sf;
 int m,found,parity,check,ng,check1,check2,crit;
 ng=int(_n/_nefx);
// dx=PI/1000;
 dx=PI/2;
// crit=int(PI/2/dx);
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*ctan(x)-c2*tan(x)=(c3*x^2+c4)/c5/x
   c1=_eps/tanh(_kx*_a);
   c2=tanh(_kx*_a);
   c3=1.0/_kx/_kx/(_b-_a)/(_b-_a); // approx. DM
   c4=-_eps;
   c5=1.0/_kx/(_b-_a);
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*dx;
      x2=(i+1)*dx;
//      check1=i-int(i/crit)*crit;
//      check2=i+1-int((i+1)/crit)*crit;
//      if(check1==0)
//        x1=x1+EPS;
//      if(check2==0)
//        x2=x1-EPS;
//      f1=calcFunc10(x1,c1,c2,c3,c4,c5);
//      f2=calcFunc10(x2,c1,c2,c3,c4,c5);
//      sf=0;
//      if(f1*f2<0)
//        {
//         sf=1;
         while(fabs(x2-x1)>SMALL)
           {
            cx=(x2+x1)/2.;
            f1=calcFunc10(x1,c1,c2,c3,c4,c5);
            cf=calcFunc10(cx,c1,c2,c3,c4,c5);
            if(f1*cf<0)
              {
               x2=cx;
              }
            else
              {
               x1=cx;
              }
            }
//         }
//      if(sf>0)
//        {
         _freqsSym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasSym[l*ng+found]=CSPEED/_freqsSym[l*ng+found]*CONV2GHZ;
//         cout<<"DM (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<" i: "<<i<<endl;
         ++found;
//        }
     } // end "i"
  } // end "l";
// cout<<_freqsSym[0]<<" "<<_freqsSym[1]<<endl;
}

void tremaine_roz::calcEigenFreqs3DAsymDM()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac;
 int m,found,parity,check,ng;
 ng=int(_n/_nefx);
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*ctan(x)-c2*tan(x)=(c3*x^2+c4)/c5
   c1=_eps*tanh(_kx*_a);
   c2=1.0/tanh(_kx*_a);
   c3=1.0/_kx/_kx/(_b-_a)/(_b-_a); // approx. DM
   c4=-_eps;
   c5=1.0/_kx/(_b-_a);
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*PI/2;
      x2=(i+1)*PI/2;
      while(fabs(x2-x1)>SMALL)
        {
         cx=(x2+x1)/2.;
         f1=calcFunc10(x1,c1,c2,c3,c4,c5);
         cf=calcFunc10(cx,c1,c2,c3,c4,c5);
         if(f1*cf<0)
           {
            x2=cx;
           }
         else
           {
            x1=cx;
           }
         }
      _freqsAsym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                       _kx*_kx)/2./PI*CONV2GHZ;
      _lambdasAsym[l*ng+found]=CSPEED/_freqsAsym[l*ng+found]*CONV2GHZ;
//      cout<<"DM (A): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<endl;
      ++found;
     } // end "i"
  } // end "l";
// cout<<_freqsSym[0]<<" "<<_freqsSym[1]<<endl;
}

void tremaine_roz::calcEigenFreqs3DSymDM2()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac,dx,sf;
 int m,found,parity,check,ng,check1,check2,crit;
 ng=int(_n/_nefx);
 dx=PI/2;
 parity=0;
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*cot(x)=c2*x (LSM)
// solve eqn.: c1*cot(x)=-c3/x (LSE)
   c1=1.0/tanh(_kx*_a);
   c2=1.0/_kx/(_b-_a)/_eps; // approx. DM
   c3=_kx*(_b-_a);
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*dx;
      x2=(i+1)*dx;
      x1=x1+EPS;
      x2=x2-EPS;
      parity=int(i/2);
      if(parity*2-i==0)
       {
        parity=0;
       }
      else
       {
        parity=1;
       }
         while(fabs(x2-x1)>SMALL)
           {
            cx=(x2+x1)/2.;
            if(parity==0) // LSM
             {
              f1=calcFunc11(x1,c1,c2);
              cf=calcFunc11(cx,c1,c2);
             }
            else // LSE
             {
              f1=calcFunc12(x1,c1,c3);
              cf=calcFunc12(cx,c1,c3);
             }
            if(f1*cf<0)
              {
               x2=cx;
              }
            else
              {
               x1=cx;
              }
            }
         _freqsSym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasSym[l*ng+found]=CSPEED/_freqsSym[l*ng+found]*CONV2GHZ;
//         cout<<"DM2 (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<" i: "<<i<<endl;
//         cout<<"f_mn [GHz]: "<<_freqsSym[l*ng+found]<<" parity: "<<parity<<" m: "<<l<<" kx: "<<_kx<<" n: "<<i<<endl;
         ++found;
//        }
     } // end "i"
  } // end "l";
// cout<<_freqsSym[0]<<" "<<_freqsSym[1]<<endl;
}

void tremaine_roz::calcEigenFreqs3DAsymDM2()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac,dx,sf;
 int m,found,parity,check,ng,check1,check2,crit;
 ng=int(_n/_nefx);
 dx=PI/2;
 parity=0;
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*cot(x)=c2*x (LSM)
// solve eqn.: c1*cot(x)=-c3/x (LSE)
   c1=tanh(_kx*_a);
   c2=1.0/_kx/(_b-_a)/_eps; // approx. DM
   c3=_kx*(_b-_a);
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*dx;
      x2=(i+1)*dx;
      parity=int(i/2);
      if(parity*2-i==0)
       {
        parity=0;
       }
      else
       {
        parity=1;
       }
         while(fabs(x2-x1)>SMALL)
           {
            cx=(x2+x1)/2.;
            if(parity==0) // LSM
             {
              f1=calcFunc11(x1,c1,c2);
              cf=calcFunc11(cx,c1,c2);
             }
            else // LSE
             {
              f1=calcFunc12(x1,c1,c3);
              cf=calcFunc12(cx,c1,c3);
             }
            if(f1*cf<0)
              {
               x2=cx;
              }
            else
              {
               x1=cx;
              }
            }
         _freqsAsym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasAsym[l*ng+found]=CSPEED/_freqsAsym[l*ng+found]*CONV2GHZ;
//         cout<<"DM (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<endl;
         ++found;
//        }
     } // end "i"
  } // end "l";
// cout<<_freqsSym[0]<<" "<<_freqsSym[1]<<endl;
}

void tremaine_roz::calcEigenFreqs3DSymDM2VpOut(double flaser,char fname[80])
{
 fstream fout;
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac,dx,sf;
 int m,found,parity,check,ng,check1,check2,crit;
 double omega,fcut,vph,ktrans2;
 fout.open(fname,ios::out);
 omega=2.*PI*flaser*1e12;
 ng=int(_n/_nefx);
 dx=PI/2;
 cout<<" f laser [THz] = "<<flaser<<endl;
 cout<<"------------------"<<endl;
 parity=0;
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*cot(x)=c2*x (LSM)
// solve eqn.: c1*cot(x)=-c3/x (LSE)
   c1=1.0/tanh(_kx*_a);
   c2=1.0/_kx/(_b-_a)/_eps; // approx. DM
   c3=_kx*(_b-_a);
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*dx;
      x2=(i+1)*dx;
      x1=x1+EPS;
      x2=x2-EPS;
      parity=int(i/2);
      if(parity*2-i==0)
       {
        parity=0;
       }
      else
       {
        parity=1;
       }
         while(fabs(x2-x1)>SMALL)
           {
            cx=(x2+x1)/2.;
            if(parity==0) // LSM
             {
              f1=calcFunc11(x1,c1,c2);
              cf=calcFunc11(cx,c1,c2);
             }
            else // LSE
             {
              f1=calcFunc12(x1,c1,c3);
              cf=calcFunc12(cx,c1,c3);
             }
            if(f1*cf<0)
              {
               x2=cx;
              }
            else
              {
               x1=cx;
              }
            }
         _freqsSym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasSym[l*ng+found]=CSPEED/_freqsSym[l*ng+found]*CONV2GHZ;
//         cout<<"DM2 (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<" i: "<<i<<endl;
         ktrans2=cx*cx/(_b-_a)/(_b-_a)+_kx*_kx;
         fcut=CSPEED*sqrt(ktrans2)/sqrt(_eps)/2./PI/1e12; // THz
         vph=_eps-ktrans2*CSPEED*CSPEED/omega/omega;
         if(vph>=0)
           {
            vph=1./sqrt(vph);
           }
         else
           {
            vph=-1;
           }
         if(vph>0)
           {
            if(parity==0)
              {
               cout<<"f_mn [THz]: "<<_freqsSym[l*ng+found]/1e3<<" type: LSM"<<" m: "<<l<<" kx: "<<_kx<<" n: "<<i<<" fcut: "<<fcut<<" vp: "<<vph<<endl;
              }
            else
              {
               cout<<"f_mn [THz]: "<<_freqsSym[l*ng+found]/1e3<<" type: LSE"<<" m: "<<l<<" kx: "<<_kx<<" n: "<<i<<" fcut: "<<fcut<<" vp: "<<vph<<endl;
              }
// parity=0 => LSM parity=1 => LSE
            fout<<_freqsSym[l*ng+found]/1e3<<" "<<parity<<" "<<l<<" "<<i<<" "<<fcut<<" "<<vph<<endl;
           }
         ++found;
//        }
     } // end "i"
  } // end "l";
 fout.close();
}


void tremaine_roz::calcEigenFreqs3DSymLSM()
{
 double f1,f2,cf,x1,x2,cx,approx,c1,c2,c3,c4,c5,_kx,fac,dx,sf;
 int m,found,parity,check,ng,check1,check2,crit;
 ng=int(_n/_nefx);
 dx=PI/2;
 parity=0;
 for(int l=0; l<_nefx; ++l)
  {
   _kx=_kxs[l];
// solve eqn.: c1*cot(x)=c2*x (LSM)
   c1=1.0/tanh(_kx*_a);
   c2=1.0/_kx/(_b-_a)/_eps; // approx. DM
   found=0;
   for(int i=0; found<ng; ++i)
     {
      x1=i*dx;
      x2=(i+1)*dx;
      x1=x1+EPS;
      x2=x2-EPS;
      parity=int(i/2);
      if(parity*2-i==0)
       {
        parity=0;
       }
      else
       {
        parity=1;
       }
         while(fabs(x2-x1)>SMALL)
           {
            cx=(x2+x1)/2.;
            if(parity==0) // LSM
             {
              f1=calcFunc11(x1,c1,c2);
              cf=calcFunc11(cx,c1,c2);
             }
            if(f1*cf<0)
              {
               x2=cx;
              }
            else
              {
               x1=cx;
              }
            }
         _freqsSym[l*ng+found]=CSPEED/sqrt(_eps-1)*sqrt(cx*cx/(_b-_a)/(_b-_a)+
                          _kx*_kx)/2./PI*CONV2GHZ;
         _lambdasSym[l*ng+found]=CSPEED/_freqsSym[l*ng+found]*CONV2GHZ;
//         cout<<"DM2 (S): "<<_freqsSym[l*ng+found]<<" cx: "<<cx<<" kx: "<<_kx<<" i: "<<i<<endl;
         ++found;
//        }
     } // end "i"
  } // end "l";
}

double tremaine_roz::calcAmplitudes3DSym(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitSym;
 int mi,ng;
 _kx=_kxs[nkx]; 
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=4*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsSym[mi];
    sn=kn*sqrt(_eps-1);
    t1=sinh(2*_kx*_a)/2/_kx*((_kx/kn)*(_kx/kn)+1.0);
    t2=_eps*cosh(_kx*_a)*cosh(_kx*_a)/sin(sn*(_b-_a))/sin(sn*(_b-_a));
    t2=t2+sinh(_kx*_a)*sinh(_kx*_a)/cos(sn*(_b-_a))/cos(sn*(_b-_a));
    t2=t2*(_kx*_kx/sn/sn+2.0)*(_b-_a)/2.0;
    t3=sinh(_kx*_a)*sinh(_kx*_a)/cos(sn*(_b-_a))/cos(sn*(_b-_a));
    t3=t3-_eps*cosh(_kx*_a)*cosh(_kx*_a)/sin(sn*(_b-_a))/sin(sn*(_b-_a));
    t3=t3*_kx*_kx/sn/sn;
    t3=t3+4*_eps*_kx*sn*cosh(_kx*_a)*sinh(_kx*_a)/kn/kn/(_eps-1)/
       sin(sn*(_b-_a))/cos(sn*(_b-_a));
    t3=t3*sin(2*sn*(_b-_a))/4/sn;
    amplitSym=num/(t1+t2+t3)*_lam[nkx];
    _amplitSym[nkx*ng+nm]=amplitSym;
    return amplitSym;
}

double tremaine_roz::calcAmplitudes3DAsym(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitAsym;
 int mi,ng;
 _kx=_kxs[nkx];
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=4*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsAsym[mi];
    sn=kn*sqrt(_eps-1);
    t1=sinh(2*_kx*_a)/2/_kx*((_kx/kn)*(_kx/kn)+1.0);
    t2=_eps*sinh(_kx*_a)*sinh(_kx*_a)/sin(sn*(_b-_a))/sin(sn*(_b-_a));
    t2=t2+cosh(_kx*_a)*cosh(_kx*_a)/cos(sn*(_b-_a))/cos(sn*(_b-_a));
    t2=t2*(_kx*_kx/sn/sn+2.0)*(_b-_a)/2.0;
    t3=cosh(_kx*_a)*cosh(_kx*_a)/cos(sn*(_b-_a))/cos(sn*(_b-_a));
    t3=t3-_eps*sinh(_kx*_a)*sinh(_kx*_a)/sin(sn*(_b-_a))/sin(sn*(_b-_a));
    t3=t3*_kx*_kx/sn/sn;
    t3=t3+4*_eps*_kx*sn*cosh(_kx*_a)*sinh(_kx*_a)/kn/kn/(_eps-1)/sin(sn*(_b-_a))/cos(sn*(_b-_a));
    t3=t3*sin(2*sn*(_b-_a))/4/sn;
    amplitAsym=num/(t1+t2+t3)*_lam[nkx];
    _amplitAsym[nkx*ng+nm]=amplitAsym;
    return amplitAsym;
}

double tremaine_roz::calcAmplitudes3DSymLSM(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitSym,ky;
 int mi,ng;
 _kx=_kxs[nkx]; 
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=2*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsSym[mi];
    sn=kn*sqrt(_eps-1);
    ky=sqrt(kn*kn*(_eps-1)-_kx*_kx);
    t1=sinh(2*_kx*_a)/2/_kx;
    t2=_eps*cosh(_kx*_a)*cosh(_kx*_a)/sin(ky*(_b-_a))/sin(ky*(_b-_a));
    t3=(_b-_a)*(1.0+_eps*_kx*_kx/ky/ky)/2.;
    t3=t3-sin(2*ky*(_b-_a))/4./ky*(1.0-_eps*_kx*_kx/ky/ky);
    t2=t2*t3;
    amplitSym=num/(t1+t2)*_lam[nkx];
    _amplitSym[nkx*ng+nm]=amplitSym;
    return amplitSym;
}

double tremaine_roz::calcAmplitudes3DAsymLSM(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitAsym,ky;
 int mi,ng;
 _kx=_kxs[nkx]; 
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=2*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsSym[mi];
    sn=kn*sqrt(_eps-1);
    ky=sqrt(kn*kn*(_eps-1)-_kx*_kx);
    t1=cosh(2*_kx*_a)/2/_kx;
    t2=_eps*sinh(_kx*_a)*sinh(_kx*_a)/sin(ky*(_b-_a))/sin(ky*(_b-_a));
    t3=(_b-_a)*(1.0+_eps*_kx*_kx/ky/ky)/2.;
    t3=t3-sin(2*ky*(_b-_a))/4./ky*(1.0-_eps*_kx*_kx/ky/ky);
    t2=t2*t3;
    amplitAsym=num/(t1+t2)*_lam[nkx];
    _amplitAsym[nkx*ng+nm]=amplitAsym;
    return amplitAsym;
}

double tremaine_roz::calcAmplitudes3DSymLSE(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitSym,ky;
 int mi,ng;
 _kx=_kxs[nkx]; 
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=2*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsSym[mi];
    sn=kn*sqrt(_eps-1);
    ky=sqrt(kn*kn*(_eps-1)-_kx*_kx);
    t1=sinh(2*_kx*_a)/2/_kx;
    t2=cosh(_kx*_a)*cosh(_kx*_a)/sin(ky*(_b-_a))/sin(ky*(_b-_a));
    t3=(_b-_a)/2.*(_eps+ky*ky/_kx/_kx);
    t3=t3-sin(2*ky*(_b-_a))/4./ky*(_eps-ky*ky/_kx/_kx);
    t2=t2*t3;
    amplitSym=num/(t1+t2)*_lam[nkx];
    _amplitSym[nkx*ng+nm]=amplitSym;
    return amplitSym;
}

double tremaine_roz::calcAmplitudes3DAsymLSE(double x, double y, int nkx, int nm)
{
 double A,num,kn,sn,t1,t2,t3,_kx,amplitAsym,ky;
 int mi,ng;
 _kx=_kxs[nkx]; 
 ng=_n/_nefx;
 mi=ng*nkx+nm;
 num=2*PI;
    kn=2*PI/CONV2GHZ/CSPEED*_freqsSym[mi];
    sn=kn*sqrt(_eps-1);
    ky=sqrt(kn*kn*(_eps-1)-_kx*_kx);
    t1=cosh(2*_kx*_a)/2/_kx;
    t2=sinh(_kx*_a)*sinh(_kx*_a)/sin(ky*(_b-_a))/sin(ky*(_b-_a));
    t3=(_b-_a)/2.*(_eps+ky*ky/_kx/_kx);
    t3=t3-sin(2*ky*(_b-_a))/4./ky*(_eps-ky*ky/_kx/_kx);
    t2=t2*t3;
    amplitAsym=num/(t1+t2)*_lam[nkx];
    _amplitAsym[nkx*ng+nm]=amplitAsym;
    return amplitAsym;
}

void tremaine_roz::outputFreqsSym(char* fname)
{
 fstream fout;
 fout.open(fname,ios::out);
 int nkx,nm,ng;
 ng=int(_n/_nefx);
 for(int i=0; i<_n; ++i)
    {
     nkx=int(i/ng);
     fout<<_freqsSym[i]<<" "<<_amplitSym[i]<<" "<<_lam[nkx]<<" "<<_kxs[nkx]<<endl;
//     fout<<2*PI*_freqsSym[i]<<" "<<_amplitSym[i]<<" "<<_lam[nkx]<<" "<<_kxs[nkx]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputFreqsAsym(char* fname)
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_n; ++i)
    {
     fout<<_freqsAsym[i]<<" "<<_amplitAsym[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::outputFreqsInf(char* fname)
{
 fstream fout;
 fout.open(fname,ios::out);
  for(int i=0; i<_n; ++i)
    {
     fout<<_freqs[i]<<" "<<_amplit[i]<<endl;
    }
 fout.close();
}

void tremaine_roz::calcTrRatio(double z0, double sigz)
{
 double min,max;
 int in,di;
 min=1e20;
 max=-1e20;
 for(int i=0; i<NDIV; ++i)
   {
    if(_z[i]<z0-sigz)
      {
       if(_q3d<0)
        {
         if(_wakeZ[i]<min)
           min=_wakeZ[i];
        }
       else
        {
         if(_wakeZ[i]>max)
           max=_wakeZ[i];
        }
      }
    else
      {
       if(_q3d<0)
        {
         if(_wakeZ[i]>max)
           max=_wakeZ[i];
        }
       else
        {
         if(_wakeZ[i]<min)
           min=_wakeZ[i];
        }
      }
   }

 if(_q3d<0)
  _trRatio=fabs(min/max);
 else
  _trRatio=fabs(max/min);
 cout<<"Transformer Ratio: "<<_trRatio<<endl;

// cout<<"Wmin: "<<endl;
// for(int i=in-di; i<in+di+1; ++i)
//   {
//    cout<<"z: "<<_z[i]<<" Ez: "<<_wakeZ[i]<<endl;
   
//   }
// cout<<"zmin: "<<_z[0]<<" zmax: "<<_z[NDIV]<<endl;
}

double tremaine_roz::getTrRatio()
{
 return _trRatio;
}

void tremaine_roz::calcWakefields3DConv(double zmin, double zmax, 
                   double x, double y)
{
// zmax should not be >Lz
// zmin should not be <0 or >Lz
 if(zmax>_lz || zmin<0 || zmin>_lz)
   {
    cout<<"ERROR: choose 0 < zmin,max < Lz "<<endl;
    exit(0);
   }
 double z,dz,zp,dzp,zpmin,zmaxint,y0,weight_y,ampS,amp1S,ampA,amp1A;
 double wakeScos,wakeSsin,wakeAcos,wakeAsin,knS,knA,kx,ky;
 double wakeScos0,wakeSsin0,wakeAcos0,wakeAsin0,alpha,beta,al1,al4;
 int ni,nkx,nm,nmodes,noi,index,lmin,lmax,parity;

 initWakes();
 dz=(_zmaxCh-_zminCh)/_nzImp;
 zmaxint=zmax;
// if(zmaxint>_lz)
//   zmaxint=_lz;
// calcEigenFreqs3DSym();
 nmodes=_n/_nefx;
 noi=int((_zminCh-zmin)/dz);
 z=_zminCh-(noi+1)*dz;
 lmin=(_ny-_nyImp)/2;
 lmax=lmin+_nyImp;
 _nzOut=0;
 for(int i=0; z<zmax; ++i)
   {
    z=z+dz;
    ni=int((zmaxint-z)/dz)+1;
//    cout<<"current z: "<<z<<"  i: "<<i<<endl;
    for(int k=0; k<_n; ++k)
      {
       nkx=int(k/nmodes);
       nm=k-nkx*nmodes;
       kx=_kxs[nkx];
       knS=2*PI*_freqsSym[k]/CONV2GHZ/CSPEED; 
       knA=2*PI*_freqsAsym[k]/CONV2GHZ/CSPEED;

     if(knS*knS*(_eps-1)-kx*kx>0)
      { 
// determine whether it is a LSM or LSE mode
// k even => LSM
       if(k-int(k/2)*2==0)
         parity=0; // LSM
       else
         parity=1; // LSE
       ky=sqrt(knS*knS*(_eps-1)-kx*kx);
//       cout<<"kx: "<<kx<<" kn: "<<knS<<" ky: "<<ky<<endl;
//       cout<<"kx: "<<kx<<" kn: "<<knS<<" kx/kz: "<<kx/knS<<endl;
       if(parity == 0)
         {
          alpha=-kx*kx/knS/knS;
          beta=1.-alpha;
         }
       else
         {
          alpha=1.0;
          beta=0;
         }

// old version (T-R) paper:
//       alpha=1./2;
//       beta=1./2;

//       al1=-cosh(kx*_a)/sin(ky*(_b-_a));
//       al4=sinh(kx*_a)/cos(ky*(_b-_a));
//       alpha=-kx*(al1*kx+al4*ky)/(kx*kx+ky*ky)*sin(ky*(_b-_a))/cosh(kx*_a);
//       beta=_eps*kx*(-al1*ky+al4*kx)/(kx*kx+ky*ky)*cos(ky*(_b-_a))/sinh(kx*_a);

//       cout<<"alpha: "<<alpha<<" beta: "<<beta<<" alpha+beta: "<<alpha+beta<<" kx: "<<kx<<" ky: "<<ky<<" k: "<<k<<" al1: "<<al1<<" al4: "<<al4<<endl;
//       cout<<" knS: "<<knS<<" freq: "<<_freqsSym[k]<<endl;
//       ampS=calcAmplitudes3DSym(x,y,nkx,nm);

// new version
       if(parity==0)
        {
         ampS=calcAmplitudes3DSymLSM(x,y,nkx,nm);
         ampA=calcAmplitudes3DAsymLSM(x,y,nkx,nm);
        }
       else
        {
         ampS=calcAmplitudes3DSymLSE(x,y,nkx,nm);
         ampA=calcAmplitudes3DAsymLSE(x,y,nkx,nm);
        }

// old version T-R paper
//       ampS=calcAmplitudes3DSym(x,y,nkx,nm);
//       ampA=calcAmplitudes3DAsym(x,y,nkx,nm);

//       cout<<" kx: "<<_kxs[nkx]*1000<<" ky: "<<ky*1000<<" "<<ampS*1000<<" lam: "<<_lam[nkx]/_lam[0]<<endl;
       wakeSsin=0;
       wakeScos=0;
       wakeAsin=0;
       wakeAcos=0;
       for(int l=lmin; l<lmax; ++l)
         {
          y0=_y[l];
          amp1S=ampS*cosh(kx*y0)*_recRhoy[l];
          amp1A=ampA*sinh(kx*y0)*_recRhoy[l];
//          cout<<"amp1S: "<<amp1S<<" amp1A: "<<amp1A<<" sinh: "<<sinh(kx*y0)<<" cosh: "<<cosh(kx*y0)<<" y0: "<<y0<<" l: "<<l<<endl;
          wakeSsin0=0;
          wakeScos0=0;
          wakeAsin0=0;
          wakeAcos0=0;
          for(int j=0; j<ni; ++j)
            {
             zp=z+(j+0.5)*dz;
             if(zp>=_zminCh && zp<=_zmaxCh)
               {
                index=int((zp-_zminCh)/dz);
                wakeScos0=wakeScos0+_rhoz[index]*cos(knS*(zp-z));
                wakeAcos0=wakeAcos0+_rhoz[index]*cos(knA*(zp-z));
                wakeSsin0=wakeSsin0+_rhoz[index]*sin(knS*(zp-z));
                wakeAsin0=wakeAsin0+_rhoz[index]*sin(knA*(zp-z));
               } // end if
            } // end "j" ->z
          wakeScos=wakeScos+wakeScos0*amp1S*dz;
          wakeSsin=wakeSsin+wakeSsin0*amp1S*dz;
          wakeAcos=wakeAcos+wakeAcos0*amp1A*dz;
          wakeAsin=wakeAsin+wakeAsin0*amp1A*dz;
         } // end "l" ->y0
//        cout<<"S: "<<wakeSsin<<" A: "<<wakeAsin<<endl;
//        cout<<"S: "<<wakeScos<<" A: "<<wakeAcos<<endl;
       if(y<_a && y>-_a)
        {
       _wakeZ[i]=_wakeZ[i]+wakeScos*cosh(kx*y)*cos(kx*x);
       _wakeZ[i]=_wakeZ[i]+wakeAcos*sinh(kx*y)*cos(kx*x);
       _wakeBZ[i]=_wakeBZ[i]+wakeScos*sinh(kx*y)*sin(kx*x);
       _wakeBZ[i]=_wakeBZ[i]+wakeAcos*cosh(kx*y)*sin(kx*x);
       _wakeEX[i]=_wakeEX[i]+knS/kx*alpha*wakeSsin*cosh(kx*y)*sin(kx*x);
       _wakeEX[i]=_wakeEX[i]+knA/kx*alpha*wakeAsin*sinh(kx*y)*sin(kx*x);
       _wakeEY[i]=_wakeEY[i]+knS/kx*beta*wakeSsin*sinh(kx*y)*cos(kx*x);
       _wakeEY[i]=_wakeEY[i]+knA/kx*beta*wakeAsin*cosh(kx*y)*cos(kx*x);
       _wakeBX[i]=_wakeBX[i]-(knS/kx*beta-kx/knS)*wakeSsin*sinh(kx*y)*cos(kx*x);
       _wakeBX[i]=_wakeBX[i]-(knA/kx*beta-kx/knA)*wakeAsin*cosh(kx*y)*cos(kx*x);
       _wakeBY[i]=_wakeBY[i]+(knS/kx*alpha+kx/knS)*wakeSsin*cosh(kx*y)*sin(kx*x);
       _wakeBY[i]=_wakeBY[i]+(knA/kx*alpha+kx/knA)*wakeAsin*sinh(kx*y)*sin(kx*x);
       _FX[i]=_FX[i]-kx/knS*wakeSsin*cosh(kx*y)*sin(kx*x);
       _FX[i]=_FX[i]-kx/knA*wakeAsin*sinh(kx*y)*sin(kx*x);
       _FY[i]=_FY[i]+kx/knS*wakeSsin*sinh(kx*y)*cos(kx*x);
       _FY[i]=_FY[i]+kx/knA*wakeAsin*cosh(kx*y)*cos(kx*x);
        }
       else
        {
         if(parity==0)
          {
          _wakeZ[i]=_wakeZ[i]+wakeScos*sin(ky*(_b-fabs(y)))*cos(kx*x)*cosh(kx*_a)
                     /sin(ky*(_b-_a));
           _wakeEX[i]=_wakeEX[i]-kx/knS*wakeSsin*sin(ky*(_b-fabs(y)))*sin(kx*x)*
                      cosh(kx*_a)/sin(ky*(_b-_a));
           _wakeEY[i]=_wakeEY[i]-1./knS/ky*(knS*knS+kx*kx)*wakeSsin*
                      cos(kx*x)*cos(ky*(_b-fabs(y)))*cosh(kx*_a)/sin(ky*(_b-_a))
                      *fabs(y)/(-y);
                      
          }
         else
          {
          _wakeZ[i]=_wakeZ[i]+wakeScos*sin(ky*(_b-fabs(y)))*cos(kx*x)*cosh(kx*_a)
                     /sin(ky*(_b-_a));
           _wakeEX[i]=_wakeEX[i]+knS/kx*wakeSsin*sin(ky*(_b-fabs(y)))*sin(kx*x)*
                      cosh(kx*_a)/sin(ky*(_b-_a));
          }
        }
      } // end if for mode consistency 
      } // end "k" ->modes
    // "-" sign shows that wakefields brake the driving charge
    _wakeZ[i]=-_wakeZ[i]*CSI/_dz*_q3d; //convert to MV/m; convert q to rho; 
    _wakeEX[i]=_wakeEX[i]*CSI/_dz*_q3d; 
    _wakeEY[i]=_wakeEY[i]*CSI/_dz*_q3d; 
    _wakeBZ[i]=-_wakeBZ[i]*CONVB/_dz*_q3d; // T
    _wakeBX[i]=_wakeBX[i]*CONVB/_dz*_q3d; 
    _wakeBY[i]=_wakeBY[i]*CONVB/_dz*_q3d; 
    _FX[i]=_FX[i]/_dz*_q3d; // MV/m
    _FY[i]=_FY[i]/_dz*_q3d; // MV/m
//    cout<<_wakeZ[i]<<endl;
    _z[i]=z;
    ++_nzOut;
   }
}

void tremaine_roz::initWakes() 
{
 for(int i=0; i<NDIVZ; ++i)
   {
    _wakeZ[i]=0;
    _wakeEX[i]=0;
    _wakeEY[i]=0;
    _wakeBX[i]=0;
    _wakeBY[i]=0;
    _wakeBZ[i]=0;
    _FX[i]=0;
    _FY[i]=0;
   }
}

int tremaine_roz::getNzout() 
{
 return _nzOut;
}

double tremaine_roz::getEx(int in) 
{
 return _wakeEX[in];
}

double tremaine_roz::getEy(int in) 
{
 return _wakeEY[in];
}

double tremaine_roz::getEz(int in) 
{
 return _wakeZ[in];
}


double tremaine_roz::getBx(int in) 
{
 return _wakeBX[in];
}

double tremaine_roz::getBy(int in) 
{
 return _wakeBY[in];
}

double tremaine_roz::getBz(int in) 
{
 return _wakeBZ[in];
}

/*
void tremaine_roz::overlapWakefields3DConvKx(double zmin, double zmax, double sigz, double z0,
                                             int nkx)
{
// zmax should not be >Lz
// zmin should not be <0 or >Lz
 if(zmax>_lz || zmin<0 || zmin>_lz)
   {
    cout<<"ERROR: choose 0 < zmin,max < Lz "<<endl;
    exit(0);
   }
 double z,dz,wake,zp,dzp,zpmin,zmaxint,dxp,x,cx,_kx;
 dz=(zmax-zmin)/NDIV;
// zmaxint=z0+NSIGZ*sigz;
// if(zmaxint<zmax)
//   zmaxint=zmax;
 for(int i=0; i<NDIV; ++i)
   {
    z=zmin+i*dz;
    wake=0;
    cout<<"current z: "<<z<<"  i: "<<i<<endl;
    for(int p=0; p<nkx; ++p)
    {
     _kx=_kxs[nkx];
    for(int k=0; k<_n; ++k)
      {
       dzp=(_lz-z)/NDIV;
       for(int j=0; j<NDIV; ++j)
         {
          zp=z+j*dzp;
          dxp=_lx/NDIV;
//          calcAmplitudes3DSym(0,0.001,0.001);
//wake=wake+calcFunc2(zp,sigz,z0)*cos(2*PI*_freqsSym[k]/CONV2GHZ/CSPEED*(zp-z))*
//               _amplitSym[k]*dzp;


          for(int m=0; m<NDIV; ++m)
            {
             x=-_lx/2+m*dxp;
             calcAmplitudes3DSym(x,0.5,0.0);
             wake=wake+calcFunc2(zp,sigz,z0)*cos(2*PI*_freqsSym[k]/CONV2GHZ/CSPEED*(zp-z))*
                  cos(_kx*x)*_amplitSym[k]*dzp*dxp*_kx/2;
             calcAmplitudes3DAsym(x,0.5,0.0);
             wake=wake+calcFunc2(zp,sigz,z0)*cos(2*PI*_freqsAsym[k]/CONV2GHZ/CSPEED*(zp-z))*
                  cos(_kx*x)*_amplitAsym[k]*dzp*dxp*_kx/2;
            }

         }
     }
     for(int i1=0; i1<NDIV; ++i1)
      {
       cx=-_lx/2+i1*_lx/NDIV;
       _rhox[i1]=_rhox[i1]+cos(_kx*cx)*_kx/2;
      }
    } // end "p" loop
    _wakeZ[i]=wake/nkx;
    _z[i]=z;
   }
}
*/

void tremaine_roz::generateRhox(double xcutoff, int nxImp, double sigx)
{
 double dx,x,lnorm,parity,alx,lnormt,ll,ul;
 int nx;
 lnorm=0;
 dx=2*xcutoff/(nxImp-1);
 _dx=dx;
 _nxImp=nxImp;
 ll=dx/2;
 ul=(nxImp/2-1)*dx;
 for(int j=0; j<_nefx; ++j)
   {
//    _kxs[j]=(2*j+1)*PI/2/xcutoff;
    _kxs[j]=(2*j+1)*PI/_lx;
    parity=-2*(j-2*int(j/2))+1;
    _lam[j]=rectRhox(ll,ul,nxImp/2,sigx,xcutoff,j);
    lnorm=lnorm+_lam[j]/(2*j+1)*parity;
   }
////???
 lnorm=lnorm*2*_lx/PI;  
///????
// cout<<"lnorm (lam): "<<lnorm<<endl;
 for(int j=0; j<_nefx; ++j)
   {
    _lam[j]=_lam[j]/lnorm; // must have units of nC/mm
    //cout<<"Ratios: "<<_lam[j]/_lam[0]<<endl;
   }
 _nx=2+2*int((_lx-dx)/2/dx);
 alx=(_nx-1)*dx;
// cout<<"nx: "<<_nx<<" dx: "<<dx<<" alx: "<<alx<<endl;
// lnorm=0;
 lnormt=0;
 for(int i=0; i<_nx; ++i)
   {
    _x[i]=-alx/2+i*dx;
    _rhox[i]=calcFuncRhox(_x[i],sigx,xcutoff);
    lnormt=lnormt+_rhox[i];
//    if(_x[i]>xcutoff || _x[i]<-xcutoff)
//      {
//       _recRhox[i]=0;
//      }
//    else 
//      {
       _recRhox[i]=0;
       for(int j=0; j<_nefx; ++j)
         {
          _recRhox[i]=_recRhox[i]+_lam[j]*cos((2*j+1)*PI/_lx*_x[i])*dx;
         }
//      }
////????
//    lnorm=lnorm+_recRhox[i];
///???
   }
 for(int i=0; i<_nx; ++i)
   {
    _rhox[i]=_rhox[i]/lnormt;
   }
// cout<<"normx: "<<lnorm<<" norm (true): "<<lnormt<<endl;
// for(int j=0; j<_nefx; ++j)
//   {
//    _lam[j]=_lam[j]*lnorm;
//   }
}

void tremaine_roz::generateRhoy(double ycutoff, int nyImp, double sigy)
{
 double dy,aly,lnorm;
 dy=2*ycutoff/(nyImp-1);
 _dy=dy;
 _nyImp=nyImp;
 _ny=2+2*int((2*_a-dy)/2/dy);
 aly=(_ny-1)*dy;
// cout<<"ny: "<<_ny<<" dy: "<<dy<<" aly: "<<aly<<endl;
 lnorm=0;
 for(int i=0; i<_ny; ++i)
   {
    _y[i]=-aly/2+i*dy;
    _rhoy[i]=calcFuncRhoy(_y[i],sigy);
    if(_y[i]>ycutoff+EPS || _y[i]<-ycutoff-EPS)
      {
       _recRhoy[i]=0;
      }
    else 
      {
       _recRhoy[i]=calcFuncRhoy(_y[i],sigy);
      }
    lnorm=lnorm+_recRhoy[i];
   }
// lnorm=lnorm*dy;
// cout<<"normy: "<<lnorm<<endl;
 for(int i=0; i<_ny; ++i)
   {
    _recRhoy[i]=_recRhoy[i]/lnorm;
   }
}

void tremaine_roz::generateRhoz(double zmin, double zmax, int nzImp, double a,
                                double sigz1, double z01, double sigz2, double z02)
{
 double dz,lnorm;
 _zminCh=zmin;
 _zmaxCh=zmax;
 _nzImp=nzImp;
 dz=(zmax-zmin)/(nzImp-1);
 _nz=nzImp;
 _dz=dz;
 lnorm=0;
 for(int i=0; i<_nz; ++i)
   {
    _zb[i]=zmin+i*dz;
    _rhoz[i]=calcFuncRhoz(_zb[i],a,sigz1,z01,sigz2,z02);
//    _rhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    _recRhoz[i]=calcFuncRhoz(_zb[i],a,sigz1,z01,sigz2,z02);
//    _recRhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    lnorm=lnorm+_recRhoz[i];
   }
// lnorm=lnorm*dz;
// cout<<"normz: "<<lnorm<<endl;
 for(int i=0; i<_nz; ++i)
   {
    _recRhoz[i]=_recRhoz[i]/lnorm;
    _rhoz[i]=_rhoz[i]/lnorm;
   }
}

void tremaine_roz::generateRhozRamp(double zmin, double zmax, int nzImp, double a,
                                double sigz1, double z01, double sigz2, double z02)
{
 double dz,lnorm;
 _zminCh=zmin;
 _zmaxCh=zmax;
 _nzImp=nzImp;
 dz=(zmax-zmin)/(nzImp-1);
 _nz=nzImp;
 _dz=dz;
 lnorm=0;
 for(int i=0; i<_nz; ++i)
   {
    _zb[i]=zmin+i*dz;
    _rhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
//    _rhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    _recRhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
//    _recRhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    lnorm=lnorm+_recRhoz[i];
   }
// lnorm=lnorm*dz;
 cout<<"normz: "<<lnorm<<" _zb[0]: "<<_zb[0]<<endl;
 for(int i=0; i<_nz; ++i)
   {
    _recRhoz[i]=_recRhoz[i]/lnorm;
    _rhoz[i]=_rhoz[i]/lnorm;
   }
}

void tremaine_roz::generateRhozRampG(double zmin, double zmax, int nzImp, 
                   double z01, double d, double sigz, int ng)
{
 double dz,lnorm;
 _zminCh=zmin;
 _zmaxCh=zmax;
 _nzImp=nzImp;
 dz=(zmax-zmin)/(nzImp-1);
 _nz=nzImp;
 _dz=dz;
 lnorm=0;
 for(int i=0; i<_nz; ++i)
   {
    _zb[i]=zmin+i*dz;
    _rhoz[i]=calcFuncRhozRampG(_zb[i],z01,d,sigz,ng);
//    _rhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    _recRhoz[i]=calcFuncRhozRampG(_zb[i],z01,d,sigz,ng);
//    _recRhoz[i]=calcFuncRhozRamp(_zb[i],a,sigz1,z01,sigz2,z02);
    lnorm=lnorm+_recRhoz[i];
   }
// lnorm=lnorm*dz;
// cout<<"normz: "<<lnorm<<endl;
 for(int i=0; i<_nz; ++i)
   {
    _recRhoz[i]=_recRhoz[i]/lnorm;
    _rhoz[i]=_rhoz[i]/lnorm;
   }
}

void tremaine_roz::generateRhoXY()
{
 for(int i=0; i<_nx; ++i)
   {
    for(int j=0; j<_ny; ++j)
      {
       _rhoxy[i][j]=_rhox[i]*_rhoy[j];
       _recRhoXY[i][j]=_recRhox[i]*_recRhoy[j];
      }
   }
}

void tremaine_roz::checkCharge()
{
 double sum;
 sum=0;
 for(int i=0; i<_nx; ++i)
   {
    for(int j=0; j<_ny; ++j)
      {
       for(int k=0; k<_nz; ++k)
         {
          sum=sum+_recRhoz[k]*_recRhoy[j]*_recRhox[i];
         }
      }
   }
 cout<<"nx: "<<_nx<<" ny: "<<_ny<<" nz: "<<_nz<<endl;
 cout<<"Total Charge (nC): "<<sum*_q3d*1.0e9<<endl;
}

double tremaine_roz::simpRhox(double ll, double ul, int n, double sigx,
                              double xcutoff, int j)
{
 double x,dx,f0,fn,os,es;
 dx=(ul-ll)/(n-1);
 f0=calcFuncRhox(ll,sigx,xcutoff)*cos((2*j+1)*PI/_lx*ll);
 fn=calcFuncRhox(ul,sigx,xcutoff)*cos((2*j+1)*PI/_lx*ul);
 os=0;
 for(int i=1; i<n; i=i+2)
   {
    x=ll+i*dx;
    os=os+calcFuncRhox(x,sigx,xcutoff)*cos((2*j+1)*PI/_lx*x);
   }
 os=4*os;
 es=0;
 for(int i=2; i<n; i=i+2)
   {
    x=ll+i*dx;
    es=es+calcFuncRhox(x,sigx,xcutoff)*cos((2*j+1)*PI/_lx*x);
   }
 es=2*es;
 return 1.0/3*(f0+es+os+fn);
}

double tremaine_roz::trapezRhox(double ll, double ul, int n, double sigx,
                              double xcutoff, int j)
{
 double x,dx,f0,fn,s,f1;
 dx=(ul-ll)/(n-1);
 f0=calcFuncRhox(ll,sigx,xcutoff)*cos((2*j+1)*PI/_lx*ll);
 fn=calcFuncRhox(ul,sigx,xcutoff)*cos((2*j+1)*PI/_lx*ul);
 f1=calcFuncRhox(ul,sigx,xcutoff)*cos((2*j+1)*PI/_lx*(ll+dx));
 s=0;
 for(int i=1; i<n-1; ++i)
   {
    x=ll+i*dx;
    s=s+calcFuncRhox(x,sigx,xcutoff)*cos((2*j+1)*PI/_lx*x);
   }
 s=2*s;
 s=dx/2*(f0+s+fn);
 s=s+dx*(f0+f1/4.0);
 return s;
}

double tremaine_roz::rectRhox(double ll, double ul, int n, double sigx,
                              double xcutoff, int j)
{
 double x,dx,f0,fn,s;
 dx=(ul-ll)/(n-1);
 s=0;
 for(int i=0; i<n; ++i)
   {
    x=ll+i*dx;
    s=s+calcFuncRhox(x,sigx,xcutoff)*cos((2*j+1)*PI/_lx*x);
   }
 return s*dx;
}

double tremaine_roz::calcFuncRhox(double x, double sigx, double xcutoff)
{
 double f;
 if(x<-xcutoff || x>xcutoff)
   {
    f=0;
   }
 else
   {
    f=CGAUSS*exp(-x*x/2/sigx/sigx)/sigx;
   }
 return f;
}

double tremaine_roz::calcFuncRhoy(double y, double sigy)
{
 return CGAUSS*exp(-y*y/2/sigy/sigy)/sigy;
}

double tremaine_roz::calcFuncRhoz(double z, double a, double sigz1, double z01,
                                  double sigz2, double z02)
{
 double s;
 s=CGAUSS*a*exp(-(z-z01)*(z-z01)/2/sigz1/sigz1)/sigz1;
 s=s+CGAUSS*(1.0-a)*exp(-(z-z02)*(z-z02)/2/sigz2/sigz2)/sigz2;
 return s;
}

double tremaine_roz::calcFuncRhozRamp(double z, double a, double sigz1, double z01,
                                  double sigz2, double z02)
{
 double s,m,n;
 s=0;
 if(z>z01-sigz1/3. && z<z01+sigz1*2./3.)
  s=-2.0/sigz1/sigz1*a*(z-z01-2.0*sigz1/3.0);
 
 if(z>z02-sigz2/3. && z<z02+sigz2*2./3.)
  s=s-2.0/sigz2/sigz2*(1.0-a)*(z-z02-2.0*sigz2/3.0);
 return s;
}

double tremaine_roz::calcFuncRhozRampG(double z, double z01, double d,
                                  double sigz, int ng)
{
 double s,q[ng],sq;
 sq=0;
 for(int i=0; i<ng; ++i)
   {
    double sl;
    sl=0;
    q[i]=sl*i+1;
    sq=sq+q[i];
   }
 s=0;
 for(int i=0; i<ng; ++i)
   {
    s=s+q[i]*exp(-(z-z01+i*d)*(z-z01+i*d)/2/sigz/sigz);
   }
 s=s*CGAUSS/sigz/sq;
 return s;
}

double tremaine_roz::calcFunc1(double x)
{
 return x*tan(x)-(_b-_a)*_eps/_a;
}

double tremaine_roz::calcFunc2(double z, double sigz, double z0)
{
 return CGAUSS/sigz*exp(-(z-z0)*(z-z0)/2/sigz/sigz);
}

double tremaine_roz::calcFunc3(double z, double sigz, double z0, double delz, 
                               int nb)
{
 double val;
 val=0;
 for(int i=0; i<nb; ++i)
   {
    val=val+1./sqrt(2*PI)/sigz*exp(-(z-z0+i*delz)*(z-z0+i*delz)/2/sigz/sigz);
   }
 return val;
}

double tremaine_roz::calcFunc4(double x, double c1, double c2, double c3, double c4)
{
 return tan(c1*x)-c2*x/(c3+c4*x*x);
}

double tremaine_roz::calcFunc5(double x, double c1, double c2, double c3)
{
 return tan(x)-c1*x/(c3*x*x-c2);
}

double tremaine_roz::calcFunc10(double x, double c1, double c2, double c3, double c4,
                                double c5)
{
 return c1/tan(x)-c2*tan(x)-(c3*x*x+c4)/c5/x;
}

double tremaine_roz::calcFunc11(double x, double c1, double c2)
{ 
 return c1/tan(x)-c2*x;
}

double tremaine_roz::calcFunc12(double x, double c1, double c3)
{ 
 return c1/tan(x)+c3/x;
}

