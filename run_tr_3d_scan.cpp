
 #include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "tremaine_roz.h"
#include "tremaine_roz.cpp"

using namespace std;
int main()
{
 int nxImp,nefx,nyImp,nzImp,nmodes,ng;
 double x, a,b,xcutoff,ycutoff,sigx,sigy,zmin,zmax,zminCh,zmaxCh;
 double z01,z02,sigz1,sigz2,ratio,eps,q,lz,y,lx;
 double dx, dz, dy, zreal, ztarget;
 int nxsample,nysample, in, ix, iy;
 char fname[80];
// x info: 
 lx=10.0; // mm
 xcutoff=3.0; // charge extent in x
 sigx=1e5;  
 nxImp=16; // grid points in x
 nefx=10; // eigenfuntions in x
 x=0; // x-coordinate where wakefields are evaluated

// y info: 
 a=2.5/2.; // half-gap (mm)
 b=0.150; // dielectric thickness
 ycutoff=0.3;
 sigy=1e5;
 nyImp=8;
 y=0.0; // y-coordinate where wakefields are evaluated

// z info:  
 lz=100.0; // mm
 zmin=0.0; 
 zmax=100.0; // zmax-zmin=lz
 zminCh=95.0; // zmin-cutoff for charge
 zmaxCh=100.0; // zmax-cutoff for charge
 z01=97.5; // z-position of charge centroid
 sigz1=1.0; // bunch length
 nzImp=32; // grid points in z-direction

// dielectric constant
 eps=3.75;
// beam charge (drive+test)
 q=-0.5; // nC
// total # of modes: #modes/eigenfunction * nefx
 nmodes=100; // (nefx*ny)

// The next 3 lines reffer to witness charge
// If ratio=1 there is no witness charge 
 ratio=1.0;
 z02=1.0;
 sigz2=1.0;

 tremaine_roz st;
 st.setStructure(a,b,lx,lz,eps);// tremaine paper (a,b,lx,lz,eps)
 st.setNModes(nmodes);
 st.setNEigenX(nefx);
 st.setQ3D(q); // nC

 st.generateRhox(xcutoff,nxImp,sigx); // xcutoff(mm), nxImp, sigx(mm)
 st.generateRhoy(ycutoff,nyImp,sigy); // ycutoff(mm), nyImp, sigy(mm)
 st.generateRhoz(zminCh,zmaxCh,nzImp,ratio,sigz1,z01,sigz2,z02); 
 st.calcEigenFreqs3DSymDM2();
 st.calcEigenFreqs3DAsymDM2();

/// dump 1-D longitudinal wakefield on axis
 st.outputWakefields1D("wake_z_wz.dat");

// dump fields on a given surface  
 fstream fout;
 fout.open("wake_2Dsurf.dat",ios::out);
 
 // here consider the (x,y) plane at a given z location "ztarget" (for now no interpolation)
 ztarget=4; //
 dz=(zmaxCh-zminCh)/nzImp; // mm
 in=int(ztarget/dz+0.5)-1;
 zreal=(in+1)*dz;
 
// cout<<"ztarget: "<<ztarget<<" zreal: "<<zreal<<" in: "<<in<<endl;
 nxsample=21; // # of x bins
 nysample=21; // # of y bins 
 dx=lx/(nxsample-1);
 dy=20*b/(nysample-1);
// for loop for an integer "iy", where "iy" must be less than nysample (see 3 lines up), and "iy" is increased by 1
 for(int iy=0; iy<nysample; ++iy)
    for(int ix=0; ix<nxsample; ++ix)
      {
       x=-lx/2.+ix*dx;
       y=-10*b+iy*dy; //Calculating y, so "b" is in the 'yinfo' section and dy is right above
       st.calcWakefields3DConv(zmin,zmax,x,y); // zmin,zmax,sigz,z0,x,y
       fout<<x<<" "<<y<<" "<<zreal<<" "<<st.getEx(in)<<" "<<st.getEy(in)<<" "<<st.getBx(in)<<" "<<st.getBy(in)<<endl; 
       cout<<ix <<" "<<iy<<" x: "<<x<<" y: "<<y<<" z: "<<zreal<<" Ex: "<<st.getEx(in)<<" Ey: "<<st.getEy(in)<<" Bx: "<<st.getBx(in)<<" By: "<<st.getBy(in)<<endl; 
      }
 fout.close();
 st.checkCharge();
// st.calcTrRatio(27.0);
 return 0;
}
