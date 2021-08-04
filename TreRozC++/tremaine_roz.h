#ifndef tremaine_roz_H
#define tremaine_roz_H

const int NMODEMAX=2000; // max # of modes
const int NDIV=1000;
const int NDIVX=500;
const int NDIVY=500;
const int NDIVZ=1000;
const int NSIGZ=3;
const int NEIGENX=100;
const int NEIGENY=100;
const double PI=3.141592654;
const double CSPEED=3.0e11; // mm/s
const double CSPEEDSI=3.0e8; // m/s
const double CSI=9.0e9; // 1/(4*PI*EPS0) = 9*10^9
const double CSIM=1.0e-7; // mu0/(4*pi)
const double CONVB=3.0e7; // convert B to SI (9*10^9*10^6(because of mm)/c)
const double CONV2GHZ=1.0e-9;
const double SMALL=1e-10;
const double EPS=1.0e-10;
const double CGAUSS=0.39894228; // 1.0/sqrt(2*pi)

class tremaine_roz{
  public:
    tremaine_roz();
    ~tremaine_roz();
    void setStructure(double a, double b, double lz, double eps); // a,b mm
    void setStructure(double a, double b, double lx, double lz, double eps); // a,b mm
    void setNModes(int n);
    void setNEigenX(int nefx);
    void calcEigenFreqsInf();
    void calcAmplitudesInf();
    void calcWakefieldsDelta(double zmin, double zmax, double z0);
    void calcWakefieldsConv(double zmin, double zmax, double sigz, double z0);
    void overlapWakefields(double zmin, double zmax, double sigz, double z0,
                          double delz, int nb);
    void normWakefieldsToGeV();
    void convZtoMicrons();
    void outputWakefields1D(char* fname);
    void outputWakefields1Dall(char* fname);
    void outputRhox(char* fname);
    void outputRhoy(char* fname);
    void outputRhoz(char* fname);
    void outputRhoxy(char* fname);
    void outputFreqsSym(char* fname);
    void outputFreqsAsym(char* fname);
    void outputFreqsInf(char* fname);
    void setQInf(double q); // q is in nC/mm
    void setQ3D(double q); // q is in nC
    void calcEigenFreqs3DSym();
    void calcEigenFreqs3DSymDM();
    void calcEigenFreqs3DSymDM2();
    void calcTrRatio(double z0, double z01);
    void calcEigenFreqs3DAsym();
    void calcEigenFreqs3DAsymDM();
    void calcEigenFreqs3DAsymDM2();
    void calcEigenFreqs3DSymDM2VpOut(double flaser, char* fname);
    void calcEigenFreqs3DSymLSM();
    double calcAmplitudes3DSym(double x, double y, int nkx, int nm);
    double calcAmplitudes3DAsym(double x, double y, int nkx, int nm);
    double calcAmplitudes3DSymLSM(double x, double y, int nkx, int nm);
    double calcAmplitudes3DAsymLSM(double x, double y, int nkx, int nm);
    double calcAmplitudes3DSymLSE(double x, double y, int nkx, int nm);
    double calcAmplitudes3DAsymLSE(double x, double y, int nkx, int nm);
    void calcWakefields3DConv(double zmin, double zmax, double x, double y);
//    void overlapWakefields3DConvKx(double zmin, double zmax, double sigz, double z0, int nkx);
    void initWakes();
    int getNzout();
    double getEx(int in);
    double getEy(int in);
    double getEz(int in);
    double getBx(int in);
    double getBy(int in);
    double getBz(int in);
    double getKx(int nkx);
    void generateRhox(double xcutoff, int nxImp, double sigx);
    void generateRhoy(double ycutoff, int nyImp, double sigy);
    void generateRhoz(double zmin, double zmax, int nyImp, double a,double sigz1,
                      double z01, double sigz2, double z02);
    void generateRhozRamp(double zmin, double zmax, int nyImp, double a,double sigz1,
                      double z01, double sigz2, double z02);
    void generateRhozRampG(double zmin, double zmax, int nyImp, 
                      double z01, double d, double sigz, int ng);
    void generateRhoXY();
    void checkCharge();
    double getTrRatio();
    double simpRhox(double ll, double ul, int n, double sigz, double xcutoff,
                    int j);
    double trapezRhox(double ll, double ul, int n, double sigz, double xcutoff,
                    int j);
    double rectRhox(double ll, double ul, int n, double sigz, double xcutoff,
                    int j);
    double calcFuncRhox(double x, double sigx, double xcutoff);
    double calcFuncRhoy(double y, double sigy);
    double calcFuncRhoz(double z, double a, double sigz1, double z01,
                        double sigz2, double z02);
    double calcFuncRhozRamp(double z, double a, double sigz1, double z01,
                        double sigz2, double z02);
    double calcFuncRhozRampG(double z, double z01, double d, 
                        double sigz, int ng);
    double calcFunc1(double x);
    double calcFunc2(double z, double sigz, double z0);
    double calcFunc3(double z, double sigz, double z0, double delz, int nb);
    double calcFunc4(double x, double c1, double c2, double c3, double c4);
    double calcFunc5(double x, double c1, double c2, double c3);
    double calcFunc10(double x, double c1, double c2, double c3, double c4, double c5);
    double calcFunc11(double x, double c1, double c2);
    double calcFunc12(double x, double c1, double c3);
  private:
    int _n,_nefx,_nx,_ny,_nz,_nzOut;
    int _nxImp,_nyImp,_nzImp;
    double _a,_b,_l,_eps;
    double _q3d,_dx,_dy,_dz;
    double _freqs[NMODEMAX],_lambdas[NMODEMAX],_amplit[NMODEMAX];
    double _z[NDIV],_wakeZ[NDIVZ],_zb[NDIVZ];
    double _wakeBZ[NDIVZ],_wakeEX[NDIVZ],_wakeBX[NDIVZ],_wakeEY[NDIVZ];
    double _wakeBY[NDIVZ],_FX[NDIVZ],_FY[NDIVZ];
    double _rhox[NDIVX],_kxs[NEIGENX],_x[NDIVX],_lam[NEIGENX],_recRhox[NDIVX];
    double _rhoy[NDIVY],_y[NDIVY],_recRhoy[NDIVY];
    double _rhoxy[NDIVX][NDIVY],_recRhoXY[NDIVX][NDIVY];
    double _rhoz[NDIVZ],_recRhoz[NDIVZ];
    double _freqsSym[NMODEMAX],_lambdasSym[NMODEMAX],_freqsAsym[NMODEMAX];
    double _lambdasAsym[NMODEMAX];
    double _amplitSym[NMODEMAX],_amplitAsym[NMODEMAX];
    double _lx,_lz,_trRatio;
    double _zminCh,_zmaxCh;
};

#endif
