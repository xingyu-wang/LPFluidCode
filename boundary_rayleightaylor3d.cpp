#include "boundary_rayleightaylor3d.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;
static const double onesq2 = 1.0/std::sqrt(2.0);

////////////////////////////////////////////////////////////////////////////////////////
// Start of Shocktube3DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////

RayleighTaylor3DBoundary::RayleighTaylor3DBoundary():lengthX(1), lengthY(1), lengthZ(1), thickness(0.1) {
        rb = lengthX/2.-thickness;
        lb = -lengthX/2.+thickness;
        nb = lengthY/2.-thickness;
        sb = -lengthY/2.+thickness;
        tb = lengthZ/2.-thickness;
        bb = -lengthZ/2.+thickness;
        rbo = lengthX/2.;
        lbo = -lengthX/2.;
        nbo = lengthY/2.;
        sbo = -lengthY/2.;
        tbo = lengthZ/2;
        bbo = -lengthZ/2;

        epsilon=1e-3;
}

int RayleighTaylor3DBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

        if(x<rb && x>lb && y<nb && y>sb && z<tb && z>bb) {
//              cout<<"inside"<<endl;
                //cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
                //cout<<"rb="<<rb<<" lb="<<lb<<" nb="<<nb<<" sb="<<sb<<endl;
                return 0; // inside
        }
        if(x>rbo || x<lbo || y>nbo || y<sbo || z>tbo || z<bbo) {
//      cout<<"outside"<<endl;
        return 0;
        } // outside       

        int bx, by, bz;
        int count = 0;
        bx=0;
        by=0;
        bz=0;

        if (x>=rb) bx = 1;
        if (x<=lb) bx = -1;
        if (y>=nb) by = 1;
        if (y<=sb) by = -1;
        if (z>=tb) bz = 1;
        if (z<=bb) bz = -1;
//      cout<<bx<<" "<<by<<" "<<bz<<" "<<endl;
        for (int ix = 0; ix <= abs(bx); ix++)
                for (int iy = 0; iy <= abs(by); iy++)
                        for (int iz = 0; iz <= abs(bz); iz++)
                                if(ix+iy+iz)
                                {
                                        reflect(x, y, z, pressure, vx, vy, vz, ix*bx, iy*by, iz*bz, xb, yb, zb, pressureb, vxb, vyb, vzb);
                                        count=count+1;
                                }
//      cout<<count<<endl;      
        return count;
}

void RayleighTaylor3DBoundary::reflect(double x, double y, double z, double pressure, double vx, double vy, double vz, int bx, int by, int bz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb, vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

        if(bx==1)
                xb.push_back(2*rbo-x);
        else if(bx==-1)
                xb.push_back(2*lbo-x);
        else
                xb.push_back(x);

        if(by==1)
                yb.push_back(2*nbo-y);
        else if(by==-1)
                yb.push_back(2*sbo-y);
        else
                yb.push_back(y);

        if(bz==1)
                zb.push_back(2*tbo-z);
        else if(bz==-1)
                zb.push_back(2*bbo-z);
        else
                zb.push_back(z);

        pressureb.push_back(pressure);

        double mag=sqrt(bx*bx+by*by+bz*bz);
        double ip=(vx*bx+vy*by+vz*bz)/mag;
        if((ip+epsilon)<=0) { //approaching solid boundary
                vxb.push_back(0);
                vyb.push_back(0);
                vzb.push_back(0);
        }
        else {
                vxb.push_back(vx-2.*ip*bx);
                vyb.push_back(vy-2.*ip*by);
                vzb.push_back(vz-2.*ip*bz);
        }
}
