#include "boundary_kelvinhelmholtz.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;
static const double onesq2 = 1.0/std::sqrt(2.0);

////////////////////////////////////////////////////////////////////////////////////////
// Start of RayleighTaylor2DBoundary
////////////////////////////////////////////////////////////////////////////////////////

KelvinHelmholtz2DBoundary::KelvinHelmholtz2DBoundary():lengthX(1), lengthY(1), thickness(0.2) {
	rb = lengthX-thickness;
	lb = thickness;
	nb = lengthY-thickness;
	sb = thickness;
	
	rbo = lengthX;
	lbo = 0;
	nbo = lengthY;
	sbo = 0;

        rbl = lengthX+thickness;
        lbl = -thickness;
        nbl = lengthY+thickness;
        sbl = -thickness;
}

int KelvinHelmholtz2DBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	
	
	if(x<rb && x>lb && y<nb && y>sb) {
		return 0; // inside
	}
	if(x>rbo || x<lbo || y>nbo || y<sbo) {
		if(x>rbo)
			xb.push_back(x-ceil((x-rbo)/lengthX)*lengthX);
		else if(x<lbo)
                        xb.push_back(x+ceil((lbo-x)/lengthX)*lengthX);
		else
			xb.push_back(x);
                if(y>nbo)
                        yb.push_back(y-ceil((y-nbo)/lengthY)*lengthY);
                else if(y<sbo)
                        yb.push_back(y+ceil((sbo-y)/lengthY)*lengthY);
                else
                        yb.push_back(y);
		return -1; // outside relocate particle	
	}

	int count=0;
	for(int ix=-1;ix<2;ix++)
		for(int iy=-1;iy<2;iy++)
		{
			if (ix==0 && iy==0) continue;
			double tx=x+ix*lengthX;
			double ty=y+iy*lengthY;
			if (tx<rbl && tx>lbl && ty<nbl &&ty>sbl)
			{
				xb.push_back(tx);
				yb.push_back(ty);
				pressureb.push_back(pressure);
				vxb.push_back(vx);
				vyb.push_back(vy);
				count++;
			}		
		}

	return count;
}
