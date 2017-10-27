#include "geometry_jet.h"
#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet3D
////////////////////////////////////////////////////////////////////////////////////////

Jet3D::Jet3D():radius(0.5), length(16), xCen(0), yCen(0), zCen(0), magnitude(0.3) {}

bool Jet3D::operator()(double x, double y, double z) const {	

	double distance = sqrt((y-yCen)*(y-yCen)+(z-zCen)*(z-zCen));
	
	/*
	//surface purturbation
	double tr = atan(y/z);	 
	double tmp = double(rand()) / double(RAND_MAX);
	double ratio = double(rand()) / double(RAND_MAX);	
	if(tmp<=0.5) ratio*=-1; 	
	//std::cout<<"tmp="<<tmp<<" ratio="<<ratio<<std::endl;
	double r = radius/(1.0+(0.08*cos(10.0*tr)+0.06*cos(20.0*tr)+0.07*sin(60.0*x/4.0)))+radius*magnitude*ratio;
	*/

	//surface purturbation
	double tr = atan(y/z);	 	
	double r = radius/(1.0+(0.08*cos(10.0*tr)+0.06*cos(20.0*tr)+0.07*sin(60.0*x/4.0)));

	return ( (x<=xCen+length/2. && x>=xCen-length/2.) && (distance <= r) ); 	
}

void Jet3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-length/2.;
	xmax = xCen+length/2.;
	ymin = yCen-(1+3*magnitude)*radius;
	ymax = yCen+(1+3*magnitude)*radius;
	zmin = zCen-(1+3*magnitude)*radius;
	zmax = zCen+(1+3*magnitude)*radius;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet3D
////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2D
////////////////////////////////////////////////////////////////////////////////////////

Jet2D::Jet2D():radius(0.5), length(16), xCen(0), yCen(0), magnitude(0.3) {}

bool Jet2D::operator()(double x, double y, double z) const {	

	double distance = sqrt((y-yCen)*(y-yCen));	

	//surface purturbation
	//double tr = atan(y/z);	 	
	double tr = atan(0);
	double r = radius/(1.0+(0.08*cos(10.0*tr)+0.06*cos(20.0*tr)+0.07*sin(60.0*x/4.0)));

	return ( (x<=xCen+length/2. && x>=xCen-length/2.) && (distance <= r) ); 	
}

void Jet2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-length/2.;
	xmax = xCen+length/2.;
	ymin = yCen-(1+3*magnitude)*radius;
	ymax = yCen+(1+3*magnitude)*radius;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2D
////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet1D
////////////////////////////////////////////////////////////////////////////////////////

Jet1D::Jet1D():length(16), xCen(0) {}

bool Jet1D::operator()(double x, double y, double z) const {		

	return (x<=xCen+length/2. && x>=xCen-length/2.); 	
}

void Jet1D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-length/2.;
	xmax = xCen+length/2.;
	ymin = 0;
	ymax = 0;
	zmin = 0;
	zmax = 0;
}

Jet1DLeft::Jet1DLeft():length(0.08), xCen(-0.25) {}

bool Jet1DLeft::operator()(double x, double y, double z) const {

        return (x<=xCen+length/2. && x>=xCen-length/2.);
}

void Jet1DLeft::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = xCen-length/2.;
        xmax = xCen+length/2.;
        ymin = 0;
        ymax = 0;
        zmin = 0;
        zmax = 0;
}

Jet1DRight::Jet1DRight():length(0.08), xCen(0.25) {}

bool Jet1DRight::operator()(double x, double y, double z) const {

        return (x<=xCen+length/2. && x>=xCen-length/2.);
}

void Jet1DRight::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = xCen-length/2.;
        xmax = xCen+length/2.;
        ymin = 0;
        ymax = 0;
        zmin = 0;
        zmax = 0;
}

Jet1DCenter::Jet1DCenter():length(0.16), xCen(0.0) {}

bool Jet1DCenter::operator()(double x, double y, double z) const {

        return (x<=xCen+length/2. && x>=xCen-length/2.);
}

void Jet1DCenter::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = xCen-length/2.;
        xmax = xCen+length/2.;
        ymin = 0;
        ymax = 0;
        zmin = 0;
        zmax = 0;
}
////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2D
////////////////////////////////////////////////////////////////////////////////////////


















////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DExpansion
////////////////////////////////////////////////////////////////////////////////////////

Jet2DExpansion::Jet2DExpansion():radius(2.5), length(20), xCen(0), yCen(0) {}

bool Jet2DExpansion::operator()(double x, double y, double z) const {	
	
	double l = length/2., yCenU=yCen+l, yCenD=yCen-l;

	if(y<=yCenU && y>=yCenD && x<=xCen+radius && x>=xCen-radius) return 1;
	else if(y>yCenU && pow(x-xCen,2)+pow(y-yCenU,2) <= pow(radius,2)) return 1;
	else if(y<yCenD && pow(x-xCen,2)+pow(y-yCenD,2) <= pow(radius,2)) return 1;
	else return 0;

}

void Jet2DExpansion::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-radius;
	xmax = xCen+radius;
	ymin = yCen-length/2.-radius;
	ymax = yCen+length/2.+radius;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DExpansion
////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet3DExpansion
////////////////////////////////////////////////////////////////////////////////////////

Jet3DExpansion::Jet3DExpansion():radius(2.5), length(10), xCen(0), yCen(0), zCen(0) {}

bool Jet3DExpansion::operator()(double x, double y, double z) const {	
	
	double l = length/2., yCenU=yCen+l, yCenD=yCen-l;
	double xs = (x-xCen)*(x-xCen), zs = (z-zCen)*(z-zCen);
	double ysu = (y-yCenU)*(y-yCenU), ysd = (y-yCenD)*(y-yCenD);

	if(y<=yCenU && y>=yCenD && xs+zs <= radius*radius) return 1;
	else if(y>yCenU && xs+ysu+zs<=radius*radius) return 1;
	else if(y<yCenD && xs+ysd+zs<=radius*radius) return 1;
	else return 0;

}

void Jet3DExpansion::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-radius;
	xmax = xCen+radius;
	ymin = yCen-length/2.-radius;
	ymax = yCen+length/2.+radius;
	zmin = zCen-radius;
	zmax = zCen+radius;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet3DExpansion
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMergeUpper
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMergeUpper::Jet2DMergeUpper():radius(2.5), length(20), xCen(0), yCen(23), alpha(-12.0/180.0*M_PI) {}

bool Jet2DMergeUpper::operator()(double x, double y, double z) const {	

  double l = length/2., pCen=xCen*sin(alpha)-yCen*cos(alpha), qCen=xCen*cos(alpha)+yCen*sin(alpha);
  double qCenU=qCen+l, qCenD=qCen-l;	
  double p=x*sin(alpha)-y*cos(alpha), q=x*cos(alpha)+y*sin(alpha);
	if(q<=qCenU && q>=qCenD && p<=pCen+radius && p>=pCen-radius) return 1;
	else if(q>qCenU && pow(p-pCen,2)+pow(q-qCenU,2) <= pow(radius,2)) return 1;
	else if(q<qCenD && pow(p-pCen,2)+pow(q-qCenD,2) <= pow(radius,2)) return 1;
	else return 0;

}

void Jet2DMergeUpper::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-fabs(cos(alpha))*length/2.-radius;
	xmax = xCen+fabs(cos(alpha))*length/2.+radius;
	ymin = yCen-fabs(sin(alpha))*length/2.-radius;
	ymax = yCen+fabs(sin(alpha))*length/2.+radius;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMergeUpper
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMergeLower
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMergeLower::Jet2DMergeLower():radius(2.5), length(20), xCen(0), yCen(-23), alpha(12.0/180.0*M_PI) {}

bool Jet2DMergeLower::operator()(double x, double y, double z) const {	

  double l = length/2., pCen=xCen*sin(alpha)-yCen*cos(alpha), qCen=xCen*cos(alpha)+yCen*sin(alpha);
  double qCenU=qCen+l, qCenD=qCen-l;	
  double p=x*sin(alpha)-y*cos(alpha), q=x*cos(alpha)+y*sin(alpha);
	if(q<=qCenU && q>=qCenD && p<=pCen+radius && p>=pCen-radius) return 1;
	else if(q>qCenU && pow(p-pCen,2)+pow(q-qCenU,2) <= pow(radius,2)) return 1;
	else if(q<qCenD && pow(p-pCen,2)+pow(q-qCenD,2) <= pow(radius,2)) return 1;
	else return 0;

}

void Jet2DMergeLower::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-fabs(cos(alpha))*length/2.-radius;
	xmax = xCen+fabs(cos(alpha))*length/2.+radius;
	ymin = yCen-fabs(sin(alpha))*length/2.-radius;
	ymax = yCen+fabs(sin(alpha))*length/2.+radius;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMergeLower
////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DCollision
////////////////////////////////////////////////////////////////////////////////////////

Jet2DCollision::Jet2DCollision():xLength(5), yLength(20), xLeftCen(-5), xRightCen(5), yCen(0) {}

bool Jet2DCollision::operator()(double x, double y, double z) const {	
	
	if(y<=(yCen+yLength/2.) && y>=(yCen-yLength/2.) && x<=(xLeftCen+xLength/2.) && x>=(xLeftCen-xLength/2.)) return 1;
	else if(y<=(yCen+yLength/2.) && y>=(yCen-yLength/2.) && x<=(xRightCen+xLength/2.) && x>=(xRightCen-xLength/2.)) return 1;
	else return 0;

}

void Jet2DCollision::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xLeftCen-xLength/2.;
	xmax = xRightCen+xLength/2.;
	ymin = yCen-yLength/2.;
	ymax = yCen+yLength/2.;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DCollision
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMerge
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMerge::Jet2DMerge():radius(2.5), length(20), xCen(0), yCen(23), alpha(-12.0/180.0*M_PI) {}

bool Jet2DMerge::operator()(double x, double y, double z) const {	

  double l = length/2., pCen=xCen*sin(alpha)-yCen*cos(alpha), qCen=xCen*cos(alpha)+yCen*sin(alpha);
  double qCenU=qCen+l, qCenD=qCen-l;	
  double p=x*sin(alpha)-fabs(y)*cos(alpha), q=x*cos(alpha)+fabs(y)*sin(alpha);
	if(q<=qCenU && q>=qCenD && p<=pCen+radius && p>=pCen-radius) return 1;
	else if(q>qCenU && pow(p-pCen,2)+pow(q-qCenU,2) <= pow(radius,2)) return 1;
	else if(q<qCenD && pow(p-pCen,2)+pow(q-qCenD,2) <= pow(radius,2)) return 1;
	else return 0;

}

void Jet2DMerge::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-fabs(cos(alpha))*length/2.-radius;
	xmax = xCen+fabs(cos(alpha))*length/2.+radius;
	ymin = -yCen-fabs(sin(alpha))*length/2.-radius;
	ymax = yCen+fabs(sin(alpha))*length/2.+radius;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMerge
////////////////////////////////////////////////////////////////////////////////////////
