#include "geometry_shocktube.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of Shocktube2D
////////////////////////////////////////////////////////////////////////////////////////

Shocktube2D::Shocktube2D():lengthX(5), lengthY(2) {}

bool Shocktube2D::operator()(double x, double y, double z) const {	

	return (x<=lengthX/2. && x>=-lengthX/2. && y<=lengthY/2. && y>=-lengthY/2.); 	
}

void Shocktube2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = -lengthX/2.;
	xmax = lengthX/2.;
	ymin = -lengthY/2.;
	ymax = lengthY/2.;
	zmin = 0;
	zmax = 0;
}

BigShocktube2D::BigShocktube2D():lengthX(10), lengthY(5) {}

bool BigShocktube2D::operator()(double x, double y, double z) const {

        return (x<=lengthX/2. && x>=-lengthX/2. && y<=lengthY/2. && y>=-lengthY/2.);
}

void BigShocktube2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = -lengthX/2.;
        xmax = lengthX/2.;
        ymin = -lengthY/2.;
        ymax = lengthY/2.;
        zmin = 0;
        zmax = 0;
}

TPShocktube2D::TPShocktube2D():lengthX(7), lengthY(3) {}

bool TPShocktube2D::operator()(double x, double y, double z) const {

        return (x<=lengthX && x>=0 && y<=lengthY && y>=0);
}

void TPShocktube2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = 0;
        xmax = lengthX;
        ymin = 0;
        ymax = lengthY;
        zmin = 0;
        zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Shocktube2D
////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
// Start of Shocktube2DLeft
////////////////////////////////////////////////////////////////////////////////////////

Shocktube2DLeft::Shocktube2DLeft():lengthX(10), lengthY(5) {}

bool Shocktube2DLeft::operator()(double x, double y, double z) const {

        return (x<=0.00001 && x>=-lengthX/2. && y<=lengthY/2. && y>=-lengthY/2.);
}

void Shocktube2DLeft::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = -lengthX/2.;
        xmax = 0.00001;
        ymin = -lengthY/2.;
        ymax = lengthY/2.;
        zmin = 0;
        zmax = 0;
}

Shocktube2DRight::Shocktube2DRight():lengthX(5), lengthY(2) {}

bool Shocktube2DRight::operator()(double x, double y, double z) const {

        return (x<=lengthX/2. && x>=0.00001 && y<=lengthY/2. && y>=-lengthY/2.);
}

void Shocktube2DRight::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = 0.00001;
        xmax = lengthX/2.;
        ymin = -lengthY/2.;
        ymax = lengthY/2.;
        zmin = 0;
        zmax = 0;
}

RayleighTaylor2D::RayleighTaylor2D():lengthX(1), lengthY(1) {}

bool RayleighTaylor2D::operator()(double x, double y, double z) const {

        return (x<=lengthX/2. && x>=-lengthX/2. && y<=lengthY/2. && y>=-lengthY/2.);
}

void RayleighTaylor2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = -lengthX/2.;
        xmax = lengthX/2.;
        ymin = -lengthY/2.;
        ymax = lengthY/2.;
        zmin = 0;
        zmax = 0;
}

DamBreak2D::DamBreak2D():lengthX(0.5), lengthY(0.5) {}

bool DamBreak2D::operator()(double x, double y, double z) const {

        return (x<=lengthX && x>=0 && y<=lengthY && y>=0);
}

void DamBreak2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = 0.;
        xmax = lengthX;
        ymin = 0.;
        ymax = lengthY;
        zmin = 0;
        zmax = 0;
}


BoundaryTest2D::BoundaryTest2D():lengthX(1), lengthY(1) {}

bool BoundaryTest2D::operator()(double x, double y, double z) const {

        return (x<=lengthX && x>=0 && y<=lengthY && y>=0);
}

void BoundaryTest2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = 0.;
        xmax = lengthX;
        ymin = 0.;
        ymax = lengthY;
        zmin = 0;
        zmax = 0;
}

KelvinHelmholtz2D::KelvinHelmholtz2D():lengthX(1), lengthY(1) {}

bool KelvinHelmholtz2D::operator()(double x, double y, double z) const {

        return (x<=lengthX && x>=0 && y<=lengthY && y>=0);
}

void KelvinHelmholtz2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = 0;
        xmax = lengthX;
        ymin = 0;
        ymax = lengthY;
        zmin = 0;
        zmax = 0;
}
