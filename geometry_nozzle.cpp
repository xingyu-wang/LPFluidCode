#include "geometry_nozzle.h"
#include <iostream>
#include <cmath>

Nozzle2DSimple::Nozzle2DSimple():left(0), radius(0.75e-3), length(0.2e-3){}

bool Nozzle2DSimple::operator()(double x, double y, double z) const {

	return x>left && x<left+length && y>-radius && y<radius;
}

void Nozzle2DSimple::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = left;
        xmax = left+length;
        ymin = -radius;
        ymax = radius;
        zmin = 0;
        zmax = 0;
}

Nozzle2D::Nozzle2D():left(0), radius(0.75e-3), length(0.2e-3), k(0.125){}

bool Nozzle2D::operator()(double x, double y, double z) const {

        return x>left && x<left+length && y>-radius+k*x && y<radius-k*x;
}

void Nozzle2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = left;
        xmax = left+length;
        ymin = -radius;
        ymax = radius;
        zmin = 0;
        zmax = 0;
}

Nozzle2DComplete::Nozzle2DComplete():x0(0),r0(0.75e-3),x1(4e-3),r1(0.25e-3),x2(7e-3),r2(0.5e-3),x3(8e-3),r3(1.65e-3){}

bool Nozzle2DComplete::operator()(double x, double y, double z) const {
	if(x<x0) return 0;
	if(x<x1) return fabs(y)<(r0+(x-x0)*(r1-r0)/(x1-x0));
        if(x<x2) return fabs(y)<(r1+(x-x1)*(r2-r1)/(x2-x1));
        if(x<x3) return fabs(y)<(r2+(x-x2)*(r3-r2)/(x3-x2));
	return 0;
}

void Nozzle2DComplete::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = x0;
        xmax = x3;
        ymin = -r3;
        ymax = r3;
        zmin = 0;
        zmax = 0;
}

Nozzle3D::Nozzle3D():left(0), radius(0.75e-3), length(0.4e-3), k(0.125){}

bool Nozzle3D::operator()(double x, double y, double z) const {

        return x>left && x<left+length && sqrt(y*y+z*z)<radius-1*k*x;
}

void Nozzle3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = left;
        xmax = left+length;
        ymin = -radius;
        ymax = radius;
        zmin = -radius;
        zmax = radius;
}

Nozzle2DRothe::Nozzle2DRothe():left(0), radius(8.3e-3), length(5e-3), k(0.125){}

bool Nozzle2DRothe::operator()(double x, double y, double z) const {

        return x>left && x<left+length && y>-radius+0*k*x && y<radius-0*k*x;
}

void Nozzle2DRothe::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = left;
        xmax = left+length;
        ymin = -radius;
        ymax = radius;
        zmin = 0;
        zmax = 0;
}

Nozzle3DRothe::Nozzle3DRothe():left(0), radius(8.3e-3), length(5e-3), k(0.125){}

bool Nozzle3DRothe::operator()(double x, double y, double z) const {

        return x>left && x<left+length && sqrt(y*y+z*z)<radius-0*k*x;
}

void Nozzle3DRothe::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
        xmin = left;
        xmax = left+length;
        ymin = -radius;
        ymax = radius;
        zmin = -radius;
        zmax = radius;
}
