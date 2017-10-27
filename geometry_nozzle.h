#ifndef __GEOMETRY_NOZZLE_H__
#define __GEOMETRY_NOZZLE_H__

#include "geometry.h"

#include <math.h>

class Nozzle2DSimple: public Geometry {
public:
        Nozzle2DSimple();
        virtual ~Nozzle2DSimple() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
	double left;
        double radius;
        double length;
};

class Nozzle2D: public Geometry {
public:
        Nozzle2D();
        virtual ~Nozzle2D() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double left;
        double radius;
        double length;
	double k;
};

class Nozzle2DComplete: public Geometry {
public:
        Nozzle2DComplete();
        virtual ~Nozzle2DComplete() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double x0;
	double r0;
	double x1;
	double r1;
	double x2;
	double r2;
	double x3;
	double r3;
};

class Nozzle3D: public Geometry {
public:
        Nozzle3D();
        virtual ~Nozzle3D() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double left;
        double radius;
        double length;
        double k;
};

class Nozzle2DRothe: public Geometry {
public:
        Nozzle2DRothe();
        virtual ~Nozzle2DRothe() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double left;
        double radius;
        double length;
        double k;
};

class Nozzle3DRothe: public Geometry {
public:
        Nozzle3DRothe();
        virtual ~Nozzle3DRothe() {}
        virtual bool operator()(double x, double y, double z) const;
        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double left;
        double radius;
        double length;
        double k;
};
#endif
