#ifndef __BOUNDARY_NOZZLE_H__
#define __BOUNDARY_NOZZLE_H__

#include "boundary.h"
#include "eos.h"
#include "particle_data.h"
#include <vector>

class NozzleInflowBoundary: public Boundary {
public:
	NozzleInflowBoundary();
	virtual ~NozzleInflowBoundary() {};
	virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
	double left;
	double right;
	double radius;
        double Uinflow;
	double Pinflow;
	double Vinflow;
};

class NozzleInflowFixPressureBoundary: public Boundary {
public:
        NozzleInflowFixPressureBoundary();
        virtual ~NozzleInflowFixPressureBoundary() {};
        virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
        double left;
        double right;
	double average_rightlimit;
        double radius;
        double Uinflow;
	double OldUinflow;
        double Pinflow;
        double Vinflow;
	double Pinitial;
	double Vinitial;
};

class Nozzle3DInflowBoundary: public Boundary {
public:
        Nozzle3DInflowBoundary();
        virtual ~Nozzle3DInflowBoundary() {};
        virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
        double left;
        double right;
        double radius;
        double Uinflow;
        double Pinflow;
        double Vinflow;
};

class Nozzle3DInflowFixPressureBoundary: public Boundary {
public:
        Nozzle3DInflowFixPressureBoundary();
        virtual ~Nozzle3DInflowFixPressureBoundary() {};
        virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
        double left;
        double right;
        double average_rightlimit;
        double radius;
        double Uinflow;
        double OldUinflow;
        double Pinflow;
        double Vinflow;
        double Pinitial;
        double Vinitial;
};

class Nozzle2DSimpleSolidBoundary: public Boundary {
public:
	Nozzle2DSimpleSolidBoundary();
	virtual ~Nozzle2DSimpleSolidBoundary() {};
	virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
	double radius;
	double thickness;
	double bo;
};

class Nozzle2DSolidBoundary: public Boundary {
public:
        Nozzle2DSolidBoundary();
        virtual ~Nozzle2DSolidBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double x0,r0,x1,r1,x2,r2,x3,r3,x4,r4;
        double thickness;
	double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class Nozzle2DBNLSolidBoundary: public Boundary {
public:
        Nozzle2DBNLSolidBoundary();
        virtual ~Nozzle2DBNLSolidBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double x0,r0,x1,r1,x2,r2,x3,r3,x4,r4;
        double thickness;
        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class Nozzle3DSimpleSolidBoundary: public Boundary {
public:
        Nozzle3DSimpleSolidBoundary();
        virtual ~Nozzle3DSimpleSolidBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double radius;
        double thickness;
        double bo;
	double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d
};

class Nozzle3DSolidRightBoundary: public Boundary {
public:
        Nozzle3DSolidRightBoundary();
        virtual ~Nozzle3DSolidRightBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double x0,r0;
        double thickness;
        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class Nozzle3DSolidBoundary: public Boundary {
public:
        Nozzle3DSolidBoundary();
        virtual ~Nozzle3DSolidBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double x0,r0,x1,r1,x2,r2,x3,r3,x4,r4;
        double thickness;
        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class Nozzle3DBNLSolidBoundary: public Boundary {
public:
        Nozzle3DBNLSolidBoundary();
        virtual ~Nozzle3DBNLSolidBoundary() {};
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
private:
        double x0,r0,x1,r1,x2,r2,x3,r3,x4,r4;
        double thickness;
        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

class NozzleOutflowBoundary: public Boundary {
public:
        NozzleOutflowBoundary();
        virtual ~NozzleOutflowBoundary() {};
        virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
	double xmin,xmax,ymin,ymax;
};

class Nozzle3DOutflowBoundary: public Boundary {
public:
        Nozzle3DOutflowBoundary();
        virtual ~Nozzle3DOutflowBoundary() {};
        virtual int UpdateInflowBoundary(ParticleData *ParticleData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
        virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb,
        std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb){return 0;};
private:
        double xmin,xmax,rmax;
};
#endif
