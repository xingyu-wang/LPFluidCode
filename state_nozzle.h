#ifndef __STATE_NOZZLE_H__
#define __STATE_NOZZLE_H__

#include "state.h"

#include <math.h>

class NozzleState: public State {
public:
        NozzleState();

        virtual ~NozzleState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
	double m_fVelY;
	double m_fVelZ;
        double m_fPressure;
};

class NozzleRotheState: public State {
public:
        NozzleRotheState();

        virtual ~NozzleRotheState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
        double m_fVelY;
        double m_fVelZ;
        double m_fPressure;
};
#endif
