#include "state_nozzle.h"
#include <iostream>

NozzleState::NozzleState():
m_fDen(0.08), m_fVelX(0), m_fVelY(0), m_fVelZ(0), m_fPressure(1e5) {}


double NozzleState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double NozzleState::density(double x, double y, double z) {
        return m_fDen;
}

void NozzleState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = m_fVelX;
        vY = m_fVelY;
        vZ = m_fVelZ;
}

NozzleRotheState::NozzleRotheState():
m_fDen(0.0053198), m_fVelX(0), m_fVelY(0), m_fVelZ(0), m_fPressure(473.86) {}


double NozzleRotheState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double NozzleRotheState::density(double x, double y, double z) {
        return m_fDen;
}

void NozzleRotheState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = m_fVelX;
        vY = m_fVelY;
        vZ = m_fVelZ;
}
