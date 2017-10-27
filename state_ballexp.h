/**
 * \file   state_ballexp.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects in the 3D ballexp simulation
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/08/11
 *
 * Created on: 2015/08/11 
 *
 */

#ifndef __STATE_BALLEXP_H__
#define __STATE_BALLEXP_H__

#include "state.h"


/**
 * \class Ballexp3DState
 * 
 * \brief A class that implements the initial state of 3D ballexp 
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/08/11
 *
 * Created on: 2015/08/11
 * 
 *
 */
class Ballexp3DState: public State {
public:
	/// constructor
	Ballexp3DState();

	/// destructor
	virtual ~Ballexp3DState() {};
	
	/**
	 * \brief         Specifies constant pressure as specified in construtor implementation 
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constnat pressure value 
	 */
	virtual double pressure(double x, double y, double z);
	
	/**
	 * \brief         Specifies a constant value as specified in constructor implementation
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return        A constant density value
	 */
	virtual double density(double x, double y, double z);
	
	/**
	 * \brief           Specifies constant velocity to the right-hand-side along the x-coordinate 
	 *                  with magnitude specified in constructor implementation
	 * \param [in]  x   The x-coordinate
	 * \param [in]  y   The y-coordinate
	 * \param [in]  z   The z-coordinate
	 * \param [out] vX  A constant value (>0)   
	 * \param [out] vY  zero
	 * \param [out] vZ  zero
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDen;
	double m_fPressure;
};

class Ballrotate3DState: public State {
public:
        /// constructor
        Ballrotate3DState();

        /// destructor
        virtual ~Ballrotate3DState() {};

        /**
         * \brief         Specifies constant pressure as specified in construtor implementation 
         * \param [in] x  The x-coordinate
         * \param [in] y  The y-coordinate
         * \param [in] z  The z-coordinate
         * \return        A constnat pressure value 
         */
        virtual double pressure(double x, double y, double z);

        /**
         * \brief         Specifies a constant value as specified in constructor implementation
         * \param [in] x  The x-coordinate
         * \param [in] y  The y-coordinate
         * \param [in] z  The z-coordinate
         * \return        A constant density value
         */
        virtual double density(double x, double y, double z);
        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fPressure;
};

class Ballpressurewave3DState: public State {
public:
        /// constructor
        Ballpressurewave3DState();

        /// destructor
        virtual ~Ballpressurewave3DState() {};

        /**
         * \brief         Specifies constant pressure as specified in construtor implementation 
         * \param [in] x  The x-coordinate
         * \param [in] y  The y-coordinate
         * \param [in] z  The z-coordinate
         * \return        A constnat pressure value 
         */
        virtual double pressure(double x, double y, double z);

        /**
         * \brief         Specifies a constant value as specified in constructor implementation
         * \param [in] x  The x-coordinate
         * \param [in] y  The y-coordinate
         * \param [in] z  The z-coordinate
         * \return        A constant density value
         */
        virtual double density(double x, double y, double z);
        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fPressure;
};
#endif //__STATE_BALLEXP_H__
