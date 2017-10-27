/**
 * \file   state_gresho.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects in the gresho simulation
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/08/06
 *
 * Created on: 2015/08/06 
 *
 */

#ifndef __STATE_GRESHO_H__
#define __STATE_GRESHO_H__

#include "state.h"


/**
 * \class Gresho2DState
 * 
 * \brief A class that implements the initial state of 2D gresho splash
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/13 
 *
 * Created on: 2015/06/13
 * 
 *
 */
class Gresho2DState: public State {
public:
	/// constructor
	Gresho2DState();

	/// destructor
	virtual ~Gresho2DState() {};
	
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
	double m_fDensity;	
};

class Yee2DState: public State {
public:
        /// constructor
        Yee2DState();

        /// destructor
        virtual ~Yee2DState() {};

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
};

class Sedov2DState: public State {
public:
        /// constructor
        Sedov2DState();

        /// destructor
        virtual ~Sedov2DState() {};

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
};


class Noh2DState: public State {
public:
        /// constructor
        Noh2DState();

        /// destructor
        virtual ~Noh2DState() {};

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
};


class ConvergentShock2DState: public State {
public:
        /// constructor
        ConvergentShock2DState();

        /// destructor
	virtual ~ConvergentShock2DState() {};
	/**
 * 	 * \brief         Specifies constant pressure as specified in construtor implementation 
 * 	 * \param [in] x  The x-coordinate
 * 	 * \param [in] y  The y-coordinate
 * 	 * \param [in] z  The z-coordinate
 * 	 * \return        A constnat pressure value 
 * 	 */
	virtual double pressure(double x, double y, double z);
	
	/**
 * 	 * \brief         Specifies a constant value as specified in constructor implementation
 * 	 * \param [in] x  The x-coordinate
 * 	 * \param [in] y  The y-coordinate
 * 	 * \param [in] z  The z-coordinate
 * 	 * \return        A constant density value
 * 	 */
	virtual double density(double x, double y, double z);
	
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDensity;
	double m_fPressure;
	double m_fVelocity;
};
#endif //__STATE_GRESHO_H__
