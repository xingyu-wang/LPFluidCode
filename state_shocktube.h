/**
 * \file   state_shocktube.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects in the shocktube simulation
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

#ifndef __STATE_SHOCKTUBE_H__
#define __STATE_SHOCKTUBE_H__

#include "state.h"


/**
 * \class Shocktube2DState
 * 
 * \brief A class that implements the initial state of 2D shocktube splash
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
class Shocktube2DState: public State {
public:
	/// constructor
	Shocktube2DState();

	/// destructor
	virtual ~Shocktube2DState() {};
	
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
	double m_fDenL, m_fDenR;	
	double m_fPressureL, m_fPressureR;
};

class TPShocktube2DState: public State {
public:
	/// constructor
	TPShocktube2DState();

	/// destructor
	virtual ~TPShocktube2DState() {};
	
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
	double m_fDen1, m_fDen2, m_fDen3;	
	double m_fPressure1, m_fPressure2, m_fPressure3;
};

class SodShocktube2DLaterState: public State {
public:
        /// constructor
        SodShocktube2DLaterState();

        /// destructor
        virtual ~SodShocktube2DLaterState() {};

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
        double m_fDenL, m_fDenR;
        double m_fPressureL, m_fPressureR;
	double time;
};

/**
 * \class NormalShock2DState
 * 
 * \brief A class that implements the initial state of 2D normal shock
 *
 * \author Xingyu Wang (xingyu.wang@stonybrook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/12/09 
 *
 * Created on: 2015/12/09
 * 
 *
 */
class NormalShock2DState: public State {
public:
	/// constructor
	NormalShock2DState();

	/// destructor
	virtual ~NormalShock2DState() {};
	
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
	 * \param [out] vX  A velocity value depends on the sign of x 
	 * \param [out] vY  zero
	 * \param [out] vZ  zero
	 * \return None
	 */
	virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
	double m_fDenL, m_fDenR;	
	double m_fPressureL, m_fPressureR;
	double m_fVelocityL, m_fVelocityR;
};


/**
 * \class SimpleWave2DState
 * 
 * \brief A class that implements the initial state of 2D simple wave
 *
 * \author Xingyu Wang (xingyu.wang@stonybrook.edu)
 *
 * \version 1.0 
 *
 * \date 2016/06/16 
 *
 * Created on: 2016/06/16
 * 
 *
 */
class SimpleWave2DState: public State {
public:
        /// constructor
        SimpleWave2DState();

        /// destructor
        virtual ~SimpleWave2DState() {};

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
         * \param [out] vX  A velocity value depends on the sign of x 
         * \param [out] vY  zero
         * \param [out] vZ  zero
         * \return None
         */
        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDenL;
        double m_fPressureL;
        double m_fVPost;
	double m_fGamma;
	double m_fITime;
};

class SodShocktube2DState: public State {
public:
        /// constructor
        SodShocktube2DState();

        /// destructor
        virtual ~SodShocktube2DState() {};

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
        double m_fDenL, m_fDenR;
        double m_fPressureL, m_fPressureR;
};

class RayleighTaylor2DState: public State {
public:
        /// constructor
        RayleighTaylor2DState();

        /// destructor
        virtual ~RayleighTaylor2DState() {};

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
        double m_fDenU, m_fDenL;
        double m_fPressure, m_fWaveLength, m_fMagnitude, m_fG;
};

class RayleighTaylor3DState: public State {
public:
        /// constructor
        RayleighTaylor3DState();

        /// destructor
        virtual ~RayleighTaylor3DState() {};

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
        double m_fDenU, m_fDenL;
        double m_fPressure, m_fWaveLengthX, m_fWaveLengthY, m_fMagnitude, m_fG;
};

class DamBreak2DState: public State {
public:
        /// constructor
        DamBreak2DState();

        /// destructor
        virtual ~DamBreak2DState() {};

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
        double m_fPressure, m_fYLength, m_fG;
};

class BoundaryTest2DState: public State {
public:
        /// constructor
        BoundaryTest2DState();

        /// destructor
        virtual ~BoundaryTest2DState() {};
        virtual double pressure(double x, double y, double z);
        virtual double density(double x, double y, double z);
        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fPressure, m_fYLength, m_fG;
};

class KelvinHelmholtz2DState: public State {
public:
        /// constructor
        KelvinHelmholtz2DState();

        /// destructor
        virtual ~KelvinHelmholtz2DState() {};

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
        double m_fDen1, m_fDen2, m_fVel1, m_fVel2, m_fDelta, m_fDeltay;
        double m_fPressure, m_fWaveLength;
};

#endif //__STATE_SHOCKTUBE_H__
