/**
 * \file   state_jet.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects in the 3D jet simulation
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/06/13
 *
 * Created on: 2015/06/13 
 *
 */

#ifndef __STATE_JET_H__
#define __STATE_JET_H__

#include "state.h"

#include <math.h>


/**
 * \class Jet3DState
 * 
 * \brief A class that implements the initial state of 3D jet splash
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
class Jet3DState: public State {
public:
	/// constructor
	Jet3DState();

	/// destructor
	virtual ~Jet3DState() {};
	
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
	double m_fVelX, m_fVelY, m_fVelZ;
	double m_fPressure;
	double m_fPMax; // g/(cm*ms^2)
	double m_fMagnitude; // magnitude of perturbation on initial pressure
	double m_fAngle; // angle of proton beam
};






/**
 * \class Jet1DState
 * 
 * \brief A class that implements the initial state of 1D jet splash
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/25 
 *
 * Created on: 2015/06/25
 * 
 *
 */
class Jet1DState: public State {
public:
	/// constructor
	Jet1DState();

	/// destructor
	virtual ~Jet1DState() {};
	
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
	double m_fVelX;
	double m_fPressure;
	double m_fPCenX; ///< center of the Gaussian pressure profile
	double m_fPPeak; ///< peak of pressure value
	double m_fPCoeff; ///< Gaussian coefficient	
};

class Jet1DLeftState: public State {
public:
        Jet1DLeftState();

        virtual ~Jet1DLeftState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
        double m_fPressure;
};

class Jet1DRightState: public State {
public:
        Jet1DRightState();

        virtual ~Jet1DRightState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
        double m_fPressure;
};

class Jet1DCenterState: public State {
public:
        Jet1DCenterState();

        virtual ~Jet1DCenterState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
        double m_fPressure;
};

class Jet1DLaterState: public State {
public:
        Jet1DLaterState();

        virtual ~Jet1DLaterState() {};

        virtual double pressure(double x, double y, double z);

        virtual double density(double x, double y, double z);

        virtual void velocity(double x, double y, double z, double& vX, double& vY, double& vZ);
private:
        double m_fDen;
        double m_fVelX;
        double m_fPressure;
	double m_fDenPost;
	double m_fPressurePost;
	double m_fShockSpeed;
	double m_fTime;
};


/**
 * \class Jet2DExpansionState
 * 
 * \brief A class that implements the initial state of 2D jet expansion
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
class Jet2DExpansionState: public State {
public:
	/// constructor
	Jet2DExpansionState();

	/// destructor
	virtual ~Jet2DExpansionState() {};
	
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
	double m_fVelocity; 
};


/**
 * \class Jet3DExpansionState
 * 
 * \brief A class that implements the initial state of 3D jet expansion
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
class Jet3DExpansionState: public State {
public:
	/// constructor
	Jet3DExpansionState();

	/// destructor
	virtual ~Jet3DExpansionState() {};
	
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
	double m_fVelocity; 
};


/**
 * \class Jet2DMergeUpperState
 * 
 * \brief A class that implements the initial state of 2D upper jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/04 
 *
 * Created on: 2015/11/04
 * 
 *
 */
class Jet2DMergeUpperState: public State {
public:
	/// constructor
	Jet2DMergeUpperState();

	/// destructor
	virtual ~Jet2DMergeUpperState() {};
	
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
	double m_fVelocity; 
  double alpha;
};

/**
 * \class Jet2DMergeLowerState
 * 
 * \brief A class that implements the initial state of 2D lower jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/04 
 *
 * Created on: 2015/11/04
 * 
 *
 */
class Jet2DMergeLowerState: public State {
public:
	/// constructor
	Jet2DMergeLowerState();

	/// destructor
	virtual ~Jet2DMergeLowerState() {};
	
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
	double m_fVelocity; 
  double alpha;
};


/**
 * \class Jet2DCollisionState
 * 
 * \brief A class that implements the initial state of 2D jet collision
 *
 * \author Xingyu Wang(xingyu.wang@stonybrook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/09 
 *
 * Created on: 2015/11/09
 * 
 *
 */
class Jet2DCollisionState: public State {
public:
	/// constructor
	Jet2DCollisionState();

	/// destructor
	virtual ~Jet2DCollisionState() {};
	
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
	double m_fVelocity; 
};

/**
 * \class Jet2DMergeState
 * 
 * \brief A class that implements the initial state of 2D upper and lower jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/10
 *
 * Created on: 2015/11/10
 * 
 *
 */
class Jet2DMergeState: public State {
public:
	/// constructor
	Jet2DMergeState();

	/// destructor
	virtual ~Jet2DMergeState() {};
	
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
	double m_fVelocity; 
  double alpha;
};

#endif //__STATE_JET_H__
