/**
 * \file   state_powder_target.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects for simulations of the CERN powder target
 *
 * \author Roman Samulyak
 *
 * \version 1.0 
 *
 * \date 2015/08/06
 *
 * Created on: 2015/08/06 
 *
 */

#ifndef __STATE_POWDER_TARGET_H__
#define __STATE_POWDER_TARGET_H__

#include "state.h"


/**
 * \class PowderTarget2DState
 * 
 * \brief A class that implements the initial state of 2D CERN powder target
 *
 * \author R.S.
 *
 * \version 1.0 
 *
 * \date 2015/06/13 
 *
 * Created on: 2015/06/13
 * 
 *
 */
class PowderTarget2DState: public State {
public:
	/// constructor
	PowderTarget2DState();

	/// destructor
	virtual ~PowderTarget2DState() {};
	
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
	double rho_0, p0, p_deposition, sigma_x, sigma_y;	
};

#endif //__POWDER_TARGET_H__
