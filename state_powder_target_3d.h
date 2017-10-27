/**
 * \file   state_powder_target_3d.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the state of fluid objects for simulations of the CERN powder target
 *
 * \author sijia
 *
 * \version 1.0 
 *
 * \date 2016/03/15
 *
 * Created on: 2016/03/15
 *
 */

#ifndef __STATE_POWDER_TARGET_3D_H__
#define __STATE_POWDER_TARGET_3D_H__

#include "state.h"


/**
 * \class PowderTarget2DState
 * 
 * \brief A class that implements the initial state of 3D CERN powder target
 *
 * \author sijia
 *
 * \version 1.0 
 *
 * \date 2016/03/15
 *
 * Created on: 2016/03/15
 * 
 *
 */
class PowderTarget3DState: public State {
public:
	/// constructor
	PowderTarget3DState();

	/// destructor
	virtual ~PowderTarget3DState() {};
	
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

#endif //__POWDER_TARGET_3D_H__
