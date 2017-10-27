/**
 * \file   boundary_solid_shocktube.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the solid boundary in the shocktube simulation
 *
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/08/08
 *
 * Created on: 2015/08/08 
 *
 */

#ifndef __BOUNDARY_KELVINHELMHOLTZ_H__
#define __BOUNDARY_KELVINHELMHOLTZ_H__

#include "boundary.h"
#include <vector>

/**
 * \class Shocktube2DSolidBoundary
 * 
 * \brief Function object (level set function) for generating a 2D shocktube solid boundary
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/08/08 
 *
 * Created on: 2015/08/08
 *
 */
class KelvinHelmholtz2DBoundary: public Boundary {
public:
	/// constructor
	KelvinHelmholtz2DBoundary();

	/// destructor
	virtual ~KelvinHelmholtz2DBoundary() {}
	
	/**
	 * \brief Get a boundary particle based on a fluid particle      
	 * \param [in] x  The x-coordinate of fluid particle
	 * \param [in] y  The y-coordinate of fluid particle
	 * \param [in] z  The z-coordinate of fluid particle
	 * \param [out] xb  The x-coordinate of boundary particle
	 * \param [out] yb  The y-coordinate of boundary particle
	 * \param [out] zb  The z-coordinate of boundary particle		
	 */
	virtual int operator()(double x, double y, double z, double pressure, double vx, double vy, double vz, 
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, 
	std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
		
		
private:
	double lengthX;
	double lengthY;	
	double thickness;

	double rb; ///< right boundary
	double lb; ///< left boundary
	double nb; ///< north boundary
	double sb; ///< south boundary
	
	double rbo; ///< outer right boundary
	double lbo; ///< outer left boundary
	double nbo; ///< outer north boundary
	double sbo; ///< outer south boundary


	double rbl; ///< right buffer limit
	double lbl; ///< left buffer limit
	double nbl; ///< north buffer limit
	double sbl; ///< south buffer limit
};

#endif		
