/**
 * \file   boundary_powder_target.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the solid boundary in the power target simulation
 *
 *
 * \author R.S.
 *
 * \version 1.0 
 *
 * \date 2015/08/08
 *
 * Created on: 2015/08/08 
 *
 */

#ifndef __BOUNDARY_POWDER_TARGET_H__
#define __BOUNDARY_POWDER_TARGET_H__

#include "boundary.h"
#include <vector>

/**
 * \class PowderTarget2DSolidBoundary
 * 
 * \brief Function object (level set function) for generating a 2D powder target boundary
 *
 * \author R.S.
 *
 * \version 1.0 
 *
 * \date 2015/08/08 
 *
 * Created on: 2015/08/08
 *
 */
class PowderTarget2DSolidBoundary: public Boundary {
public:
	/// constructor
	PowderTarget2DSolidBoundary();

	/// destructor
	virtual ~PowderTarget2DSolidBoundary() {}
	
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
	double R; // container radius	
	double L; // container thickness
	double H; // container height

	double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d
	void reflectRight(double x, double y, double pressure, double vx, double vy,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& pressureb, 
	std::vector<double>& vxb,	std::vector<double>& vyb);

	void reflectLeft(double x, double y, double pressure, double vx, double vy,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& pressureb, 
	std::vector<double>& vxb,	std::vector<double>& vyb);	
};


#endif // __BOUNDARY_POWDER_TARGET_H__
