/**
 * \file   boundary_powder_target_3d.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the solid boundary in the 3D power target simulation
 *
 *
 * \author Sijia
 *
 * \version 1.0 
 *
 * \date 2016/03/15
 *
 * Created on: 2016/03/15
 *
 */

#ifndef __BOUNDARY_POWDER_TARGET_3D_H__
#define __BOUNDARY_POWDER_TARGET_3D_H__

#include "boundary.h"
#include <vector>

/**
 * \class PowderTarget2DSolidBoundary
 * 
 * \brief Function object (level set function) for generating a 3D powder target boundary
 *
 * \author Sijia
 *
 * \version 1.0 
 *
 * \date 2015/08/08 
 *
 * Created on: 2015/08/08
 *
 */
class PowderTarget3DSolidBoundary: public Boundary {
public:
	/// constructor
	PowderTarget3DSolidBoundary();

	/// destructor
	virtual ~PowderTarget3DSolidBoundary() {}
	
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
	double length; // container length

	double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

	double ip3d(double x1, double y1, double z1, double x2, double y2, double z2){ return x1*x2 + y1*y2 + z1*z2;} // inner product 3d


	void reflectRight(double x, double y, double z, double pressure, double vx, double vy, double vz,
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectLeft(double x, double y, double z, double pressure, double vx, double vy, double vz,
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);    

	void reflectSemiCircle(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectFront(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);
	void reflectBack(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectLeftFront(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);
    
	void reflectRightFront(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectSemiCircleFront(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);


	void reflectLeftBack(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectRightBack(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);

	void reflectSemiCircleBack(double x, double y, double z, double pressure, double vx, double vy,  double vz,  
	std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb, 
        std::vector<double>& vxb,	std::vector<double>& vyb, std::vector<double>& vzb);



};


#endif // __BOUNDARY_POWDER_TARGET_3D_H__
