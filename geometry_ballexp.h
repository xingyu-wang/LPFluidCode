/**
 * \file   geometry_ballexp.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the 3D ballexp splash simulation
 *
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

#ifndef __GEOMETRY_BALLEXP_H__
#define __GEOMETRY_BALLEXP_H__

#include "geometry.h"



/**
 * \class Ballexp3D
 * 
 * \brief Function object (level set function) for generating a 3D ballexp geometry
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
class Ballexp3D: public Geometry {
public:
	/// constructor
	Ballexp3D();

	/// destructor
	virtual ~Ballexp3D() {}
	
	/**
	 * \brief Level set function of a 3D ballexp 
	 *        along the x-axis with small perturbations on the surface 
	 *
	 * The level set function is 
	 * \f$
	 *      
	 * \f$ 
	 *
	 * \param [in] x  The x-coordinate
	 * \param [in] y  The y-coordinate
	 * \param [in] z  The z-coordinate
	 * \return  \c true if (x,y,z) is inside the level set function; \c false otherwise
	 */
	virtual bool operator()(double x, double y, double z) const;	
	
	/**
	 * \brief       Calculates the bounding box  
	 *
	 * \param [out] xmin  The minimum in x-coordinate
	 * \param [out] xmax  The maximum in x-coordinate
	 * \param [out] ymin  The minimum in y-coordinate
	 * \param [out] ymax  The maximum in y-coordinate
	 * \param [out] zmin  The minimum in z-coordinate
	 * \param [out] zmax  The maximum in z-coordinate
	 * \return      None
	 */
	virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
	double radius;	
	double xCen;
	double yCen;
	double zCen;
};


#endif // __GEOMETRY_BALLEXP_H__
