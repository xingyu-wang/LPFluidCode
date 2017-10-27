/**
 * \file   geometry_powder_target.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the powder target simulation
 *
 *
 * \author Roman Samulyak
 *
 * \version 1.0 
 *
 * \date 2015/11/02
 *
 *
 */

#ifndef __GEOMETRY_POWDER_TARGET_H__
#define __GEOMETRY_POWDER_TARGET_H__

#include "geometry.h"


/**
 * \class PowderTarget2D
 * 
 * \brief Function object (level set function) for generating a 2D powder target geometry
 *
 * \author R.S.
 *
 * \version 1.0 
 *
 * \date 2015/08/06 
 *
 * Created on: 2015/08/06
 *
 */
class PowderTarget2D: public Geometry {
public:
	/// constructor
	PowderTarget2D();

	/// destructor
	virtual ~PowderTarget2D() {}
	
	/**
	 * \brief Level set function of a 2D powder target along the x-axis  
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
	double height;
};


#endif // __GEOMETRY_POWDER_TARGET_H__
