/**
 * \file   geometry_gresho.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the gresho simulation
 *
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

#ifndef __GEOMETRY_GRESHO_H__
#define __GEOMETRY_GRESHO_H__

#include "geometry.h"


/**
 * \class Gresho2D
 * 
 * \brief Function object (level set function) for generating a 2D gresho geometry
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
class Gresho2D: public Geometry {
public:
	/// constructor
	Gresho2D();

	/// destructor
	virtual ~Gresho2D() {}
	
	/**
	 * \brief Level set function of a 2D gresho along the x-axis  
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
};

class Yee2D: public Geometry {
public:
        /// constructor
        Yee2D();

        /// destructor
        virtual ~Yee2D() {}

        /**
         * \brief Level set function of a 2D gresho along the x-axis  
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
};

class Sedov2D: public Geometry {
public:
        /// constructor
        Sedov2D();

        /// destructor
        virtual ~Sedov2D() {}

        /**
         * \brief Level set function of a 2D gresho along the x-axis  
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
};

#endif // __GEOMETRY_GRESHO_H__
