/**
 * \file   boundary_solid_gresho.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the solid boundary in the gresho simulation
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

#ifndef __BOUNDARY_SOLID_GRESHO_H__
#define __BOUNDARY_SOLID_GRESHO_H__

#include "boundary.h"
#include <vector>

/**
 * \class Gresho2DSolidBoundary
 * 
 * \brief Function object (level set function) for generating a 2D gresho solid boundary
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
class Gresho2DSolidBoundary: public Boundary {
public:
	/// constructor
	Gresho2DSolidBoundary();

	/// destructor
	virtual ~Gresho2DSolidBoundary() {}
	
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
	double radius;	
	double thickness;
	double bo; 

	double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d
	
};

class Yee2DSolidBoundary: public Boundary {
public:
        /// constructor
        Yee2DSolidBoundary();

        /// destructor
        virtual ~Yee2DSolidBoundary() {}

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
        double radius;
        double thickness;
        double bo;

        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};


class Sedov2DSolidBoundary: public Boundary {
public:
        /// constructor
        Sedov2DSolidBoundary();

        /// destructor
        virtual ~Sedov2DSolidBoundary() {}

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
        double radius;
        double thickness;
        double bo;

        double ip2d(double x1, double y1, double x2, double y2) {return x1*x2+y1*y2;}///< inner product 2d

};

#endif // __BOUNDARY_SOLID_GRESHO_H__
