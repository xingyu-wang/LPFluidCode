/**
 * \file   geometry_shocktube.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the shocktube simulation
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

#ifndef __GEOMETRY_SHOCKTUBE_H__
#define __GEOMETRY_SHOCKTUBE_H__

#include "geometry.h"


/**
 * \class Shocktube2D
 * 
 * \brief Function object (level set function) for generating a 2D shocktube geometry
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
class Shocktube2D: public Geometry {
public:
	/// constructor
	Shocktube2D();

	/// destructor
	virtual ~Shocktube2D() {}
	
	/**
	 * \brief Level set function of a 2D shocktube along the x-axis  
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
	double lengthX;
	double lengthY;
};


class TPShocktube2D: public Geometry {
public:
	/// constructor
	TPShocktube2D();

	/// destructor
	virtual ~TPShocktube2D() {}
	
	/**
	 * \brief Level set function of a 2D shocktube along the x-axis  
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
	double lengthX;
	double lengthY;
};


class BigShocktube2D: public Geometry {
public:
        /// constructor
        BigShocktube2D();

        /// destructor
        virtual ~BigShocktube2D() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};


class Shocktube2DLeft: public Geometry {
public:
        /// constructor
        Shocktube2DLeft();

        /// destructor
        virtual ~Shocktube2DLeft() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};

class Shocktube2DRight: public Geometry {
public:
        /// constructor
        Shocktube2DRight();

        /// destructor
        virtual ~Shocktube2DRight() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};

class RayleighTaylor2D: public Geometry {
public:
        /// constructor
        RayleighTaylor2D();

        /// destructor
        virtual ~RayleighTaylor2D() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};

class DamBreak2D: public Geometry {
public:
        /// constructor
        DamBreak2D();

        /// destructor
        virtual ~DamBreak2D() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};

class BoundaryTest2D: public Geometry {
public:
        /// constructor
        BoundaryTest2D();

        /// destructor
        virtual ~BoundaryTest2D() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};


class KelvinHelmholtz2D: public Geometry {
public:
        /// constructor
        KelvinHelmholtz2D();

        /// destructor
        virtual ~KelvinHelmholtz2D() {}

        /**
         * \brief Level set function of a 2D shocktube along the x-axis  
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
        double lengthX;
        double lengthY;
};

#endif // __GEOMETRY_SHOCKTUBE_H__
