#ifndef __GEOMETRY_SHOCKTUBE_THREED_H__
#define __GEOMETRY_SHOCKTUBE_THREED_H__

#include "geometry.h"

class Shocktube3D: public Geometry {
public:
        /// constructor
        Shocktube3D();

        /// destructor
        virtual ~Shocktube3D() {}

        /**
         * \brief Level set function of a 3D shocktube along the x-axis  
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
	double lengthZ;
};

class RayleighTaylor3D: public Geometry {
public:
        /// constructor
        RayleighTaylor3D();

        /// destructor
        virtual ~RayleighTaylor3D() {}

        /**
         * \brief Level set function of a 3D shocktube along the x-axis  
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
        double lengthZ;
};

#endif
