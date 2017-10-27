#ifndef __BOUNDARY_RAYLEIGHTAYLOR_THREED_H__
#define __BOUNDARY_RAYLEIGHTAYLOR_THREED_H__

#include "boundary.h"
#include <vector>


class RayleighTaylor3DBoundary: public Boundary {
public:
        /// constructor
        RayleighTaylor3DBoundary();

        /// destructor
        virtual ~RayleighTaylor3DBoundary() {}

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
        double lengthZ;
        double thickness;

        double rb; ///< right boundary
        double lb; ///< left boundary
        double nb; ///< north boundary
        double sb; ///< south boundary
        double tb; ///< top boundary
        double bb; ///< bottom boundary


        double rbo; ///< outer right boundary
        double lbo; ///< outer left boundary
        double nbo; ///< outer north boundary
        double sbo; ///< outer south boundary
        double tbo; ///< outer top boundary
        double bbo; ///< outer bottom boundary
        double epsilon;

        void reflect(double x, double y, double z, double pressure, double vx, double vy, double vz, int bx, int by, int bz,
        std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<double>& pressureb,
        std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb);
};

#endif
