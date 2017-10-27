/**
 * \file   geometry_jet.h
 *
 * \brief  This header file contains classes for the 
 *         initialization of the geometry of fluid objects in the 3D jet splash simulation
 *
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * \version 1.0 
 *
 * \date 2015/06/13
 *
 * Created on: 2015/06/13 
 *
 */

#ifndef __GEOMETRY_JET_H__
#define __GEOMETRY_JET_H__

#include "geometry.h"

#include <math.h>



/**
 * \class Jet3D
 * 
 * \brief Function object (level set function) for generating a 3D jet (cylinder) geometry
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/13 
 *
 * Created on: 2015/06/13
 *
 */
class Jet3D: public Geometry {
public:
	/// constructor
	Jet3D();

	/// destructor
	virtual ~Jet3D() {}
	
	/**
	 * \brief Level set function of a 3D jet (cylinder) 
	 *        along the x-axis with small perturbations on the surface 
	 *
	 * The level set function is 
	 * \f$
	 *   (y-\mbox{ycen})^2+(z-\mbox{zcen})^2 \leq r^2   
	 * \f$
	 * \f$
	 *  |x-\mbox{xcen}| \leq \mbox{length}/2   
	 * \f$
	 * , where (xcen,ycen,zcen), length, and \e r 
	 * refer to the center, the length, and the radius of the jet 
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
	double length;
	double xCen;
	double yCen;
	double zCen;
	double magnitude; ///< magnitude of the surface perturbation
};



/**
 * \class Jet2D
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/25 
 *
 * Created on: 2015/06/25
 *
 */
class Jet2D: public Geometry {
public:
	/// constructor
	Jet2D();

	/// destructor
	virtual ~Jet2D() {}
	
	/**
	 * \brief Level set function of a 2D jet
	 *        along the x-axis with small perturbations on the surface 
	 *
	 * The level set function is 
	 * \f$
	 *   (y-\mbox{ycen})^2 \leq r^2   
	 * \f$
	 * \f$
	 *  |x-\mbox{xcen}| \leq \mbox{length}/2   
	 * \f$
	 * , where (xcen,ycen,zcen), length, and \e r 
	 * refer to the center, the length, and the radius of the jet 
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
	double length;
	double xCen;
	double yCen;
	double magnitude; ///< magnitude of the surface perturbation
};


/**
 * \class Jet1D
 * 
 * \brief Function object (level set function) for generating a 1D jet geometry
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/25 
 *
 * Created on: 2015/06/25
 *
 */
class Jet1D: public Geometry {
public:
	/// constructor
	Jet1D();

	/// destructor
	virtual ~Jet1D() {}
	
	/**
	 * \brief Level set function of a 1D jet 
	 *
	 * The level set function is 
	 * \f$
	 *  |x-\mbox{xcen}| \leq \mbox{length}/2   
	 * \f$
	 * , where xcen and length refer to the center and  length of the jet 
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
	double length;
	double xCen;
};

class Jet1DLeft: public Geometry {
public:
        Jet1DLeft();

        virtual ~Jet1DLeft() {}

        virtual bool operator()(double x, double y, double z) const;

        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double length;
        double xCen;
};

class Jet1DRight: public Geometry {
public:
        Jet1DRight();

        virtual ~Jet1DRight() {}

        virtual bool operator()(double x, double y, double z) const;

        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double length;
        double xCen;
};

class Jet1DCenter: public Geometry {
public:
        Jet1DCenter();

        virtual ~Jet1DCenter() {}

        virtual bool operator()(double x, double y, double z) const;

        virtual void getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
private:
        double length;
        double xCen;
};

/**
 * \class Jet2DExpansion
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/25 
 *
 * Created on: 2015/06/25
 *
 */
class Jet2DExpansion: public Geometry {
public:
	/// constructor
	Jet2DExpansion();

	/// destructor
	virtual ~Jet2DExpansion() {}
	
	/**
	 * \brief Level set function of a 2D jet along the x-axis  
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
	double length;
	double xCen;
	double yCen;
};






/**
 * \class Jet3DExpansion
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * \version 1.0 
 *
 * \date 2015/06/25 
 *
 * Created on: 2015/06/25
 *
 */
class Jet3DExpansion: public Geometry {
public:
	/// constructor
	Jet3DExpansion();

	/// destructor
	virtual ~Jet3DExpansion() {}
	
	/**
	 * \brief Level set function of a 3D jet along the x-axis  
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
	double length;
	double xCen;
	double yCen;
	double zCen;
};



/**
 * \class Jet2DMergeUpper
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry for the upper jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/04 
 *
 * Created on: 2015/11/04
 *
 */
class Jet2DMergeUpper: public Geometry {
public:
	/// constructor
	Jet2DMergeUpper();

	/// destructor
	virtual ~Jet2DMergeUpper() {}
	
	/**
	 * \brief Level set function of a 2D jet along the x-axis  
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
	double length;
	double xCen;
	double yCen;
  double alpha;
};

/**
 * \class Jet2DMergeLower
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry for the lower jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/04 
 *
 * Created on: 2015/11/04
 *
 */
class Jet2DMergeLower: public Geometry {
public:
	/// constructor
	Jet2DMergeLower();

	/// destructor
	virtual ~Jet2DMergeLower() {}
	
	/**
	 * \brief Level set function of a 2D jet along the x-axis  
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
	double length;
	double xCen;
	double yCen;
  double alpha;
};

/**
 * \class Jet2DCollision
 * 
 * \brief Function object (level set function) for generating two 2D jet geometry
 *
 * \author Xingyu Wang(xingyu.wang@stonybrook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/09 
 *
 * Created on: 2015/11/09
 *
 */
class Jet2DCollision: public Geometry {
public:
	/// constructor
	Jet2DCollision();

	/// destructor
	virtual ~Jet2DCollision() {}
	
	/**
	 * \brief Level set function of a 2D jet along the x-axis  
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
	double xLength;
	double yLength;
	double xLeftCen;
	double xRightCen;
	double yCen;
};


/**
 * \class Jet2DMerge
 * 
 * \brief Function object (level set function) for generating a 2D jet geometry for both the upper and lower jet in jet merge simulation
 *
 * \author Xingyu Wang (xingyu.wang@stonyrbook.edu)
 *
 * \version 1.0 
 *
 * \date 2015/11/10 
 *
 * Created on: 2015/11/10
 *
 */
class Jet2DMerge: public Geometry {
public:
	/// constructor
	Jet2DMerge();

	/// destructor
	virtual ~Jet2DMerge() {}
	
	/**
	 * \brief Level set function of a 2D jet along the x-axis  
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
	double length;
	double xCen;
	double yCen;
  double alpha;
};

#endif // __GEOMETRY_JET_H__
