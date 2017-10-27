/**
 * \file   boundary.h
 *
 * \brief  This header file contains classes for the initialization of boundary objects
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

#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include <unordered_map>
#include <string>
#include <vector>
#include "eos.h"
#include "particle_data.h"
/**
 * \class Boundary
 * 
 * \brief An abstract class for the initialization of the boundary of fluid objects
 *
 * This class can be used as a function object: the operator() is overloaded
 * as how to get a boundary particle based on a fluid particle 
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
class Boundary {
public:
	/// Destructor
	virtual ~Boundary() {}
	
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
	std::vector<double>& pressureb, std::vector<double>& vxb, std::vector<double>& vyb, std::vector<double>& vzb)=0; 	
	virtual int UpdateInflowBoundary(ParticleData* pData, EOS* m_pEOS, double dt, double m_fInitParticleSpacing){return 0;};
};


/**
 * \class BoundaryFactory
 * 
 * \brief The abstract factory class for creating objects in the Boundary family 
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
 *
 */
class BoundaryFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the Boundary family 
	 */
	typedef Boundary* (*GeoCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class BoundaryFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class BoundaryFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          BoundaryFactory& factory = BoundaryFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one BoundaryFactory instance is allowed in each program
	 */
	static BoundaryFactory& instance(); 
	
	/**
     * \brief      Registers (links) the boundary name \e name 
	 *		       with the function \e func for creating objects in the Boundary family
	 *			   
	 *             After registration, \e name can be used as an argument in the createBoundary member function
	 *             for creating objects of the linked type in the Boundary family
	 *	           
	 *  
	 * \param [in] name the boundary name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the Boundary family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the BoundaryRegistrar class for
	 *		       the purpose of linking a boundary name and a specific class in the Boundary family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerBoundary(std::string name, GeoCreateFunc func);
	
	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the Boundary family via the 
	 *				    registrerBoundary member function
	 *  
	 * \return     A Boundary * pointer pointing to an object of a specific class in the Boundary family 
	 * 
	 * Example usage: 
	 * \code
	 *            BoundaryFactory& factory = BoundaryFactory::instance();
	 *            Boundary* newBoundary = factory.createBoundary(name);
	 * \endcode
	 *
	 */
	Boundary* createBoundary(std::string name);
private:
	std::unordered_map<std::string,GeoCreateFunc> bTable; ///< hash table for the (name,creatFunction) pair
	BoundaryFactory() {} ///< for singleton design pattern
	BoundaryFactory(const BoundaryFactory& other); ///< for singleton design pattern (Don't implement)
	BoundaryFactory& operator=(const BoundaryFactory& other); ///< for singleton design pattern (Don't implement)
};


#endif // __BOUNDARY_H__
