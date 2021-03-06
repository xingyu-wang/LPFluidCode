/**
 * \file   lp_solver.h
 *
 * \brief  This header file contains classes of the main Lagrangian Particle solvers such as 
 *         the hyperbolic solvers and the elliptic solvers
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com) 
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design 
 *            and the design of data pointer swaping algorithms in the Strang splitting method               
 *
 *
 * \version 1.0 
 *
 * \date 2014/10/09
 *
 * Created on: 2014/9/20 
 *
 */


#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__

#include <cstddef>
#include <vector>
#include <string>
#include <fstream>

class Initializer;
class ParticleData;
class NeighbourSearcher;
class EOS;


/**
 * \class LPSolver
 * 
 * \brief An abstract class for the family of Lagrangian Particle solvers
 *
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *
 * \version 1.0 
 *
 * \date 2014/10/09 
 *
 * Created on: 2014/09/20 
 *
 */
class LPSolver {

public:
	/// Destructor
	virtual ~LPSolver() {}

	/**
	 * \brief         The black box main Lagrangian particle solver for one iteration step
	 * 
	 * The method should be called by TimeController repeated at every time step
	 *
	 * \param [in] dt The length of physical time for this iteration
	 * \return        0 if the iteration is success 
	 * \warning       The function should always return 0 because all exceptions should be handled inside this class
	 */ 
	virtual int solve(double dt) = 0;
	
	/**
	 * \brief   Getter function of the minimum inter-particle distance among all fluid particles 
	 * \param   None
	 * \return  The minimum inter-particle distance among all fluid particles
	 */		
	virtual double getMinParticleSpacing() const {return m_fMinParticleSpacing;}
	
	/**
	 * \brief   Getter function of the maximum sound speed among all fluid particles 
	 * \param   None
	 * \return  The maximum sound speed among all fluid particles
	 */	
	virtual double getMaxSoundSpeed() const {return m_fMaxSoundSpeed;}

	/**
	 * \brief   Getter function of the maximum absolute value velocity among all fluid particles 
	 * \param   None
	 * \return  The maximum absolute value velocity among all fluid particles
	 */	
	virtual double getMaxFluidVelocity() const {return m_fMaxFluidVelocity;}

protected:

	double m_fMinParticleSpacing; ///< Minimum inter-particle spacing among fluid particles at a time step		
	double m_fMaxSoundSpeed; ///< Maximum sound speed of fluid particles at a time step	
	double m_fMaxFluidVelocity; ///< Maximum absolute value velocity of fluid particles at a time step
	bool m_iIfDebug;///< if true then print debug info
	std::ofstream debug;///< output information for debugging	
};


/**
 * \class HyperbolicLPSolver
 * 
 * \brief The default Lagrangian Particle solver for the compressible Euler's equation in 2D and 3D
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *
 * Co-author: Yu, Kwangmin (yukwangmin@gmail.com) on initial interface design
 *            and the design of data pointer swaping algorithms in the Strang splitting method
 *
 * \version 1.0 
 *
 * \date 2014/10/09 
 *
 * Created on: 2014/09/20 
 *
 */
class HyperbolicLPSolver : public LPSolver {

public:
	/**
	 * \brief       Constructor
	 * 
	 * Get and set up parameters and obtain access to objects needed for the main solver 
	 *
	 * \param [in] init   To retrieve information from \e init   
	 * \param [in] pData  To obtain access to an object of the PaticleData clas
	 * \param [in] ns     To obtain access to an object in the NeighbourSearcher class
	 */
	HyperbolicLPSolver(const Initializer& init, ParticleData* pData, NeighbourSearcher* ns);
	
	/**
	 * \brief         The Lagrangian particle solver for the compressible Euler's equations for one iteration step
	 * 
	 * The method should be called by TimeController repeated at every time step
	 *
	 * \param [in] dt The length of physical time for this iteration
	 * \return        0 if the iteration is success 
	 * \warning       The function should always return 0 because all exceptions should be handled inside this class
	 */
	virtual int solve(double dt);	
	
	virtual ~HyperbolicLPSolver();
private:

	//-----------------------------------------Data----------------------------------------

	//--------------------------Info got from input argument list---------------------------

	ParticleData* m_pParticleData; ///< Pointer to the object containing major particle data arrays 	
	NeighbourSearcher* m_pNeighbourSearcher; ///< Pointer to the object for the neighbour search task
	
	//--------------------------------------------------------------------------------------	

	//--------------------------Info get from Initializer class------------------------------
	
	EOS* m_pEOS; ///< Pointer to the object for computing eos
	double m_pGamma;
	double m_pPinf;
	double m_pEinf;
	int m_iNumThreads; ///< Number of threads	
	bool m_iIfMultiThreads;	///< true if use multithreads
	int m_iDimension; ///< dimension
	bool m_iRandomDirSplitOrder; ///< if true then the order of directional splitting is randomly set
	int m_iLPFOrder; ///< the order of Local Polynomial Fitting (LPF)
	std::size_t m_iNumRow2ndOrder; ///< the smallest number of rows of A to solve 2nd order LPF
	std::size_t m_iNumRow1stOrder; ///< the smallest number of rows of A to solve 1st order LPF
	std::size_t m_iNumCol2ndOrder; ///< the number of columns of A when solving 2nd order LPF
	std::size_t m_iNumCol1stOrder; ///< the number of columns of A when solving 1st order LPF	
	bool m_iMovingBoxForGhostParticle; ///< if true then the fluid box will be updated
	double m_fAvgParticleSpacing; ///< the average particle spacing
        double m_fInitParticleSpacing;//< the initial particle spacing for uniform density 
	double m_fGravity; ///< Gravity 
	double m_fInvalidPressure; ///< if p < invalid pressure => invalid state
	double m_fInvalidDensity; ///< volume cannot be negative: if volume < invalid volume => invalid state
	int    m_iUseCriticalPressure; ///< use critical pressure or not. 1:yes 0:no
	double m_fCriticalPressure; ///< critical pressure for the cavitation model
	int m_iVariableNeiSearchRadius; ///< is use variable neighbour search radius 1:yes 0:no
	double m_fNeiSearchRadius; ///< the radius for neighbour search
	size_t m_iNumParticleWithinSearchRadius; ///< the number of particles within the search radius at initialization time
	double m_fContactLength; ///< defined length such that for two fluid particles from different fluid object, if the distance from each other is shorter than the length the two fluid particles start to interact with each other
	double m_fTimesNeiSearchRadius;///< how many times is the neighbour search radius wrt average inter-particle spacing	
	double m_fTimesContactLength;///< how many times is the neighbour search radius wrt average inter-particle spacing
	double m_fTimesBoundingBox;///< how many times is the buffer zone wrt average inter-particle spacing
	bool m_iIfRestart;///< if a restart run
	//--------------------------------------------------------------------------------------
	
	//---------------------------------Other parameters-------------------------------------

	int m_iNumPhase; ///< number of phases in directional splitting
	/**
	 *\brief 2D: A 2X3 table which maps direction split order and split phase to 0(x) or 1(y)\n
		     3D: A 6X5 table which maps direction split order and split phase to 0(x), 1(y), or 2(z)
	*/
	std::vector<std::vector<int> > m_vDirSplitTable; 

	int m_iDirSplitOrder;///< In 3D: 0=xyzyx, 1=xzyzx, 2=yxzxy, 3=yzxzy, 4=zxyxz, 5=zyx. In 2D: 0=xyx, 1=yxy	
	
	double m_fDt; ///< the time length of this iteration 

	int m_iContactAlert; ///< If true two or more fluid objects may soon in contact with each other; the value is determined by whether the bounding boxes of these fluid objects are overlap or not
	
	//std::ofstream debug;///< output information for debugging
	//bool m_iIfDebug;///< if true then print debug info
	
	bool m_iFreeBoundary; ///< if there is free boundary condition
	bool m_iPeriodicBoundary; ///< if there is periodic boundary condition
	bool m_iSolidBoundary; ///< if there is solid boundary condition
	bool m_iUseLimiter; ///< if use limiter
	bool m_iIfLaxWendroff; ///< if use Lax Wendroff scheme for second order
        bool m_iIfNoSplit; ///< if use non directional splitting Lax Wendroff scheme
	int m_iIfSPH; //< if use SPH density estimator

	double m_fTotalTime;
	double m_fSolverTime;
	double m_fSPHTime;
	double m_fOctreeTime;
	double m_fNeighbourTime;
	double m_fBoundaryTime;
	int m_iCount;

	std::vector<bool> m_vDoNotUpdate;
	std::vector<bool> m_vFillGhost;
	std::vector<bool> m_vGamma;
	std::vector<std::size_t> m_vMirrorIndex; ///< The index of the corresponding fluid particle of a mirror particle
	//-------------------------------------------------------------------------------------




	//-------------------------------------Methods-----------------------------------------
	
	/**
	 * \brief
	 *
	 */
	void checkInvalid();

	/**
	 * \brief
	 *
	 */
	void identifyFluidBoundary();

	/**  
	 * \brief compute average velocity in local neighbourhood of each fluid particle for
	 *        the moving of particles
	 * 
	 */
	void computeLocalWeightedAverageVelocity(size_t index, std::vector<double>& result);

	/**  
	 * \brief set ghost velocities as the corresponding fluid particle 
	 *
	 
	 */
	void setGhostVelocity(int phase);
	
	
	/**  
	 * \brief set the pressure and velocities of a mirror particle
	 *
	 */
	void setMirrorPressureAndVelocity(int phase);


	/**  
	 * \brief A composite function that calls a bunch of other methods in order to set up
	 *        the environment for the next iteration based on the update of this iteration
	 * 
	 */	
	void computeSetupsForNextIteration(); 
	
		
	/**  
	 * \brief Update the bounding boxes of fluid objects based on the updated particle location  
	 *
	 *
	 *
	 */
	void updateFluidBoundingBox();
	
	/**  
	 * \brief Checks if two fluid objects in the simulation is in contact or not 
	 *
	 * The check is based on if the bounding box of the two fluid objects overlap or not 
	 *
	 *
	 */
	void checkForContactAlert();
	
	/**  
	 * \brief Generates ghost particles  
	 *  
	 *
	 */
	bool generateGhostParticle();
	
	/**  
	 * \brief Determines if a ghost candidate particle is added as a ghost particle based on some criteria
	 *
	 *  
	 *
	 */
	bool isValidGhostParticle(double x, double y, double z, const int* neiList, const double* neiListDist, 
	size_t numNei, int objectTag); 
	
	/**  
	 * \brief Determines if a ghost candidate particle is added as a ghost particle based on some criteria
	 *
	 *  
	 *
	 */
	bool isValidGhostParticle(double x, double y, double z, const int* neiList, size_t numNei, int objectTag);


	/*
	 * \brief compute the average number of particles within the search radius 
	 *
	 */
	void computeAvgNumParticleWithinSearchRadius();

	/**  
	 * \brief Searches neighbours for fluid particles based on octree neighobur search  
	 */
	void searchNeighbourForFluidParticle();
        void searchNeighbourForFluidParticle(int choice);
        void SPHDensityEstimatorForFluidParticle(int choice);
	
	/**  
	 * \brief Searches neighbours for ghost particles based on octree neighobur search  
	 */
	void searchNeighbourForGhostParticle();

	/**  
	 * \brief Searches neighbours for all particles including fluid, boundary, and ghost particles 
	 *        based on octree neighobur search 
	 *
	 *  
	 *
	 */
	void searchNeighbourForAllParticle(); 	
	
	/**  
	 * \brief Searches neighbours for all particles including fluid, boundary, and ghost particles 
	 *        based on octree neighobur search, with variable search radius for fluid particles 
	 *
	 *  
	 *
	 */
	void searchNeighbourForAllParticleVariableSearchRadius();

	/**  
	 * \brief Generate ghost particles for free boundary; locally generatye ghost particles for each fluid particle 
	 *         
	 */
	bool generateGhostParticleByFillingVacancy();
	
	#ifdef _OPENMP
	void fillGhostParticle2D(int dir, int count[], size_t index, 
	std::vector<double>& gX, std::vector<double>& gY, std::vector<size_t>& fIndex);	
	
	void fillGhostParticle3D(int dir, int count[], size_t index, 
	std::vector<double>& gX, std::vector<double>& gY, std::vector<double>& gZ, std::vector<size_t>& fIndex);
	#else
	bool fillGhostParticle2D(int dir, int count[], size_t index, size_t* ghostIndex);
	
	bool fillGhostParticle3D(int dir, int count[], size_t index, size_t* ghostIndex);
	#endif
	
	

	/**  
	 * \brief Generate particles for solid boundary by reflecting each fluid particle 
	 *        across the solid boundary 
	 */
	bool generateSolidBoundaryByMirrorParticles();
        bool generatePeriodicBoundaryByMirrorParticles();


	/**  
	 * \brief Set up one-sided neighbour list based on the entire neighbour list (for fluid particles only) 
	 *
	 *  
	 *
	 */
	void setUpwindNeighbourList();
	
	/**  
	* \brief Reset the order of local polynomial fitting to the pre-set value \e m_iLPFOrder  
	*
	* The order of local polynomial fitting are multiple arrays, each represent a direction 
	* (eg. two arrays for right and left in the x-coordinate) 
	*
	*/	
	void resetLPFOrder();
	
	/**
	* \brief cavitation modelling
	*
	*/
	void cavitation();

	/**  
	 * \brief Compute the minimum inter-particle spacing among fluid particles
	 *
	 * For the computation of next dt  
	 *
	 */
	void computeMinParticleSpacing();
	
	/**  
	 * \brief Compute the average inter-particle spacing among fluid particles
	 *
	 * For variable neighbour search radius  
	 *
	 */
	void computeAvgParticleSpacing();
	
	/**  
	 * \brief Update the local inter-particle spacing
	 *
	 * For variable neighbour search radius  
	 *
	 */
	void updateLocalParSpacing();
	
	/**  
	 * \brief Update the local inter-particle spacing
	 *
	 * For variable neighbour search radius  
	 *
	 */
	void updateLocalParSpacingByVolume();


	/**  
	 * \brief Update the neighbour search radius to reflect the change in the average inter-particle spacing
	 *
	 * 
	 *
	 */
	void updateNeighbourSearchRadius();
	
	/**  
	 * \brief Update the contact length to reflect the change in the average inter-particle spacing
	 *
	 * 
	 *
	 */
	void updateContactLength();


	/**  
	 * \brief Computes the maximum sound speed of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 */
	void computeMaxSoundSpeed();
	
	/**  
	 * \brief Computes the maximum absolute value velocity of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 */
	void computeMaxFluidVelocity();
	
	
	/**  
	 * \brief Changes m_iObjectTag of a fluid object if it is in contact with particles in other fluid object 
	 *
	 * If a fluid particle have a fluid neighbour from other fluid object within m_fContactLength
     * then set the m_vObjectTag of this particle to be the negative of original
	 *
	 */
	void changeFluidObjectTagIfInContact(int index, size_t numNeiFound, const double* neiListDist);
	
	/**  
	 * \brief Change the neighbour list of a fluid particle 
	 * 
	 * Change the neighbour list to only include
	 * 1. non-fluid particles and\n 
	 * 2. fluid particles in the same fluid object  
	 *
	 * \note This is done only when this fluid particle is \e not in contact with other fluid particles
	 *
	 */
	void changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(int index, size_t numNeiFound);
	
	/**  
	 * \brief Changes the neighbour list to include only fluid particles 
	 *
	 * This function is now only used by boundary and ghost particles
	 *
	 */
	void changeNeighbourhoodToIncludeOnlyFluidNei(int index, size_t numNeiFound);

	
	/**  
	 * \brief Calculates the spatial derivatives and performs time integration by the Strang splitting method
	 *
	 * Calls the following functions (names only for clarity):\n
	 * setNeighbourListPointers();\n
     * setInAndOutDataPointers();\n
	 * setLPFOrderPointers();\n
     * computeSpatialDer();\n 
	 * timeIntegration();\n 
	 * printInvalidState();\n
	 * lowerLPFOrder();
	 *
	 */
	bool directionalSplitting(int phase);
	
        bool density_derivative();

	bool directionalSplitting_van_leer(int phase);

        bool directionalSplitting_lax_wendroff(int phase);

	bool directionalSplitting_upwind(int phase);

        bool nodirectionalSplitting();

        bool nodirectionalSplitting_noextreme();


	/**  
	 * \brief Alias two one-sided neighbour lists based on if this round is on x, y, or z 
	 *
	 *
	 */	
	void setNeighbourListPointers(int dir, // input
		const int **neighbourList0, const int **neighbourList1, // output
		const int **neighbourListSize0, const int **neighbourListSize1);
	
	/**  
	 * \brief Alias input and output data pointers for this round 
	 *
	 * \note There are 3 rounds in total for 2D and 5 rounds in toal for 3D
	 */
	void setInAndOutDataPointers(int phase, int dir,
		const double** inVelocity, const double** inPressure, const double** inVolume, const double** inSoundSpeed, 
		double** outVelocity, double** outPressure, double** outVolume, double** outSoundSpeed);
		
	/**  
	 * \brief Alias the pointers to the arrays of the order of local polynomial fitting
	 *
	 *
	 */
	void setLPFOrderPointers(int dir, // input
		int** LPFOrder0, int** LPFOrder1, std::vector<int*>& LPFOrderOther); // output

	
	/**  
	 * \brief Alias the pointers to the arrays of the order of local polynomial fitting
	 *
	 *
	 */
	void setLPFOrderPointers(int dir, // input
		int** LPFOrder0, int** LPFOrder1, int** LPFFirstOrder0, int** LPFFirstOrder1); // output


	/**  
	 * \brief Compute the spatial derivatives by solving least squares problem for this round
	 *
	 *
	 */
	void computeSpatialDer(int dir, size_t index, // input
		int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),
		const double* inPressure, const double* inVelocity,
		const int *neighbourList, const int *neighbourListSize,int additional,  
		int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd); // output

	void computeSpatialDer(size_t index, size_t offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),
                                                  const double* inVolume, const int *neighbourList, const int *neighbourListSize,
                                                  int *LPFOrder, double* volume_x, double* volume_y, double* volume_z);

	
        /**  
         * \brief Compute the spatial derivatives by solving least squares problem for this round, and output turncation error in GFD approximation
         *
         *
         */
        void computeSpatialDer(int dir, size_t index, // input
                int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*,double*),
                const double* inPressure, const double* inVelocity,
                const int *neighbourList, const int *neighbourListSize,int additional,
                int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd, double* p_e, double *vel_e); // output


        /**  
         * \brief Compute the spatial derivatives by solving least squares problem for non directional splitting Lax Wendroff method
         *
         *
         */
	void computeSpatialDer(size_t index,  size_t offset, //input
void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),  const double* inPressure, const double* inVelocityU, const double* inVelocityV, const double* inVelocityW, const double* inVolume, const int *neiList, const int *neiListSize,
                                                          int *LPFOrder, double* Pd, double* Ud, double* Vd, double* Wd, double* Volumed, int number_of_derivative);// output	
	/**  
	 * \brief Performs time integration for this round
	 *
	 * \note There are 3 rounds in total for 2D and 5 rounds in toal for 3D
	 */
	void timeIntegration(double real_dt, double multiplier1st, double multiplier2nd, 
		double gravity, double inVolume, double inVelocity, double inPressure,
		double inSoundSpeed, 
		double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
		double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
		double* outVolume, double* outVelocity, double* outPressure); // output

        void timeIntegration(double real_dt, double multiplier1st, double multiplier2nd,
                double gravity, double inVolume, double inVelocity, double inPressure,
                double inSoundSpeed,
                double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
                double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
                double* outVolume, double* outVelocity, double* outPressure, double* volumet, double* velocityt, double* pressuret); // output

	/**  
	 * \brief Performs time integration for this round
	 *
	 * \note There are 3 rounds in total for 2D and 5 rounds in toal for 3D
	 */
	void timeIntegration_van_leer(double realDt, double multiplier1st, double multiplier2nd, 
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double phi,
	double vel_d_0_first, double vel_d_1_first, double p_d_0_first, double p_d_1_first,
	double vel_d_0_second, double vel_d_1_second, double p_d_0_second, double p_d_1_second,
	double vel_dd_0_second, double vel_dd_1_second, double p_dd_0_second, double p_dd_1_second,
	double* outVolume, double* outVelocity, double* outPressure); // output

        /**  
         * \brief Performs time integration for non directional splitting Lax Wendroff method
         *
         * \note There is only 1 round in total
         */

	void timeIntegration(double Dt, double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure);//output
        void timeIntegration(int index, double Dt, double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,//rotate
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure);//output
	
	/**  
	 * \brief Print out info about the particle which evolve into invalid states 
	 *
	 *
	 */
	void printInvalidState(int phase, int dir,  int index, double positionX, double positionY, double positionZ,
		double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0, 
		double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1);
	
	
	
	/**  
	 * \brief Lower the order of local polynomial fitting for a particle in a specific direction (x,y, or z)
	 *
	 * \return \c true if at least one of the LPFOrder array for this particle in a direction is not zero;
	 *         \c false otherwise
	 *
	 * \note   If return \c true then we will redo this particle with the lowered order of local polynomial fitting;
	 *         if return \c false we will go back to phase for all particles, with this particle basically not 
	 *         updated for the following one entire iteration
	 *
	 */
	bool lowerLPFOrder(int index, const std::vector<int*>& LPFOrderOther, // input
		int* LPFOrder0, int* LPFOrder1); // output	
	
	
	/**  
	 * \brief Computes the number of rows and columns used for a particle at this round 
	 *
	 * \note A helper function of computeSpatialDer()
	 *
	 */
	void computeNumRowAndNumColAndLPFOrder(size_t index, // input
		const int *neighbourList, const int *neighbourListSize, size_t numRow2nd, size_t numRow1st,
		int* LPFOrder, size_t *numRow, size_t *numCol); // output


	/**  
	 * \brief Computes the matrix A in the least squares problem Ax~b in the 2D context 
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeA2D(size_t index, const int *neighbourList, const int* LPFOrder, size_t numRow, size_t numCol, // input
					double *A, double *distance); // output 
	
	/**  
	 * \brief Computes the matrix A in the least squares problem Ax~b in the 3D context
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeA3D(size_t index, const int *neighbourList, const int* LPFOrder, size_t numRow, size_t numCol, // input
					double *A, double *distance); // output
	
	/**  
	 * \brief Computes the vector b in the least squares problem Ax~b
	 *
	 * \note A helper function of computeSpatialDer()
	 */
	void computeB(size_t index, const int *neighbourList, size_t numRow, const double* inData, 
				  double *b); // output
	
        /**  
         * \brief Computes the error in the least squares problem Ax~b
         *
         * \note A helper function of computeSpatialDer()
         */
        void compute_GFD_error(size_t numRow, size_t numCol, double *A, double *b, double *x,
				 double *e); // output
	
	/**  
	 * \brief Compute the pressure and velocity by 0-th order local polynomial fitting for boundary particles
	 *
	 * \note calls the function computeOthOrderWeightedLPF()
	 */
	void setBoundaryPressureAndVelocity(int phase);
	
	/**  
	 * \brief Compute the pressure and velocity by 0-th order local polynomial fitting for ghost particles
	 *
	 * \note calls the function computeOthOrderWeightedLPF()
	 */
	void setGhostPressureAndVelocity(int phase); 
	
	/**  
	 * \brief Computes the 0-th order local polynomial fitting
	 *
	 *  
	 */
	void compute0thOrderWeightedLPF(std::vector<const double*>& position, 
									std::size_t startIndex, std::size_t numParticle,
									std::vector<double*>& data);
	
	/**  
	 * \brief Updates the states of fluid particles at the end of one iteration by swapping pointers 
	 *
	 *
	 */
	void updateFluidState();
	
	/**  
	 * \brief Updates the location of fluid particles based on velocities at the end of one iteration
	 *
	 * Based on a combination of forward and backward Euler's method
	 */
	void moveFluidParticle();
	
	/**  
	 * \brief Updates the location of fluid particles based on velocities at the end of one iteration
	 *
	 * Based on a combination of forward and backward Euler's method plus adjustment relative velocities
	 * calculated from the local neighbourhood
	 */
	void moveFluidParticleAdjusted();

	/**  
	 * \brief Update the velocities of fluid particles at the end of one iteration by swapping pointers
	 *
	 * 
	 */
	void updateFluidVelocity();


	/**  
	 * \brief Tests on the neighbour search methods 
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void testNeighbourSearch();
	
	/**  
	 * \brief Use brute force neighbour search to search neighbours for fluid particles
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourBruteForce(int, int* neighbourList, int* neighbourListSize); 

	/**  
	 * \brief Use brute force neighbour search to search neighbours for fluid particles
	 *
	 * \note This is a function for testing the validity of methods in this class
	 *       This overloaded function has the additional neighbourListDist output
	 */
	void searchNeighbourBruteForce(int, int* neighbourList, double* neighoburListDist, int* neighbourListSize);


	/**  
	 * \brief Use brute force neighbour search to search neighbours for ghost particles
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourBruteForce(double x, double y, double z, int* neighbourList, size_t& numNeiFound); 
	
	/**  
	 * \brief Generate ghost particles by the brute force neighbour search
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void generateGhostParticleByBruteForceNeighbourSearch();
	
	/**  
	 * \brief Sum the distance of ghost particles to the origin to test consistency of ghost particle
	 *        search result of the serial and multithreaded version
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void verifyGhostParticleBySumOfDist(std::ofstream& ofs);

	/**  
	 * \brief A wrapper function which uses brute force neighbour search to search neighbours for all particles 
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	void searchNeighbourForAllParticleByBruteForceNeighbourSearch();
	
	/**  
	 * \brief Check the validity of all one-sided neighbour lists
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	bool checkUpwindNeighbourList();
	
	/**  
	 * \brief helper function of writeResult()
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	std::string rightFlush(size_t writeStep, size_t numDigits);	
	
	/**  
	 * \brief Writes results for debugging purposes
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
	int writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle);
};





/**
 * \class HyperbolicLPSolver1D
 * 
 * \brief The default Lagrangian Particle solver for the compressible Euler's equation in 1D
 *
 * \author Chen, Hsin-Chiang (morrischen2008@gmail.com)
 *          
 *
 * \version 1.0 
 *
 * \date 2015/03/13 
 *
 * Created on: 2015/03/13 
 *
 */
class HyperbolicLPSolver1D : public LPSolver {

public:
	/**
	 * \brief Constructor
	 * 
	 * Get and set up parameters and obtain access to objects needed for the main solver 
	 *
	 * \param [in] init   To retrieve information from \e init   
	 * \param [in] pData  To obtain access to an object of the PaticleData clas
	 * 
	 */
	HyperbolicLPSolver1D(Initializer& init, ParticleData* pData);
	
	/**
	 * \brief The 1D Lagrangian particle solver for the compressible Euler's equations for one iteration step
	 * 
	 * The method should be called by TimeController repeated at every time step
	 *
	 * \param [in] dt The length of physical time for this iteration
	 * \return        0 if the iteration is success 
	 * \warning       The function should always return 0 because all exceptions should be handled inside this class
	 */
	virtual int solve(double dt);	

private:

	//-----------------------------------------Data----------------------------------------

	//--------------------------Info got from input argument list---------------------------

	ParticleData* m_pParticleData; ///< Pointer to the object containing major particle data arrays 		
	
	//--------------------------------------------------------------------------------------	

	//--------------------------Info get from Initializer class------------------------------
	
	EOS* m_pEOS; ///< Pointer to the object for computing eos
	int m_iDimension; ///< dimension
	int m_iLPFOrder; ///< the order of Local Polynomial Fitting (LPF)
	std::size_t m_iNumRow2ndOrder; ///< the smallest number of rows of A to solve 2nd order LPF
	std::size_t m_iNumRow1stOrder; ///< the smallest number of rows of A to solve 1st order LPF
	std::size_t m_iNumCol2ndOrder; ///< the number of columns of A when solving 2nd order LPF
	std::size_t m_iNumCol1stOrder; ///< the number of columns of A when solving 1st order LPF		
	double m_fAvgParticleSpacing; ///< the average particle spacing
	double m_fInitParticleSpacing;//< the initial particle spacing for uniform density 
	double m_fInvalidPressure; ///< if p < invalid pressure => invalid state
	double m_fInvalidDensity; ///< volume cannot be negative: if volume < invalid volume => invalid state
	
	int    m_iUseCriticalPressure; ///< use critical pressure or not. 1:yes 0:no
	double m_fCriticalPressure; ///< critical pressure for the cavitation model
	int m_iVariableNeiSearchRadius; ///< is use variable neighbour search radius 1:yes 0:no


	std::string m_sBoundaryType; ///< ONLY 1D 
	bool m_iUseLimiter;///< If use limiter or not; ONLY 1D
	double m_fThresholdP;///< Threshold value on pressure if limiter is used 
	//--------------------------------------------------------------------------------------
	
	//---------------------------------Other parameters-------------------------------------
	
	double m_fPeriodicLeftBoundary; ///< left periodic boundary; ONLY 1D
	double m_fPeriodicRightBoundary; ///< right periodic boundary; ONLY 1D
	double m_fDisBetweenPeriodicBoundaries; ///< distance between left and right periodic boundaries; ONLY 1D
        int m_iIfSPH; //< if use SPH density estimator
	
	double m_fDt; ///< the time length of this iteration
	//bool m_iIfDebug;///< if true then print debug info
	//std::ofstream debug;///< output information for debugging		
	//-------------------------------------------------------------------------------------



	//-------------------------------------Methods-----------------------------------------
	
	/**  
	 * \brief A composite function that calls a bunch of other methods in order to set up
	 *        the environment for the next iteration based on the update of this iteration
	 *
	 * This function calls the following methods (names only for clarity)\n
	 * 
	 * updateFreeBoundaryLocation();\n
	 * updateFreeBoundaryPressureAndVelocity();\n
	 * updateSolidBoundaryPressureAndVelocity();\n
	 * updateLimiter();\n 
	 * computeMinParticleSpacing();\n
	 * computeMaxSoundSpeed();\n
	 * computeMaxFluidVelocity();\n 
	 */	
	void computeSetupsForNextIteration(); 


	/**  
	 * \brief Update the location of free boundary particle (the two with index 0 and total_num_particle - 1) 
	 *        by reconstruction based on desnity 
	 */		
	void updateFreeBoundaryLocation();


	/**  
	 * \brief Assign the pressure and velocity of free boundary particle (index 0 and totalNum-1)
	 *
	 * pressure = 1e-9 (pressure of vacuum)
	 * velocity = velocity of the particle next to it
	 *
	 */		
	void updateFreeBoundaryPressureAndVelocity();


	/**  
	 * \brief Assign the pressure and velocity of solid boundary particle (index 0 and totalNum-1) 
	 *  
	 * pressure = pressure of the particle next to it
	 * velocity = negate the velocity of the particle next to it   
	 *
	 */		
	void updateSolidBoundaryPressureAndVelocity();
	
	
	/**  
	 * \brief Compute the minimum inter-particle spacing among fluid particles
	 *
	 * For the computation of next dt  
	 *
	 * \note Computation includes the distance of fluid particle to boundary particles (solid, free boundary particles)
	 *
	 */
	void computeMinParticleSpacing();


	/**  
	 * \brief Computes the maximum sound speed of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 * \note For periodic boundary all particles are used; For free and solid boundary, the leftmost and rightmost
	 *       particles are not used
	 *
	 */
	void computeMaxSoundSpeed();


	/**  
	 * \brief Computes the maximum absolute value velocity of fluid particles
	 *
	 * For the computation of next dt 
	 *
	 */
	void computeMaxFluidVelocity();

	
	/**  
	 * \brief Update the values of divided difference for limiter
	 * 
	 *
	 */
	void updateLimiter();

	/**
	 * \brief Computes the order of LPF based on 
	 *        m_iLPFOrder (default), limiter (if used), and location (if near boundary)
	 *
	 *
	 */
	void computeLPFOrder(std::size_t index, int& left_order, int& right_order); 

	
	/**  
	 * \brief Compute the spatial derivatives by solving least squares problem 
	 *
	 *
	 */
	void computeSpatialDer(int order, int direction, int i, 
                           const double *local_u, const double *local_x,
                           double& dudx, double& dudxx); //output 
        void computeSpatialDerUpwind(int order, int direction, int i,
                           const double *local_u, const double *local_p, const double *local_x,
                           double& dudx, double& dpdx); //output	
        void computeSpatialDerCenter(int order, int i,
                           const double *local_u, const double *local_p, const double *local_x,
                           double& dudx, double& dudxx, double& dpdx, double& dpdxx); //output
	/**  
	 * \brief Solving least squares problem by QR decomposition 
	 *
	 * \note The matrix A is always recomputed even though px_left and ux_left use same neighbour particles.
	 *       This is fine in 1D but still can be optimised. 
	 * \warning This function will exit the program if either\n 
	 *          1. QR do not have sufficient rank, or\n
	 *          2. QR get nan derivatives
	 */
	void solveByQR(int order, int num_nei, const double *local_u, const double *local_x, 
                  double &dudx, double &dudxx);	//output

	
	/**  
	 * \brief Solving \b weighted least squares problem by using Cramer's rule 
	 *
	 * \note The matrix A is always recomputed even though px_left and ux_left use same neighbour particles.
	 *       This is fine in 1D but still can be optimised.\n
	 *       A const weight function is called and one could use other weight functions if necessary 
	 *       (but needs to write a weight function in this class and make the corresponding function pointer)
	 * \warning This function will exit the program if get nan derivatives
	 */
	void solveByCramer(int order, int num_nei, const double *local_u, const double *local_x, 
					   double &dudx, double &dudxx);	//output
	

	/*
	 * \brief returns a constant weight independent of the input
	 *
	 * \note A helper function for solveByCramer(), which always returns 1 
	 */
	double constantWeight(double h);

	/**  
	 * \brief Performs time integration for this time step
	 *
	 * 
	 */
	void timeIntegration(int i, const double* V_old, const double* up_old,
						 const double* p_old, const double* cs_old,
						 double ux_left,  double uxx_left,
						 double px_left,  double pxx_left,
						 double ux_right, double uxx_right,
						 double px_right, double pxx_right,
						 double* V, double* up, double* p); // output
	
        void timeIntegration(double dt, int i, const double* V_old, const double* up_old,
                                                 const double* p_old, const double* cs_old,
                                                 double ux_left,  double uxx_left,
                                                 double px_left,  double pxx_left,
                                                 double ux_right, double uxx_right,
                                                 double px_right, double pxx_right,
                                                 double* V, double* up, double* p); // output
        void timeIntegrationCenter(int i, const double* V_old, const double* up_old,
                                                 const double* p_old, const double* cs_old, double mass,
                                                 double ux,  double uxx,
                                                 double px,  double pxx,
                                                 double* V, double* up, double* p); // output
	
	/**  
	 * \brief Print out info about the particle which evolve into invalid states 
	 *
	 *
	 */
	void printInvalidState(int i, 
						   double ux_left, double uxx_left, double px_left, double pxx_left, 
						   double ux_right, double uxx_right, double px_right, double pxx_right);
	
			
	/**  
	 * \brief Updates the states of fluid particles at the end of one iteration by swapping pointers 
	 *
	 *
	 */
	void updateFluidState(); 
        void updateFluidState_pc1();
        void updateFluidState_pc2();
	
	
	/**  
	 * \brief Updates the location of fluid particles based on velocities at the end of one iteration
	 *
	 * Based on a combination of forward and backward Euler's method
	 */
	void moveFluidParticle(); 
        void moveFluidParticle_pc1();
        void moveFluidParticle_pc2();
	void reorderFluidParticle();	
	
        int solveDefault(double dt);
	int solve_van_leer(double dt);
        int solveUpwind(double dt);
        int solveUpwindPredictorCorrector(double dt);
        int solveNewUpwind(double dt);

	void SPHDensityEstimatorForFluidParticle();
	
	/**  
	 * \brief Time integration by Van Leer limiter
	 *
	 */
	void timeIntegration_van_leer(
		int i, const double* V_old, const double* up_old, const double* p_old, const double* cs_old, double phi,	
		double ux_left1, double px_left1, double ux_right1, double px_right1, 
		double ux_left2,  double uxx_left, double px_left2,  double pxx_left,
		double ux_right2, double uxx_right, double px_right2, double pxx_right,
		double* V, double* up, double* p);// output
	


	/**
	* \brief cavitation modelling
	*
	*/
	void cavitation();




	/**  
	 * \brief Update the velocities of fluid particles at the end of one iteration by swapping pointers
	 *
	 * 
	 */
	//void updateFluidVelocity(double*& up, double*& up_old);
	

	/**  
	 * \brief Lower the order of local polynomial fitting for a particle in a specific direction (x,y, or z)
	 *
	 * \return \c true if at least one of the LPFOrder array for this particle in a direction is not zero;
	 *         \c false otherwise
	 *
	 * \note   If return \c then we will redo this particle with the lowered order of local polynomial fitting;
	 *         if return \c false we will go back to phase for all particles, with this particle basically not 
	 *         updated for the following one entirew iteration
	 *
	 */
//	bool lowerLPFOrder(int index, const std::vector<int*>& LPFOrderOther, // input
//		int* LPFOrder0, int* LPFOrder1); // output
	
	/**  
	 * \brief helper function of writeResult()
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
//	std::string rightFlush(size_t writeStep, size_t numDigits);	
	
	/**  
	 * \brief Writes results for debugging purposes
	 *
	 * \note This is a function for testing the validity of methods in this class
	 */
//	int writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle);


};





#endif // __LP_SOLVER_H__
