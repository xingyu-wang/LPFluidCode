#include "lp_solver.h"
#include "boundary.h"
#include "neighbour_searcher.h"
#include "eos.h"
#include "particle_data.h"
#include "initializer.h"
#include "ls_solver.h"
#include "hexagonal_packing.h"
#include "omp.h"
//#include "voronoi_area_estimator.h"
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <memory> // shared_ptr
#include <iomanip> // setw
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of HyperbolicLPSolver
////////////////////////////////////////////////////////////////////////////////////////


HyperbolicLPSolver::HyperbolicLPSolver(const Initializer& init, ParticleData* pData, NeighbourSearcher* ns) {
	
	srand(time(0));

	m_pParticleData = pData; 
	m_pNeighbourSearcher = ns;
	m_pEOS = init.getEOS();
	std::vector<double> eos_parameters;
	m_pEOS->getParameters(eos_parameters);

	if(eos_parameters.size()==1)
	{
		m_pGamma=eos_parameters[0];
		m_pPinf=0.;
		m_pEinf=0.;
		cout<<"Polytropic EOS, gamma = "<<m_pGamma<<endl;	
	}
	else if(eos_parameters.size()==3)
        {
                m_pGamma=eos_parameters[0];
                m_pPinf=eos_parameters[1];
                m_pEinf=eos_parameters[2];
                cout<<"Stiff polytropic EOS, gamma = "<<m_pGamma<<", P_inf = "<<m_pPinf<<", E_inf = "<<m_pEinf<<endl;
        }
	else
	{
		m_pGamma=0.;
		m_pPinf=0.;
		m_pEinf=0.;
		cout<<"Warning: Cannot recognize EOS."<<endl;
	}
	
	// get parameters from init
	m_iNumThreads = init.getNumThreads();
	m_iDimension = init.getDimension();
	m_iNumPhase = m_iDimension==3? 5:3;
	m_iRandomDirSplitOrder = init.getRandomDirSplitOrder();
	m_iLPFOrder = init.getLPFOrder(); 
	m_iNumRow1stOrder = init.getNumRow1stOrder();
	m_iNumRow2ndOrder = init.getNumRow2ndOrder();
	m_iNumCol1stOrder = init.getNumCol1stOrder(); 	
	m_iNumCol2ndOrder = init.getNumCol2ndOrder();
	m_iMovingBoxForGhostParticle = init.getMovingBoxForGhostParticle(); 
	m_fAvgParticleSpacing = init.getInitParticleSpacing();
        m_fInitParticleSpacing = m_fAvgParticleSpacing;
	m_fGravity = init.getGravity(); 
	m_fInvalidPressure = init.getInvalidPressure(); 
	m_fInvalidDensity = init.getInvalidDensity();
	m_iUseCriticalPressure = init.getUseCriticalPressure();
	m_fCriticalPressure = init.getCriticalPressure();
	m_iVariableNeiSearchRadius = init.getVariableNeiSearchRadius(); 	
	//cout<<"m_iUseCriticalPressure="<<m_iUseCriticalPressure<<endl;
	//cout<<"m_fCriticalPressure="<<m_fCriticalPressure<<endl;
	m_fNeiSearchRadius = init.getNeiSearchRadius(); 
	m_iNumParticleWithinSearchRadius = init.getNumParticleWithinSearchRadius(); 
	m_fContactLength = init.getContactLength();
	m_fTimesNeiSearchRadius = init.getTimesNeiSearchRadius();
	m_fTimesContactLength = init.getTimesContactLength();
	m_fTimesBoundingBox = init.getTimesBoundingBox();
	m_iIfRestart = init.getIfRestart();
	m_iUseLimiter = init.getUseLimiter();
	m_iIfLaxWendroff = init.getIfLaxWendroff();
	m_iIfNoSplit=init.getIfNoSplit();
	m_iIfSPH=init.getIfSPH();
	if(m_iIfLaxWendroff)
		printf("Lax Wendroff\n");
	else 
		printf("Beam Warming\n");
	if(m_iIfNoSplit)
	{
		printf("No directional splitting\n");
		if((m_iLPFOrder==1)||(m_iUseLimiter==1)||(m_iIfLaxWendroff==0))
			printf("Warning: non directional splitting only support Lax Wendroff scheme without limiter.\n");
	}
	if(m_iIfSPH)
		printf("SPH density estimator\n");
	else
		printf("PDE density updator\n");

	//m_iIfDebug = init.getIfDebug();
	//debug.open(init.getDebugfileName(), std::ofstream::out | std::ofstream::app);	
	
	// all fluid objects should be distant at initialization
	// this variable should always be false in the single fluid object case
	//m_iContactAlert = false; 

	// set OpenMP environment
	m_iIfMultiThreads = false;
	if(m_iNumThreads > 0) {//> 1
		// set the number of threads
		omp_set_num_threads(min(omp_get_max_threads(), m_iNumThreads));	
		m_iIfMultiThreads = true;	
		
		cout<<"-------HyperbolicLPSolver::HyperbolicLPSolver()-------"<<endl;
		cout<<"m_iNumThreads = "<<m_iNumThreads<<endl;
		cout<<"omp_get_num_procs() = "<<omp_get_num_procs()<<endl;
		cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<" (after omp_set_num_threads())"<<endl;
		cout<<"------------------------------------------------------"<<endl;
	}

	if(m_iDimension==2) {
		m_vDirSplitTable = vector<vector<int> >({{0,1,0},{1,0,1}});
//                m_vDirSplitTable = vector<vector<int> >({{1,0,1},{0,1,0}});

//		cout<<m_vDirSplitTable[0][0]<<" "<<m_vDirSplitTable[0][1]<<" "<<m_vDirSplitTable[0][2]<<endl;
//		cout<<m_vDirSplitTable[1][0]<<" "<<m_vDirSplitTable[1][1]<<" "<<m_vDirSplitTable[1][2]<<endl;
//		assert(false);
	}
	else if(m_iDimension==3)
		m_vDirSplitTable = vector<vector<int> >
		({{0,1,2,1,0},
		  {0,2,1,2,0},
		  {1,0,2,0,1},
		  {1,2,0,2,1},
		  {2,0,1,0,2},
		  {2,1,0,1,2}});
	
	// for completeness initialize to zero
	m_iDirSplitOrder = 0;
	m_fDt = 0;
		
	//m_vDoNotUpdate = vector<bool>(m_pParticleData->getFluidNum(),false); 
	
	m_iFreeBoundary = false;
	m_iPeriodicBoundary = false;
	m_iSolidBoundary = false;
	for(auto s:m_pParticleData->m_vBoundaryObjTypes) {
		if(s=="free") {
			if(m_iFreeBoundary) continue; // avoid re-initialize memory
			m_iFreeBoundary = true;
			m_vFillGhost = vector<bool>(m_pParticleData->m_iCapacity,false);
		}
		else if(s=="periodic") {
			m_iPeriodicBoundary = true;
		}
		else if(s=="solid" || s=="inflow" || s=="outflow") {
			m_iSolidBoundary = true;
		}
	}	

	m_fTotalTime=0;
	m_fSolverTime=0;
        m_fSPHTime=0;
        m_fOctreeTime=0;
        m_fNeighbourTime=0;
        m_fBoundaryTime=0;

/*
	if(m_iIfRestart) {
		assert(false);
		if(m_iVariableNeiSearchRadius==1 || m_iVariableNeiSearchRadius==2) {// global and local aveage particle spacing
			searchNeighbourForAllParticle();
			computeAvgParticleSpacing();
			updateNeighbourSearchRadius();
		}
		if(m_iVariableNeiSearchRadius==2) // local average particle spacing
			updateLocalParSpacing();
	}	
*/
	
//	searchNeighbourForFluidParticle();
//	identifyFluidBoundary();
        searchNeighbourForFluidParticle(0);

//threepoint
//	m_vGamma= vector<bool>(m_pParticleData->m_iCapacity,false);
	
//        for(size_t index=m_pParticleData->m_iFluidStartIndex;
//        index<m_pParticleData->m_iFluidStartIndex+m_pParticleData->m_iFluidNum; index++) {
//		if(m_pParticleData->m_vPressure[index]*m_pParticleData->m_vVolume[index]>0.5)
//	                m_vGamma[index] = true;
//        }


//threepoint end
	computeSetupsForNextIteration();
		
        m_fTotalTime=0;
        m_fSolverTime=0;
        m_fSPHTime=0;
        m_fOctreeTime=0;
        m_fNeighbourTime=0;
        m_fBoundaryTime=0;
	m_iCount=0;
}

HyperbolicLPSolver::~HyperbolicLPSolver() {
	delete m_pEOS;
}


void HyperbolicLPSolver::computeSetupsForNextIteration() {
		
	//testNeighbourSearch();
	//assert(false);
/*	 
	if(m_iMovingBoxForGhostParticle) 
		updateFluidBoundingBox();// update the fluid bounding box
*/
	// will check for contact alert when having multiple fluid objects and when it is false
	// when m_iContactAlert becomes true it remains true for good
	//if(m_pParticleData->m_vFluidBoundingBox.size() > 1 && !m_iContactAlert)
	//	checkForContactAlert();
/*		
	while(!generateGhostParticle()) {// generate ghost particles
		cout<<"Not enough space for ghost particles!!!"<<endl;
		cout<<"m_pParticleData->m_iCapacity changed from "<<m_pParticleData->m_iCapacity;
		size_t newCapacity = 1.5*m_pParticleData->m_iCapacity;
		m_pParticleData->augmentAllDataArrays(newCapacity);
		m_pNeighbourSearcher->setMaxParticleNum(newCapacity);
		cout<<" to "<<m_pParticleData->m_iCapacity<<endl;
	}
*/
	//ofstream ofs1("comp_ghost_parallel", ios::out | ios::app);
	//verifyGhostParticleBySumOfDist(ofs1);
	//ofstream ofs2("comp_ghost_serial", ios::out | ios::app);
	//generateGhostParticleByBruteForceNeighbourSearch();
	//verifyGhostParticleBySumOfDist(ofs2);
	//writeResult(0,0,0,m_pParticleData->getTotalNum());

	double startTime;

	startTime = omp_get_wtime();
	if(m_iSolidBoundary) generateSolidBoundaryByMirrorParticles();
	if(m_iPeriodicBoundary) generatePeriodicBoundaryByMirrorParticles();
	if(m_iSolidBoundary || m_iPeriodicBoundary) 
	{
                m_fBoundaryTime+=omp_get_wtime() - startTime;
//		printf("Create boundary particles takes %.16g seconds\n", omp_get_wtime() - startTime);
	}
//Octree: build octree and use it to search neighbours
	startTime = omp_get_wtime();
	searchNeighbourForFluidParticle(0);
        m_fOctreeTime+=omp_get_wtime() - startTime;
//	printf("Search neighbours for all fluid particles takes %.16g seconds\n", omp_get_wtime() - startTime);

        updateLocalParSpacingByVolume();

        checkInvalid();

	startTime = omp_get_wtime();
	if(m_iFreeBoundary) generateGhostParticleByFillingVacancy();
	if(m_iFreeBoundary)
	{ 
                m_fBoundaryTime+=omp_get_wtime() - startTime;
//		printf("Fill ghost particles takes %.16g seconds\n", omp_get_wtime() - startTime);
	}
//Stencil: set upwind and central GFD stencils
	startTime = omp_get_wtime();
setUpwindNeighbourList();
m_fNeighbourTime+=omp_get_wtime() - startTime;
//	printf("Set upwind neighbour list takes %.16g seconds\n", omp_get_wtime() - startTime);

startTime = omp_get_wtime();
if(m_iIfSPH) density_derivative();
SPHDensityEstimatorForFluidParticle(m_iIfSPH);
if(m_iIfSPH) 
{
	m_fSPHTime+=omp_get_wtime() - startTime;
//		printf("SPH density estimator for all fluid particles takes %.16g seconds\n", omp_get_wtime() - startTime);
}
// initialize the LPF order (1 or 2) in directions right, left, north, and south
resetLPFOrder();


// to determine the dt for next step
computeMinParticleSpacing();
computeMaxSoundSpeed();
computeMaxFluidVelocity();
/*	
if(m_iVariableNeiSearchRadius==1 || m_iVariableNeiSearchRadius==2) {
	computeAvgParticleSpacing();
	updateNeighbourSearchRadius();
}
if(m_iVariableNeiSearchRadius==2)
	updateLocalParSpacing();
*/

//updateContactLength();	


}

int HyperbolicLPSolver::solve(double dt) {	
//cout<<"--------------HyperbolicLPSolver::solve()--------------"<<endl;
//	cout<<"---------------------------------------------------------"<<endl;
double currentstepstartTime;
currentstepstartTime = omp_get_wtime();
// dt for this time step 
m_fDt = dt;
if(m_iIfNoSplit==0)//the previous LP algorithm
{
double startTime;
startTime = omp_get_wtime();
if(m_iRandomDirSplitOrder) {// get a random order of split
	if(m_iDimension==2) m_iDirSplitOrder = rand()%2; // 0, 1
	else if(m_iDimension==3) m_iDirSplitOrder = rand()%6; // 0, 1, 2, 3, 4, 5
}
cout<<"m_iDirSplitOrder="<<m_iDirSplitOrder<<endl;  

//fill_n(m_vDoNotUpdate.begin(),m_vDoNotUpdate.size(),false);//all fluid particles should be updated	
//size_t numParNotUpdate = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
for(size_t index=m_pParticleData->m_iFluidStartIndex; 
index<m_pParticleData->m_iFluidStartIndex+m_pParticleData->m_iFluidNum; index++) {
	m_pParticleData->m_vVolumeOld[index] = m_pParticleData->m_vVolume[index];
}

//	double startTime = omp_get_wtime();
for(int phase=0; phase<m_iNumPhase; ) {
	
	cout<<"phase="<<phase<<endl;

	bool phase_success; 
	
	if(m_iUseLimiter)
	{
		if(m_iIfLaxWendroff)
			phase_success = directionalSplitting_lax_wendroff(phase);
		else
			phase_success = directionalSplitting_van_leer(phase);	
	}
	else
		phase_success = directionalSplitting(phase);
		
	if(!phase_success) {
		phase=0; 
		//numParNotUpdate++;
		cout<<"GO BACK TO PHASE 0!!!!!!!"<<endl;
		continue;
	} 
	
	if(m_iSolidBoundary) setMirrorPressureAndVelocity(phase);
	if(m_iPeriodicBoundary) setMirrorPressureAndVelocity(phase);
	if(m_iFreeBoundary) setGhostVelocity(phase);
	
	//setGhostPressureAndVelocity(phase);
	//setBoundaryPressureAndVelocity(phase);

	phase++;
}
m_fSolverTime+=omp_get_wtime() - startTime;
//cout<<"numParNotUpdate="<<numParNotUpdate<<endl;
//	printf("Strang splitting takes %.16g seconds\n", omp_get_wtime() - startTime);
}
else//The current non directional splitting Lax Wendroff scheme
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
for(size_t index=m_pParticleData->m_iFluidStartIndex;
index<m_pParticleData->m_iFluidStartIndex+m_pParticleData->m_iFluidNum; index++) {
	m_pParticleData->m_vVolumeOld[index] = m_pParticleData->m_vVolume[index];
}

double startTime;
startTime = omp_get_wtime();
//Solver: upwind scheme
for(int phase=0; phase<m_iNumPhase; ) {

	cout<<"upwind phase="<<phase<<endl;

	bool phase_success;

	phase_success = directionalSplitting_upwind(phase);
	
	if(!phase_success) {
		phase=0;
		//numParNotUpdate++;
		cout<<"GO BACK TO PHASE 0!!!!!!!"<<endl;
		continue;
	}

	if(m_iSolidBoundary) setMirrorPressureAndVelocity(phase);
	if(m_iPeriodicBoundary) setMirrorPressureAndVelocity(phase);
	if(m_iFreeBoundary) setGhostVelocity(phase);

	//setGhostPressureAndVelocity(phase);
	//setBoundaryPressureAndVelocity(phase);

	phase++;
}

//        if(m_iSolidBoundary) generateSolidBoundaryByMirrorParticles();
//        if(m_iFreeBoundary) generateGhostParticleByFillingVacancy();

//Solver: Lax-Wendroff scheme
bool phase_success;
phase_success = nodirectionalSplitting_noextreme();
m_fSolverTime+=omp_get_wtime() - startTime;
//        printf("No directioinal Splitting takes %.16g seconds\n", omp_get_wtime() - startTime);
if(!phase_success) {
	cout<<"Error in nodirectionalSplitting"<<endl;
	return 1;
}
}
updateFluidState();
moveFluidParticle();
//moveFluidParticleAdjusted();
updateFluidVelocity();	

computeSetupsForNextIteration();
m_fTotalTime+=omp_get_wtime() - currentstepstartTime;
m_iCount++;
cout<<"****************************************************************************"<<endl;
cout<< setw(60) <<"Number of time steps = "<<m_iCount<<endl;
cout<< setw(60) <<"Total running time = "<<m_fTotalTime<<endl;
cout<< setw(60) <<"Time to solve GFD derivative and update states = "<<m_fSolverTime<<endl;
cout<< setw(60) <<"Time to estimate SPH density = "<<m_fSPHTime<<endl;
cout<< setw(60) <<"Time to bluid octree and search neighbours = "<<m_fOctreeTime<<endl;
cout<< setw(60) <<"Time to bluid GFD stencil = "<<m_fNeighbourTime<<endl;	
cout<< setw(60) <<"Time to generate boundary particles = "<<m_fBoundaryTime<<endl;
cout<<"****************************************************************************"<<endl;
//cout<<"-------------------------------------------------------"<<endl;

return 0;
}


void HyperbolicLPSolver::checkInvalid() {
size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
#ifdef _OPENMP
#pragma omp parallel for
#endif
for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
	if(std::isnan(m_pParticleData->m_vPositionX[index]) || std::isinf(m_pParticleData->m_vPositionX[index])) {
		cout<<"invalid x"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vPositionY[index]) || std::isinf(m_pParticleData->m_vPositionY[index])) {
		cout<<"invalid y"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vVelocityU[index]) || std::isinf(m_pParticleData->m_vVelocityU[index])) {
		cout<<"invalid u"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vVelocityV[index]) || std::isinf(m_pParticleData->m_vVelocityV[index])) {
		cout<<"invalid v"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vPressure[index]) || std::isinf(m_pParticleData->m_vPressure[index])) {
		cout<<"invalid p"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vVolume[index]) || std::isinf(m_pParticleData->m_vVolume[index])) {
		cout<<"invalid vol"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vSoundSpeed[index]) || std::isinf(m_pParticleData->m_vSoundSpeed[index])) {
		cout<<"invalid cs"<<endl;
		assert(false);
	}
	if(std::isnan(m_pParticleData->m_vLocalParSpacing[index]) || std::isinf(m_pParticleData->m_vLocalParSpacing[index]            || m_pParticleData->m_vLocalParSpacing[index] <=0)) {
		cout<<"invalid spacing, spacing="<<m_pParticleData->m_vLocalParSpacing[index]<<endl;
		assert(false);
	}
}
}

void HyperbolicLPSolver::computeLocalWeightedAverageVelocity(size_t index, vector<double>& result) {

//cout<<"-------HyperbolicLPSolver::computeLocalWeightedAverageVelocity()-------"<<endl;	
	
size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;

if(m_iDimension==3) {
	double w2_total = 0, du_w_sum = 0, dv_w_sum = 0, dw_w_sum = 0;
	for(int i=0; i<m_pParticleData->m_vNeighbourListSize[index]; i++) {
				
		size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
		
		if(neiIndex>=fluidStartIndex && neiIndex<fluidEndIndex) {
			double dx = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double dy = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];
			double dz = m_pParticleData->m_vPositionZ[neiIndex] - m_pParticleData->m_vPositionZ[index];
			double dis2 = dx*dx+dy*dy+dz*dz;
			double weight = exp(-dis2);
			double w2 = weight*weight;
			w2_total += w2;
			
			double du = m_pParticleData->m_vTemp1VelocityU[neiIndex]-m_pParticleData->m_vTemp1VelocityU[index];
			double dv = m_pParticleData->m_vTemp1VelocityV[neiIndex]-m_pParticleData->m_vTemp1VelocityV[index];
			double dw = m_pParticleData->m_vTemp1VelocityW[neiIndex]-m_pParticleData->m_vTemp1VelocityW[index];		
			
			du_w_sum += du*w2;
			dv_w_sum += dv*w2;
			dw_w_sum += dw*w2;
		}		
	
	}
	if(w2_total != 0) {
		result[0] = du_w_sum/w2_total;
		result[1] = dv_w_sum/w2_total;
		result[2] = dw_w_sum/w2_total;
		
	}
	else {
		result[0] = 0;
		result[1] = 0;
		result[2] = 0;
	}	
}
else if(m_iDimension==2) {
	double w2_total = 0, du_w_sum = 0, dv_w_sum = 0;
	for(int i=0; i<m_pParticleData->m_vNeighbourListSize[index]; i++) {
					
		size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
		
		if(neiIndex>=fluidStartIndex && neiIndex<fluidEndIndex) {
			double dx = m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index];
			double dy = m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index];
			
			double dis2 = dx*dx+dy*dy;
			double weight = exp(-dis2);
			double w2 = weight*weight;
			w2_total += w2;
			
			double du = m_pParticleData->m_vTemp1VelocityU[neiIndex]-m_pParticleData->m_vTemp1VelocityU[index];
			double dv = m_pParticleData->m_vTemp1VelocityV[neiIndex]-m_pParticleData->m_vTemp1VelocityV[index];	
			
			du_w_sum += du*w2;
			dv_w_sum += dv*w2;
		}		
	
	}
	if(w2_total != 0) {
		result[0] = du_w_sum/w2_total;
		result[1] = dv_w_sum/w2_total;
		
	}
	else {
		result[0] = 0;
		result[1] = 0;
	}	
}

for(size_t i=0; i<result.size(); i++) assert( !(std::isnan(result[i]) || std::isinf(result[i])) );

//cout<<"-------HyperbolicLPSolver::computeLocalWeightedAverageVelocity()-------"<<endl;

}


void HyperbolicLPSolver::updateFluidBoundingBox() {

//cout<<"-------HyperbolicLPSolver::updateFluidBoundingBox()-------"<<endl;

double *positionX = m_pParticleData->m_vPositionX;
double *positionY = m_pParticleData->m_vPositionY;
double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;

double sp = m_fTimesBoundingBox*m_fAvgParticleSpacing;
cout<<"m_fTimesBoundingBox="<<m_fTimesBoundingBox<<endl;
cout<<"m_fAvgParticleSpacing="<<m_fAvgParticleSpacing<<endl;

for(size_t p=0; p<fluidBoxes.size(); p++) {
	size_t startIndex = fluidBoxes[p]->getStartIndex();
	size_t number = fluidBoxes[p]->getNumber();
		
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"Before update:"<<endl;
	cout<<"fluidBoxes[p]->getObjectTag()="<<fluidBoxes[p]->getObjectTag()<<endl;
	cout<<"fluidBoxes[p]->getStartIndex()="<<fluidBoxes[p]->getStartIndex()<<endl;	
	cout<<"fluidBoxes[p]->getNumber()="<<fluidBoxes[p]->getNumber()<<endl;
	cout<<"fluidBoxes[p]->getXmin()="<<fluidBoxes[p]->getXmin()<<endl;
	cout<<"fluidBoxes[p]->getXmax()="<<fluidBoxes[p]->getXmax()<<endl;
	cout<<"fluidBoxes[p]->getYmin()="<<fluidBoxes[p]->getYmin()<<endl;
	cout<<"fluidBoxes[p]->getYmax()="<<fluidBoxes[p]->getYmax()<<endl;
	cout<<"fluidBoxes[p]->getZmin()="<<fluidBoxes[p]->getZmin()<<endl;
	cout<<"fluidBoxes[p]->getZmax()="<<fluidBoxes[p]->getZmax()<<endl;
	cout<<"----------------------------------------------------------"<<endl;	

	// TODO THIS IS A PLACE TO USE OPENMP LOOP LEVEL PARALLELIZATION
	fluidBoxes[p]->setXmin(
	*min_element(positionX+startIndex, positionX+startIndex+number) - sp);
	
	fluidBoxes[p]->setXmax(
	*max_element(positionX+startIndex, positionX+startIndex+number) + sp);
	
	fluidBoxes[p]->setYmin(
	*min_element(positionY+startIndex, positionY+startIndex+number) - sp);
	
	fluidBoxes[p]->setYmax(
	*max_element(positionY+startIndex, positionY+startIndex+number) + sp);
	
	if(m_iDimension==3) {
		fluidBoxes[p]->setZmin(
		*min_element(positionZ+startIndex, positionZ+startIndex+number) - sp);
		
		fluidBoxes[p]->setZmax(
		*max_element(positionZ+startIndex, positionZ+startIndex+number) + sp);
	}

	cout<<"----------------------------------------------------------"<<endl;
	cout<<"After update:"<<endl;
	cout<<"fluidBoxes[p]->getObjectTag()="<<fluidBoxes[p]->getObjectTag()<<endl;
	cout<<"fluidBoxes[p]->getStartIndex()="<<fluidBoxes[p]->getStartIndex()<<endl;
	cout<<"fluidBoxes[p]->getNumber()="<<fluidBoxes[p]->getNumber()<<endl;
	cout<<"fluidBoxes[p]->getXmin()="<<fluidBoxes[p]->getXmin()<<endl;
	cout<<"fluidBoxes[p]->getXmax()="<<fluidBoxes[p]->getXmax()<<endl;
	cout<<"fluidBoxes[p]->getYmin()="<<fluidBoxes[p]->getYmin()<<endl;
	cout<<"fluidBoxes[p]->getYmax()="<<fluidBoxes[p]->getYmax()<<endl;
	cout<<"fluidBoxes[p]->getZmin()="<<fluidBoxes[p]->getZmin()<<endl;
	cout<<"fluidBoxes[p]->getZmax()="<<fluidBoxes[p]->getZmax()<<endl;
	cout<<"----------------------------------------------------------"<<endl;	

}	

}


void HyperbolicLPSolver::checkForContactAlert() {
cout<<"-------HyperbolicLPSolver::checkForContactAlert()-------"<<endl;	

vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;

if(fluidBoxes.size() <= 1) return;

for(size_t p=0; p<fluidBoxes.size(); p++) {
	double xmin1 = fluidBoxes[p]->getXmin();
	double xmax1 = fluidBoxes[p]->getXmax();
	double ymin1 = fluidBoxes[p]->getYmin();
	double ymax1 = fluidBoxes[p]->getYmax();
	double zmin1 = fluidBoxes[p]->getZmin();
	double zmax1 = fluidBoxes[p]->getZmax();
	for(size_t q=p+1; q<fluidBoxes.size(); q++) {
		double xmin2 = fluidBoxes[q]->getXmin();
		double xmax2 = fluidBoxes[q]->getXmax();
		double ymin2 = fluidBoxes[q]->getYmin();
		double ymax2 = fluidBoxes[q]->getYmax();
		double zmin2 = fluidBoxes[q]->getZmin();
		double zmax2 = fluidBoxes[q]->getZmax();
		
		bool x_intersect = (xmax1>=xmax2 && xmin1<=xmax2) || (xmax2>=xmax1 && xmin2<=xmax1);
		bool y_intersect = (ymax1>=ymax2 && ymin1<=ymax2) || (ymax2>=ymax1 && ymin2<=ymax1);
		bool z_intersect = (zmax1>=zmax2 && zmin1<=zmax2) || (zmax2>=zmax1 && zmin2<=zmax1);

		if(x_intersect && y_intersect && z_intersect) { 
			cout<<"fluid object "<<p+1<<" is approaching fluid object "<<q+1<<endl;
			m_iContactAlert = true;
			return;
		}
	}	
}
cout<<"All fluid objects are distant so far"<<endl;

}

bool HyperbolicLPSolver::generateGhostParticle() {
//cout<<"-------HyperbolicLPSolver::generateGhostParticle()-------"<<endl;

double *positionX = m_pParticleData->m_vPositionX;
double *positionY = m_pParticleData->m_vPositionY;
double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
const size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
const size_t fluidNum = m_pParticleData->getFluidNum();	
	
// build the search structure on "fluid" particles	
//double startTime = omp_get_wtime();	
m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, fluidStartIndex, fluidNum);
//double elipsedTime = omp_get_wtime() - startTime;
//printf("Building search structure takes %.16g seconds\n", elipsedTime);

const vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;

//cout<<"m_fAvgParticleSpacing="<<m_fAvgParticleSpacing<<endl;
const double h_r = 0.5*m_fAvgParticleSpacing;//ghost density is based on "global" average inter-particle spacing
const size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	

#ifdef _OPENMP
size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	cout<<"numThreads="<<numThreads<<endl;
//	cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<endl;
vector<vector<double>> ghostX(numThreads);
vector<vector<double>> ghostY(numThreads);
vector<vector<double>> ghostZ;
if(m_pParticleData->m_iDimension == 3) ghostZ = vector<vector<double>>(numThreads);
vector<size_t> numGhost(numThreads,0);
#endif
	
//for(size_t p=0; p<fluidBoxes.size(); p++) {
//	cout<<"m_vFluidBoundingBox["<<p<<"]:"<<endl;
//	cout<<"m_iObjectTag="<<fluidBoxes[p]->getObjectTag()<<endl;
//	cout<<"m_iNumber="<<fluidBoxes[p]->getNumber()<<endl;
//	cout<<"m_fXmin="<<fluidBoxes[p]->getXmin()<<endl;
//	cout<<"m_fXmax="<<fluidBoxes[p]->getXmax()<<endl;
//	cout<<"m_fYmin="<<fluidBoxes[p]->getYmin()<<endl;
//	cout<<"m_fYmax="<<fluidBoxes[p]->getYmax()<<endl;
//	cout<<"m_fZmin="<<fluidBoxes[p]->getZmin()<<endl;
//	cout<<"m_fZmax="<<fluidBoxes[p]->getZmax()<<endl;
//}

size_t ghostIndex = m_pParticleData->getGhostStartIndex();
if(m_iDimension==2) {
	
	for(size_t p=0; p<fluidBoxes.size(); p++) {
		
		int tag = fluidBoxes[p]->getObjectTag();
		//cout<<"p="<<p<<" tag="<<tag<<endl;
		double xmin = fluidBoxes[p]->getXmin();	
		double xmax = fluidBoxes[p]->getXmax();
		double ymin = fluidBoxes[p]->getYmin();	
		double ymax = fluidBoxes[p]->getYmax();
			
		HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
		// get parameters of hexagonal packing
		size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
		hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);		

		// compute the location of particles	
		#ifdef _OPENMP
		#pragma omp parallel  
		{
		
		int tid = omp_get_thread_num();		
		#endif

		int neiListTemp[maxNeiNum]; // the index of the neighbours
		double neiListDistTemp[maxNeiNum]; // the distance of the ghsot to its neighbours //TODO NOT USED
		size_t numNeiFound;
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t j=m0; j<=m1; j++) { 
			
			if((j+1)%2 != 0) { // odd-numbered rows 
				for(size_t k=n0_odd; k<=n1_odd; k++) { 
					double x = hex2D.computeX(0,k);
					double y = hex2D.computeY(j);
				
				
					#ifdef _OPENMP				
					m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
														  neiListTemp, neiListDistTemp, numNeiFound, tid); //output
					if(isValidGhostParticle(x,y,0,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
						ghostX[tid].push_back(x);	
						ghostY[tid].push_back(y);
						numGhost[tid]++;
					}
					#else
					m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
														  neiListTemp, neiListDistTemp, numNeiFound); //output
					if(isValidGhostParticle(x,y,0,neiListTemp,neiListDistTemp,numNeiFound,tag)) {
						if(ghostIndex>=m_pParticleData->m_iCapacity) return false; // exceed array size capacity
						positionX[ghostIndex] = x;
						positionY[ghostIndex] = y;	
						ghostIndex++;
					}
					#endif
					
				}
			} 
			else{ // even-numbered rows
				for(size_t k=n0_even; k<=n1_even; k++) {
					double x = hex2D.computeX(1,k);
					double y = hex2D.computeY(j);
					
					//size_t numNeiFound;
					#ifdef _OPENMP				
					m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
														  neiListTemp, neiListDistTemp, numNeiFound, tid); //output	
					if(isValidGhostParticle(x,y,0,neiListTemp,neiListDistTemp,numNeiFound,tag)) {
						ghostX[tid].push_back(x);	
						ghostY[tid].push_back(y);
						numGhost[tid]++;
					}
					#else
					m_pNeighbourSearcher->searchNeighbour(x, y, 0, m_fNeiSearchRadius, 
														  neiListTemp, neiListDistTemp, numNeiFound); // output
					if(isValidGhostParticle(x,y,0,neiListTemp,neiListDistTemp,numNeiFound,tag)) {
						if(ghostIndex>=m_pParticleData->m_iCapacity) return false; // exceed array size capacity
						positionX[ghostIndex] = x;
						positionY[ghostIndex] = y;	
						ghostIndex++;
					}
					#endif		
						
				}
			}
		} 	
		
		#ifdef _OPENMP	
		} // omp parallel section
		#endif
	}	
}
else if(m_iDimension==3) {

	for(size_t p=0; p<fluidBoxes.size(); p++) {
		
		int tag = fluidBoxes[p]->getObjectTag();

		double xmin = fluidBoxes[p]->getXmin();	
		double xmax = fluidBoxes[p]->getXmax();
		double ymin = fluidBoxes[p]->getYmin();	
		double ymax = fluidBoxes[p]->getYmax();
		double zmin = fluidBoxes[p]->getZmin();	
		double zmax = fluidBoxes[p]->getZmax();

		HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
		//get parameters of hexagonal packing
		size_t l0,l1;
		size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
		size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
		hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
							n0_odd, n1_odd, n0_even, n1_even, 
							nn0_odd, nn1_odd, nn0_even, nn1_even);	
		
		#ifdef _OPENMP
		#pragma omp parallel 
		{
		int tid = omp_get_thread_num();			
		#endif
		
		int neiListTemp[maxNeiNum]; // the index of the neighbours
		double neiListDistTemp[maxNeiNum]; // the distance of the ghsot to its neighbours //TODO NOT USED
		size_t numNeiFound;
		
		#ifdef _OPENMP
		#pragma omp for 	
		#endif
		for(size_t i=l0; i<=l1; i++) { // compute the location of particles
			if((i+1)%2 != 0) { //odd-numbered layers
				for(size_t j=m0_odd; j<=m1_odd; j++) { 
					if((j+1)%2 != 0) { //odd-numbered rows 
						for(size_t k=n0_odd; k<=n1_odd; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);	
							//cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;	
							//size_t numNeiFound;
							#ifdef _OPENMP
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound,tid); 
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {	
								ghostX[tid].push_back(x);	
								ghostY[tid].push_back(y);
								ghostZ[tid].push_back(z);
								numGhost[tid]++;
							}
							#else
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound); // output
							
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {	
								if(ghostIndex>=m_pParticleData->m_iCapacity) return false;	
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;
							}
							#endif
							
						}
					} 
					else{ //even-numbered rows
						for(size_t k=n0_even; k<=n1_even; k++) {
							double x = hex3D.computeX(1,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);
							//cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
							//size_t numNeiFound;
							#ifdef _OPENMP
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound,tid); 
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								ghostX[tid].push_back(x);	
								ghostY[tid].push_back(y);
								ghostZ[tid].push_back(z);
								numGhost[tid]++;
							}
							#else
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound); // output
							
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								if(ghostIndex>=m_pParticleData->m_iCapacity) return false;	
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;
							}
							#endif		

						}
					}
				}
					
			} 
			else { //even-numbered layers
				for(size_t j=m0_even; j<=m1_even; j++) { 
					if((j+1)%2 != 0) { //odd-numbered rows
						for(size_t k=nn0_odd; k<=nn1_odd; k++) { 
							double x = hex3D.computeX(1,k);
							double y = hex3D.computeY(1,j);
							double z = hex3D.computeZ(i);
							//cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
							//size_t numNeiFound;
							#ifdef _OPENMP
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound,tid);
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								ghostX[tid].push_back(x);	
								ghostY[tid].push_back(y);
								ghostZ[tid].push_back(z);
								numGhost[tid]++;
							}
							#else
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound); // output
							
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								if(ghostIndex>=m_pParticleData->m_iCapacity) return false;	
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;
							}
							#endif		

						}
					} 
					else { //even-numbered rows
						for(size_t k=nn0_even; k<=nn1_even; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(1,j);
							double z = hex3D.computeZ(i);
							//cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
							//size_t numNeiFound;
							#ifdef _OPENMP
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound,tid); 
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								ghostX[tid].push_back(x);	
								ghostY[tid].push_back(y);
								ghostZ[tid].push_back(z);
								numGhost[tid]++;
							}
							#else
							m_pNeighbourSearcher->searchNeighbour(x, y, z, m_fNeiSearchRadius, 
																  neiListTemp,neiListDistTemp,numNeiFound); // output
		
							if(isValidGhostParticle(x,y,z,neiListTemp,neiListDistTemp,numNeiFound,tag)) {		
								if(ghostIndex>=m_pParticleData->m_iCapacity) return false;	
								positionX[ghostIndex] = x;
								positionY[ghostIndex] = y;
								positionZ[ghostIndex] = z;	
								ghostIndex++;
							}
							#endif		
						}
					}
				}	
			}    
		}
		
		#ifdef _OPENMP	
		} //omp parallel section
		#endif

	}	

}

#ifdef _OPENMP
size_t sum = 0;
for(size_t i=0; i<numThreads; i++) 
	sum += numGhost[i];
if(ghostIndex+sum > m_pParticleData->m_iCapacity) return false; // not enough space -> augment data array size!

if(m_iDimension==3) {
	for(size_t i=0; i<numThreads; i++) {
		for(size_t j=0; j<numGhost[i]; j++) {
			positionX[ghostIndex] = ghostX[i][j];
			positionY[ghostIndex] = ghostY[i][j];
			positionZ[ghostIndex] = ghostZ[i][j];
			ghostIndex++;
		}
	}
}
else if(m_iDimension==2) {
	for(size_t i=0; i<numThreads; i++) {
		for(size_t j=0; j<numGhost[i]; j++) {
			positionX[ghostIndex] = ghostX[i][j];
			positionY[ghostIndex] = ghostY[i][j];
			ghostIndex++;
		}
	}
}
#endif


m_pParticleData->m_iGhostNum = ghostIndex - m_pParticleData->m_iGhostStartIndex;
m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
							   m_pParticleData->m_iBoundaryNum + 
							   m_pParticleData->m_iGhostNum;


//cout<<"m_pParticleData->m_iGhostStartIndex="<<m_pParticleData->m_iGhostStartIndex<<endl;
cout<<"m_pParticleData->m_iFluidNum="<<m_pParticleData->m_iFluidNum<<endl;
cout<<"m_pParticleData->m_iGhostNum="<<m_pParticleData->m_iGhostNum<<endl;
//cout<<"m_pParticleData->m_iTotalNum="<<m_pParticleData->m_iTotalNum<<endl;
//cout<<"---------------------------------------------------------"<<endl;

return true;
}


void HyperbolicLPSolver::computeAvgNumParticleWithinSearchRadius() {

int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();

double sum = 0;
for(size_t i=fluidStartIndex; i<fluidEndIndex; i++) {
	sum += neighbourListSize[i];
}

cout<<"m_iNumParticleWithinSearchRadius changed from "<<m_iNumParticleWithinSearchRadius;
m_iNumParticleWithinSearchRadius = sum/(double)m_pParticleData->getFluidNum();
cout<<" to "<<m_iNumParticleWithinSearchRadius<<endl;
}


bool HyperbolicLPSolver::isValidGhostParticle(double x, double y, double z, 
const int* neiList, size_t numNei, int objectTag) {

//cout<<"objectTag="<<objectTag<<endl;

// 1. have at leat one fluid neighbour
if(numNei < 1) return false; 

// 2. the distance from the nearest neighbour is not too small
double xDiff = m_pParticleData->m_vPositionX[neiList[0]] - x;
double yDiff = m_pParticleData->m_vPositionY[neiList[0]] - y;
double zDiff = m_pParticleData->m_vPositionZ[neiList[0]] - z;
double dist = sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff);
if(dist < 1.1*m_fAvgParticleSpacing) return false;
//cout<<"m_fContactLength="<<m_fContactLength<<"  neiListDist[0]="<<neiListDist[0]<<endl;
//if(neiListDist[0] < m_fContactLength) return false;

// 3. no neighbour comes from other fluid object 
if(m_pParticleData->m_vFluidBoundingBox.size() > 1) {
//if(m_iContactAlert) {
	for(size_t l=0; l<numNei; l++) {
		
		if(m_pParticleData->m_vObjectTag[neiList[l]] != objectTag &&
		   m_pParticleData->m_vObjectTag[neiList[l]] != -objectTag) { 
			cout<<"objectTag="<<objectTag<<endl;
			cout<<"(x,y,z)="<<"("<<x<<","<<y<<","<<z<<")"<<endl;
			cout<<"m_pParticleData->m_vObjectTag[neiList["<<l<<"]]="
			<<m_pParticleData->m_vObjectTag[neiList[l]]<<endl;		
			return false;
		}
	}
}

// 3. the percentage of fluid is not too large
//cout<<"numNei="<<numNei<<"  m_iNumParticleWithinSearchRadius="<<m_iNumParticleWithinSearchRadius<<endl;
//cout<<(double)numNei/(double)m_iNumParticleWithinSearchRadius<<endl;
//if((double)numNei/(double)m_iNumParticleWithinSearchRadius > 0.5) return false;
			

//TODO the ObjectTag may be replaced by info in bounding box
// 4. cannot have neighbours from different fluid objects
// only do when m_iContactAlert is true
//if(m_iContactAlert) { 
//	int firstTag = abs(m_pParticleData->m_vObjectTag[neiList[0]]); // first tag, +-num are from same object
//	for(size_t l=1; l<numNei; l++) 
//		if(abs(m_pParticleData->m_vObjectTag[neiList[l]]) != firstTag) return false; 	
//}

// passed all 4 criteria
return true;
}

bool HyperbolicLPSolver::isValidGhostParticle(double x, double y, double z, 
const int* neiList, const double* neiListDist, size_t numNei, int objectTag) {

//cout<<"objectTag="<<objectTag<<endl;

// 1. have at leat one fluid neighbour
if(numNei < 1) return false; 

// 2. the distance from the nearest neighbour is not too small
//double xDiff = m_pParticleData->m_vPositionX[neiList[0]] - x;
//double yDiff = m_pParticleData->m_vPositionY[neiList[0]] - y;
//double zDiff = m_pParticleData->m_vPositionZ[neiList[0]] - z;
//double dist = sqrt(xDiff*xDiff+yDiff*yDiff+zDiff*zDiff);
//if(dist < 1.1*m_fAvgParticleSpacing) return false;
//cout<<"m_fContactLength="<<m_fContactLength<<"  neiListDist[0]="<<neiListDist[0]<<endl;
if(neiListDist[0] < 1.1*m_fAvgParticleSpacing) return false;

// 3. no neighbour comes from other fluid object 
if(m_pParticleData->m_vFluidBoundingBox.size() > 1) {
//if(m_iContactAlert) {
	for(size_t l=0; l<numNei; l++) {
		
		if(m_pParticleData->m_vObjectTag[neiList[l]] != objectTag &&
		   m_pParticleData->m_vObjectTag[neiList[l]] != -objectTag) { 
			//cout<<"objectTag="<<objectTag<<endl;
			//cout<<"(x,y,z)="<<"("<<x<<","<<y<<","<<z<<")"<<endl;
			//cout<<"m_pParticleData->m_vObjectTag[neiList["<<l<<"]]="
			//<<m_pParticleData->m_vObjectTag[neiList[l]]<<endl;		
			return false;
		}
	}
}

// 3. the percentage of fluid is not too large
//cout<<"numNei="<<numNei<<"  m_iNumParticleWithinSearchRadius="<<m_iNumParticleWithinSearchRadius<<endl;
//cout<<(double)numNei/(double)m_iNumParticleWithinSearchRadius<<endl;
//if((double)numNei/(double)m_iNumParticleWithinSearchRadius > 0.5) return false;
			

//TODO the ObjectTag may be replaced by info in bounding box
// 4. cannot have neighbours from different fluid objects
// only do when m_iContactAlert is true
//if(m_iContactAlert) { 
//	int firstTag = abs(m_pParticleData->m_vObjectTag[neiList[0]]); // first tag, +-num are from same object
//	for(size_t l=1; l<numNei; l++) 
//		if(abs(m_pParticleData->m_vObjectTag[neiList[l]]) != firstTag) return false; 	
//}

// passed all 4 criteria
return true;
}


void HyperbolicLPSolver::searchNeighbourForFluidParticle() {
	searchNeighbourForFluidParticle(0);
}


void HyperbolicLPSolver::searchNeighbourForFluidParticle(int choice) {
	
	cout<<"-------HyperbolicLPSolver::searchNeighbourForFluidParticle()-------"<<endl;

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case	
	const double *mass = m_pParticleData->m_vMass;
	double *VolumeVoronoi= m_pParticleData->m_vVolumeVoronoi;

	if(choice==2)
	{
		std::cout<<"Voronoi density estimator has been removed from LP code. SPH density estimator is used.\n"<<std::endl;
		choice=1;
/*		std::cout<<"Calculating Voronoi Area"<<std::endl;
	        VoronoiAreaEstimator voronoi(m_iDimension, m_pParticleData->m_iFluidNum + m_pParticleData->m_iBoundaryNum, positionX, positionY, positionZ, mass, VolumeVoronoi);
		int voronoi_error=voronoi.ComputeVoronoiArea();
		if (voronoi_error){
			std::cout<<"Error in voronoi area estimator"<<std::endl;
		}*/
	}

//	double startTime = omp_get_wtime();
//	cout<<"Start to build octree"<<endl;
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, mass, VolumeVoronoi,
	m_pParticleData->m_iFluidStartIndex ,m_pParticleData->m_iFluidNum + m_pParticleData->m_iBoundaryNum);		
//	printf("Building tree takes %.16g seconds\n", omp_get_wtime() - startTime);
//	cout<<"end building octree"<<endl;

	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;		
	 
//	startTime = omp_get_wtime();
//	cout<<"Search neighbours and calculate density"<<endl;
	#ifdef _OPENMP
	#pragma omp parallel
	{
	
	int tid = omp_get_thread_num();
	
	#endif	
	
	double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
	size_t numNeiFound;	

	// fluid
	//if(m_iContactAlert) { // TODO multiple fluid objects 
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1) { // multiple fluid objects
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			/*	
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			*/
				
			double radius = m_fTimesNeiSearchRadius * sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbourQuadrant(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbourQuadrant(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif

//sph density estimator, the option is keeped but not used. The current version call the SPH density estimator using a seperate function, SPHDensityEstimatorForFluidParticle 
			if(choice==1){
                        double radius2multiplier=3.0;
                        double radius2 = radius2multiplier*sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
                        double count_density=1.0/m_pParticleData->m_vVolumeOld[index];
                        m_pNeighbourSearcher->densityEstimator(
                                positionX[index],positionY[index],positionZ[index],radius2,
                                &count_density,m_pParticleData->m_vVolume_x[index],m_pParticleData->m_vVolume_y[index],m_pParticleData->m_vVolume_z[index]);
                        if(m_iDimension==2){
                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                        }
			else{//3d TODO
			}
			}
//end of sph density estimator

			if(numNeiFound > maxNeiNum) {
				cout<<"numNeiFound="<<numNeiFound<<" > maxNeiNum="<<maxNeiNum<<endl;
				assert(false);
			}

			neighbourListSize[index] = numNeiFound;		    
			changeFluidObjectTagIfInContact(index,numNeiFound,neiListDist);
//disk collision
//			changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(index,numNeiFound);
		}

	}
	else { // single fluid object
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;	
			/*
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif	
			*/
				
			double radius;
                        if(m_iDimension==2) radius = m_fTimesNeiSearchRadius * sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
			if(m_iDimension==3) radius = m_fTimesNeiSearchRadius * cbrt(1.4142*mass[index]*m_pParticleData->m_vVolumeOld[index]);	

//Using octree to search neighbours. The neighbourhood should contains similar number of particles from each direction.
			if(m_iDimension==3) {
	                        #ifdef _OPENMP
        	                m_pNeighbourSearcher->searchNeighbour(//Direction
                	                positionX[index],positionY[index],positionZ[index],radius, 
                        	        neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output   
                        	#else
	                        m_pNeighbourSearcher->searchNeighbour(//Direction
        	                        positionX[index],positionY[index],positionZ[index],radius, 
                	                neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
                        	#endif
			}
                        if(m_iDimension==2) {
                                #ifdef _OPENMP
                                m_pNeighbourSearcher->searchNeighbour(//Quadrant
                                        positionX[index],positionY[index],positionZ[index],radius,
                                        neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output   
                                #else
                                m_pNeighbourSearcher->searchNeighbour(//Quadrant
                                        positionX[index],positionY[index],positionZ[index],radius,
                                        neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
                                #endif
                        }

			if(0)//BNL nozzle leaking solution
			{
				int bnlj=neiListStartIndex;
				for(int bnli=bnlj;bnli<neiListStartIndex+numNeiFound;bnli++)
				{
					if(neighbourList[bnli]>=m_pParticleData->m_iBoundaryStartIndex + m_pParticleData->m_iInflowNum && neighbourList[bnli]<m_pParticleData->m_iBoundaryStartIndex + m_pParticleData->m_iBoundaryNum)
					{
//						if((positionX[index]-0.004)*(positionX[m_vMirrorIndex[neighbourList[bnli]-m_pParticleData->m_iBoundaryStartIndex - m_pParticleData->m_iInflowNum]]-0.004)<0)
						if((positionX[index]<0.004)&&(positionX[m_vMirrorIndex[neighbourList[bnli]-m_pParticleData->m_iBoundaryStartIndex - m_pParticleData->m_iInflowNum]]>0.004))
						{
//							std::cout<<positionX[index]<<" "<<positionX[m_vMirrorIndex[neighbourList[bnli]-m_pParticleData->m_iBoundaryStartIndex - m_pParticleData->m_iInflowNum]]<<std::endl;
							continue;
						}
					}
					neighbourList[bnlj++]=neighbourList[bnli];
				}
				numNeiFound=bnlj-neiListStartIndex;

			}

//SPH density estimator option is keeped but not used. The current LP algorithm call sthe SPH density estimator using a seperate function, SPHDensityEstimatorForFluidParticle 
			if (choice==1){
                        double radius2multiplier;
			if (m_iDimension==2)
				radius2multiplier=5.0;
			if(m_iDimension==3)
				radius2multiplier=5.0;
			double radius2;
			if(m_iDimension==2) radius2 = radius2multiplier*sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
                        if(m_iDimension==3) radius2 = radius2multiplier*cbrt(1.4142*mass[index]*m_pParticleData->m_vVolumeOld[index]);
                        double count_density=1.0/m_pParticleData->m_vVolumeOld[index];
			if(m_iDimension==3) count_density=5.0*radius2/7.0*count_density;
                        m_pNeighbourSearcher->densityEstimator(
                                index, positionX[index],positionY[index],positionZ[index],radius2,
                                &count_density,m_pParticleData->m_vVolume_x[index],m_pParticleData->m_vVolume_y[index],m_pParticleData->m_vVolume_z[index]);
                        if(m_iDimension==2){
                                m_pParticleData->m_vVolume[index]=1.0/count_density;  
                        }
                        else{//3d
				m_pParticleData->m_vVolume[index]=5.0*radius2/7.0/count_density;
                        }
			}
//end of sph density estimator
//Voronoi density estimator option is keeped but not used.
                        if (choice==2){
                        double radius2multiplier=5.0;
                        double radius2 = radius2multiplier*sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
			double count_density=0;
                        m_pNeighbourSearcher->VoronoiDensityEstimator(index,
                                positionX[index],positionY[index],positionZ[index],radius2,
                                &count_density);
                        if(m_iDimension==2){
                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                        }
                        else{//3d TODO
                        }
                        }
//end of voronoi density estimator


			if(numNeiFound > maxNeiNum) {
				cout<<"numNeiFound="<<numNeiFound<<" > maxNeiNum="<<maxNeiNum<<endl;
				assert(false);
			}

			neighbourListSize[index] = numNeiFound;
			//cout<<"numNeiFound"<<numNeiFound<<endl;
			//cout<<"(x,y)=("<<positionX[index]<<","<<positionY[index]<<endl;
		}
	}
	
    #ifdef _OPENMP
	}
	#endif 
	
//	printf("Searching neighbours takes %.16g seconds\n", omp_get_wtime() - startTime);

	cout<<"-------END HyperbolicLPSolver::searchNeighbourForFluidParticle()-------"<<endl;

}

void HyperbolicLPSolver::SPHDensityEstimatorForFluidParticle(int choice) {
	if(choice)
	        cout<<"-------HyperbolicLPSolver::SPHDensityEstimatorForFluidParticle()-------"<<endl;

        const double *positionX = m_pParticleData->m_vPositionX;
        const double *positionY = m_pParticleData->m_vPositionY;
        const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case 
        const double *mass = m_pParticleData->m_vMass;
        double *VolumeVoronoi= m_pParticleData->m_vVolumeVoronoi;

        if(choice==2)
        {
                std::cout<<"Voronoi density estimator has been removed from LP code. SPH density estimator is used.\n"<<std::endl;
                choice=1;
        }
	if(choice==1||choice==3){
        m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, mass, VolumeVoronoi,
        m_pParticleData->m_iFluidStartIndex ,m_pParticleData->m_iFluidNum + m_pParticleData->m_iBoundaryNum);

//        double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
        size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
        size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();

        cout<<"Calculate density"<<endl;
        #ifdef _OPENMP
        #pragma omp parallel
        {

//        int tid = omp_get_thread_num();

        #endif

//        double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
//        size_t numNeiFound;

        if(m_pParticleData->m_vFluidBoundingBox.size() > 1) { // multiple fluid objects

                #ifdef _OPENMP
                #pragma omp for
                #endif
                for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
                        if(choice==1){
                        double radius2multiplier=5.0;
//                        double radius2 = sqrt(m_pParticleData->m_vVolume[index])*m_fInitParticleSpacing*radius2multiplier;
                        double radius2 = radius2multiplier*sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
                        double count_density=1.0/m_pParticleData->m_vVolumeOld[index];
                        m_pNeighbourSearcher->densityEstimator(
                                positionX[index],positionY[index],positionZ[index],radius2,
                                &count_density,m_pParticleData->m_vVolume_x[index],m_pParticleData->m_vVolume_y[index],m_pParticleData->m_vVolume_z[index]);
                        if(m_iDimension==2){
                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                        }
                        else{//3d TODO
                        }
                        }
//end of sph density estimator
		}
	}
        else { // single fluid object
                #ifdef _OPENMP
                #pragma omp for
                #endif
//                for(size_t index=fluidStartIndex; index<fluidEndIndex+m_pParticleData->getInflowNum(); index++) {
                for(size_t index=fluidStartIndex; index<fluidEndIndex+m_pParticleData->getInflowNum(); index++) {
//sph density estimator
//			double tolerance=1;
//			double gra_mag=sqrt(m_pParticleData->m_vVolume_x[index]*m_pParticleData->m_vVolume_x[index]+m_pParticleData->m_vVolume_y[index]*m_pParticleData->m_vVolume_y[index]+m_pParticleData->m_vVolume_z[index]*m_pParticleData->m_vVolume_z[index]);
//			m_pParticleData->m_vIfSPHDensity[index]=0;
//                        if (choice==1||((choice==3)&&(gra_mag>tolerance))){
			double alpha=3.0;
			double max_gradient=3.0;
			double vmin=1e10;
			double vmax=0;
			double ratio=1.5;
			if(choice==1||choice==3){
//                        m_pParticleData->m_vIfSPHDensity[index]+=1;
                        double radius2multiplier;
                        if (m_iDimension==2)
                                radius2multiplier=5.0;
                        if(m_iDimension==3)
                                radius2multiplier=5.0;
			if(choice==3)
				radius2multiplier=8.0;
                        double gra_mag=sqrt(m_pParticleData->m_vVolume_x[index]*m_pParticleData->m_vVolume_x[index]+m_pParticleData->m_vVolume_y[index]*m_pParticleData->m_vVolume_y[index]+m_pParticleData->m_vVolume_z[index]*m_pParticleData->m_vVolume_z[index]);
                        double radius2;
                        if(m_iDimension==2) radius2 = radius2multiplier*sqrt(2*mass[index]*m_pParticleData->m_vVolumeOld[index]/1.7321);
                        if(m_iDimension==3) radius2 = radius2multiplier*cbrt(1.4142*mass[index]*m_pParticleData->m_vVolumeOld[index]);
//                      cout<<radius2<<endl;
                        double count_density=1.0/m_pParticleData->m_vVolumeOld[index];
                        if(m_iDimension==3) count_density=5.0*radius2/7.0*count_density;
			if(0)
	                        m_pNeighbourSearcher->densityEstimator(
        	                        index, positionX[index],positionY[index],positionZ[index],radius2,
                	                &count_density,m_pParticleData->m_vVolume_x[index],m_pParticleData->m_vVolume_y[index],m_pParticleData->m_vVolume_z[index]);
			if(choice==3 || choice==1)
                                m_pNeighbourSearcher->densityEstimator(
                                        index, positionX[index],positionY[index],positionZ[index],radius2,
                                        &count_density,m_pParticleData->m_vVolumeOld,&vmin,&vmax);
//                                        &count_density,m_pParticleData->m_vVolume_x[index],m_pParticleData->m_vVolume_y[index],m_pParticleData->m_vVolume_z[index]);
                        if(m_iDimension==2){
				if(choice==1)
	                                m_pParticleData->m_vVolume[index]=1.0/count_density;
				if(choice==3)
				{
/*					if(vmax/vmin>ratio)
					{
                                                m_pParticleData->m_vIfSPHDensity[index]+=1;
                                                m_pParticleData->m_vVolume[index]=1.0/count_density;
					}*/
					if((m_pParticleData->m_vVolume[index]*count_density)>alpha||(m_pParticleData->m_vVolume[index]*count_density)<1.0/alpha)
					{
						m_pParticleData->m_vIfSPHDensity[index]+=1;
						m_pParticleData->m_vVolume[index]=1.0/count_density;
					}
/*					if(gra_mag>max_gradient)
                                        {
                                                m_pParticleData->m_vIfSPHDensity[index]+=1;
                                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                                        }*/
				}	
                        }
                        else{//3d
				if(choice==1)
                                	m_pParticleData->m_vVolume[index]=5.0*radius2/7.0/count_density;
                                if(choice==3)
                                {
					count_density=count_density*7.0/5.0/radius2;
/*                                        if(vmax/vmin>ratio)
                                        {
                                                m_pParticleData->m_vIfSPHDensity[index]+=1;
                                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                                        }*/
                                        if((m_pParticleData->m_vVolume[index]*count_density)>alpha||(m_pParticleData->m_vVolume[index]*count_density)<1.0/alpha)
                                        {
                                                m_pParticleData->m_vIfSPHDensity[index]+=1;
                                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                                        }
/*                                        if(gra_mag>max_gradient)
                                        {
                                                m_pParticleData->m_vIfSPHDensity[index]+=1;
                                                m_pParticleData->m_vVolume[index]=1.0/count_density;
                                        }*/
                                }

                        }
                        }
//end of sph density estimator
                }
        }

        #ifdef _OPENMP
        }
        #endif
	}
	if(choice)
	        cout<<"-------END HyperbolicLPSolver::SPHDensityEstimatorForFluidParticle()-------"<<endl;

}







void HyperbolicLPSolver::searchNeighbourForGhostParticle() {
	
	cout<<"-------HyperbolicLPSolver::searchNeighbourForGhostParticle()-------"<<endl;

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	
	//double startTime = omp_get_wtime();
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, 
	m_pParticleData->m_iFluidStartIndex ,m_pParticleData->m_iFluidNum);	
	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Building search structure takes %.16g seconds\n", elipsedTime);
	
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

//	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
//	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
//	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
//	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;		
	 

	#ifdef _OPENMP
	#pragma omp parallel
	{
	
	int tid = omp_get_thread_num();
	
	#endif	
	
	double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
	size_t numNeiFound;	

	// ghost
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;	
		double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		
		//cout<<"index="<<index<<"  numNeiFound="<<numNeiFound<<endl;
		/*
		//ghost particle take only fluid particles as neighbours
		size_t incr = 0;
		for(size_t k=0; k<numNeiFound; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;
		*/

		neighbourListSize[index] = numNeiFound;	
	//	changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	
	}
	
    #ifdef _OPENMP
	}
	#endif 
	

	cout<<"-------END HyperbolicLPSolver::searchNeighbourForGhostParticle()-------"<<endl;
}




void HyperbolicLPSolver::searchNeighbourForAllParticle() {
	
	cout<<"-------HyperbolicLPSolver::searchNeighbourForAllParticle()-------"<<endl;

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	
	const size_t totalNumParticle = m_pParticleData->getTotalNum();

	//double startTime = omp_get_wtime();
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, totalNumParticle);	
	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Building search structure takes %.16g seconds\n", elipsedTime);
	
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	//double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;		
	 

	#ifdef _OPENMP
	#pragma omp parallel
	{
	
	int tid = omp_get_thread_num();
	
	#endif	
	
	double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
	size_t numNeiFound;	

	// fluid
	//if(m_iContactAlert) { // multiple fluid objects 
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1) { // multiple fluid objects
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			
			/*	
			double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			*/
			neighbourListSize[index] = numNeiFound;		    
			changeFluidObjectTagIfInContact(index,numNeiFound,neiListDist);
			changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(index,numNeiFound);
		}

	}
	else { // single fluid object
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif	
			
			/*	
			assert(localParSpacing[index] >= m_fAvgParticleSpacing);
			double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			*/
			neighbourListSize[index] = numNeiFound;	
			//cout<<"numNeiFound"<<numNeiFound<<endl;
			//cout<<"(x,y)=("<<positionX[index]<<","<<positionY[index]<<endl;
		}
	}
	

	// boundary
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;		
		
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		
		/*	
		double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		*/
		neighbourListSize[index] = numNeiFound;
		changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	}

		
	// ghost
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;	
		
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
	
		//ghost particle take only fluid particles as neighbours
		size_t incr = 0;
		for(size_t k=0; k<numNeiFound; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

	//	neighbourListSize[index] = numNeiFound;	
	//	changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	
	}
	
	

	#ifdef _OPENMP
	}
	#endif
	
	
	


//	cout<<"totalNumParticle="<<totalNumParticle<<endl;
//	cout<<"maxNeiNum="<<maxNeiNum<<endl;
//	cout<<"Fluid particles:"<<endl;
//	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
//	cout<<"Boundary particles:"<<endl;
//	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
//	cout<<"Ghost particles:"<<endl;
//	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
//		
//		size_t neiListStartIndex = index*maxNeiNum;			
//		cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//		cout<<"-----------------------------"<<endl;
//		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
//			cout<<"neiIndex="<<neighbourList[k]<<endl;
//		cout<<"-----------------------------"<<endl;
//	}
	cout<<"-----------------------------------------------------------------"<<endl;
	

}



void HyperbolicLPSolver::searchNeighbourForAllParticleVariableSearchRadius() {
	
	cout<<"-------HyperbolicLPSolver::searchNeighbourForAllParticleVariableSearchRadius()-------"<<endl;

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	
	const size_t totalNumParticle = m_pParticleData->getTotalNum();

	//double startTime = omp_get_wtime();
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, totalNumParticle);	
	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Building search structure takes %.16g seconds\n", elipsedTime);
	
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;		
	 

	#ifdef _OPENMP
	#pragma omp parallel
	{
	
	int tid = omp_get_thread_num();
	
	#endif	
	
	double neiListDist[maxNeiNum]; // a temp array for dist between a particle and its neighbours
	size_t numNeiFound;	

	// fluid
	//if(m_iContactAlert) { // multiple fluid objects 
	if(m_pParticleData->m_vFluidBoundingBox.size() > 1) { // multiple fluid objects
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			
			/*
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			*/
				
			double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			
			neighbourListSize[index] = numNeiFound;		    
			changeFluidObjectTagIfInContact(index,numNeiFound,neiListDist);
			changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(index,numNeiFound);
		}

	}
	else { // single fluid object
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
			size_t neiListStartIndex = index*maxNeiNum;		
			/*
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif	
			*/
				
			assert(localParSpacing[index] >= m_fAvgParticleSpacing);
			double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
			#ifdef _OPENMP
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
			#else
			m_pNeighbourSearcher->searchNeighbour(
				positionX[index],positionY[index],positionZ[index],radius, 
				neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
			#endif
			
			neighbourListSize[index] = numNeiFound;	
			//cout<<"numNeiFound"<<numNeiFound<<endl;
			//cout<<"(x,y)=("<<positionX[index]<<","<<positionY[index]<<endl;
		}
	}
	

	// boundary
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;		
		/*	
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		*/
			
		double radius = m_fTimesNeiSearchRadius * localParSpacing[index];
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],radius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
		
		neighbourListSize[index] = numNeiFound;
		changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	}

		
	// ghost
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;	
		
		#ifdef _OPENMP
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,tid,index); // output	
		#else
		m_pNeighbourSearcher->searchNeighbour(
			positionX[index],positionY[index],positionZ[index],m_fNeiSearchRadius, 
			neighbourList+neiListStartIndex,neiListDist,numNeiFound,index); // output
		#endif
	
		//ghost particle take only fluid particles as neighbours
		size_t incr = 0;
		for(size_t k=0; k<numNeiFound; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

	//	neighbourListSize[index] = numNeiFound;	
	//	changeNeighbourhoodToIncludeOnlyFluidNei(index,numNeiFound);		    
	
	}
	
	

	#ifdef _OPENMP
	}
	#endif
	

	cout<<"-----------------------------------------------------------------"<<endl;
	

}


void HyperbolicLPSolver::identifyFluidBoundary() {
	
	double *x = m_pParticleData->m_vPositionX;
	double *y = m_pParticleData->m_vPositionY;
	double *z = m_pParticleData->m_vPositionZ;
	//int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
	//size_t ghostIndex = m_pParticleData->getGhostStartIndex();
	//cout<<"ghostStartIndex="<<ghostIndex<<endl;

	if(m_iDimension==2) {
	
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=m_pParticleData->getFluidStartIndex(); 
		index<m_pParticleData->getFluidStartIndex()+m_pParticleData->getFluidNum(); index++)  {
			
			double x0 = x[index], y0 = y[index];
			size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
			
			int count[8] = {0};
			
			//cout<<"index="<<index<<endl;
			//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
			size_t other = 0;
			m_vFillGhost[index]=false;
			for(int i=0; i<neighbourListSize[index]; i++) {
				
				size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
				
				double x1 = x[neiIndex], y1 = y[neiIndex];
				double dx = x1-x0, dy = y1-y0;

				if(dy>0 && dy<dx)       count[0]++;
				else if(dx>0 && dy>dx)  count[1]++;
				else if(dx<0 && dy>-dx) count[2]++;
				else if(dy>0 && dy<-dx) count[3]++;
				else if(dy<0 && dy>dx)  count[4]++;
				else if(dx<0 && dy<dx)  count[5]++;
				else if(dx>0 && dy<-dx) count[6]++;
				else if(dy<0 && dy>-dx) count[7]++;
				else other++;	
			}
			for(int i=0; i<8; i++) {
				//cout<<"count["<<i<<"]="<<count[i]<<endl;
				if(count[i]==0) {
					m_vFillGhost[index] = true;
					break;	
				}
			}
			//cout<<"other="<<other<<endl;
		}	
	}
	else if(m_iDimension==3) {
		
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=m_pParticleData->getFluidStartIndex(); 
		index<m_pParticleData->getFluidStartIndex()+m_pParticleData->getFluidNum(); index++)  {
			
			double x0 = x[index], y0 = y[index], z0 = z[index];
			size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
			
			int count[16] = {0};
			
			//cout<<"index="<<index<<endl;
			//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
			size_t other = 0;
			m_vFillGhost[index]=false;
			for(int i=0; i<neighbourListSize[index]; i++) {
				
				size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
				
				double x1 = x[neiIndex], y1 = y[neiIndex], z1 = z[neiIndex];
				double dx = x1-x0, dy = y1-y0, dz = z1-z0;
				if(dz>0) {
					if(dy>0 && dy<dx)       count[0]++;
					else if(dx>0 && dy>dx)  count[1]++;
					else if(dx<0 && dy>-dx) count[2]++;
					else if(dy>0 && dy<-dx) count[3]++;
					else if(dy<0 && dy>dx)  count[4]++;
					else if(dx<0 && dy<dx)  count[5]++;
					else if(dx>0 && dy<-dx) count[6]++;
					else if(dy<0 && dy>-dx) count[7]++;
					else other++;
				}
				if(dz<0) {
					if(dy>0 && dy<dx)       count[8]++;
					else if(dx>0 && dy>dx)  count[9]++;
					else if(dx<0 && dy>-dx) count[10]++;
					else if(dy>0 && dy<-dx) count[11]++;
					else if(dy<0 && dy>dx)  count[12]++;
					else if(dx<0 && dy<dx)  count[13]++;
					else if(dx>0 && dy<-dx) count[14]++;
					else if(dy<0 && dy>-dx) count[15]++;
					else other++;
				}
				else other++;
					
			}
			for(int i=0; i<16; i++) {
				//cout<<"count["<<i<<"]="<<count[i]<<endl;
				if(count[i]==0) {
					m_vFillGhost[index] = true;
					break;	
				}
			}
			//cout<<"other="<<other<<endl;
		}	
	}
}


bool HyperbolicLPSolver::generateSolidBoundaryByMirrorParticles() {
	
	cout<<"-------HyperbolicLPSolver::generateSolidBoundaryByMirrorParticles()-------"<<endl;	

	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;


	#ifdef _OPENMP
	size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
	vector<vector<double>> bX(numThreads);
	vector<vector<double>> bY(numThreads);
	vector<vector<double>> bZ(numThreads);
	vector<vector<double>> bPressure(numThreads);
	vector<vector<double>> bVx(numThreads);
	vector<vector<double>> bVy(numThreads);
	vector<vector<double>> bVz(numThreads);
	vector<vector<size_t>> fIndex(numThreads);// the corresponding fluid index of a mirror particle
/*        for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {

                if(m_pParticleData->m_vBoundaryObjTypes[p]!="solid") continue;
                cout<<"m_pParticleData->m_vBoundaryObjTypes["<<p<<"]="<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;
	}*/

        for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {
                if(m_pParticleData->m_vBoundaryObjTypes[p]!="outflow") continue;
//		std::cout<<"Outflow begin"<<p<<std::endl;
                m_pParticleData->m_vBoundaryObj[p]->UpdateInflowBoundary(m_pParticleData,m_pEOS,m_fDt,m_fInitParticleSpacing);
//                std::cout<<"Outflow end"<<std::endl;
        }

	for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {
//		cout<<p<<" "<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;	
		if(m_pParticleData->m_vBoundaryObjTypes[p]!="inflow") continue;
		if(m_pParticleData->m_vBoundaryObj[p]->UpdateInflowBoundary(m_pParticleData,m_pEOS,m_fDt,m_fInitParticleSpacing))
		{
			std::cout<<"Error: too many inflow particles."<<std::endl;
			assert(0);
			return false;
		}
	}
	fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;

	#pragma omp parallel  
	{
	
	int tid = omp_get_thread_num();
	#else
	vector<double> bX, bY, bZ, bPressure, bVx, bVy, bVz;
	vector<size_t> fIndex;
	#endif

	for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {
//		cout<<p<<" "<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;	
		if(m_pParticleData->m_vBoundaryObjTypes[p]!="solid") continue;
//		cout<<"m_pParticleData->m_vBoundaryObjTypes["<<p<<"]="<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;
		if(m_iDimension==2) {
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for(size_t index=fluidStartIndex; index<fluidEndIndex+m_pParticleData->m_iInflowNum; index++)  {
				#ifdef _OPENMP
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				0,
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				0,
				bX[tid],bY[tid],bZ[tid],bPressure[tid],bVx[tid],bVy[tid],bVz[tid]); 
				
				for(int k=0; k<num; k++) fIndex[tid].push_back(index);
				#else
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				0,
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				0,	
				bX,bY,bZ,bPressure,bVx,bVy,bVz); 
				for(int k=0; k<num; k++) fIndex.push_back(index);
				#endif
			}
		}
		else if(m_iDimension==3) {
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for(size_t index=fluidStartIndex; index<fluidEndIndex+m_pParticleData->m_iInflowNum; index++)  {
				#ifdef _OPENMP
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				m_pParticleData->m_vPositionZ[index],
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				m_pParticleData->m_vVelocityW[index],
				bX[tid],bY[tid],bZ[tid],bPressure[tid],bVx[tid],bVy[tid],bVz[tid]);
				for(int k=0; k<num; k++) fIndex[tid].push_back(index);
				#else
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				m_pParticleData->m_vPositionZ[index],
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				m_pParticleData->m_vVelocityW[index],
				bX,bY,bZ,bPressure,bVx,bVy,bVz); 
				for(int k=0; k<num; k++) fIndex.push_back(index);
				#endif
			}
		}
	}

	#ifdef _OPENMP
	}
	#endif

//	cout<<"Start to put solid boundary particles into main arrays"<<endl;
	size_t boundaryIndex = m_pParticleData->getBoundaryStartIndex()+m_pParticleData->m_iInflowNum;
//	cout<<"boundaryIndex="<<boundaryIndex<<endl;
	

	#ifdef _OPENMP
	size_t sum = 0;
	for(size_t tid=0; tid<numThreads; tid++) 
		sum += bX[tid].size();
	if(boundaryIndex+sum > m_pParticleData->m_iCapacity) {
		cout<<m_pParticleData->m_iCapacity<<endl;
		cout<<fluidStartIndex<<" "<<fluidEndIndex<<endl;
		cout<<boundaryIndex<<" "<<sum<<endl;
		cout<<"Error: Not enough memory for solid boundary particles!!!"<<endl;
		assert(0);
		return false; // not enough space -> augment data array size!
	}	
	m_vMirrorIndex.resize(sum);

	if(m_iDimension==3) {
		size_t count = 0;
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<bX[tid].size(); j++) {
				m_pParticleData->m_vPositionX[boundaryIndex] = bX[tid][j];
				m_pParticleData->m_vPositionY[boundaryIndex] = bY[tid][j];
				m_pParticleData->m_vPositionZ[boundaryIndex] = bZ[tid][j];
				m_pParticleData->m_vPressure[boundaryIndex] = bPressure[tid][j];
				m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[tid][j];
				m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[tid][j];
				m_pParticleData->m_vVelocityW[boundaryIndex] = bVz[tid][j];	

				m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[tid][j]];
                                m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[tid][j]];
                                m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[tid][j]];

				boundaryIndex++;
				m_vMirrorIndex[count++] = fIndex[tid][j]; 
			}
		}
	}
	else if(m_iDimension==2) {
		size_t count = 0;
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<bX[tid].size(); j++) {
				m_pParticleData->m_vPositionX[boundaryIndex] = bX[tid][j];
				m_pParticleData->m_vPositionY[boundaryIndex] = bY[tid][j];
				m_pParticleData->m_vPressure[boundaryIndex] = bPressure[tid][j];
				m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[tid][j];
				m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[tid][j];

                                m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[tid][j]];
                                m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[tid][j]];
                                m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[tid][j]];

				boundaryIndex++;
				m_vMirrorIndex[count++] = fIndex[tid][j];
			}
		}
	}	
	#else	
	if(boundaryIndex+bX.size() > m_pParticleData->m_iCapacity) {
		cout<<"Not enough memory for solid boundary particles!!!"<<endl;
		return false; // not enough space -> augment data array size!
	}
	m_vMirrorIndex.resize(bX.size());
	assert(bX.size()==fIndex.size());

	if(m_iDimension==3) {
		
		for(size_t j=0; j<bX.size(); j++) {
			m_pParticleData->m_vPositionX[boundaryIndex] = bX[j];
			m_pParticleData->m_vPositionY[boundaryIndex] = bY[j];
			m_pParticleData->m_vPositionZ[boundaryIndex] = bZ[j];
			m_pParticleData->m_vPressure[boundaryIndex] = bPressure[j];
			m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[j];
			m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[j];
			m_pParticleData->m_vVelocityW[boundaryIndex] = bVz[j];

                        m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[j]];
                        m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[j]];
                        m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[j]];

			boundaryIndex++;
			m_vMirrorIndex[j] = fIndex[j];
		}
		
	}
	else if(m_iDimension==2) {
		
		for(size_t j=0; j<bX.size(); j++) {
			m_pParticleData->m_vPositionX[boundaryIndex] = bX[j];
			m_pParticleData->m_vPositionY[boundaryIndex] = bY[j];	
			m_pParticleData->m_vPressure[boundaryIndex] = bPressure[j];
			m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[j];
			m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[j];
			
                        m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[j]];
                        m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[j]];
                        m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[j]];

			boundaryIndex++;
			m_vMirrorIndex[j] = fIndex[j];
		}
		
	}
	#endif
	
	m_pParticleData->m_iBoundaryNum = boundaryIndex - m_pParticleData->m_iBoundaryStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum;

	cout<<"Generated "<<m_pParticleData->m_iBoundaryNum<<" solid boundary particles"<<endl;
	m_pParticleData->setGhostStartIndex(m_pParticleData->m_iTotalNum);	

	cout<<"-------END HyperbolicLPSolver::generateSolidBoundaryByMirrorParticles()-------"<<endl;

	return true;
}


bool HyperbolicLPSolver::generatePeriodicBoundaryByMirrorParticles() {
	
	cout<<"-------HyperbolicLPSolver::generatePeriodicBoundaryByMirrorParticles()-------"<<endl;	

	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;


	#ifdef _OPENMP
	size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
	vector<vector<double>> bX(numThreads);
	vector<vector<double>> bY(numThreads);
	vector<vector<double>> bZ(numThreads);
	vector<vector<double>> bPressure(numThreads);
	vector<vector<double>> bVx(numThreads);
	vector<vector<double>> bVy(numThreads);
	vector<vector<double>> bVz(numThreads);
	vector<vector<size_t>> fIndex(numThreads);// the corresponding fluid index of a mirror particle
/*        for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {

                if(m_pParticleData->m_vBoundaryObjTypes[p]!="solid") continue;
                cout<<"m_pParticleData->m_vBoundaryObjTypes["<<p<<"]="<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;
	}*/
	#pragma omp parallel  
	{
	
	int tid = omp_get_thread_num();
	#else
	vector<double> bX, bY, bZ, bPressure, bVx, bVy, bVz;
	vector<size_t> fIndex;
	#endif

	for(size_t p=0; p<m_pParticleData->m_vBoundaryObjTypes.size(); p++) {
//		cout<<p<<" "<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;	
		if(m_pParticleData->m_vBoundaryObjTypes[p]!="periodic") continue;
//		cout<<"m_pParticleData->m_vBoundaryObjTypes["<<p<<"]="<<m_pParticleData->m_vBoundaryObjTypes[p]<<endl;
		if(m_iDimension==2) {
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for(size_t index=fluidStartIndex; index<fluidEndIndex; index++)  {
				
				#ifdef _OPENMP
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				0,
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				0,
				bX[tid],bY[tid],bZ[tid],bPressure[tid],bVx[tid],bVy[tid],bVz[tid]); 

				if(num==-1)
				{
					m_pParticleData->m_vPositionX[index]=bX[tid].back();
					m_pParticleData->m_vPositionY[index]=bY[tid].back();
					bX[tid].pop_back();
					bY[tid].pop_back();
					num = m_pParticleData->m_vBoundaryObj[p]->operator()(
                                		m_pParticleData->m_vPositionX[index],
                                		m_pParticleData->m_vPositionY[index],
                                		0,
                                		m_pParticleData->m_vPressure[index],
                                		m_pParticleData->m_vVelocityU[index],
                                		m_pParticleData->m_vVelocityV[index],
                                		0,
                                		bX[tid],bY[tid],bZ[tid],bPressure[tid],bVx[tid],bVy[tid],bVz[tid]);
				}
				assert(num>=0);				
				for(int k=0; k<num; k++) fIndex[tid].push_back(index);
				#else
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				0,
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				0,	
				bX,bY,bZ,bPressure,bVx,bVy,bVz);

				if(num==-1)
				{
					m_pParticleData->m_vPositionX[index]=bX.back();
					m_pParticleData->m_vPositionY[index]=bY.back();
					bX.pop_back();
					bY.pop_back();
                                	num = m_pParticleData->m_vBoundaryObj[p]->operator()(
                                		m_pParticleData->m_vPositionX[index],
                                		m_pParticleData->m_vPositionY[index],
                                		0,
                                		m_pParticleData->m_vPressure[index],
                                		m_pParticleData->m_vVelocityU[index],
                                		m_pParticleData->m_vVelocityV[index],
                                		0,
                                		bX,bY,bZ,bPressure,bVx,bVy,bVz);
				}
				assert(num>=0);
				for(int k=0; k<num; k++) fIndex.push_back(index);
				#endif
			}
		}
		else if(m_iDimension==3) {
//TODO: implement 3D periodic boundary condition
			printf("3D periodic boundary condition has not beem implemented yet!\n");
			assert(false);
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for(size_t index=fluidStartIndex; index<fluidEndIndex; index++)  {
				#ifdef _OPENMP
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				m_pParticleData->m_vPositionZ[index],
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				m_pParticleData->m_vVelocityW[index],
				bX[tid],bY[tid],bZ[tid],bPressure[tid],bVx[tid],bVy[tid],bVz[tid]);
				for(int k=0; k<num; k++) fIndex[tid].push_back(index);
				#else
				int num = m_pParticleData->m_vBoundaryObj[p]->operator()(
				m_pParticleData->m_vPositionX[index],
				m_pParticleData->m_vPositionY[index],
				m_pParticleData->m_vPositionZ[index],
				m_pParticleData->m_vPressure[index],
				m_pParticleData->m_vVelocityU[index],
				m_pParticleData->m_vVelocityV[index],
				m_pParticleData->m_vVelocityW[index],
				bX,bY,bZ,bPressure,bVx,bVy,bVz); 
				for(int k=0; k<num; k++) fIndex.push_back(index);
				#endif
			}
		}
	}

	#ifdef _OPENMP
	}
	#endif

//	cout<<"Start to put solid boundary particles into main arrays"<<endl;
	size_t boundaryIndex = m_pParticleData->getBoundaryStartIndex();
//	cout<<"boundaryIndex="<<boundaryIndex<<endl;
	

	#ifdef _OPENMP
	size_t sum = 0;
	for(size_t tid=0; tid<numThreads; tid++) 
		sum += bX[tid].size();
	if(boundaryIndex+sum > m_pParticleData->m_iCapacity) {
		cout<<m_pParticleData->m_iCapacity<<endl;
		cout<<fluidStartIndex<<" "<<fluidEndIndex<<endl;
		cout<<boundaryIndex<<" "<<sum<<endl;
		cout<<"Not enough memory for periodic boundary particles!!!"<<endl;
		return false; // not enough space -> augment data array size!
	}	
	m_vMirrorIndex.resize(sum);

	if(m_iDimension==3) {
		size_t count = 0;
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<bX[tid].size(); j++) {
				m_pParticleData->m_vPositionX[boundaryIndex] = bX[tid][j];
				m_pParticleData->m_vPositionY[boundaryIndex] = bY[tid][j];
				m_pParticleData->m_vPositionZ[boundaryIndex] = bZ[tid][j];
				m_pParticleData->m_vPressure[boundaryIndex] = bPressure[tid][j];
				m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[tid][j];
				m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[tid][j];
				m_pParticleData->m_vVelocityW[boundaryIndex] = bVz[tid][j];	

				m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[tid][j]];
                                m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[tid][j]];
                                m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[tid][j]];

				boundaryIndex++;
				m_vMirrorIndex[count++] = fIndex[tid][j]; 
			}
		}
	}
	else if(m_iDimension==2) {
		size_t count = 0;
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<bX[tid].size(); j++) {
				m_pParticleData->m_vPositionX[boundaryIndex] = bX[tid][j];
				m_pParticleData->m_vPositionY[boundaryIndex] = bY[tid][j];
				m_pParticleData->m_vPressure[boundaryIndex] = bPressure[tid][j];
				m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[tid][j];
				m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[tid][j];

                                m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[tid][j]];
                                m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[tid][j]];
                                m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[tid][j]];

				boundaryIndex++;
				m_vMirrorIndex[count++] = fIndex[tid][j];
			}
		}
	}	
	#else	
	if(boundaryIndex+bX.size() > m_pParticleData->m_iCapacity) {
		cout<<"Not enough memory for periodic boundary particles!!!"<<endl;
		return false; // not enough space -> augment data array size!
	}
	m_vMirrorIndex.resize(bX.size());
	assert(bX.size()==fIndex.size());

	if(m_iDimension==3) {
		
		for(size_t j=0; j<bX.size(); j++) {
			m_pParticleData->m_vPositionX[boundaryIndex] = bX[j];
			m_pParticleData->m_vPositionY[boundaryIndex] = bY[j];
			m_pParticleData->m_vPositionZ[boundaryIndex] = bZ[j];
			m_pParticleData->m_vPressure[boundaryIndex] = bPressure[j];
			m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[j];
			m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[j];
			m_pParticleData->m_vVelocityW[boundaryIndex] = bVz[j];

                        m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[j]];
                        m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[j]];
                        m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[j]];

			boundaryIndex++;
			m_vMirrorIndex[j] = fIndex[j];
		}
		
	}
	else if(m_iDimension==2) {
		
		for(size_t j=0; j<bX.size(); j++) {
			m_pParticleData->m_vPositionX[boundaryIndex] = bX[j];
			m_pParticleData->m_vPositionY[boundaryIndex] = bY[j];	
			m_pParticleData->m_vPressure[boundaryIndex] = bPressure[j];
			m_pParticleData->m_vVelocityU[boundaryIndex] = bVx[j];
			m_pParticleData->m_vVelocityV[boundaryIndex] = bVy[j];
			
                        m_pParticleData->m_vVolume[boundaryIndex] = m_pParticleData->m_vVolume[fIndex[j]];
                        m_pParticleData->m_vVolumeOld[boundaryIndex] = m_pParticleData->m_vVolumeOld[fIndex[j]];
                        m_pParticleData->m_vMass[boundaryIndex] = m_pParticleData->m_vMass[fIndex[j]];

			boundaryIndex++;
			m_vMirrorIndex[j] = fIndex[j];
		}
		
	}
	#endif
	
	m_pParticleData->m_iBoundaryNum = boundaryIndex - m_pParticleData->m_iBoundaryStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum;

	cout<<"Generated "<<m_pParticleData->m_iBoundaryNum<<" periodic boundary particles"<<endl;
	m_pParticleData->setGhostStartIndex(m_pParticleData->m_iTotalNum);	

	cout<<"-------END HyperbolicLPSolver::generatePeriodicBoundaryByMirrorParticles()-------"<<endl;

	return true;
}

bool HyperbolicLPSolver::generateGhostParticleByFillingVacancy() {
	
	cout<<"-------HyperbolicLPSolver::generateGhostParticleByFillingVacancy()-------"<<endl;

	double *x = m_pParticleData->m_vPositionX;
	double *y = m_pParticleData->m_vPositionY;
	double *z = m_pParticleData->m_vPositionZ;
	//int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
	size_t ghostIndex = m_pParticleData->getGhostStartIndex();
//	cout<<"ghostStartIndex="<<ghostIndex<<endl;
	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;

	#ifdef _OPENMP
	size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	cout<<"numThreads="<<numThreads<<endl;
//	cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<endl;
	vector<vector<size_t>> fIndex(numThreads);
	vector<vector<double>> gX(numThreads);
	vector<vector<double>> gY(numThreads);
	vector<vector<double>> gZ;
	if(m_iDimension == 3) gZ = vector<vector<double>>(numThreads);
	
	#pragma omp parallel  
	{
	
	int tid = omp_get_thread_num();
	
	#endif


//	#ifdef _OPENMP
//	#pragma omp parallel  
//	{	
//	int tid = omp_get_thread_num();		
//	#endif

	if(m_iDimension==2) {
	
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++)  {
//			cout<<index<<endl;			
//			if(!m_vFillGhost[index]) continue;
			m_vFillGhost[index]=false;
			if(m_pParticleData->m_vObjectTag[index]<0) continue; // fluid particle in contact region

			double x0 = x[index], y0 = y[index];
			size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
			
			int count[8] = {0};
			
//			cout<<"index="<<index<<endl;
//			cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
			size_t other = 0;	
			for(int i=0; i<neighbourListSize[index]; i++) {
				
				size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
				
				double x1 = x[neiIndex], y1 = y[neiIndex];
				double dx = x1-x0, dy = y1-y0;

				if(dy>0 && dy<dx)       count[0]++;
				else if(dx>0 && dy>dx)  count[1]++;
				else if(dx<0 && dy>-dx) count[2]++;
				else if(dy>0 && dy<-dx) count[3]++;
				else if(dy<0 && dy>dx)  count[4]++;
				else if(dx<0 && dy<dx)  count[5]++;
				else if(dx>0 && dy<-dx) count[6]++;
				else if(dy<0 && dy>-dx) count[7]++;
				else other++;	
			}
			for(int i=0; i<8; i++) {
//				cout<<"count["<<i<<"]="<<count[i]<<endl;
				if(count[i]<8) {
					m_vFillGhost[index] = true;
					#ifdef _OPENMP
					fillGhostParticle2D(i,count,index,gX[tid],gY[tid],fIndex[tid]);
					#else
					bool b = fillGhostParticle2D(i,count,index,&ghostIndex);
					if(!b) return false;
					#endif
				}
			}
//			cout<<"other="<<other<<endl;

		}
		
		
	}
	else if(m_iDimension==3) {
		
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
//			if(!m_vFillGhost[index]) continue;
			m_vFillGhost[index]=false;
			if(m_pParticleData->m_vObjectTag[index]<0) continue; // fluid particle in contact region

			double x0 = x[index], y0 = y[index], z0 = z[index];
			size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
			
			int count[6] = {0};
			
//			cout<<"index="<<index<<endl;
//			cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
			size_t other = 0;	
			for(int i=0; i<neighbourListSize[index]; i++) {
				
				size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];
				
				double x1 = x[neiIndex], y1 = y[neiIndex], z1 = z[neiIndex];
				double dx = x1-x0, dy = y1-y0, dz = z1-z0;
/*				if(dz>0) {
					if(dy>0 && dy<dx)       count[0]++;
					else if(dx>0 && dy>dx)  count[1]++;
					else if(dx<0 && dy>-dx) count[2]++;
					else if(dy>0 && dy<-dx) count[3]++;
					else if(dy<0 && dy>dx)  count[4]++;
					else if(dx<0 && dy<dx)  count[5]++;
					else if(dx>0 && dy<-dx) count[6]++;
					else if(dy<0 && dy>-dx) count[7]++;
					else other++;
				}
				if(dz<0) {
					if(dy>0 && dy<dx)       count[8]++;
					else if(dx>0 && dy>dx)  count[9]++;
					else if(dx<0 && dy>-dx) count[10]++;
					else if(dy>0 && dy<-dx) count[11]++;
					else if(dy<0 && dy>dx)  count[12]++;
					else if(dx<0 && dy<dx)  count[13]++;
					else if(dx>0 && dy<-dx) count[14]++;
					else if(dy<0 && dy>-dx) count[15]++;
					else other++;
				}
				else other++;*/
				if(dx>fabs(dy) && dx>fabs(dz)) count[0]++;
				else if(-dx>fabs(dy) && -dx>fabs(dz)) count[1]++;
				else if(dy>fabs(dx) && dy>fabs(dz)) count[2]++;
                                else if(-dy>fabs(dx) && -dy>fabs(dz)) count[3]++;
                                else if(dz>fabs(dx) && dz>fabs(dy)) count[4]++;
                                else if(-dz>fabs(dx) && -dz>fabs(dy)) count[5]++;
				else other++;
			}
			for(int i=0; i<6; i++) {
				//cout<<"count["<<i<<"]="<<count[i]<<endl;
				if(count[i]==0) {//=8 FOR POWDER TARGET
					m_vFillGhost[index] = true;
					#ifdef _OPENMP
					fillGhostParticle3D(i,count,index,gX[tid],gY[tid],gZ[tid],fIndex[tid]);
					#else
					bool b = fillGhostParticle3D(i,count,index,&ghostIndex);
					printf("Error: signle thread code for 3D ghost particles are old version and may produce inaccurate result.\n");assert(false);//TODO
					if(!b) return false;
					#endif
				}
			}
			//cout<<"other="<<other<<endl;

		}	
	
	
	}
	

	#ifdef _OPENMP
	}
	#endif
	
	cout<<"Start to fill ghost particles"<<endl;

	#ifdef _OPENMP
	size_t sum = 0;
	for(size_t tid=0; tid<numThreads; tid++) 
		sum += gX[tid].size();
	if(ghostIndex+sum > m_pParticleData->m_iCapacity) return false; // not enough space -> augment data array size!
	
	double *pressure = m_pParticleData->m_vPressure;
	double *velocityU = m_pParticleData->m_vVelocityU;
	double *velocityV = m_pParticleData->m_vVelocityV;
	double *velocityW;
	if(m_iDimension==3) velocityW = m_pParticleData->m_vVelocityW;
	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

	if(m_iDimension==3) {
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<gX[tid].size(); j++) {
				x[ghostIndex] = gX[tid][j];
				y[ghostIndex] = gY[tid][j];
				z[ghostIndex] = gZ[tid][j];
				size_t index = fIndex[tid][j];

				pressure[ghostIndex] = 0;
				velocityU[ghostIndex] = velocityU[index];
				velocityV[ghostIndex] = velocityV[index];
				velocityW[ghostIndex] = velocityW[index];
				localParSpacing[ghostIndex] = localParSpacing[index];
				m_pParticleData->m_vMass[ghostIndex] = 0;
				m_pParticleData->m_vVolume[ghostIndex] = 1.0e5;
                                m_pParticleData->m_vVolumeOld[ghostIndex] = 1.0e5;
				size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
				size_t sz = m_pParticleData->m_vNeighbourListSize[index];
				if(sz<m_pParticleData->m_iMaxNeighbourNum){
					m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=ghostIndex;
					m_pParticleData->m_vNeighbourListSize[index]++;
					ghostIndex++;
				}
				else{
					std::cout<<"Too many ghost neighbours for particle ("<<x[index]<<" "<<y[index]<<" "<<z[index]<<")."<<std::endl;
				}
			}
		}
	}
	else if(m_iDimension==2) {
		for(size_t tid=0; tid<numThreads; tid++) {
			for(size_t j=0; j<gX[tid].size(); j++) {
				x[ghostIndex] = gX[tid][j];
				y[ghostIndex] = gY[tid][j];
				size_t index = fIndex[tid][j];
				
				pressure[ghostIndex] = 0;
				velocityU[ghostIndex] = velocityU[index];
				velocityV[ghostIndex] = velocityV[index];	
				localParSpacing[ghostIndex] = localParSpacing[index];
                                m_pParticleData->m_vMass[ghostIndex] = 0;
                                m_pParticleData->m_vVolume[ghostIndex] = 1.0e5;
                                m_pParticleData->m_vVolumeOld[ghostIndex] = 1.0e5;
				size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
				size_t sz = m_pParticleData->m_vNeighbourListSize[index];
				if(sz<m_pParticleData->m_iMaxNeighbourNum){
					m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=ghostIndex;
					m_pParticleData->m_vNeighbourListSize[index]++;
					ghostIndex++;
				}
				else{
					std::cout<<"Too many ghost neighbours for particle ("<<x[index]<<" "<<y[index]<<" "<<z[index]<<")."<<std::endl;
				}
			}
		}
	}	
	#endif
	
	m_pParticleData->m_iGhostNum = ghostIndex - m_pParticleData->m_iGhostStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum + 
								   m_pParticleData->m_iGhostNum;

	cout<<"Generated "<<m_pParticleData->m_iGhostNum<<" ghost particles"<<endl;
        if(m_pParticleData->m_iTotalNum > m_pParticleData->m_iCapacity) {
                cout<<m_pParticleData->m_iTotalNum<<" "<<m_pParticleData->m_iCapacity<<endl;
                cout<<"Error: Not enough memory for vacuum particles!!!"<<endl;
		assert(0);
                return false; // not enough space -> augment data array size!
        }
	cout<<"-------END HyperbolicLPSolver::generateGhostParticleByFillingVacancy()-------"<<endl;
	return true;	
}

#ifdef _OPENMP
void HyperbolicLPSolver::fillGhostParticle2D(int dir, int count[], size_t index, 
vector<double>& gX, vector<double>& gY, vector<size_t>& fIndex) {
	
	double delta = m_pParticleData->m_vLocalParSpacing[index];
	double r = m_fTimesNeiSearchRadius * delta;
	int k = m_fTimesNeiSearchRadius+1;
	
	if(dir==0 || dir==3 || dir==4 || dir==7) {
		bool repeat = false;
		if(dir==0 && count[1]==0) repeat = true;
		else if(dir==3 && count[2]==0)repeat = true;
		else if(dir==4 && count[5]==0) repeat = true;
		else if(dir==7 && count[6]==0) repeat = true;
		for(int i=1; i<=k; i++) {
			int limit = (repeat==true)? i-1:i;
			for(int j=1; j<=limit; j++) {
				if(((i*i+j*j)*delta*delta <= r*r)&&((i*i+j*j)>0)) {
					double x=0, y=0;
					if(dir==0 || dir==7)
						x = m_pParticleData->m_vPositionX[index] + i*delta;
					else
						x = m_pParticleData->m_vPositionX[index] - i*delta;
					if(dir==0 || dir==3)
						y = m_pParticleData->m_vPositionY[index] + j*delta;
					else
						y = m_pParticleData->m_vPositionY[index] - j*delta;
					
					gX.push_back(x);
					gY.push_back(y);
					fIndex.push_back(index);	
				}
			}
		}
	}
	else {
		for(int j=1; j<=k; j++) {
			for(int i=1; i<=j; i++) {
				if((i*i+j*j)*delta*delta <= r*r) {
					double x=0, y=0;
					if(dir==1 || dir==6)
						x = m_pParticleData->m_vPositionX[index] + i*delta;
					else
						x = m_pParticleData->m_vPositionX[index] - i*delta;
					if(dir==1 || dir==2)
						y = m_pParticleData->m_vPositionY[index] + j*delta;
					else
						y = m_pParticleData->m_vPositionY[index] - j*delta;

					gX.push_back(x);
					gY.push_back(y);
					fIndex.push_back(index);	
				}
			}
		}		
	}
	
	//cout<<"index="<<index<<" num ghost="<<*ghostIndex-saveg<<endl;
	//cout<<"ghostIndex="<<*ghostIndex<<endl;

}
#else
bool HyperbolicLPSolver::fillGhostParticle2D(int dir, int count[], size_t index, size_t* ghostIndex) {
	//size_t saveg = *ghostIndex;
	/*
	double r = m_fNeiSearchRadius;
	double delta = m_fAvgParticleSpacing;
	int k = r/delta+1;
	*/

	double delta = m_pParticleData->m_vLocalParSpacing[index];
	double r = m_fTimesNeiSearchRadius * delta;
	int k = m_fTimesNeiSearchRadius;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
	if(dir==0 || dir==3 || dir==4 || dir==7) {
		bool repeat = false;
		if(dir==0 && count[1]==0) repeat = true;
		else if(dir==3 && count[2]==0)repeat = true;
		else if(dir==4 && count[5]==0) repeat = true;
		else if(dir==7 && count[6]==0) repeat = true;
		for(int i=1; i<=k; i++) {
			int limit = (repeat==true)? i-1:i;
			for(int j=1; j<=limit; j++) {
				if((i*i+j*j)*delta*delta <= r*r) {
					double x=0, y=0;
					if(dir==0 || dir==7)
						x = m_pParticleData->m_vPositionX[index] + i*delta;
					else
						x = m_pParticleData->m_vPositionX[index] - i*delta;
					if(dir==0 || dir==3)
						y = m_pParticleData->m_vPositionY[index] + j*delta;
					else
						y = m_pParticleData->m_vPositionY[index] - j*delta;

					if(*ghostIndex>=m_pParticleData->m_iCapacity) return false;
					
					m_pParticleData->m_vPositionX[*ghostIndex] = x;
					m_pParticleData->m_vPositionY[*ghostIndex] = y;
					
					m_pParticleData->m_vPressure[*ghostIndex] = 0;
					m_pParticleData->m_vVelocityU[*ghostIndex] = m_pParticleData->m_vVelocityU[index];
					m_pParticleData->m_vVelocityV[*ghostIndex] = m_pParticleData->m_vVelocityV[index];	
					m_pParticleData->m_vLocalParSpacing[*ghostIndex] = m_pParticleData->m_vLocalParSpacing[index];

					int sz = m_pParticleData->m_vNeighbourListSize[index];
					m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=(*ghostIndex);
					m_pParticleData->m_vNeighbourListSize[index]++;
					(*ghostIndex)++;
					
				}
			}
		}
	}
	else {
		for(int j=1; j<=k; j++) {
			for(int i=1; i<=j; i++) {
				if((i*i+j*j)*delta*delta <= r*r) {
					double x=0, y=0;
					if(dir==1 || dir==6)
						x = m_pParticleData->m_vPositionX[index] + i*delta;
					else
						x = m_pParticleData->m_vPositionX[index] - i*delta;
					if(dir==1 || dir==2)
						y = m_pParticleData->m_vPositionY[index] + j*delta;
					else
						y = m_pParticleData->m_vPositionY[index] - j*delta;
					
					if(*ghostIndex>=m_pParticleData->m_iCapacity) return false;

					m_pParticleData->m_vPositionX[*ghostIndex] = x;
					m_pParticleData->m_vPositionY[*ghostIndex] = y;
					
					m_pParticleData->m_vPressure[*ghostIndex] = 0;
					m_pParticleData->m_vVelocityU[*ghostIndex] = m_pParticleData->m_vVelocityU[index];
					m_pParticleData->m_vVelocityV[*ghostIndex] = m_pParticleData->m_vVelocityV[index];	
					m_pParticleData->m_vLocalParSpacing[*ghostIndex] = m_pParticleData->m_vLocalParSpacing[index];

					int sz = m_pParticleData->m_vNeighbourListSize[index];
					m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=(*ghostIndex);
					m_pParticleData->m_vNeighbourListSize[index]++;
					(*ghostIndex)++;
					
				}
			}
		}		
	}
	
	return true;
	//cout<<"index="<<index<<" num ghost="<<*ghostIndex-saveg<<endl;
	//cout<<"ghostIndex="<<*ghostIndex<<endl;

}
#endif



#ifdef _OPENMP
void HyperbolicLPSolver::fillGhostParticle3D(int dir, int count[], size_t index, 
vector<double>& gX, vector<double>& gY, vector<double>& gZ, vector<size_t>& fIndex) {
	
	double delta = m_pParticleData->m_vLocalParSpacing[index];
	double r = (m_fTimesNeiSearchRadius+1) * delta;
//	int k = m_fTimesNeiSearchRadius+1;
	int k=3;
	int min=0;//=8 for powder target

	if(dir==0 || dir==1){
                for(int i=1; i<=k; i++) {
                        for(int j=-i+(count[3]>=min); j<=i-(count[2]>=min); j++) {
                                for(int l=-i+(count[5]>=min); l<=i-(count[4]>=min); l++) {
                                        if((i*i+j*j+l*l)*delta*delta <= r*r) {
                                                double x=0, y=0, z=0;
                                                if(dir==0)
                                                        x = m_pParticleData->m_vPositionX[index] + i*delta;
                                                else
                                                        x = m_pParticleData->m_vPositionX[index] - i*delta;
                                                y = m_pParticleData->m_vPositionY[index] + j*delta;
                                                z = m_pParticleData->m_vPositionZ[index] + l*delta;
                                                gX.push_back(x);
                                                gY.push_back(y);
                                                gZ.push_back(z);
                                                fIndex.push_back(index);
                                        }
                                }
                        }
                }
	}

        if(dir==2 || dir==3){
                for(int i=1; i<=k; i++) {
                        for(int j=-i+(count[1]>=min); j<=i-(count[0]>=min); j++) {
                                for(int l=-i+(count[5]>=min); l<=i-(count[4]>=min); l++) {
                                        if((i*i+j*j+l*l)*delta*delta <= r*r) {
                                                double x=0, y=0, z=0;
                                                if(dir==2)
                                                        y = m_pParticleData->m_vPositionY[index] + i*delta;
                                                else    
                                                        y = m_pParticleData->m_vPositionY[index] - i*delta;
                                                x = m_pParticleData->m_vPositionX[index] + j*delta;
                                                z = m_pParticleData->m_vPositionZ[index] + l*delta;
                                                gX.push_back(x);
                                                gY.push_back(y);
                                                gZ.push_back(z);
                                                fIndex.push_back(index);
                                        }
                                }
                        }
                }
        }

        if(dir==4 || dir==5){
                for(int i=1; i<=k; i++) {
                        for(int j=-i+(count[1]>=min); j<=i-(count[0]>=min); j++) {
                                for(int l=-i+(count[3]>=min); l<=i-(count[2]>=min); l++) {
                                        if((i*i+j*j+l*l)*delta*delta <= r*r) {
                                                double x=0, y=0, z=0;
                                                if(dir==4)
                                                        z = m_pParticleData->m_vPositionZ[index] + i*delta;
                                                else    
                                                        z = m_pParticleData->m_vPositionZ[index] - i*delta;
                                                x = m_pParticleData->m_vPositionX[index] + j*delta;
                                                y = m_pParticleData->m_vPositionY[index] + l*delta;
                                                gX.push_back(x);
                                                gY.push_back(y);
                                                gZ.push_back(z);
                                                fIndex.push_back(index);
                                        }
                                }
                        }
                }
        }

/*	if(dir==0 || dir==3 || dir==4 || dir==7 || dir==8 || dir==11 || dir==12 || dir==15) {
		bool repeat = false;
		if(dir==0 && count[1]<min) repeat = true;
		else if(dir==3 && count[2]<min)   repeat = true;
		else if(dir==4 && count[5]<min)   repeat = true;
		else if(dir==7 && count[6]<min)   repeat = true;
		else if(dir==8 && count[9]<min)   repeat = true;
		else if(dir==11 && count[10]<min) repeat = true;
		else if(dir==12 && count[13]<min) repeat = true;
		else if(dir==15 && count[14]<min) repeat = true;

		bool repeat_a = false;
                if(dir==0 && count[7]<min) repeat_a = true;
                else if(dir==1 && count[2]<min)   repeat_a = true;
                else if(dir==3 && count[4]<min)   repeat_a = true;
                else if(dir==5 && count[6]<min)   repeat_a = true;
                else if(dir==8 && count[15]<min)   repeat_a = true;
                else if(dir==9 && count[10]<min)   repeat_a = true;
                else if(dir==11 && count[12]<min)   repeat_a = true;
                else if(dir==13 && count[14]<min)   repeat_a = true;
		
		for(int i=1; i<=k; i++) {
			int limit = (repeat==true)? i-1:i;
			int limit_a = (repeat_a=true)? 1:0;
			for(int j=limit_a; j<=limit; j++) {
				for(int l=0; l<=k; l++) {
					if((i*i+j*j+l*l)*delta*delta <= r*r) {
						double x=0, y=0, z=0;
						if(dir==0 || dir==7 || dir==8 || dir==15)
							x = m_pParticleData->m_vPositionX[index] + i*delta;
						else
							x = m_pParticleData->m_vPositionX[index] - i*delta;
						if(dir==0 || dir==3 || dir==8 || dir==11)
							y = m_pParticleData->m_vPositionY[index] + j*delta;
						else
							y = m_pParticleData->m_vPositionY[index] - j*delta;
						if(dir==0 || dir==3 || dir==4 || dir==7)
							z = m_pParticleData->m_vPositionZ[index] + l*delta;
						else
							z = m_pParticleData->m_vPositionZ[index] - l*delta;
						gX.push_back(x);
						gY.push_back(y);
						gZ.push_back(z);
						fIndex.push_back(index);	
					}
				}
			}
		}
	}
	else {
                bool repeat_a = false;
                if(dir==0 && count[7]<min) repeat_a = true;
                else if(dir==1 && count[2]<min)   repeat_a = true;
                else if(dir==3 && count[4]<min)   repeat_a = true;
                else if(dir==5 && count[6]<min)   repeat_a = true;
                else if(dir==8 && count[15]<min)   repeat_a = true;
                else if(dir==9 && count[10]<min)   repeat_a = true;
                else if(dir==11 && count[12]<min)   repeat_a = true;
                else if(dir==13 && count[14]<min)   repeat_a = true;

		for(int j=1; j<=k; j++) {
                        int limit_a = (repeat_a=true)? 1:0;
			for(int i=limit_a; i<=j; i++) {
				for(int l=1; l<=k; l++) {
					if((i*i+j*j+l*l)*delta*delta <= r*r) {
						double x=0, y=0, z=0;
						if(dir==1 || dir==6 || dir==9 || dir==14)
							x = m_pParticleData->m_vPositionX[index] + i*delta;
						else
							x = m_pParticleData->m_vPositionX[index] - i*delta;
						if(dir==1 || dir==2 || dir==9 || dir==10)
							y = m_pParticleData->m_vPositionY[index] + j*delta;
						else
							y = m_pParticleData->m_vPositionY[index] - j*delta;
						if(dir==1 || dir==2 || dir==5 || dir==6)
							z = m_pParticleData->m_vPositionZ[index] + l*delta;
						else
							z = m_pParticleData->m_vPositionZ[index] - l*delta;
						gX.push_back(x);
						gY.push_back(y);
						gZ.push_back(z);
						fIndex.push_back(index);	
					}
				}
			}
		}		
	}*/
	
	//cout<<"index="<<index<<" num ghost="<<*ghostIndex-saveg<<endl;
	//cout<<"ghostIndex="<<*ghostIndex<<endl;

}
#else
bool HyperbolicLPSolver::fillGhostParticle3D(int dir, int count[], size_t index, size_t* ghostIndex) {
//	cout<<"ghostIndex="<<*ghostIndex<<endl;	

	double delta = m_pParticleData->m_vLocalParSpacing[index];
	double r = m_fTimesNeiSearchRadius * delta;
//	int k = m_fTimesNeiSearchRadius+1;
	int k = 2;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
	if(dir==0 || dir==3 || dir==4 || dir==7 || dir==8 || dir==11 || dir==12 || dir==15) {
		bool repeat = false;
		if(dir==0 && count[1]==0) repeat = true;
		else if(dir==3 && count[2]==0)   repeat = true;
		else if(dir==4 && count[5]==0)   repeat = true;
		else if(dir==7 && count[6]==0)   repeat = true;
		else if(dir==8 && count[9]==0)   repeat = true;
		else if(dir==11 && count[10]==0) repeat = true;
		else if(dir==12 && count[13]==0) repeat = true;
		else if(dir==15 && count[14]==0) repeat = true;
		for(int i=1; i<=k; i++) {
			int limit = (repeat==true)? i-1:i;
			for(int j=1; j<=limit; j++) {
				for(int l=1; l<=k; l++) {
					if((i*i+j*j+l*l)*delta*delta <= r*r) {
						double x=0, y=0, z=0;
						if(dir==0 || dir==7 || dir==8 || dir==15)
							x = m_pParticleData->m_vPositionX[index] + i*delta;
						else
							x = m_pParticleData->m_vPositionX[index] - i*delta;
						if(dir==0 || dir==3 || dir==8 || dir==11)
							y = m_pParticleData->m_vPositionY[index] + j*delta;
						else
							y = m_pParticleData->m_vPositionY[index] - j*delta;
						if(dir==0 || dir==3 || dir==4 || dir==7)
							z = m_pParticleData->m_vPositionZ[index] + l*delta;
						else
							z = m_pParticleData->m_vPositionZ[index] - l*delta;

						if(*ghostIndex>=m_pParticleData->m_iCapacity) return false;
						
						m_pParticleData->m_vPositionX[*ghostIndex] = x;
						m_pParticleData->m_vPositionY[*ghostIndex] = y;
						m_pParticleData->m_vPositionZ[*ghostIndex] = z;

						m_pParticleData->m_vPressure[*ghostIndex] = m_pParticleData->m_vPressure[index];
						m_pParticleData->m_vVelocityU[*ghostIndex] = m_pParticleData->m_vVelocityU[index];
						m_pParticleData->m_vVelocityV[*ghostIndex] = m_pParticleData->m_vVelocityV[index];	
						m_pParticleData->m_vVelocityW[*ghostIndex] = m_pParticleData->m_vVelocityW[index];
						m_pParticleData->m_vLocalParSpacing[*ghostIndex] = m_pParticleData->m_vLocalParSpacing[index];

						int sz = m_pParticleData->m_vNeighbourListSize[index];
						m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=(*ghostIndex);
						m_pParticleData->m_vNeighbourListSize[index]++;
						(*ghostIndex)++;
					}	
				}
			}
		}
	}
	else {
		for(int j=1; j<=k; j++) {
			for(int i=1; i<=j; i++) {
				for(int l=1; l<=k; l++) {
					if((i*i+j*j+l*l)*delta*delta <= r*r) {
						double x=0, y=0, z=0;
						if(dir==1 || dir==6 || dir==9 || dir==14)
							x = m_pParticleData->m_vPositionX[index] + i*delta;
						else
							x = m_pParticleData->m_vPositionX[index] - i*delta;
						if(dir==1 || dir==2 || dir==9 || dir==10)
							y = m_pParticleData->m_vPositionY[index] + j*delta;
						else
							y = m_pParticleData->m_vPositionY[index] - j*delta;
						if(dir==1 || dir==2 || dir==5 || dir==6)
							z = m_pParticleData->m_vPositionZ[index] + l*delta;
						else
							z = m_pParticleData->m_vPositionZ[index] - l*delta;
					
						if(*ghostIndex>=m_pParticleData->m_iCapacity) return false;

						m_pParticleData->m_vPositionX[*ghostIndex] = x;
						m_pParticleData->m_vPositionY[*ghostIndex] = y;
						m_pParticleData->m_vPositionZ[*ghostIndex] = z;

						m_pParticleData->m_vPressure[*ghostIndex] = 0;
						m_pParticleData->m_vVelocityU[*ghostIndex] = m_pParticleData->m_vVelocityU[index];
						m_pParticleData->m_vVelocityV[*ghostIndex] = m_pParticleData->m_vVelocityV[index];	
						m_pParticleData->m_vVelocityW[*ghostIndex] = m_pParticleData->m_vVelocityW[index];
						m_pParticleData->m_vLocalParSpacing[*ghostIndex] = m_pParticleData->m_vLocalParSpacing[index];

						int sz = m_pParticleData->m_vNeighbourListSize[index];
						m_pParticleData->m_vNeighbourList[neiListStartIndex+sz]=(*ghostIndex);
						m_pParticleData->m_vNeighbourListSize[index]++;
						(*ghostIndex)++;
					}	
				}
			}
		}		
	}
	
	return true;
	//cout<<"index="<<index<<" num ghost="<<*ghostIndex-saveg<<endl;
	//cout<<"ghostIndex="<<*ghostIndex<<endl;

}
#endif







//
void HyperbolicLPSolver::changeFluidObjectTagIfInContact(int index, size_t numNeiFound, const double* neiListDist) {
	
	if(numNeiFound==0) return;	

	int* objectTag = m_pParticleData->m_vObjectTag;
//	const double* localParSpacing = m_pParticleData->m_vLocalParSpacing;
	
	// tag < 0  means have already contacted with other fluid object	
	if(objectTag[index] < 0) return;

	// tag == 0 means this is not a fluid particle (This function should not be called by non-fluid particle!)
	assert(objectTag[index] != 0);

	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	 	
	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	double len = 2. * m_pParticleData->m_vLocalParSpacing[index];
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
	
		//if(neiDist > m_fContactLength) continue; 
		if(neiListDist[k] > len) break;

		if(objectTag[neiIndex] == 0) continue; // neighbour is not a fluid particle
		if(objectTag[neiIndex] == objectTag[index] || 
		   objectTag[neiIndex] == -objectTag[index]) // neighbour is the same fluid object
			continue; 	

		objectTag[index] = -objectTag[index];
		break; 	
	}
}

 
void HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNeiFromSameObject(int index, size_t numNeiFound) {
	
	if(numNeiFound==0) return;
	
	int* objectTag = m_pParticleData->m_vObjectTag;
		
	if(objectTag[index] < 0) return; // already in contact -> can take all neighbours
	
	assert(objectTag[index] != 0); // should be a fluid particle

	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	int count = 0;
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
		// only admit neighbours that are 1: non-fluid 2: fluid from the same fluid object
		if(objectTag[neiIndex]==0 || 
		   objectTag[neiIndex] == objectTag[index] ||
		   objectTag[neiIndex] == -objectTag[index]) {  
			neighbourList[neiListStartIndex+count] = neiIndex;
			count++;
		}	
	}	
	// modify neighbourListSize
	neighbourListSize[index] = count;
}

void HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNei(int index, size_t numNeiFound) {
	int* objectTag = m_pParticleData->m_vObjectTag;

	if(objectTag[index] != 0) return; // This function should be called by non-fluid particles
	
	int *neighbourList = m_pParticleData->m_vNeighbourList;	
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

	size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;	
	int count = 0;
	for(size_t k=0; k<numNeiFound; k++) {
		int neiIndex = neighbourList[neiListStartIndex+k];
		// only admit neighbours that are fluid particles
		if(objectTag[neiIndex] != 0) {  
			neighbourList[neiListStartIndex+count] = neiIndex;
			count++;
		}	
	}
//	cout<<"-------HyperbolicLPSolver::changeNeighbourhoodToIncludeOnlyFluidNei()-------"<<endl;
//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
//	cout<<"fluid neighbour count ="<<count<<endl;
//	for(int k=0; k<count; k++) {
//		int neiIndex = neighbourList[neiListStartIndex+k];
//		cout<<"neiIndex="<<neiIndex<<endl;	
//	}
//	cout<<"----------------------------------------------------------------------------"<<endl;
	
	// modify neighbourListSize
	neighbourListSize[index] = count;


}


// a1>a0  1 2 
// a1<a0  3 4
// b1>=b0 1 3 
// b1<b0  2 4 
/*void setListInOneDir2D(size_t neiIndex, double a0, double a1, double b0, double b1, 
					   int* arr1, size_t& n1, int* arr2, size_t& n2,
					   int* arr3, size_t& n3, int* arr4, size_t& n4) {
	
	if(a1 > a0) {
		//printf("a1 > a0\n");
		//printf("a1=%.16g\n",a1);
		//printf("a0=%.16g\n",a0);
		if(b1 >= b0) arr1[n1++]=neiIndex;
		else arr2[n2++]=neiIndex;
	}
	else if(a1 < a0) { 
		//printf("a1 < a0\n");
		//printf("a1=%.16g\n",a1);
		//printf("a0=%.16g\n",a0);
		if(b1 >= b0) arr3[n3++]=neiIndex;
		else arr4[n4++]=neiIndex;
	}

}*/

void setListInOneDir2D(size_t neiIndex, double a0, double a1, double b0, double b1,
                                           int* arr1, size_t& n1, int* arr2, size_t& n2,
                                           int* arr3, size_t& n3, int* arr4, size_t& n4) {

if((5*fabs(a1-a0))>=fabs(b1-b0)){
        if(a1 > a0) {
                //printf("a1 > a0\n");
                //printf("a1=%.16g\n",a1);
                //printf("a0=%.16g\n",a0);
                if(b1 >= b0) arr1[n1++]=neiIndex;
                else arr2[n2++]=neiIndex;
        }
        else if(a1 < a0) {
                //printf("a1 < a0\n");
                //printf("a1=%.16g\n",a1);
                //printf("a0=%.16g\n",a0);
                if(b1 >= b0) arr3[n3++]=neiIndex;
                else arr4[n4++]=neiIndex;
        }
}
}

// a1>a0  1 2 3 4
// a1<a0  5 6 7 8
// b1>=b0 1 2 5 6
// b1<b0  3 4 7 8
// c1>=c0 1 3 5 7
// c1<c0  2 4 6 8
void setListInOneDir3D(size_t neiIndex, double a0, double a1, double b0, double b1, double c0, double c1,
					   int* arr1, size_t& n1, int* arr2, size_t& n2,
					   int* arr3, size_t& n3, int* arr4, size_t& n4,
					   int* arr5, size_t& n5, int* arr6, size_t& n6,
					   int* arr7, size_t& n7, int* arr8, size_t& n8) {
if(((5*fabs(a1-a0))>=fabs(b1-b0))&&((5*fabs(a1-a0))>=fabs(c1-c0))){
	
	if(a1 > a0) { 
		if(b1 >= b0) {
			if(c1 >= c0) arr1[n1++]=neiIndex;
			else arr2[n2++]=neiIndex;
		}
		else {
			if(c1 >= c0) arr3[n3++]=neiIndex;
			else arr4[n4++]=neiIndex;
		}
	
	}
	else if(a1 < a0) { 
		if(b1 >= b0) {
			if(c1 >= c0) arr5[n5++]=neiIndex;
			else arr6[n6++]=neiIndex;
		}
		else {
			if(c1 >= c0) arr7[n7++]=neiIndex;
			else arr8[n8++]=neiIndex;
		}	
	}
}
}


// helper function of HyperbolicLPSolver::setUpwindNeighbourList()
void setListInOneDir2D(size_t index, size_t maxNeiNumInOneDir, 
					   const int* list1, size_t sz1, const int* list2, size_t sz2,
					   int* upwindNeighbourList, int* upwindNeighbourListSize) { // output
	
	size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
	size_t numInOneDir = 0, n1 = 0, n2 = 0;

	while(numInOneDir < maxNeiNumInOneDir) {
		if(n1 < sz1) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n2 < sz2) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
			numInOneDir++;
		}
		if(n1==sz1 && n2==sz2) break;
	}
	upwindNeighbourListSize[index] = numInOneDir;
}

void setListInOneDir2D(size_t index, size_t maxNeiNumInOneDir,
                                           const int* list1, size_t sz1, const int* list2, size_t sz2,
                                           int* upwindNeighbourList, int* upwindNeighbourListSize, const double* positionX, const double* positionY) { // output

        size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
        size_t numInOneDir = 0, n1 = 0, n2 = 0;

//        const double *positionX = m_pParticleData->m_vPositionX;
//        const double *positionY = m_pParticleData->m_vPositionY;
//        const double *positionZ = m_pParticleData->m_vPositionZ;
        int nei1,nei2;
        double d1,d2;
        while((numInOneDir < (maxNeiNumInOneDir-1))&&(n1<sz1)&&(n2<sz2)) {
                nei1=list1[n1];
                nei2=list2[n2];
                d1=(positionX[nei1]-positionX[index])*(positionX[nei1]-positionX[index])+(positionY[nei1]-positionY[index])*(positionY[nei1]-positionY[index]);
                d2=(positionX[nei2]-positionX[index])*(positionX[nei2]-positionX[index])+(positionY[nei2]-positionY[index])*(positionY[nei2]-positionY[index]);
                if (d1<d2) {
                        n1=n1+1;
                        n2=n2+1;
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = nei1;
                        numInOneDir++;
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = nei2;
                        numInOneDir++;
                }
                else {
                        n1=n1+1;
                        n2=n2+1;
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = nei2;
                        numInOneDir++;
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = nei1;
                        numInOneDir++;
                }
//	        if(n1==sz1 || n2==sz2) break;
        }

        while(numInOneDir < maxNeiNumInOneDir) {
                if(n1 < sz1) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
                        numInOneDir++;
                }
                if(numInOneDir >= maxNeiNumInOneDir) break;
                if(n2 < sz2) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
                        numInOneDir++;
                }
                if(n1==sz1 && n2==sz2) break;
        }
        upwindNeighbourListSize[index] = numInOneDir;
}


void setListInFourDir2D(size_t index, size_t neiIndex, size_t maxNeiNumInOneDir, double dx, double dy,
                                                                int* neighbourListRight, int* neighbourListRightSize, int* neighbourListLeft, int* neighbourListLeftSize, int* neighbourListNorth, int* neighbourListNorthSize, int* neighbourListSouth, int* neighbourListSouthSize)//output
{
	size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;

	if((fabs(3*dx)>fabs(dy))&&(dx>0))//right
	{
		neighbourListRight[neiListInOneDirStartIndex+neighbourListRightSize[index]] = neiIndex;
		neighbourListRightSize[index]++;
	}
        if((fabs(3*dx)>fabs(dy))&&(dx<0))//left
        {
                neighbourListLeft[neiListInOneDirStartIndex+neighbourListLeftSize[index]] = neiIndex;
                neighbourListLeftSize[index]++;
        }
        if((fabs(dx)<fabs(3*dy))&&(dy>0))//north
        {
                neighbourListNorth[neiListInOneDirStartIndex+neighbourListNorthSize[index]] = neiIndex;
                neighbourListNorthSize[index]++;
        }
        if((fabs(dx)<fabs(3*dy))&&(dy<0))//south
        {
                neighbourListSouth[neiListInOneDirStartIndex+neighbourListSouthSize[index]] = neiIndex;
                neighbourListSouthSize[index]++;
        }
}

// helper function of HyperbolicLPSolver::setUpwindNeighbourList()
void setListInOneDir3D(size_t index, size_t maxNeiNumInOneDir, 
					   const int* list1, size_t sz1, const int* list2, size_t sz2,
					   const int* list3, size_t sz3, const int* list4, size_t sz4,
					   int* upwindNeighbourList, int* upwindNeighbourListSize) { // output
	
	size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
	size_t numInOneDir = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	while(numInOneDir < maxNeiNumInOneDir) {
		if(n1 < sz1) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n2 < sz2) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n3 < sz3) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list3[n3++];
			numInOneDir++;
		}
		if(numInOneDir >= maxNeiNumInOneDir) break;
		if(n4 < sz4) {
			upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list4[n4++];
			numInOneDir++;
		}	
		
		if(n1==sz1 && n2==sz2 && n3==sz3 && n4==sz4) break;
	}
	upwindNeighbourListSize[index] = numInOneDir;
}

bool ascByDist(const pair<int,double>& l, const pair<int,double>& r ) {
        return l.second < r.second;
}

void setListInOneDir3D(size_t index, size_t maxNeiNumInOneDir,
                                           const int* list1, size_t sz1, const int* list2, size_t sz2,
                                           const int* list3, size_t sz3, const int* list4, size_t sz4,
                                           int* upwindNeighbourList, int* upwindNeighbourListSize, const double* positionX, const double* positionY, const double* positionZ) { // output

        size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
        size_t numInOneDir = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	double dx, dy, dz, dis;
        vector<pair<int,double>> index_dis(4);
	while((numInOneDir < (maxNeiNumInOneDir-3)) && (n1<sz1) && (n2<sz2) && (n3<sz3) && (n4<sz4)) {
		dx=positionX[list1[n1]]-positionX[index];
                dy=positionY[list1[n1]]-positionY[index];
                dz=positionZ[list1[n1]]-positionZ[index];
		dis=dx*dx+dy*dy+dz*dz;
		index_dis[0]={list1[n1],dis};
                dx=positionX[list2[n2]]-positionX[index];
                dy=positionY[list2[n2]]-positionY[index];
                dz=positionZ[list2[n2]]-positionZ[index];
                dis=dx*dx+dy*dy+dz*dz;
                index_dis[1]={list2[n2],dis};
                dx=positionX[list3[n3]]-positionX[index];
                dy=positionY[list3[n3]]-positionY[index];
                dz=positionZ[list3[n3]]-positionZ[index];
                dis=dx*dx+dy*dy+dz*dz;
                index_dis[2]={list3[n3],dis};
                dx=positionX[list4[n4]]-positionX[index];
                dy=positionY[list4[n4]]-positionY[index];
                dz=positionZ[list4[n4]]-positionZ[index];
                dis=dx*dx+dy*dy+dz*dz;
                index_dis[3]={list4[n4],dis};

                std::sort(index_dis.begin(), index_dis.end(), ascByDist);
		for(int add=0; add<4; add++)
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir+add] = index_dis[add].first;
		n1++;
		n2++;
		n3++;
		n4++;
		numInOneDir=numInOneDir+4;
	}
        while(numInOneDir < maxNeiNumInOneDir) {
                if(n1 < sz1) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list1[n1++];
                        numInOneDir++;
                }
                if(numInOneDir >= maxNeiNumInOneDir) break;
                if(n2 < sz2) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list2[n2++];
                        numInOneDir++;
                }
                if(numInOneDir >= maxNeiNumInOneDir) break;
                if(n3 < sz3) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list3[n3++];
                        numInOneDir++;
                }
                if(numInOneDir >= maxNeiNumInOneDir) break;
                if(n4 < sz4) {
                        upwindNeighbourList[neiListInOneDirStartIndex+numInOneDir] = list4[n4++];
                        numInOneDir++;
                }

                if(n1==sz1 && n2==sz2 && n3==sz3 && n4==sz4) break;
        }
        upwindNeighbourListSize[index] = numInOneDir;
}

void reorderNeighbour2D(size_t index, size_t neiListStartIndex, int* neighbourList, int* neighbourListSize, const double* positionX, const double* positionY) {
	vector<int> rnlist(neighbourListSize[index]),rslist(neighbourListSize[index]),lnlist(neighbourListSize[index]),lslist(neighbourListSize[index]);
        vector<pair<int,double>> index_dis(4);
	int rnsize=0,rssize=0,lnsize=0,lssize=0;
	for(size_t i=neiListStartIndex;i<neiListStartIndex+neighbourListSize[index];i++)
	{
		int k=neighbourList[i];
		double dx=positionX[k]-positionX[index];
		double dy=positionY[k]-positionY[index];
		if((dx+0.001*dy)>0&&(dy-0.001*dx>0))	rnlist[rnsize++]=k;
                if((dx+0.001*dy)<0&&(dy-0.001*dx>0))    lnlist[lnsize++]=k;
                if((dx+0.001*dy)>0&&(dy-0.001*dx<0))    rslist[rssize++]=k;
                if((dx+0.001*dy)<0&&(dy-0.001*dx<0))    lslist[lssize++]=k;
	}
	int rnindex=0,rsindex=0,lnindex=0,lsindex=0;
	int neighbourIndex=0;
	double dx,dy,dis;

        while(rnindex<rnsize && lnindex<lnsize && rsindex<rssize && lsindex<lssize) {
                dx=positionX[rnlist[rnindex]]-positionX[index];
                dy=positionY[rnlist[rnindex]]-positionY[index];
                dis=dx*dx+dy*dy;
                index_dis[0]={rnlist[rnindex++],dis};
                dx=positionX[rslist[rsindex]]-positionX[index];
                dy=positionY[rslist[rsindex]]-positionY[index];
                dis=dx*dx+dy*dy;
                index_dis[1]={rslist[rsindex++],dis};
                dx=positionX[lnlist[lnindex]]-positionX[index];
                dy=positionY[lnlist[lnindex]]-positionY[index];
                dis=dx*dx+dy*dy;
                index_dis[2]={lnlist[lnindex++],dis};
                dx=positionX[lslist[lsindex]]-positionX[index];
                dy=positionY[lslist[lsindex]]-positionY[index];
                dis=dx*dx+dy*dy;
                index_dis[3]={lslist[lsindex++],dis};

                std::sort(index_dis.begin(), index_dis.end(), ascByDist);
                for(int add=0; add<4; add++)
                        neighbourList[neiListStartIndex+neighbourIndex+add] = index_dis[add].first;
                neighbourIndex=neighbourIndex+4;
        }
        while(rnindex<rnsize || lsindex<lssize || lnindex<lnsize || rsindex<rssize) {
                if(rnindex < rnsize) {
                        neighbourList[neiListStartIndex+neighbourIndex] = rnlist[rnindex++];
                        neighbourIndex++;
                }
                if(lsindex < lssize) {
                        neighbourList[neiListStartIndex+neighbourIndex] = lslist[lsindex++];
                        neighbourIndex++;
                }
                if(lnindex < lnsize) {
                        neighbourList[neiListStartIndex+neighbourIndex] = lnlist[lnindex++];
                        neighbourIndex++;
                }
                if(rsindex < rssize) {
                        neighbourList[neiListStartIndex+neighbourIndex] = rslist[rsindex++];
                        neighbourIndex++;
                }
        }

}


// only fluid particles have upwind neighbour list
void HyperbolicLPSolver::setUpwindNeighbourList() {
	
	cout<<"-------HyperbolicLPSolver::setUpwindNeighbourList()-------"<<endl;
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;		
	const double *positionZ = m_pParticleData->m_vPositionZ;
	if(m_iDimension==2) positionZ = nullptr;
 
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListRight = m_pParticleData->m_vNeighbourListRight;
	int *neighbourListLeft = m_pParticleData->m_vNeighbourListLeft;
	int *neighbourListNorth = m_pParticleData->m_vNeighbourListNorth;
	int *neighbourListSouth = m_pParticleData->m_vNeighbourListSouth;
	int *neighbourListUp, *neighbourListDown;
	if(m_iDimension==3) {
		neighbourListUp = m_pParticleData->m_vNeighbourListUp;
		neighbourListDown = m_pParticleData->m_vNeighbourListDown;
	}	
	
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
	int *neighbourListRightSize = m_pParticleData->m_vNeighbourListRightSize;
	int *neighbourListLeftSize = m_pParticleData->m_vNeighbourListLeftSize;
	int *neighbourListNorthSize = m_pParticleData->m_vNeighbourListNorthSize;
	int *neighbourListSouthSize = m_pParticleData->m_vNeighbourListSouthSize;
	int *neighbourListUpSize, *neighbourListDownSize;
	if(m_iDimension==3) {
		neighbourListUpSize = m_pParticleData->m_vNeighbourListUpSize;
		neighbourListDownSize = m_pParticleData->m_vNeighbourListDownSize;
	}
		
	// get the maximum size of neighbour lists	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
	size_t maxNeiNumInOneDir = m_pParticleData->m_iMaxNeighbourNumInOneDir;

	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	
	#ifdef _OPENMP
	#pragma omp parallel
	{
	#endif
	if(m_iDimension==2) {
		
		int rn[maxNeiNum]; int rs[maxNeiNum]; 
		int ln[maxNeiNum]; int ls[maxNeiNum]; 
		int nr[maxNeiNum]; int nl[maxNeiNum]; 
		int sr[maxNeiNum]; int sl[maxNeiNum];

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		size_t rnSz = 0, rsSz = 0, lnSz = 0, lsSz = 0;
		size_t nrSz = 0, nlSz = 0, srSz = 0, slSz = 0;
		
		size_t neiListStartIndex = index*maxNeiNum;
		size_t neiListEndIndex = neiListStartIndex + neighbourListSize[index];

		double x0 = positionX[index], y0 = positionY[index];
//Sort the neighbours by distance for free boundary. In other cases the neighbours are already sorted.
		if(m_iFreeBoundary) {
			if(m_vFillGhost[index]) {
				vector<pair<int,double>> index_dis(neighbourListSize[index]);
				for(int i=0; i<neighbourListSize[index]; i++) {
					size_t  neiIndex = neighbourList[neiListStartIndex+i];
					double dx = positionX[neiIndex] - x0, dy = positionY[neiIndex] - y0;
					double dis = dx*dx+dy*dy;
					index_dis[i]={neiIndex,dis};
				}
				std::sort(index_dis.begin(), index_dis.end(), ascByDist);
				for(size_t i=0; i<index_dis.size(); i++) {
					neighbourList[neiListStartIndex+i] = index_dis[i].first;
				}	
			}
		}
//Put neighbours in each direction in a list
		for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {	
			size_t  neiIndex = neighbourList[i];
			double x1 = positionX[neiIndex], y1 = positionY[neiIndex];
			// right or left
			setListInOneDir2D(neiIndex, x0, x1, y0, y1,
							  rn, rnSz, rs, rsSz, ln, lnSz, ls, lsSz);
			
			// north or south
			setListInOneDir2D(neiIndex, y0, y1, x0, x1,
							  nr, nrSz, nl, nlSz, sr, srSz, sl, slSz);

		}

//Set up upwind stencil
                setListInOneDir2D(index, maxNeiNumInOneDir, rn, rnSz, rs, rsSz,
                                                  neighbourListRight, neighbourListRightSize, positionX, positionY);

                setListInOneDir2D(index, maxNeiNumInOneDir, ln, lnSz, ls, lsSz,
                                                  neighbourListLeft, neighbourListLeftSize, positionX, positionY);

                setListInOneDir2D(index, maxNeiNumInOneDir, nr, nrSz, nl, nlSz,
                                                  neighbourListNorth, neighbourListNorthSize, positionX, positionY);

                setListInOneDir2D(index, maxNeiNumInOneDir, sr, srSz, sl, slSz,
                                                  neighbourListSouth, neighbourListSouthSize, positionX, positionY);

//Set up central stencil, reorder neighbour list to balance neighbours from each quedrant begin
		reorderNeighbour2D(index, neiListStartIndex, neighbourList,neighbourListSize,positionX,positionY);

//Plot neighbourhood of one particle
/*		if(index==16583) //norandom: 16583 random: 16585
		{
	                for(size_t it=fluidStartIndex; it<fluidEndIndex; it++)
				m_pParticleData->m_vNeighOfParticle[it]=0;

			m_pParticleData->m_vNeighOfParticle[index]=2;
			for(size_t neiindex=neiListStartIndex;neiindex<neiListStartIndex+36;neiindex++)
				m_pParticleData->m_vNeighOfParticle[neighbourList[neiindex]]=1;
		}*/
	}

}
else if(m_iDimension==3) {

	int rnu[maxNeiNum]; int rsu[maxNeiNum]; int rnd[maxNeiNum]; int rsd[maxNeiNum]; // right
	int lnu[maxNeiNum]; int lsu[maxNeiNum]; int lnd[maxNeiNum]; int lsd[maxNeiNum]; // left
	int nru[maxNeiNum]; int nlu[maxNeiNum]; int nrd[maxNeiNum]; int nld[maxNeiNum]; // north
	int sru[maxNeiNum]; int slu[maxNeiNum]; int srd[maxNeiNum]; int sld[maxNeiNum]; // south
	int urn[maxNeiNum]; int uln[maxNeiNum]; int urs[maxNeiNum]; int uls[maxNeiNum]; // up
	int drn[maxNeiNum]; int dln[maxNeiNum]; int drs[maxNeiNum]; int dls[maxNeiNum]; // down 
	
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		size_t rnuSz = 0, rsuSz = 0, rndSz = 0, rsdSz = 0;
		size_t lnuSz = 0, lsuSz = 0, lndSz = 0, lsdSz = 0;
		size_t nruSz = 0, nluSz = 0, nrdSz = 0, nldSz = 0;
		size_t sruSz = 0, sluSz = 0, srdSz = 0, sldSz = 0;
		size_t urnSz = 0, ulnSz = 0, ursSz = 0, ulsSz = 0; 
		size_t drnSz = 0, dlnSz = 0, drsSz = 0, dlsSz = 0; 
		
		size_t neiListStartIndex = index*maxNeiNum;
		size_t neiListEndIndex = neiListStartIndex + neighbourListSize[index];
		
		double x0 = positionX[index], y0 = positionY[index], z0 = positionZ[index];
	
		if(m_iFreeBoundary) {
			if(m_vFillGhost[index]) {
				vector<pair<int,double>> index_dis(neighbourListSize[index]);
				for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {
					size_t  neiIndex = neighbourList[i];
					double dx = positionX[neiIndex] - x0, dy = positionY[neiIndex] - y0, dz = positionZ[neiIndex] - z0;
					double dis = dx*dx+dy*dy+dz*dz;
					index_dis[i-neiListStartIndex]={neiIndex,dis};
				}
				std::sort(index_dis.begin(), index_dis.end(),ascByDist);
				for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {
					neighbourList[i] = index_dis[i-neiListStartIndex].first;
				}
			}
		}	
		
		for(size_t i=neiListStartIndex; i<neiListEndIndex; i++) {	
			size_t  neiIndex = neighbourList[i];
			double x1 = positionX[neiIndex], y1 = positionY[neiIndex], z1 = positionZ[neiIndex];
			// right or left
			setListInOneDir3D(neiIndex, x0, x1, y0, y1, z0, z1,
							  rnu, rnuSz, rnd, rndSz, rsu, rsuSz, rsd, rsdSz,
							  lnu, lnuSz, lnd, lndSz, lsu, lsuSz, lsd, lsdSz);
			// north or south
			setListInOneDir3D(neiIndex, y0, y1, x0, x1, z0, z1,
							  nru, nruSz, nrd, nrdSz, nlu, nluSz, nld, nldSz,
							  sru, sruSz, srd, srdSz, slu, sluSz, sld, sldSz);
			// up or down
			setListInOneDir3D(neiIndex, z0, z1, x0, x1, y0, y1,
							  urn, urnSz, urs, ursSz, uln, ulnSz, uls, ulsSz,	
							  drn, drnSz, drs, drsSz, dln, dlnSz, dls, dlsSz);	
		}
	

		//DEBUG INFO
		//if(index%100==0 && index<=1000) {
		//	cout<<"-------HyperbolicLPSolver::setUpwindNeighbourList()-------"<<endl;
		//	cout<<"index="<<index<<endl;
		//	cout<<"rnuSz="<<rnuSz<<endl;
		//	cout<<"rsuSz="<<rsuSz<<endl;
		//	cout<<"rndSz="<<rndSz<<endl;
		//	cout<<"rsdSz="<<rsdSz<<endl;	
		//	cout<<"right size = "<<rnuSz+rsuSz+rndSz+rsdSz<<endl;
		//	cout<<"lnuSz="<<lnuSz<<endl;
		//	cout<<"lsuSz="<<lsuSz<<endl;
		//	cout<<"lndSz="<<lndSz<<endl;
		//	cout<<"lsdSz="<<lsdSz<<endl;
		//	cout<<"left size = "<<lnuSz+lsuSz+lndSz+lsdSz<<endl;
		//	cout<<"nruSz="<<nruSz<<endl;
		//	cout<<"nluSz="<<nluSz<<endl;
		//	cout<<"nrdSz="<<nrdSz<<endl;
		//	cout<<"nldSz="<<nldSz<<endl;
		//	cout<<"north size = "<<nruSz+nluSz+nrdSz+nldSz<<endl;
		//	cout<<"sruSz="<<sruSz<<endl;
		//	cout<<"sluSz="<<sluSz<<endl;
		//	cout<<"srdSz="<<srdSz<<endl;
		//	cout<<"sldSz="<<sldSz<<endl;
		//	cout<<"south size = "<<sruSz+sluSz+srdSz+sldSz<<endl;
		//	cout<<"urnSz="<<urnSz<<endl;
		//	cout<<"ulnSz="<<ulnSz<<endl;
		//	cout<<"ursSz="<<ursSz<<endl;
		//	cout<<"ulsSz="<<ulsSz<<endl;
		//	cout<<"up size = "<<urnSz+ulnSz+ursSz+ulsSz<<endl;
		//	cout<<"drnSz="<<drnSz<<endl;
		//	cout<<"dlnSz="<<dlnSz<<endl;
		//	cout<<"drsSz="<<drsSz<<endl;
		//	cout<<"dlsSz="<<dlsSz<<endl;
		//	cout<<"down size = "<<drnSz+dlnSz+drsSz+dlsSz<<endl;
		//	cout<<"----------------------------------------------------------"<<endl;
		//}

		// call helper function to set the lists
		//right
		setListInOneDir3D(index, maxNeiNumInOneDir, rnu, rnuSz, rsu, rsuSz, rnd, rndSz, rsd, rsdSz,
						  neighbourListRight, neighbourListRightSize, positionX, positionY, positionZ); // output
		// left
		setListInOneDir3D(index, maxNeiNumInOneDir, lnu, lnuSz, lsu, lsuSz, lnd, lndSz, lsd, lsdSz,
						  neighbourListLeft, neighbourListLeftSize, positionX, positionY, positionZ); // output
		//north
		setListInOneDir3D(index, maxNeiNumInOneDir, nru, nruSz, nlu, nluSz, nrd, nrdSz, nld, nldSz,
						  neighbourListNorth, neighbourListNorthSize, positionX, positionY, positionZ); // output
		//south
		setListInOneDir3D(index, maxNeiNumInOneDir, sru, sruSz, slu, sluSz, srd, srdSz, sld, sldSz,
						  neighbourListSouth, neighbourListSouthSize, positionX, positionY, positionZ); // output
		//up
		setListInOneDir3D(index, maxNeiNumInOneDir, urn, urnSz, uln, ulnSz, urs, ursSz, uls, ulsSz,
						  neighbourListUp, neighbourListUpSize, positionX, positionY, positionZ); // output
		//down
		setListInOneDir3D(index, maxNeiNumInOneDir, drn, drnSz, dln, dlnSz, drs, drsSz, dls, dlsSz,
						  neighbourListDown, neighbourListDownSize, positionX, positionY, positionZ); // output
//reordering neighbour list to balance neighbours from each quedrant begin

/*		double alpha=0.001,beta=0.0001,gamma=0.00001;
		double A[9];
		A[0]=cos(alpha)*cos(gamma)-sin(alpha)*sin(beta)*sin(gamma);
		A[1]=-sin(alpha)*cos(beta);
		A[2]=-cos(alpha)*sin(gamma)-sin(alpha)*sin(beta)*cos(gamma);
		A[3]=cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma);
		A[4]=cos(alpha)*cos(beta);
		A[5]=cos(alpha)*sin(beta)*cos(gamma)-sin(alpha)*sin(gamma);
		A[6]=cos(beta)*sin(gamma);
		A[7]=-sin(beta);
		A[8]=cos(beta)*cos(gamma);*/

		vector<pair<int,double>> index_dis(neighbourListSize[index]);
		for(int i=0; i<neighbourListSize[index]; i++) {
			size_t  neiIndex = neighbourList[neiListStartIndex+i];
			double dx = positionX[neiIndex] - x0, dy = positionY[neiIndex] - y0, dz = positionZ[neiIndex] - z0;
			double dis = dx*dx+dy*dy+dz*dz;
			index_dis[i]={neiIndex,dis};
		}
		std::sort(index_dis.begin(), index_dis.end(), ascByDist);
		for(size_t i=0; i<index_dis.size(); i++) {
			neighbourList[neiListStartIndex+i] = index_dis[i].first;
		}
/*		double penalty[8];
		for(int i=0;i<8;i++)
			penalty[i]=1.0;
		double penalty_weight=10000.0;
		for(int i=0; i<neighbourListSize[index]; i++) {
			size_t  neiIndex = neighbourList[neiListStartIndex+i];
			double dx = positionX[neiIndex] - x0, dy = positionY[neiIndex] - y0, dz = positionZ[neiIndex] - z0;
			double dis = dx*dx+dy*dy+dz*dz;
			int region=4*((A[0]*dx+A[1]*dy+A[2]*dz)>0)+2*((A[3]*dx+A[4]*dy+A[5]*dz)>0)+((A[6]*dx+A[7]*dy+A[8]*dz)>0);
			dis=dis*penalty[region];
			penalty[region]=penalty[region]*penalty_weight;

			index_dis[i]={neiIndex,dis};

		}*/

	        double penalty[6];
                for(int i=0;i<6;i++)
                        penalty[i]=1.0;
                double penalty_weight=10000.0;
                for(int i=0; i<neighbourListSize[index]; i++) {
                        size_t  neiIndex = neighbourList[neiListStartIndex+i];
                        double dx = positionX[neiIndex] - x0, dy = positionY[neiIndex] - y0, dz = positionZ[neiIndex] - z0;
                        double dis = dx*dx+dy*dy+dz*dz;
			int region;
			if ((fabs(dx)>=fabs(dy))&&(fabs(dx)>=fabs(dz)))
			{
				if(dx>0)
					region=0;
				else
					region=1;
			}
                        else if ((fabs(dy)>=fabs(dx))&&(fabs(dy)>=fabs(dz)))
                        {
                                if(dy>0)
                                        region=2;
                                else
                                        region=3;
                        }
                        else
                        {
                                if(dz>0)
                                        region=4;
                                else
                                        region=5;
                        }
                        dis=dis*penalty[region];
                        penalty[region]=penalty[region]*penalty_weight;
                        index_dis[i]={neiIndex,dis};

                }


		std::sort(index_dis.begin(), index_dis.end(), ascByDist);
//		cout<<endl;
//		cout<<index<<" "<<x0<<" "<<y0<<" "<<z0<<endl;
		for(size_t i=0; i<index_dis.size(); i++) {
			neighbourList[neiListStartIndex+i] = index_dis[i].first;
		}

//reordering neighbour list to balance neighbours from each quedrant end

	}

}

#ifdef _OPENMP
}
#endif

// DEBUG (check if the upwind neighbour list is correct)
//checkUpwindNeighbourList();

cout<<"-------END HyperbolicLPSolver::setUpwindNeighbourList()-------"<<endl;
}



void HyperbolicLPSolver::resetLPFOrder() {
size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;	
size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
if(m_iDimension==3) {
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		m_pParticleData->m_vLPFOrderRight[index] = m_iLPFOrder;
		m_pParticleData->m_vLPFOrderLeft[index] = m_iLPFOrder-m_iIfLaxWendroff*m_iUseLimiter;
		m_pParticleData->m_vLPFOrderNorth[index] = m_iLPFOrder;
		m_pParticleData->m_vLPFOrderSouth[index] = m_iLPFOrder-m_iIfLaxWendroff*m_iUseLimiter;
		m_pParticleData->m_vLPFOrderUp[index] = m_iLPFOrder;
		m_pParticleData->m_vLPFOrderDown[index] = m_iLPFOrder-m_iIfLaxWendroff*m_iUseLimiter;	
	}
}
else if(m_iDimension==2) {
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif	
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		m_pParticleData->m_vLPFOrderRight[index] = m_iLPFOrder;
		m_pParticleData->m_vLPFOrderLeft[index] = m_iLPFOrder-m_iIfLaxWendroff*m_iUseLimiter;
		m_pParticleData->m_vLPFOrderNorth[index] = m_iLPFOrder;
		m_pParticleData->m_vLPFOrderSouth[index] = m_iLPFOrder-m_iIfLaxWendroff*m_iUseLimiter;
	}	
}

if(m_iUseLimiter) {

	if(m_iDimension==3) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vLPFFirstOrderRight[index] = 1;
			m_pParticleData->m_vLPFFirstOrderLeft[index] = 1;
			m_pParticleData->m_vLPFFirstOrderNorth[index] = 1;
			m_pParticleData->m_vLPFFirstOrderSouth[index] = 1;
			m_pParticleData->m_vLPFFirstOrderUp[index] = 1;
			m_pParticleData->m_vLPFFirstOrderDown[index] = 1;	
		}
	}
	else if(m_iDimension==2) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif	
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vLPFFirstOrderRight[index] = 1;
			m_pParticleData->m_vLPFFirstOrderLeft[index] = 1;
			m_pParticleData->m_vLPFFirstOrderNorth[index] = 1;
			m_pParticleData->m_vLPFFirstOrderSouth[index] = 1;
		}	
	}

}

}


void HyperbolicLPSolver::cavitation() {
if(m_iUseCriticalPressure==0) return; // no cavitation model

size_t startIndex = m_pParticleData->getFluidStartIndex();	
size_t endIndex = startIndex + m_pParticleData->getFluidNum();
size_t count = 0;

#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {
	if(m_pParticleData->m_vPressure[index] < m_fCriticalPressure) {
		m_pParticleData->m_vPressure[index]=0;
		m_pParticleData->m_vSoundSpeed[index] = 
		m_pEOS->getSoundSpeed(m_pParticleData->m_vPressure[index],1./m_pParticleData->m_vVolume[index]);
		count++;
	}
}

cout<<"Cavitation: number of particles: "<<count<<endl;
}


void HyperbolicLPSolver::computeMinParticleSpacing() {

// set particle position pointers (to save typing)
const double *positionX = m_pParticleData->m_vPositionX;
const double *positionY = m_pParticleData->m_vPositionY;
const double *positionZ = m_pParticleData->m_vPositionZ;

// whole neighbour list
const int *neighbourList = m_pParticleData->m_vNeighbourList;
const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;

// initial value
m_fMinParticleSpacing=numeric_limits<double>::max();	

size_t startIndex = m_pParticleData->getFluidStartIndex();
size_t endIndex = startIndex + m_pParticleData->getFluidNum();

#ifdef _OPENMP
size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	cout<<"numThreads="<<numThreads<<endl;
//	cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<endl;
vector<double> minsp(numThreads,numeric_limits<double>::max());	
#pragma omp parallel  
{
int tid = omp_get_thread_num();
#endif

if(m_iDimension==3) {
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

		size_t totalNumNei = neighbourListSize[index];	
		for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dz = positionZ[neiIndex] - positionZ[index];
				double dist  = sqrt(dx*dx + dy*dy + dz*dz);
				#ifdef _OPENMP
					minsp[tid] = min(minsp[tid],dist);
				#else
					m_fMinParticleSpacing = min(m_fMinParticleSpacing,dist);
				#endif
				
				break;
			}
		}
	}
}
else if(m_iDimension==2) {
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

		size_t totalNumNei = neighbourListSize[index];	
		for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dist  = sqrt(dx*dx + dy*dy);
				#ifdef _OPENMP
					minsp[tid] = min(minsp[tid],dist);
				#else
					m_fMinParticleSpacing = min(m_fMinParticleSpacing,dist);
				#endif	
				break;
			}
		}
	}
	#ifdef _OPENMP
	for(size_t i=0; i<numThreads; i++) {
		m_fMinParticleSpacing = min(m_fMinParticleSpacing,minsp[i]);
	}
	#endif

}

#ifdef _OPENMP
	for(size_t i=0; i<numThreads; i++) {
		m_fMinParticleSpacing = min(m_fMinParticleSpacing,minsp[i]);	
	}
}
#endif

// the corner case when no fluid particle has a fluid neighbour
assert(m_fMinParticleSpacing != numeric_limits<double>::max());

//debug<<"-------HyperbolicLPSolver::computeMinParticleSpacing()-------"<<endl;
cout<<"m_fMinParticleSpacing="<<m_fMinParticleSpacing<<endl;
//debug<<"-------------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::computeAvgParticleSpacing() {

// set particle position pointers (to save typing)
const double *positionX = m_pParticleData->m_vPositionX;
const double *positionY = m_pParticleData->m_vPositionY;
const double *positionZ = m_pParticleData->m_vPositionZ;

// whole neighbour list
const int *neighbourList = m_pParticleData->m_vNeighbourList;
const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;
	
double sumDist = 0;
size_t num = 0;
size_t startIndex = m_pParticleData->getFluidStartIndex();
size_t endIndex = startIndex + m_pParticleData->getFluidNum();
if(m_iDimension==3) {
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

		size_t totalNumNei = neighbourListSize[index];	
		for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dz = positionZ[neiIndex] - positionZ[index];
				double dist  = sqrt(dx*dx + dy*dy + dz*dz);
				sumDist += dist;
				num++;
				break;
			}
		}
	}
}
else if(m_iDimension==2) {
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle

		size_t totalNumNei = neighbourListSize[index];	
		for(size_t i=0; i<totalNumNei; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // find one fluid neighbour, compare, and break
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dist  = sqrt(dx*dx + dy*dy);
				sumDist += dist;
				num++;
				break;
			}
		}
	}
}

assert(num!=0);	

//cout<<"-------HyperbolicLPSolver::computeAvgParticleSpacing()-------"<<endl;
cout<<"m_fAvgParticleSpacing changed from "<<m_fAvgParticleSpacing;
m_fAvgParticleSpacing = sumDist/(double)num;
cout<<" to "<<m_fAvgParticleSpacing<<endl;
//cout<<"-------------------------------------------------------------"<<endl;
}


void HyperbolicLPSolver::updateNeighbourSearchRadius() {

m_fNeiSearchRadius = m_fTimesNeiSearchRadius * m_fAvgParticleSpacing;

//cout<<"-------HyperbolicLPSolver::updateNeighbourSearchRadius()-------"<<endl;
//cout<<"m_fNeiSearchRadius="<<m_fNeiSearchRadius<<endl;
//cout<<"-------------------------------------------------------------"<<endl;

}



void HyperbolicLPSolver::updateLocalParSpacing() {

// set particle position pointers (to save typing)
const double *positionX = m_pParticleData->m_vPositionX;
const double *positionY = m_pParticleData->m_vPositionY;
const double *positionZ = m_pParticleData->m_vPositionZ;

double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

// whole neighbour list
const int *neighbourList = m_pParticleData->m_vNeighbourList;
const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;

size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;
		
size_t startIndex = m_pParticleData->getFluidStartIndex();
size_t endIndex = startIndex + m_pParticleData->getFluidNum();

if(m_iDimension==3) {
	
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle
		
		double sum = 0;
		size_t num = 0;
		for(int i=0; i<neighbourListSize[index]; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // if fluid neighbour
				for(int j=0; j<neighbourListSize[neiIndex]; j++) {
					size_t neineiIndex = neighbourList[neiIndex*maxNeiNum+j];
					if(neineiIndex>=startIndex && neineiIndex<endIndex) {
						double dx = positionX[neineiIndex] - positionX[neiIndex];
						double dy = positionY[neineiIndex] - positionY[neiIndex];
						double dz = positionZ[neineiIndex] - positionZ[neiIndex];
						double dist  = sqrt(dx*dx + dy*dy + dz*dz);
						sum += dist;
						num++;
						break;	
					}
				}
			}
		}
		for(int i=0; i<neighbourListSize[index]; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // if fluid neighbour
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dz = positionZ[neiIndex] - positionZ[index];
				double dist  = sqrt(dx*dx + dy*dy + dz*dz);
				sum += dist;
				num++;
				break;		
			}
		}
		if(num>0) {
			localParSpacing[index]=max(sum/(double)num,m_fAvgParticleSpacing);	
		}
		else {
			cout<<"index="<<index<<", No Neighbour!!!"<<endl;
			cout<<"(x,y,z)=("<<positionX[index]<<","<<positionY[index]<<","<<positionZ[index]<<")"<<endl;
			localParSpacing[index]=m_fAvgParticleSpacing;
			//assert(false);
		}
	}

}
else if(m_iDimension==2) {
	for(size_t index=startIndex; index<endIndex; index++) { // for each fluid particle
		
		double sum = 0;
		size_t num = 0;
		for(int i=0; i<neighbourListSize[index]; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // if fluid neighbour
				for(int j=0; j<neighbourListSize[neiIndex]; j++) {
					size_t neineiIndex = neighbourList[neiIndex*maxNeiNum+j];
					if(neineiIndex>=startIndex && neineiIndex<endIndex) {
						double dx = positionX[neineiIndex] - positionX[neiIndex];
						double dy = positionY[neineiIndex] - positionY[neiIndex];
						double dist  = sqrt(dx*dx + dy*dy);
						sum += dist;
						num++;
						break;	
					}
				}
			}
		}
		for(int i=0; i<neighbourListSize[index]; i++) { // for all of its neighbours
			size_t neiIndex = neighbourList[index*maxNeiNum+i];
			if(neiIndex>=startIndex && neiIndex<endIndex) { // if fluid neighbour
				double dx = positionX[neiIndex] - positionX[index];
				double dy = positionY[neiIndex] - positionY[index];
				double dist  = sqrt(dx*dx + dy*dy);
				sum += dist;
				num++;
				break;		
			}
		}
		if(num>0) {
			localParSpacing[index]=max(sum/(double)num,m_fAvgParticleSpacing);	
		}
		else {
			cout<<"index="<<index<<", No Neighbour!!!"<<endl;
			cout<<"(x,y)=("<<positionX[index]<<","<<positionY[index]<<")"<<endl;
			localParSpacing[index]=m_fAvgParticleSpacing;
			//assert(false);
		}
	}
}

}


void HyperbolicLPSolver::updateLocalParSpacingByVolume() {
	
const double *volume    = m_pParticleData->m_vVolume;
const double *volumeOld = m_pParticleData->m_vVolumeOld;
		
double *localParSpacing = m_pParticleData->m_vLocalParSpacing;

size_t startIndex = m_pParticleData->getFluidStartIndex();
size_t endIndex = startIndex + m_pParticleData->getFluidNum();

#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {
	if(volume[index]==0 || volumeOld[index]==0) {
		cout<<"Cannot update local_spacing because: volume["<<index<<"]="
		<<volume[index]<<"  volumeOld["<<index<<"]="<<volumeOld[index]<<endl;
		//continue;
	}
	else {
if(m_iDimension==3) 
		localParSpacing[index] *= std::cbrt(volume[index]/volumeOld[index]);
if(m_iDimension==2)
                localParSpacing[index] *= std::sqrt(volume[index]/volumeOld[index]);

	}
}	
}


void HyperbolicLPSolver::updateContactLength() {

m_fContactLength = m_fTimesContactLength * m_fAvgParticleSpacing;

//cout<<"-------HyperbolicLPSolver::updateContactLength()-------"<<endl;
//cout<<"m_fContactLength="<<m_fContactLength<<endl;
//cout<<"-------------------------------------------------------------"<<endl;

}


void HyperbolicLPSolver::computeMaxSoundSpeed() {

const double* soundSpeed = m_pParticleData->m_vSoundSpeed;

//initial value
//m_fMaxSoundSpeed = numeric_limits<double>::min();
m_fMaxSoundSpeed = -1;

size_t startIndex = m_pParticleData->m_iFluidStartIndex;
size_t endIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;	

#ifdef _OPENMP
size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	cout<<"numThreads="<<numThreads<<endl;
//	cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<endl;
vector<double> maxcs(numThreads,-1);	
#pragma omp parallel  
{
int tid = omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp for
#endif
for(size_t index=startIndex; index<endIndex; index++) {
	#ifdef _OPENMP
//	if(m_pParticleData->m_vVolume[index]<1e+4)
		maxcs[tid] = max(maxcs[tid],soundSpeed[index]);
	#else
//        if(m_pParticleData->m_vVolume[index]<1e+4)
		m_fMaxSoundSpeed = max(m_fMaxSoundSpeed,soundSpeed[index]);
	#endif
}

#ifdef _OPENMP
	for(size_t i=0; i<numThreads; i++) {
		m_fMaxSoundSpeed = max(m_fMaxSoundSpeed,maxcs[i]);	
	}
}
#endif

//assert(m_fMaxSoundSpeed != -1);

//cout<<"-------HyperbolicLPSolver::computeMaxSoundSpeed()-------"<<endl;
cout<<"m_fMaxSoundSpeed="<<m_fMaxSoundSpeed<<endl;
//cout<<"--------------------------------------------------------"<<endl;
assert(m_fMaxSoundSpeed != -1);

}


void HyperbolicLPSolver::computeMaxFluidVelocity() {

const double* vU = m_pParticleData->m_vVelocityU;
const double* vV = m_pParticleData->m_vVelocityV;
const double* vW = (m_iDimension==3)? m_pParticleData->m_vVelocityW:nullptr;

// initial value
m_fMaxFluidVelocity = -1;

size_t startIndex = m_pParticleData->getFluidStartIndex();
size_t endIndex = startIndex + m_pParticleData->getFluidNum();	
#ifdef _OPENMP
size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	cout<<"numThreads="<<numThreads<<endl;
//	cout<<"omp_get_max_threads()="<<omp_get_max_threads()<<endl;
vector<double> maxs(numThreads,-1);	
#pragma omp parallel  
{
int tid = omp_get_thread_num();
#endif
	
if(m_iDimension==3) {
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=startIndex; index<endIndex; index++) { 
		double speed = vU[index]*vU[index]+vV[index]*vV[index]+vW[index]*vW[index];
		#ifdef _OPENMP
//	        if(m_pParticleData->m_vVolume[index]<1e+4)
			maxs[tid] = max(maxs[tid],speed);
		#else
//	        if(m_pParticleData->m_vVolume[index]<1e+4)
			m_fMaxFluidVelocity = max(m_fMaxFluidVelocity,speed);
		#endif
	}
}
else if(m_iDimension==2) {
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=startIndex; index<endIndex; index++) {
		double speed = vU[index]*vU[index]+vV[index]*vV[index];
		#ifdef _OPENMP
//	        if(m_pParticleData->m_vVolume[index]<1e+4)
			maxs[tid] = max(maxs[tid],speed);
		#else
//	        if(m_pParticleData->m_vVolume[index]<1e+4)
			m_fMaxFluidVelocity = max(m_fMaxFluidVelocity,speed);	
		#endif
	}
}

#ifdef _OPENMP
	for(size_t i=0; i<numThreads; i++) {
		m_fMaxFluidVelocity = max(m_fMaxFluidVelocity,maxs[i]);	
	}
}
#endif

assert(m_fMaxFluidVelocity != -1);

m_fMaxFluidVelocity = sqrt(m_fMaxFluidVelocity);

//cout<<"-------HyperbolicLPSolver::computeMaxFluidVelocity()-------"<<endl;
cout<<"m_fMaxFluidVelocity="<<m_fMaxFluidVelocity<<endl;
//cout<<"--------------------------------------------------------"<<endl;

}


bool HyperbolicLPSolver::directionalSplitting_van_leer(int phase) {	
cout<<"--------------HyperbolicLPSolver::directionalSplitting_van_leer()--------------"<<endl;

// determine dir: x(0), y(1), or z(2)
const int dir = m_vDirSplitTable[m_iDirSplitOrder][phase];
//cout<<"dir="<<dir<<endl;	
	
// set neighbour list pointers by dir (dir=0->right/left, dir=1->north/south, dir=2->up/down)
const int *neiList0=nullptr, *neiList1=nullptr;
const int *neiListSize0=nullptr, *neiListSize1=nullptr;
setNeighbourListPointers(dir, &neiList0, &neiList1, &neiListSize0, &neiListSize1);

// input data pointers
const double *inVelocity=nullptr, *inPressure=nullptr, *inVolume=nullptr, *inSoundSpeed=nullptr;
// output data pointers
double *outVelocity=nullptr, *outPressure=nullptr, *outVolume=nullptr, *outSoundSpeed=nullptr;
// set data pointers by phase (for temp1 or temp2) and dir(for velocity U or V) 
setInAndOutDataPointers(phase,dir,&inVelocity,&inPressure,&inVolume,&inSoundSpeed,
						&outVelocity,&outPressure,&outVolume,&outSoundSpeed);

// set local polynomail order pointers
// dir==0->right(0) & left(1), dir==1->north(0) & south(1), dir==2->up(0) & down(1)
int *LPFOrder0=nullptr, *LPFOrder1=nullptr, *LPFFirstOrder0=nullptr, *LPFFirstOrder1=nullptr;		
setLPFOrderPointers(dir,&LPFOrder0,&LPFOrder1,&LPFFirstOrder0,&LPFFirstOrder1);

bool phaseSuccess = true;
 
double gravity;
if(m_iDimension==2) gravity = dir==1? m_fGravity:0; // only in y direction
else if(m_iDimension==3) gravity = dir==2? m_fGravity:0; // only in z direction

// set real dt:  
// 2D: phase:	0		1     2 
//	 		   dt/4   dt/2   dt/4
// 3D: phase:		0		1		2		3		4
//				  dt/6    dt/6     dt/3    dt/6    dt/6
double realDt;
if(m_iDimension==2) realDt = phase==1? m_fDt/2.:m_fDt/4.;
else if(m_iDimension==3) realDt = phase==2? m_fDt/3.:m_fDt/6.;

// set the function pointer to compute the A matrix for QR solver
void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*,double*);
int offset; // offset is used to get results of computed spatial derivatives from QR solver
if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}

// the coeff before first and second order term during time integration
double multiplier1st, multiplier2nd;
if(m_iDimension==2) {multiplier1st=2; multiplier2nd=m_fDt/2.;}
else if(m_iDimension==3) {multiplier1st=3; multiplier2nd=m_fDt*3./4.;}

size_t startIndex = m_pParticleData->m_iFluidStartIndex;	
size_t endIndex = startIndex + m_pParticleData->m_iFluidNum; 

int additional = 0;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {

/*		if(index==3461)//debug 12/04
	{
		size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
		cout<<"Index = 3461, Position = ("<<m_pParticleData->m_vPositionX[index]<<", "<<m_pParticleData->m_vPositionY[index]<<"), volume = "<<inVolume[index]<<", velocity = "<< inVelocity[index]<<", pressure = "<<inPressure[index]<<", speed of sound = "<<inSoundSpeed[index]<<endl;
		if(inPressure[index]<1.0e-6)
			cout<<"Warning: pressure too small"<<endl;
			
		
		cout<<"Neighbour list 0 size"<<neiListSize0[index]<<endl;
		cout<<"Neighbour list 0: index, x, y"<<endl;
		for(int i=0;i<neiListSize0[index];i++)
		{
			int temp_index=neiList0[index*maxNeiNum+i];
			cout<<temp_index<<" ("<<m_pParticleData->m_vPositionX[temp_index]<<", "<<m_pParticleData->m_vPositionY[temp_index]<<")"<<endl;
		}
		cout<<"Neighbour list 1 size"<<neiListSize1[index]<<endl;
		cout<<"Neighbour list 1: index, x, y"<<endl;
		for(int i=0;i<neiListSize1[index];i++)
		{
			int temp_index=neiList1[index*maxNeiNum+i];
			cout<<temp_index<<" ("<<m_pParticleData->m_vPositionX[temp_index]<<", "<<m_pParticleData->m_vPositionY[temp_index]<<")"<<endl;
		}
	}	//debug 12/04*/

	if(inSoundSpeed[index]==0 || inVolume[index]==0) {
		
		cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
		<<", vol="<<inVolume[index]<<", vel="<<inVelocity[index]
		<<", cs="<<inSoundSpeed[index]<<endl;

		outVolume[index]   = inVolume[index];
		outPressure[index] = inPressure[index];
		outVelocity[index] = inVelocity[index];
		outSoundSpeed[index] = inSoundSpeed[index];
	}
	else {
		
		if(LPFOrder0[index]<=1 && LPFOrder1[index]<=1) { // high order regime drops to <=1 on both directions 
			// spatial derivatives (0:right/north/up; 1:left/south/down)
			double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output	
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList0, neiListSize0, additional,
							  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList1, neiListSize1, additional,
							  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output
				
			timeIntegration(realDt, multiplier1st, multiplier2nd, 
							gravity, inVolume[index], inVelocity[index], inPressure[index], inSoundSpeed[index], 
							vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1,
							&outVolume[index], &outVelocity[index], &outPressure[index]); // output	
		}
		else {
			// spatial derivatives (0:right/north/up; 1:left/south/down)
			double vel_d_0_first, vel_d_1_first, p_d_0_first, p_d_1_first; // first order first der
			double vel_d_0_second, vel_d_1_second, p_d_0_second, p_d_1_second; // second order first der
			double vel_dd_0_second, vel_dd_1_second, p_dd_0_second, p_dd_1_second; // second order second der
		
			double tmp1, tmp2;	
			// first prder LPF
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList0, neiListSize0, additional,
							  LPFFirstOrder0, &vel_d_0_first, &tmp1, &p_d_0_first, &tmp2); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList1, neiListSize1, additional,
							  LPFFirstOrder1, &vel_d_1_first, &tmp1, &p_d_1_first, &tmp2); // output
			//second order LPF	
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList0, neiListSize0, additional,
							  LPFOrder0, &vel_d_0_second, &vel_dd_0_second, &p_d_0_second, &p_dd_0_second); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList1, neiListSize1, additional,
							  LPFOrder1, &vel_d_1_second, &vel_dd_1_second, &p_d_1_second, &p_dd_1_second); // output
			
			// compute theta
			double theta_vel, theta_p;
			if(vel_d_1_first==0 || vel_d_0_first==0) theta_vel = -1;
//				else theta_vel = 0.5*(vel_d_1_first/vel_d_0_first+vel_d_0_first/vel_d_1_first);
			else theta_vel = 2.0/(vel_d_1_first/vel_d_0_first+vel_d_0_first/vel_d_1_first);
			
			if(p_d_1_first==0 || p_d_0_first==0) theta_p = -1;
//				else theta_p = 0.5*(p_d_1_first/p_d_0_first+p_d_0_first/p_d_1_first);
			else theta_p = 2.0/(p_d_1_first/p_d_0_first+p_d_0_first/p_d_1_first);
			
			// compute phi
			double phi_vel, phi_p;
			if(theta_vel <= 0) phi_vel = 0;
			else phi_vel = 2*theta_vel/(1+theta_vel);
			
			if(theta_p <= 0) phi_p = 0;
			else phi_p = 2*theta_p/(1+theta_p);

			double phi = min(phi_vel,phi_p);

			timeIntegration_van_leer(realDt, multiplier1st, multiplier2nd, 
									 gravity, inVolume[index], inVelocity[index], 
									 inPressure[index], inSoundSpeed[index], 
									 phi,
									 vel_d_0_first, vel_d_1_first, p_d_0_first, p_d_1_first,
									 vel_d_0_second, vel_d_1_second, p_d_0_second, p_d_1_second,
									 vel_dd_0_second, vel_dd_1_second, p_dd_0_second, p_dd_1_second,
									 &outVolume[index], &outVelocity[index], &outPressure[index]); // output

/*				if(index==3461){//12/04 debug 
				cout<<"phi="<<phi<<endl;
				cout<<"Volume= "<<outVolume[index]<<", velocity ="<<outVelocity[index]<<", pressure = "<<outPressure[index];
			}//12/04 debug*/

		}
		bool redo = false;
		if(outPressure[index]<m_fInvalidPressure) redo=true;
		else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
		else if(std::isnan(outVolume[index]) || std::isnan(outVelocity[index]) || std::isnan(outPressure[index])) {
			cout<<"nan value!!!"<<endl;
			redo=true;
		}
		else if(std::isinf(outVolume[index]) || std::isinf(outVelocity[index]) || std::isinf(outPressure[index])) {
			cout<<"inf value!!!"<<endl;
			redo=true;
		}
		if(redo) {
			if(LPFOrder0[index]>0) LPFOrder0[index]--;
			if(LPFOrder1[index]>0) LPFOrder1[index]--;
			if(LPFOrder0[index]==0 && LPFOrder1[index]==0) {
				outVolume[index]   = inVolume[index];
				outPressure[index] = inPressure[index];
				outVelocity[index] = inVelocity[index];
				outSoundSpeed[index] = inSoundSpeed[index];	
			}
			else {
				index--;	
			}
		}
		else {
			if(outVolume[index]==0)
				outSoundSpeed[index] = 0;
			else
				outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
		}

		

	}
}	

return phaseSuccess;
}


bool HyperbolicLPSolver::directionalSplitting_lax_wendroff(int phase) {
cout<<"--------------HyperbolicLPSolver::directionalSplitting_lax_wendroff()--------------"<<endl;

// determine dir: x(0), y(1), or z(2)
const int dir = m_vDirSplitTable[m_iDirSplitOrder][phase];
//cout<<"dir="<<dir<<endl;      

// set neighbour list pointers by dir (dir=0->right/left, dir=1->north/south, dir=2->up/down)
const int *neiList0=nullptr, *neiList1=nullptr;
const int *neiListSize0=nullptr, *neiListSize1=nullptr;
const int *neiListLW=m_pParticleData->m_vNeighbourList;
const int *neiListSizeLW=m_pParticleData->m_vNeighbourListSize;
setNeighbourListPointers(dir, &neiList0, &neiList1, &neiListSize0, &neiListSize1);

// input data pointers
const double *inVelocity=nullptr, *inPressure=nullptr, *inVolume=nullptr, *inSoundSpeed=nullptr;
// output data pointers
double *outVelocity=nullptr, *outPressure=nullptr, *outVolume=nullptr, *outSoundSpeed=nullptr;
// set data pointers by phase (for temp1 or temp2) and dir(for velocity U or V) 
setInAndOutDataPointers(phase,dir,&inVelocity,&inPressure,&inVolume,&inSoundSpeed,
						&outVelocity,&outPressure,&outVolume,&outSoundSpeed);

// set local polynomail order pointers
// dir==0->right(0) & left(1), dir==1->north(0) & south(1), dir==2->up(0) & down(1)
int *LPFOrder0=nullptr, *LPFOrder1=nullptr, *LPFFirstOrder0=nullptr, *LPFFirstOrder1=nullptr;
setLPFOrderPointers(dir,&LPFOrder0,&LPFOrder1,&LPFFirstOrder0,&LPFFirstOrder1);

bool phaseSuccess = true;

double gravity;
if(m_iDimension==2) gravity = dir==1? m_fGravity:0; // only in y direction
else if(m_iDimension==3) gravity = dir==2? m_fGravity:0; // only in z direction

// set real dt:  
// 2D: phase:   0               1     2 
//                         dt/4   dt/2   dt/4
// 3D: phase:           0               1               2               3               4
//                                dt/6    dt/6     dt/3    dt/6    dt/6
double realDt;
if(m_iDimension==2) realDt = phase==1? m_fDt/2.:m_fDt/4.;
else if(m_iDimension==3) realDt = phase==2? m_fDt/3.:m_fDt/6.;

// set the function pointer to compute the A matrix for QR solver
void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*);
int offset; // offset is used to get results of computed spatial derivatives from QR solver
if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}

// the coeff before first and second order term during time integration
double multiplier1st, multiplier2nd;
if(m_iDimension==2) {multiplier1st=2; multiplier2nd=m_fDt/2.;}
else if(m_iDimension==3) {multiplier1st=3; multiplier2nd=m_fDt*3./4.;}

size_t startIndex = m_pParticleData->m_iFluidStartIndex;
size_t endIndex = startIndex + m_pParticleData->m_iFluidNum;

int additional=0;
//	if(m_iDimension==2) additional=3;
//	else if (m_iDimension==3) additional=8;
int additional1=0;
if(m_iDimension==2) additional1=3;
else if (m_iDimension==3) additional1=6;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {
	double vel_error1, vel_error2;
	double p_error1, p_error2;

	if(inSoundSpeed[index]==0 || inVolume[index]==0) {

		cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
		<<", vol="<<inVolume[index]<<", vel="<<inVelocity[index]
		<<", cs="<<inSoundSpeed[index]<<endl;

		outVolume[index]   = inVolume[index];
		outPressure[index] = inPressure[index];
		outVelocity[index] = inVelocity[index];
		outSoundSpeed[index] = inSoundSpeed[index];
	}
	else {

		if(LPFOrder0[index]<=1 && LPFOrder1[index]<=1) { // high order regime drops to <=1 on both directions 
			// spatial derivatives (0:right/north/up; 1:left/south/down)
			double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output    
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList0, neiListSize0, additional,
							  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList1, neiListSize1, additional,
							  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output

			timeIntegration(realDt, multiplier1st, multiplier2nd,
							gravity, inVolume[index], inVelocity[index], inPressure[index], inSoundSpeed[index],
							vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1,
							&outVolume[index], &outVelocity[index], &outPressure[index]); // output 
		}
		else {
			// spatial derivatives (0:right/north/up; 1:left/south/down)
			double vel_d_0_first, vel_d_1_first, p_d_0_first, p_d_1_first; // first order first der
			double vel_d_0_second, vel_d_1_second, p_d_0_second, p_d_1_second; // second order first der
			double vel_dd_0_second, vel_dd_1_second, p_dd_0_second, p_dd_1_second; // second order second der

			double tmp1, tmp2;
			// first prder LPF
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList0, neiListSize0,additional,
							  LPFFirstOrder0, &vel_d_0_first, &tmp1, &p_d_0_first, &tmp2); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiList1, neiListSize1,additional,
							  LPFFirstOrder1, &vel_d_1_first, &tmp1, &p_d_1_first, &tmp2); // output
			//second order LPF      
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiListLW, neiListSizeLW,additional,
							  LPFOrder0, &vel_d_0_second, &vel_dd_0_second, &p_d_0_second, &p_dd_0_second, &vel_error2, &p_error2); // output
			computeSpatialDer(dir, index, offset, computeA,
							  inPressure, inVelocity, neiListLW, neiListSizeLW,additional1,
							  LPFOrder1, &vel_d_1_second, &vel_dd_1_second, &p_d_1_second, &p_dd_1_second, &vel_error1, &p_error1); // output
			vel_d_1_second=vel_d_0_second;
			vel_dd_1_second=vel_dd_0_second;
			p_d_1_second=p_d_0_second;
			p_dd_1_second=p_dd_0_second;

			
			double phi_vel, phi_p;
			if(vel_error1<0.00001)
				phi_vel=0.0;
			else
			{
				if(vel_error2>1.0*vel_error1)
					phi_vel=0.0;
				else
					phi_vel=(1.0-1.0*vel_error2/vel_error1);
			}

			if(p_error1<0.00001)
				phi_p=0.0;
			else
			{
				if(p_error2>1.0*p_error1)
					phi_p=0.0;
				else
					phi_p=(1.0-1.0*p_error2/p_error1);
			}				

			double phi = min(phi_vel,phi_p);
			m_pParticleData->m_vPhi[index]=phi;
			m_pParticleData->m_vPError0[index]=p_error2;
			m_pParticleData->m_vPError1[index]=p_error1;
			m_pParticleData->m_vVelError0[index]=vel_error2;
			m_pParticleData->m_vVelError1[index]=vel_error1;

			if(LPFOrder0[index]<2)
				phi=0.0;
//			        if(std::isnan(phi) || std::isinf(phi)){
//               				 printf("%d %f %f %f %f\n",index, vel_error1, vel_error2, p_error1, p_error2);
//					if(vel_error1<0.0001)
//						printf("How\n");
//         			      	 printf("%f %f %f\n",phi_vel,phi_p,phi);
//
//               				 assert(false);
//       				 }
			
			timeIntegration_van_leer(realDt, multiplier1st, multiplier2nd,
									 gravity, inVolume[index], inVelocity[index],
									 inPressure[index], inSoundSpeed[index],
									 phi,
									 vel_d_0_first, vel_d_1_first, p_d_0_first, p_d_1_first,
									 vel_d_0_second, vel_d_1_second, p_d_0_second, p_d_1_second,
									 vel_dd_0_second, vel_dd_1_second, p_dd_0_second, p_dd_1_second,
									 &outVolume[index], &outVelocity[index], &outPressure[index]); // output


		}
		bool redo = false;
		if(outPressure[index]<m_fInvalidPressure) redo=true;
		else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
		else if(std::isnan(outVolume[index]) || std::isnan(outVelocity[index]) || std::isnan(outPressure[index])) {
			cout<<"nan value!!!"<<endl;
			redo=true;
		}
		else if(std::isinf(outVolume[index]) || std::isinf(outVelocity[index]) || std::isinf(outPressure[index])) {
			cout<<"inf value!!!"<<endl;
			redo=true;
		}
		if(redo) {
			if(LPFOrder0[index]>0) LPFOrder0[index]--;
			if((LPFOrder0[index]==1)&&(LPFOrder1[index]>0)) LPFOrder1[index]--;
			if(LPFOrder0[index]==0 && LPFOrder1[index]==0) {
				outVolume[index]   = inVolume[index];
				outPressure[index] = inPressure[index];
				outVelocity[index] = inVelocity[index];
				outSoundSpeed[index] = inSoundSpeed[index];
			}
			else {
				index--;
			}
		}
		else {
			if(outVolume[index]==0)
				outSoundSpeed[index] = 0;
			else
				outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
		}



	}
}

return phaseSuccess;
}

bool HyperbolicLPSolver::density_derivative() {
const int *neiList=m_pParticleData->m_vNeighbourList;
const int *neiListSize=m_pParticleData->m_vNeighbourListSize;
const double *inVolume=m_pParticleData->m_vVolumeOld;
double *Volume_x=m_pParticleData->m_vVolume_x;
double *Volume_y=m_pParticleData->m_vVolume_y;
double *Volume_z=m_pParticleData->m_vVolume_z;
double *inDensity=m_pParticleData->m_vDensity;
int *LPFOrder0=nullptr, *LPFOrder1=nullptr;
vector<int*> LPFOrderOther;    
setLPFOrderPointers(0,&LPFOrder0,&LPFOrder1,LPFOrderOther);

bool phaseSuccess = true;

void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*);
int offset;
if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}
// iteration start index
size_t startIndex = m_pParticleData->m_iFluidStartIndex;
size_t endIndex = startIndex + m_pParticleData->m_iFluidNum;
for(size_t index=startIndex; index<(endIndex+m_pParticleData->m_iBoundaryNum+m_pParticleData->m_iGhostNum); index++)
	inDensity[index]=1.0/inVolume[index];

// iterate through fluid particles
#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {
                double volume_x, volume_y, volume_z; // output
                computeSpatialDer(index, offset, computeA,
                                                  inVolume, neiList, neiListSize,
                                                  LPFOrder0, &volume_x, &volume_y, &volume_z); // output
		Volume_x[index]=volume_x;
                Volume_y[index]=volume_y;
                Volume_z[index]=volume_z;
}
return phaseSuccess;
}





bool HyperbolicLPSolver::directionalSplitting(int phase) {	
//	cout<<"--------------HyperbolicLPSolver::directionalSplitting()--------------"<<endl;

// determine dir: x(0), y(1), or z(2)
const int dir = m_vDirSplitTable[m_iDirSplitOrder][phase];
//cout<<"dir="<<dir<<endl;	
// set particle position pointers (to save typing)
//const double *positionX = m_pParticleData->m_vPositionX;
//const double *positionY = m_pParticleData->m_vPositionY;
//const double *positionZ = m_pParticleData->m_vPositionZ;
//cout<<"1"<<endl;	
// set neighbour list pointers by dir (dir=0->right/left, dir=1->north/south, dir=2->up/down)
const int *neiList0=nullptr, *neiList1=nullptr;
const int *neiListSize0=nullptr, *neiListSize1=nullptr;
const int *neiListLW=m_pParticleData->m_vNeighbourList;
const int *neiListSizeLW=m_pParticleData->m_vNeighbourListSize;
setNeighbourListPointers(dir, &neiList0, &neiList1, &neiListSize0, &neiListSize1);
//cout<<"2"<<endl;
// input data pointers
const double *inVelocity=nullptr, *inPressure=nullptr, *inVolume=nullptr, *inSoundSpeed=nullptr;
// output data pointers
double *outVelocity=nullptr, *outPressure=nullptr, *outVolume=nullptr, *outSoundSpeed=nullptr;
// set data pointers by phase (for temp1 or temp2) and dir(for velocity U or V) 
setInAndOutDataPointers(phase,dir,&inVelocity,&inPressure,&inVolume,&inSoundSpeed,
						&outVelocity,&outPressure,&outVolume,&outSoundSpeed);
//cout<<"3"<<endl;
// set local polynomail order pointers
// dir==0->right(0) & left(1), dir==1->north(0) & south(1), dir==2->up(0) & down(1)
int *LPFOrder0=nullptr, *LPFOrder1=nullptr;
vector<int*> LPFOrderOther;	
setLPFOrderPointers(dir,&LPFOrder0,&LPFOrder1,LPFOrderOther);
//cout<<"4"<<endl;
// phase_success will be false if one particle LPF order is zero in one side
bool phaseSuccess = true;

// gravity is only on y(1) direction 
double gravity;
if(m_iDimension==2) gravity = dir==1? m_fGravity:0; // only in y direction
else if(m_iDimension==3) gravity = dir==1? m_fGravity:0; // only in z direction TODO: MODIFIED GRAVITY DIRECTION

// set real dt:  
// 2D: phase:	0		1     2 
//	 		   dt/4   dt/2   dt/4
// 3D: phase:		0		1		2		3		4
//				  dt/6    dt/6     dt/3    dt/6    dt/6
double realDt;
if(m_iDimension==2) realDt = phase==1? m_fDt/2.:m_fDt/4.;
else if(m_iDimension==3) realDt = phase==2? m_fDt/3.:m_fDt/6.;
//cout<<"5"<<endl;
// set the function pointer to compute the A matrix for QR solver
void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*);
int offset; // offset is used to get results of computed spatial derivatives from QR solver
if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}
//cout<<"6"<<endl;
// the coeff before first and second order term during time integration
double multiplier1st, multiplier2nd;
if(m_iDimension==2) {multiplier1st=2; multiplier2nd=m_fDt/2.;}
else if(m_iDimension==3) {multiplier1st=3; multiplier2nd=m_fDt*3./4.;}
/*	
if(m_iDimension==2) {
	multiplier1st=2;
	if(phase==1) multiplier2nd = m_fDt/4.;
	else multiplier2nd = m_fDt/8.;
	//if(phase==1) multiplier2nd = m_fDt/2.;
	//else multiplier2nd = m_fDt/4.;
}
else if(m_iDimension==3) {
	multiplier1st=3;
	if(phase==2) multiplier2nd = m_fDt/4.;
	else multiplier2nd = m_fDt/8.;
	//if(phase==2) multiplier2nd = m_fDt/2.;
	//else multiplier2nd = m_fDt/4.;
}
*/

// iteration start index
size_t startIndex = m_pParticleData->m_iFluidStartIndex;
//	cout<<"m_pParticleData->m_iFluidStartIndex="<<m_pParticleData->m_iFluidStartIndex<<endl;	
size_t endIndex = startIndex + m_pParticleData->m_iFluidNum; 
int additional=0;
int additional_lw=0;
//        if(m_iDimension==2) additional=3;
//        else if (m_iDimension==3) additional=8;
//cout<<"m_pParticleData->getFluidNum() = "<<m_pParticleData->getFluidNum()<<endl;
//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
//double startTime = omp_get_wtime();	

// iterate through fluid particles
#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {
			
	if(inSoundSpeed[index]==0 || inVolume[index]==0) {
		
		cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
		<<", vol="<<inVolume[index]<<", vel="<<inVelocity[index]
		<<", cs="<<inSoundSpeed[index]<<endl;

		outVolume[index]   = inVolume[index];
		outPressure[index] = inPressure[index];
		outVelocity[index] = inVelocity[index];
		outSoundSpeed[index] = inSoundSpeed[index];
	}
	else {
		bool redo = false;

//			cout<<"computeSpatialDer"<<", index="<<index<<endl;	
		// spatial derivatives (0:right/north/up; 1:left/south/down)
		double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
		if(m_iIfLaxWendroff==0)//beam warming
		{
			computeSpatialDer(dir, index, offset, computeA,//beam warming
						  inPressure, inVelocity, neiList0, neiListSize0, additional,
						  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
			computeSpatialDer(dir, index, offset, computeA,//beam warming
						  inPressure, inVelocity, neiList1, neiListSize1, additional,
						  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output
		}
		else//lax wendroff
		{
//test if it is local extreme
			int count_p=0,count_vel=0;
			int maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
			for(int i=0;i<8;i++)
			{
				if(inPressure[neiListLW[index*maxNeiNum+i]]>inPressure[index])
					count_p++;
				if(inVelocity[neiListLW[index*maxNeiNum+i]]>inVelocity[index])
					count_vel++;
			}
//			if(0)
			if((count_p==8)||(count_vel==8)||(count_p==0)||(count_vel==0))//if local exterme
			{
			LPFOrder0[index]=1;
			LPFOrder1[index]=1;
			computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList0, neiListSize0,additional,
						  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
			computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList1, neiListSize1,additional,
						  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output
//			if(m_iIfLaxWendroff&&(LPFOrder0[index]==2)&&(LPFOrder1[index]==2)){//if beam warming shows good neighbourhood, the nsue lax wendroff
			}
			else
			{
			LPFOrder0[index]=2;
			computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiListLW, neiListSizeLW,additional_lw, 
						  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
			vel_d_1=vel_d_0;
			vel_dd_1=vel_dd_0;
			p_d_1=p_d_0;
			p_dd_1=p_dd_0;
/*				if(((p_d_0*p_d_0+p_dd_0*p_dd_0+vel_d_0*vel_d_0+vel_dd_0*vel_dd_0)>0.0001)&&(fabs(m_pParticleData->m_vPositionX[index])>0.5)){
				printf("\n\n\n%d %f %f %f %f %f %f %f %f\n",index, m_pParticleData->m_vPositionX[index], m_pParticleData->m_vPositionY[index], inVelocity[index], inPressure[index], vel_d_0, vel_dd_0, p_d_0, p_dd_0);
				for(int j=m_pParticleData->m_iMaxNeighbourNum*index;j<m_pParticleData->m_iMaxNeighbourNum*index+neiListSizeLW[index];j++)
				printf("%d %f %f %f %f\n", neiListLW[j],m_pParticleData->m_vPositionX[neiListLW[j]], m_pParticleData->m_vPositionY[neiListLW[j]], inVelocity[neiListLW[j]], inPressure[neiListLW[j]]);
				printf("\n");
}*/
			if(LPFOrder0[index]<2)
			{
/*					LPFOrder0[index]=2;
				LPFOrder1[index]=2;
				redo=true;*/
				LPFOrder1[index]=LPFOrder0[index];
				computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList0, neiListSize0,additional,
						  LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
				computeSpatialDer(dir, index, offset, computeA,
						  inPressure, inVelocity, neiList1, neiListSize1,additional,
						  LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output

			}
			}
			m_pParticleData->m_vPhi[index]=LPFOrder0[index];
		}
//			cout<<"timeIntegration"<<", index="<<index<<endl;
		// update outVolume, outVelocity, outPressure
		double volumet,velocityt,pressuret;

		timeIntegration(realDt, multiplier1st, multiplier2nd, 
						gravity, inVolume[index], inVelocity[index], inPressure[index], inSoundSpeed[index], 
						vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1,
						&outVolume[index], &outVelocity[index], &outPressure[index], &volumet, &velocityt, &pressuret); // output	
/*debug 0909*/
		if(dir==0)
		{
			m_pParticleData->m_vPxl[index]=p_d_0;
			m_pParticleData->m_vPxr[index]=p_d_1;
			m_pParticleData->m_vVxl[index]=vel_d_0;
			m_pParticleData->m_vVxr[index]=vel_d_1;
			m_pParticleData->m_vPtx[index]=pressuret;
			m_pParticleData->m_vVtx[index]=velocityt;
			m_pParticleData->m_vVolumetx[index]=volumet;
		}
		if(dir==1)
		{
			m_pParticleData->m_vPyl[index]=p_d_0;
			m_pParticleData->m_vPyr[index]=p_d_1;
			m_pParticleData->m_vVyl[index]=vel_d_0;
			m_pParticleData->m_vVyr[index]=vel_d_1;
			m_pParticleData->m_vPty[index]=pressuret;
			m_pParticleData->m_vVty[index]=velocityt;
			m_pParticleData->m_vVolumety[index]=volumet;
		}
/*debug 0909*/

//			bool redo = false;
		if(outPressure[index]<m_fInvalidPressure) redo=true;
		else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
		else if(std::isnan(outVolume[index]) || std::isnan(outVelocity[index]) || std::isnan(outPressure[index])) {
			cout<<"nan value!!!"<<endl;
			redo=true;
		}
		else if(std::isinf(outVolume[index]) || std::isinf(outVelocity[index]) || std::isinf(outPressure[index])) {
			cout<<"inf value!!!"<<endl;
			redo=true;
		}

/*                        if(((m_pParticleData->m_vPositionX[index]*m_pParticleData->m_vPositionX[index]+m_pParticleData->m_vPositionY[index]*m_pParticleData->m_vPositionY[index])<0.8)&&(((LPFOrder0[index]+LPFOrder1[index])<4)||redo||(index==-1)))
		{
			printf("Warning: reduce to first order for particle %d, (x,y) = (%f %f), (V, U, P) = (%f %f %f).\n",index, m_pParticleData->m_vPositionX[index], m_pParticleData->m_vPositionY[index], inVolume[index],inVelocity[index],inPressure[index]);
			printf("dir = %d, neighboursize= (%d, %d, %d), LPFOrder= (%d %d), redo = %d\n", dir, neiListSizeLW[index], neiListSize0[index], neiListSize1[index], LPFOrder0[index],LPFOrder1[index],redo);
			printf("Output V U P = (%f %f %f)\n",outVolume[index],outVelocity[index],outPressure[index]);
			printf("vel_d_0 = %f, vel_dd_0 = %f, vel_d_1 = %f, vel_dd_1=%f, p_d_0 = %f, p_dd_0 = %f, p_d_1 = %f, p_dd_1 = %f\n",vel_d_0, vel_dd_0, vel_d_1, vel_dd_1, p_d_0, p_dd_0, p_d_1, p_dd_1);
			if(LPFOrder0[index]*LPFOrder1[index]==0)
				printf("WARNING: order = 0!\n");	
		}

		if(index==-1)
		{
			printf("Neighbour information for index = 3151.\ndx dy dv du dp\n");
			size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
			printf("List 0\n");
			 for(size_t i=0; i<neiListSize0[index]; i++) { // Note that the neighbour list does not contain the particle itself 
				int neiIndex = neiList0[index*maxNeiNum+i];
				printf("%f %f %f %f %f\n",m_pParticleData->m_vPositionX[neiIndex]-m_pParticleData->m_vPositionX[index],m_pParticleData->m_vPositionY[neiIndex]-m_pParticleData->m_vPositionY[index],inVolume[neiIndex]-inVolume[index],inVelocity[neiIndex]-inVelocity[index],inPressure[neiIndex]-inPressure[index]);
			}
			printf("List 1\n");
			 for(size_t i=0; i<neiListSize1[index]; i++) { // Note that the neighbour list does not contain the particle itself 
				int neiIndex = neiList1[index*maxNeiNum+i];
				printf("%f %f %f %f %f\n",m_pParticleData->m_vPositionX[neiIndex]-m_pParticleData->m_vPositionX[index],m_pParticleData->m_vPositionY[neiIndex]-m_pParticleData->m_vPositionY[index],inVolume[neiIndex]-inVolume[index],inVelocity[neiIndex]-inVelocity[index],inPressure[neiIndex]-inPressure[index]);
			}


		}
*/
		if(redo) {
			if(LPFOrder0[index]>0) LPFOrder0[index]--;
			if(LPFOrder1[index]>0) LPFOrder1[index]--;
			if(LPFOrder0[index]==0 && LPFOrder1[index]==0) {
				outVolume[index]   = inVolume[index];
				outPressure[index] = inPressure[index];
				outVelocity[index] = inVelocity[index];
				outSoundSpeed[index] = inSoundSpeed[index];	
			}
			else {
				index--;	
			}
		}
		else {
			if(outVolume[index]==0)
				outSoundSpeed[index] = 0;
			else
			{
				outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
//				if(m_vGamma[index])//threepoint
//				outSoundSpeed[index]=sqrt(1.5/1.4)*outSoundSpeed[index];
			}
		}
		//bool isInvalid = (outPressure[index]<m_fInvalidPressure || outVolume[index]==0 || 
		//                  1./outVolume[index]<=m_fInvalidDensity ||
		//				  std::isnan(outPressure[index]) || std::isnan(outVolume[index]) || 
		//				 std::isnan(outVelocity[index]));
	
		//if(isInvalid) {	
		//	bool atLeastOneNonzeroOrder = lowerLPFOrder(index,LPFOrderOther,LPFOrder0,LPFOrder1);	
		//	if(atLeastOneNonzeroOrder) { // then use the lowered order to recompute this partcle again
		//		index--;
		//		//continue;
		//	}
		//	else {
		//		//printInvalidState(phase,dir,index,positionX[index],positionY[index],positionZ[index], 
		//		//			  vel_d_0,vel_dd_0,p_d_0,p_dd_0,vel_d_1,vel_dd_1,p_d_1,p_dd_1); // input
		//		cout<<"invalid. phase="<<phase<<", index="<<index<<endl;
		//		phaseSuccess = false; // mark this phase as false
		//		//outVolume[index]   = inVolume[index];
		//		//outPressure[index] = inPressure[index];
		//		//outVelocity[index] = inVelocity[index];
		//	}
		//	continue;
		//}
		
		// Cavitation using critical pressure
		//if(m_iUseCriticalPressure) {
		//	if(outPressure[index] < m_fCriticalPressure) outPressure[index]=0;
		//}	

		// update outSoundSpeed
		//outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);

	}
}	

//double elipsedTime = omp_get_wtime() - startTime;
//printf("Directional splitting takes %.16g seconds\n", elipsedTime);	

//cout<<"----------------------------------------------------------------------"<<endl;

// if(phaseSuccess==true) ==> continue to next phase
// else                   ==> go back to phase 0 (all particles)
return phaseSuccess;
}

bool HyperbolicLPSolver::directionalSplitting_upwind(int phase) {
//      cout<<"--------------HyperbolicLPSolver::directionalSplitting()--------------"<<endl;

// determine dir: x(0), y(1), or z(2)
const int dir = m_vDirSplitTable[m_iDirSplitOrder][phase];
// set neighbour list pointers by dir (dir=0->right/left, dir=1->north/south, dir=2->up/down)
const int *neiList0=nullptr, *neiList1=nullptr;
const int *neiListSize0=nullptr, *neiListSize1=nullptr;
setNeighbourListPointers(dir, &neiList0, &neiList1, &neiListSize0, &neiListSize1);
// input data pointers
const double *inVelocity=nullptr, *inPressure=nullptr, *inVolume=nullptr, *inSoundSpeed=nullptr;
// output data pointers
double *outVelocity=nullptr, *outPressure=nullptr, *outVolume=nullptr, *outSoundSpeed=nullptr;
// set data pointers by phase (for temp1 or temp2) and dir(for velocity U or V) 
setInAndOutDataPointers(phase,dir,&inVelocity,&inPressure,&inVolume,&inSoundSpeed,
						&outVelocity,&outPressure,&outVolume,&outSoundSpeed);
// set local polynomail order pointers
// dir==0->right(0) & left(1), dir==1->north(0) & south(1), dir==2->up(0) & down(1)
int *LPFOrder0=nullptr, *LPFOrder1=nullptr;
vector<int*> LPFOrderOther;
setLPFOrderPointers(dir,&LPFOrder0,&LPFOrder1,LPFOrderOther);
// phase_success will be false if one particle LPF order is zero in one side
bool phaseSuccess = true;

// gravity is only on y(1) direction 
double gravity;
if(m_iDimension==2) gravity = dir==1? m_fGravity:0; // only in y direction
else if(m_iDimension==3) gravity = dir==1? m_fGravity:0; // only in z direction TODO: MODIFIED GRAVITY DIRECTION FOR POWDER TARGET SIMULATION

// set real dt:  
// 2D: phase:   0               1     2 
//                         dt/4   dt/2   dt/4
// 3D: phase:           0               1               2               3               4
//                                dt/6    dt/6     dt/3    dt/6    dt/6
double realDt;
if(m_iDimension==2) realDt = phase==1? m_fDt/2.:m_fDt/4.;
else if(m_iDimension==3) realDt = phase==2? m_fDt/3.:m_fDt/6.;
// set the function pointer to compute the A matrix for QR solver
void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*);
int offset; // offset is used to get results of computed spatial derivatives from QR solver
if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; offset=2;}
else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; offset=3;}
// the coeff before first and second order term during time integration
double multiplier1st, multiplier2nd;
if(m_iDimension==2) {multiplier1st=2; multiplier2nd=m_fDt/2.;}
else if(m_iDimension==3) {multiplier1st=3; multiplier2nd=m_fDt*3./4.;}
// iteration start index
size_t startIndex = m_pParticleData->m_iFluidStartIndex;
size_t endIndex = startIndex + m_pParticleData->m_iFluidNum;
int additional=0;
int warningcount=0;
//cout<<"m_pParticleData->getFluidNum() = "<<m_pParticleData->getFluidNum()<<endl;
//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
//double startTime = omp_get_wtime();   

// iterate through fluid particles
#ifdef _OPENMP
#pragma omp parallel for 
#endif
for(size_t index=startIndex; index<endIndex; index++) {

	if(inSoundSpeed[index]==0 || inVolume[index]==0) {

		cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
		<<", vol="<<inVolume[index]<<", vel="<<inVelocity[index]
		<<", cs="<<inSoundSpeed[index]<<endl;

		outVolume[index]   = inVolume[index];
		outPressure[index] = inPressure[index];
		outVelocity[index] = inVelocity[index];
		outSoundSpeed[index] = inSoundSpeed[index];
	}
	else {
		LPFOrder0[index]=1;
		LPFOrder1[index]=1;
		bool redo = false;

		// spatial derivatives (0:right/north/up; 1:left/south/down)
                        double vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1; // output
                        computeSpatialDer(dir, index, offset, computeA,//upwind
                                                          inPressure, inVelocity, neiList0, neiListSize0, additional,
                                                          LPFOrder0, &vel_d_0, &vel_dd_0, &p_d_0, &p_dd_0); // output
                        computeSpatialDer(dir, index, offset, computeA,//upwind
                                                          inPressure, inVelocity, neiList1, neiListSize1, additional,
                                                          LPFOrder1, &vel_d_1, &vel_dd_1, &p_d_1, &p_dd_1); // output

                        timeIntegration(realDt, multiplier1st, multiplier2nd,
                                                        gravity, inVolume[index], inVelocity[index], inPressure[index], inSoundSpeed[index],
                                                        vel_d_0, vel_dd_0, p_d_0, p_dd_0, vel_d_1, vel_dd_1, p_d_1, p_dd_1,
                                                        &outVolume[index], &outVelocity[index], &outPressure[index]); // output
			if(LPFOrder0[index]*LPFOrder1[index]==0 && warningcount++==0)
			{
				printf("Warning: order = %d %d for index = %ld, dir = %d. ",LPFOrder0[index],LPFOrder1[index],index,dir);
				printf("Neighbourlist size = %d %d.\n", neiListSize0[index], neiListSize1[index]);
			}
                        if(outPressure[index]<m_fInvalidPressure) redo=true;
                        else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
                        else if(std::isnan(outVolume[index]) || std::isnan(outVelocity[index]) || std::isnan(outPressure[index])) {
                                cout<<"nan value!!!"<<endl;
                                redo=true;
                        }
                        else if(std::isinf(outVolume[index]) || std::isinf(outVelocity[index]) || std::isinf(outPressure[index])) {
                                cout<<"inf value!!!"<<endl;
                                redo=true;
                        }
                        if(redo) {
                                if(LPFOrder0[index]>0) LPFOrder0[index]--;
                                if(LPFOrder1[index]>0) LPFOrder1[index]--;
                                if((LPFOrder0[index]==0 && LPFOrder1[index]==0)) {
                                        outVolume[index]   = inVolume[index];
                                        outPressure[index] = inPressure[index];
                                        outVelocity[index] = inVelocity[index];
                                        outSoundSpeed[index] = inSoundSpeed[index];
                                }
                                else {
                                        index--;
                                }
                        }
                        else {
                                if(outVolume[index]==0)
                                        outSoundSpeed[index] = 0;
                                else
                                        outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
                        }
                }
        }
	if(warningcount)
		printf("There are %d particles have zero order for direction = %d\n",warningcount,dir);
        return phaseSuccess;
}

bool HyperbolicLPSolver::nodirectionalSplitting() {
	cout<<"--------------HyperbolicLPSolver::nodirectionalSplitting()--------------"<<endl;
        // set neighbour list pointers	const int *neiList=m_pParticleData->m_vNeighbourList;
        const int *neiList=m_pParticleData->m_vNeighbourList;
	const int *neiListSize=m_pParticleData->m_vNeighbourListSize;
        // input data pointers
	const double *inVelocityU, *inVelocityV, *inVelocityW, *inPressure=m_pParticleData->m_vPressure, *inVolume=m_pParticleData->m_vVolume, *inSoundSpeed=m_pParticleData->m_vSoundSpeed;
	inVelocityU = m_pParticleData->m_vVelocityU;
        inVelocityV = m_pParticleData->m_vVelocityV;
	if(m_iDimension==3)
		inVelocityW = m_pParticleData->m_vVelocityW;
        // output data pointers
        double *outVelocityU, *outVelocityV, *outVelocityW, *outPressure=m_pParticleData->m_vTemp1Pressure, *outVolume=m_pParticleData->m_vTemp1Volume, *outSoundSpeed=m_pParticleData->m_vTemp1SoundSpeed;
        outVelocityU = m_pParticleData->m_vTemp1VelocityU;
        outVelocityV = m_pParticleData->m_vTemp1VelocityV;
        if(m_iDimension==3)
                outVelocityW = m_pParticleData->m_vTemp1VelocityW;
// set local polynomail order pointers
        int *LPFOrder=m_pParticleData->m_vLPFOrderRight;
// phase_success
        bool phaseSuccess = true;
// gravity is only on the last direction 
        double gravity=m_fGravity;
// set the function pointer to compute the A matrix for QR solver

        void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t, double*, double*);
	int number_of_derivative;
	size_t offset;
        if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; number_of_derivative=5; offset = 2;}
        else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; number_of_derivative=9; offset = 3;}
        // iteration start index
        size_t startIndex = m_pParticleData->m_iFluidStartIndex;
        size_t endIndex = startIndex + m_pParticleData->m_iFluidNum;
        // iterate through fluid particles
        #ifdef _OPENMP
        #pragma omp parallel for 
        #endif
        for(size_t index=startIndex; index<endIndex; index++) {
                if(inSoundSpeed[index]==0 || inVolume[index]==0) {

                        cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
                        <<", vol="<<inVolume[index]<<", velU="<<inVelocityU[index]<<", velV="<<inVelocityV[index]
                        <<", cs="<<inSoundSpeed[index]<<endl;

                        outVolume[index]   = inVolume[index];
                        outPressure[index] = inPressure[index];
                        outVelocityU[index] = inVelocityU[index];
                        outVelocityV[index] = inVelocityV[index];
                        if(m_iDimension==3) outVelocityW[index] = inVelocityW[index];
                        outSoundSpeed[index] = inSoundSpeed[index];
                }
                else {
			// spatial derivatives
			double Ud[number_of_derivative],Vd[number_of_derivative],Wd[number_of_derivative],Pd[number_of_derivative],Volumed[number_of_derivative];
                        computeSpatialDer(index, offset, computeA, inPressure, inVelocityU, inVelocityV, inVelocityW, inVolume, neiList, neiListSize,
                                                          LPFOrder, Pd, Ud, Vd, Wd, Volumed, number_of_derivative);//output
                        // update outVolume, outVelocityUVW, outPressure
/*			m_pParticleData->m_vPhi[index]=Volumed[0];
                        m_pParticleData->m_vPError0[index]=Volumed[1];
                        m_pParticleData->m_vPError1[index]=Volumed[2];
                        m_pParticleData->m_vVelError0[index]=Pd[3];
                        m_pParticleData->m_vVelError1[index]=Pd[4];*/
/*			if(m_pParticleData->m_vPositionX[index]>2)
			{
				for(int i=0;i<4;i++)
				{
					Pd[i]=0;
					Vd[i]=0;
					Ud[i]=0;
					Wd[i]=0;
					Volumed[i]=0;
				}
			}
*/
                        if(m_iDimension==3)	timeIntegration(m_fDt, gravity, inVolume[index], inVelocityU[index], inVelocityV[index], inVelocityW[index], inPressure[index], inSoundSpeed[index],
					Volumed, Ud, Vd, Wd, Pd,
                                                        &outVolume[index], &outVelocityU[index], &outVelocityV[index], &outVelocityW[index], &outPressure[index]); // output 
			else
			{
				double temp1=0,temp2;
				timeIntegration(m_fDt, gravity, inVolume[index], inVelocityU[index], inVelocityV[index], temp1, inPressure[index], inSoundSpeed[index],
                                        Volumed, Ud, Vd, Wd, Pd,
                                                        &outVolume[index], &outVelocityU[index], &outVelocityV[index], &temp2, &outPressure[index]); // output 
			}
                        if(LPFOrder[index]<2)
                        {
                                printf("Warning: reduce to first order for particle %ld, (x,y) = (%f %f), (V, U, P) = (%f %f %f).\n",index, m_pParticleData->m_vPositionX[index], m_pParticleData->m_vPositionY[index], inVolume[index],inVelocityU[index],inPressure[index]);
                                printf("Neighboursize= %d, LPFOrder= %d.\n", neiListSize[index], LPFOrder[index]);
                                printf("Output V U P = (%f %f %f)\n",outVolume[index],outVelocityU[index],outPressure[index]);
                        }

                        m_pParticleData->m_vVelError0[index]=Pd[3];
                        m_pParticleData->m_vVelError1[index]=Pd[4];
                        m_pParticleData->m_vPError0[index]=Pd[2];

                        bool redo = false;
                        if(outPressure[index]<m_fInvalidPressure) redo=true;
                        else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
                        else if(std::isnan(outVolume[index]) || std::isnan(outVelocityU[index]) || std::isnan(outVelocityV[index]) || std::isnan(outPressure[index])) {
                                cout<<"nan value!!!"<<endl;
                                redo=true;
                        }
                        else if(std::isinf(outVolume[index]) || std::isinf(outVelocityU[index]) || std::isinf(outVelocityV[index]) || std::isinf(outPressure[index])) {
                                cout<<"inf value!!!"<<endl;
                                redo=true;
                        }
                        if(redo) {
                                if(LPFOrder[index]>0) LPFOrder[index]--;
                                if(LPFOrder[index]==0) {
                                        outVolume[index]   = inVolume[index];
                                        outPressure[index] = inPressure[index];
                                        outVelocityU[index] = inVelocityU[index];
                                        outVelocityV[index] = inVelocityV[index];
					if(m_iDimension==3)	outVelocityW[index] = inVelocityW[index];
                                        outSoundSpeed[index] = inSoundSpeed[index];
                                }
                                else {
                                        index--;
                                }
                        }
                        else {
                                if(outVolume[index]==0)
                                        outSoundSpeed[index] = 0;
                                else
                                        outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
                        }
		}
	}
        cout<<"--------------END HyperbolicLPSolver::nodirectionalSplitting()--------------"<<endl;

        return phaseSuccess;
}

bool HyperbolicLPSolver::nodirectionalSplitting_noextreme() {
        cout<<"--------------HyperbolicLPSolver::nodirectionalSplitting()--------------"<<endl;
        // set neighbour list pointers  const int *neiList=m_pParticleData->m_vNeighbourList;
        const int *neiList=m_pParticleData->m_vNeighbourList;
        const int *neiListSize=m_pParticleData->m_vNeighbourListSize;
        // input data pointers
        const double *inVelocityU, *inVelocityV, *inVelocityW, *inPressure=m_pParticleData->m_vPressure, *inVolume=m_pParticleData->m_vVolume, *inSoundSpeed=m_pParticleData->m_vSoundSpeed;
        inVelocityU = m_pParticleData->m_vVelocityU;
        inVelocityV = m_pParticleData->m_vVelocityV;
        if(m_iDimension==3)
                inVelocityW = m_pParticleData->m_vVelocityW;
        // output data pointers
        double *outVelocityU, *outVelocityV, *outVelocityW, *outPressure=m_pParticleData->m_vTemp1Pressure, *outVolume=m_pParticleData->m_vTemp1Volume, *outSoundSpeed=m_pParticleData->m_vTemp1SoundSpeed;
        outVelocityU = m_pParticleData->m_vTemp1VelocityU;
        outVelocityV = m_pParticleData->m_vTemp1VelocityV;
        if(m_iDimension==3)
                outVelocityW = m_pParticleData->m_vTemp1VelocityW;
// set local polynomail order pointers
        int *LPFOrder=m_pParticleData->m_vLPFOrderRight;
// phase_success
        bool phaseSuccess = true;
// gravity is only on the last direction 
        double gravity=m_fGravity;
// set the function pointer to compute the A matrix for QR solver

        void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t, double*, double*);
        int number_of_derivative;
        size_t offset;
        if(m_iDimension==2) {computeA = &HyperbolicLPSolver::computeA2D; number_of_derivative=5; offset = 2;}
        else if(m_iDimension==3) {computeA = &HyperbolicLPSolver::computeA3D; number_of_derivative=9; offset = 3;}
        // iteration start index
        size_t startIndex = m_pParticleData->m_iFluidStartIndex;
        size_t endIndex = startIndex + m_pParticleData->m_iFluidNum;

        // iterate through fluid particles
        #ifdef _OPENMP
        #pragma omp parallel for 
        #endif
        for(size_t index=startIndex; index<endIndex; index++) {
//Test if the particle is a local extremal. If so, use upwind method instead
/*		m_pParticleData->m_vPhi[index]=1;
                int count_p=0,count_vel=0;
                int maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
                for(int i=0;i<8;i++)
                {
                	if(inPressure[neiList[index*maxNeiNum+i]]>inPressure[index])
                        	count_p++;
			double v1=inVelocityU[neiList[index*maxNeiNum+i]]*inVelocityU[neiList[index*maxNeiNum+i]]+inVelocityV[neiList[index*maxNeiNum+i]]*inVelocityV[neiList[index*maxNeiNum+i]];
			if(m_iDimension==3)
				v1=v1+inVelocityW[neiList[index*maxNeiNum+i]]*inVelocityW[neiList[index*maxNeiNum+i]];
			double v2=inVelocityU[index]*inVelocityU[index]+inVelocityV[index]*inVelocityV[index];
			if(m_iDimension==3)
				v2=v2+inVelocityW[index]*inVelocityW[index];
                        if(v1>v2)
                               	count_vel++;
                }
                if(((count_p==8)||(count_vel==8)||(count_p==0)||(count_vel==0)))//if local exterme, go to next index
                	continue;
*/
                m_pParticleData->m_vPhi[index]=1;
                int count_p=0,count_u=0,count_v=0,count_rho=0;
                int maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
                for(int i=0;i<8;i++)
                {
                        if(inPressure[neiList[index*maxNeiNum+i]]>inPressure[index])
                                count_p++;
                        if(inVelocityU[neiList[index*maxNeiNum+i]]>inVelocityU[index])
                                count_u++;
                        if(inVelocityV[neiList[index*maxNeiNum+i]]>inVelocityV[index])
                                count_v++;
                        if(inVolume[neiList[index*maxNeiNum+i]]>inVolume[index])
                                count_rho++;
                }
                if(((count_p==8)||(count_u==8)||(count_v==8)||(count_rho==8)||(count_p==0)||(count_u==0)||(count_v==0)||(count_rho==0)))//if local exterme, go to next index
                        continue;

//Test if the particle is a free boundary particle. If so, use upwind method instead.
		if(m_iFreeBoundary==true)
			if(m_vFillGhost[index]==true)
				continue;

                if(inSoundSpeed[index]==0 || inVolume[index]==0) {

                        cout<<"NOTICE!!!!!!! index="<<index<<", p="<<inPressure[index]
                        <<", vol="<<inVolume[index]<<", velU="<<inVelocityU[index]<<", velV="<<inVelocityV[index]
                        <<", cs="<<inSoundSpeed[index]<<endl;

                        outVolume[index]   = inVolume[index];
                        outPressure[index] = inPressure[index];
                        outVelocityU[index] = inVelocityU[index];
                        outVelocityV[index] = inVelocityV[index];
                        if(m_iDimension==3) outVelocityW[index] = inVelocityW[index];
                        outSoundSpeed[index] = inSoundSpeed[index];
                }
                else {
			LPFOrder[index]=2;
                        // spatial derivatives
                        double Ud[number_of_derivative],Vd[number_of_derivative],Wd[number_of_derivative],Pd[number_of_derivative],Volumed[number_of_derivative];
                        computeSpatialDer(index, offset, computeA, inPressure, inVelocityU, inVelocityV, inVelocityW, inVolume, neiList, neiListSize,
                                                          LPFOrder, Pd, Ud, Vd, Wd, Volumed, number_of_derivative);//output
			if(LPFOrder[index]<2)
			{
//				cout<<"Warning: first order used for particle "<<m_pParticleData->m_vPositionX[index]<<" "<<m_pParticleData->m_vPositionY[index]<<endl;
				continue;
			}

                        if(m_iDimension==3)     timeIntegration(index, m_fDt, gravity, inVolume[index], inVelocityU[index], inVelocityV[index], inVelocityW[index], inPressure[index], inSoundSpeed[index],
                                        Volumed, Ud, Vd, Wd, Pd,
                                                        &outVolume[index], &outVelocityU[index], &outVelocityV[index], &outVelocityW[index], &outPressure[index]); // output 
                        else
                        {
                                double temp1=0,temp2;
                                timeIntegration(index, m_fDt, gravity, inVolume[index], inVelocityU[index], inVelocityV[index], temp1, inPressure[index], inSoundSpeed[index],
                                        Volumed, Ud, Vd, Wd, Pd,
                                                        &outVolume[index], &outVelocityU[index], &outVelocityV[index], &temp2, &outPressure[index]); // output 
                        }
                        m_pParticleData->m_vPhi[index]=LPFOrder[index];

                        bool redo = false;
                        if(outPressure[index]<m_fInvalidPressure) redo=true;
                        else if(outVolume[index]!=0 && 1./outVolume[index]<m_fInvalidDensity) redo=true;
                        else if(std::isnan(outVolume[index]) || std::isnan(outVelocityU[index]) || std::isnan(outVelocityV[index]) || std::isnan(outPressure[index])) {
                                cout<<"nan value!!!"<<endl;
                                redo=true;
                        }
                        else if(std::isinf(outVolume[index]) || std::isinf(outVelocityU[index]) || std::isinf(outVelocityV[index]) || std::isinf(outPressure[index])) {
                                cout<<"inf value!!!"<<endl;
                                redo=true;
                        }
                        if(redo) {
                                outVolume[index]   = inVolume[index];
                                outPressure[index] = inPressure[index];
                                outVelocityU[index] = inVelocityU[index];
                                outVelocityV[index] = inVelocityV[index];
                                if(m_iDimension==3)     outVelocityW[index] = inVelocityW[index];
                                outSoundSpeed[index] = inSoundSpeed[index];
                        }
                        else {
                                if(outVolume[index]==0)
                                        outSoundSpeed[index] = 0;
                                else
                                        outSoundSpeed[index] = m_pEOS->getSoundSpeed(outPressure[index],1./outVolume[index]);
                        }
                }
        }
        cout<<"--------------END HyperbolicLPSolver::nodirectionalSplitting()--------------"<<endl;

        return phaseSuccess;
}





void HyperbolicLPSolver::setNeighbourListPointers(int dir, // input
	const int **neiList0, const int **neiList1, // output
	const int **neiListSize0, const int **neiListSize1) { 
	
	if(dir==0) { // x
		*neiList0 = m_pParticleData->m_vNeighbourListRight; 
		*neiList1 = m_pParticleData->m_vNeighbourListLeft;
		*neiListSize0 = m_pParticleData->m_vNeighbourListRightSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListLeftSize;
	}
	else if(dir==1) { // y
		*neiList0 = m_pParticleData->m_vNeighbourListNorth; 
		*neiList1 = m_pParticleData->m_vNeighbourListSouth;
		*neiListSize0 = m_pParticleData->m_vNeighbourListNorthSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListSouthSize;
	}
	else if(dir==2) { // z (if m_iDimension==2, dir != 2 for sure)
		*neiList0 = m_pParticleData->m_vNeighbourListUp; 
		*neiList1 = m_pParticleData->m_vNeighbourListDown;
		*neiListSize0 = m_pParticleData->m_vNeighbourListUpSize; 
		*neiListSize1 = m_pParticleData->m_vNeighbourListDownSize;	
	}
	else 
		assert(false);

}


void HyperbolicLPSolver::setInAndOutDataPointers(int phase, int dir,
	const double** inVelocity, const double** inPressure, const double** inVolume, const double** inSoundSpeed, 
	double** outVelocity, double** outPressure, double** outVolume, double** outSoundSpeed) {
	
	// assign pressure, volume, and sound_speed pointers
	if(phase==0) { // input: original output:temp1
		*inPressure   = m_pParticleData->m_vPressure;
		*inVolume     = m_pParticleData->m_vVolume;
		*inSoundSpeed = m_pParticleData->m_vSoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp1Pressure;
		*outVolume     = m_pParticleData->m_vTemp1Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;	
	}
	else if(phase==1 || phase==3) { // input:temp1 output:temp2
		*inPressure   = m_pParticleData->m_vTemp1Pressure;
		*inVolume     = m_pParticleData->m_vTemp1Volume;
		*inSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp2Pressure;
		*outVolume     = m_pParticleData->m_vTemp2Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp2SoundSpeed;	
	}
	else if(phase==2 || phase==4){ // input:temp2 output: temp1
		*inPressure   = m_pParticleData->m_vTemp2Pressure;
		*inVolume     = m_pParticleData->m_vTemp2Volume;
		*inSoundSpeed = m_pParticleData->m_vTemp2SoundSpeed;
		
		*outPressure   = m_pParticleData->m_vTemp1Pressure;
		*outVolume	  = m_pParticleData->m_vTemp1Volume;
		*outSoundSpeed = m_pParticleData->m_vTemp1SoundSpeed;	
	}	
	else assert(false);
	
	// assign velocity pointers
	if(m_iDimension==2) {
		if(phase==0 || phase==1) { // input: original output:temp2
			if(dir==0) {
				*inVelocity = m_pParticleData->m_vVelocityU;
				*outVelocity = m_pParticleData->m_vTemp2VelocityU;
			}
			else if(dir==1) {
				*inVelocity = m_pParticleData->m_vVelocityV;
				*outVelocity = m_pParticleData->m_vTemp2VelocityV;	
			}	
		}	
		else if(phase==2){ // input:temp2 output: temp1
			if(dir==0) { // u v u
				*inVelocity = m_pParticleData->m_vTemp2VelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);	
			}
			else if(dir==1) { // v u v
				*inVelocity = m_pParticleData->m_vTemp2VelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
			}		
		}	
		else assert(false);	
	}
	else if(m_iDimension==3) {
		if(phase==0 || phase==1 || phase==2) { // input: original output:temp1
			if(dir==0) {
				*inVelocity = m_pParticleData->m_vVelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
			}
			else if(dir==1) {
				*inVelocity = m_pParticleData->m_vVelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;	
			}
			else if(dir==2) {
				*inVelocity = m_pParticleData->m_vVelocityW;
				*outVelocity = m_pParticleData->m_vTemp1VelocityW;	
			}
		}	
		else if(phase==3){ // input:temp1 output: temp2
			if(dir==0) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityU;
				*outVelocity = m_pParticleData->m_vTemp2VelocityU;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==1) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityV;
				*outVelocity = m_pParticleData->m_vTemp2VelocityV;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==2) { 
				*inVelocity = m_pParticleData->m_vTemp1VelocityW;
				*outVelocity = m_pParticleData->m_vTemp2VelocityW;
				// swap pointers so that temp2 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
			}
		}
		else if(phase==4){ // input:temp2 output: temp1
			if(dir==0) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityU;
				*outVelocity = m_pParticleData->m_vTemp1VelocityU;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==1) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityV;
				*outVelocity = m_pParticleData->m_vTemp1VelocityV;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vTemp2VelocityW);
			}
			else if(dir==2) { 
				*inVelocity = m_pParticleData->m_vTemp2VelocityW;
				*outVelocity = m_pParticleData->m_vTemp1VelocityW;
				// swap pointers so that temp1 will contain new info
				swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vTemp2VelocityU);
				swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vTemp2VelocityV);
			}
		}
		else assert(false);	
	}

}


void HyperbolicLPSolver::setLPFOrderPointers(int dir, // input
	int** LPFOrder0, int** LPFOrder1, vector<int*>& LPFOrderOther) { // output
	
	if(dir==0) { // x
		*LPFOrder0 = m_pParticleData->m_vLPFOrderRight; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderLeft;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderNorth); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderSouth);
		if(m_iDimension==3) {
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderUp); 
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderDown);
		}
	}
	else if(dir==1) { // y
		*LPFOrder0 = m_pParticleData->m_vLPFOrderNorth; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderSouth;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderRight); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderLeft);
		if(m_iDimension==3) {
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderUp); 
			LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderDown);
		}
	}
	else if(dir==2) { // z
		*LPFOrder0 = m_pParticleData->m_vLPFOrderUp; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderDown;
		
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderRight); // other directions
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderLeft);
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderNorth); 
		LPFOrderOther.push_back(m_pParticleData->m_vLPFOrderSouth);
	}
	else
		assert(false);

}

void HyperbolicLPSolver::setLPFOrderPointers(int dir, // input
	int** LPFOrder0, int** LPFOrder1, int** LPFFirstOrder0, int** LPFFirstOrder1) { // output
	
	if(dir==0) { // x
		*LPFOrder0 = m_pParticleData->m_vLPFOrderRight; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderLeft;
		
		*LPFFirstOrder0 = m_pParticleData->m_vLPFFirstOrderRight; // this direction
		*LPFFirstOrder1 = m_pParticleData->m_vLPFFirstOrderLeft;
	}
	else if(dir==1) { // y
		*LPFOrder0 = m_pParticleData->m_vLPFOrderNorth; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderSouth;
		
		*LPFFirstOrder0 = m_pParticleData->m_vLPFFirstOrderNorth; // this direction
		*LPFFirstOrder1 = m_pParticleData->m_vLPFFirstOrderSouth;
	}
	else if(dir==2) { // z
		*LPFOrder0 = m_pParticleData->m_vLPFOrderUp; // this direction
		*LPFOrder1 = m_pParticleData->m_vLPFOrderDown;
		
		*LPFFirstOrder0 = m_pParticleData->m_vLPFFirstOrderUp; // this direction
		*LPFFirstOrder1 = m_pParticleData->m_vLPFFirstOrderDown;
	}
	else
		assert(false);

}

void HyperbolicLPSolver::computeSpatialDer(int dir, size_t index, // input 
	int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),
	const double* inPressure, const double* inVelocity,
	const int *neighbourList, const int *neighbourListSize,int additional,
	int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd) { // output
	
	//bool sufficientRank = false;

	// the initial value of the try-error process of finding a good A matrix with sufficient rank
	// these two variables will be increased if necessary in the while loop
	size_t numRow2nd = m_iNumRow2ndOrder+additional;
	size_t numRow1st = m_iNumRow1stOrder;
	double distance;	
	//cout<<"-------HyperbolicLPSolver::computeSpatialDer()-------"<<endl;
	//cout<<"numRow2nd="<<numRow2nd<<endl;
	//cout<<"numRow1st="<<numRow1st<<endl;

	while(true) {
	//while(!sufficientRank) {
		
		// decide row and column number based on LPFOrder (and numRow2nd, numRow1st) and the total number of neighbours
		size_t numRow, numCol;
		computeNumRowAndNumColAndLPFOrder(index, neighbourList, neighbourListSize, numRow2nd, numRow1st, // input
										  LPFOrder, &numRow, &numCol); // output	

		if(LPFOrder[index] == 0) { 
			*vel_d  = 0; *vel_dd = 0; *p_d  = 0; *p_dd = 0;
			//return;
			break;
		}
		
		// compute A
		double A[numRow*numCol];
		(this->*computeA)(index, neighbourList, LPFOrder, numRow, numCol, // input
						  A, &distance); // output
		double b[numRow]; 
		computeB(index, neighbourList, numRow, inPressure, // input: pressure first
				 b); // output

		QRSolver qrSolver(numRow,numCol,A);
		
		double result[numCol];
		int info = qrSolver.solve(result,b);

		if(info!=0) { // then need to recompute A
//                        cout<<"index="<<index<<", numRow="<<numRow<<", rank="<<info<<endl;
			if(LPFOrder[index]==2) numRow2nd++;
			else if(LPFOrder[index]==1) numRow1st++;	
		}	
		else {
			//sufficientRank = true;
			
			*p_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
			*p_dd = LPFOrder[index]==2? result[dir+offset]/distance/distance:0;

			computeB(index, neighbourList, numRow, inVelocity, // input: velocity comes second 
					 b); // output (rewrite array b)	
			
			qrSolver.solve(result,b);

			*vel_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
			*vel_dd = LPFOrder[index]==2? result[dir+offset]/distance/distance:0;
			
			if(std::isnan(*p_d) || std::isnan(*p_dd) || std::isnan(*vel_d) || std::isnan(*vel_dd) ||
			   std::isinf(*p_d) || std::isinf(*p_dd) || std::isinf(*vel_d) || std::isinf(*vel_dd)) {
				if(LPFOrder[index]==2) numRow2nd++;
				else if(LPFOrder[index]==1) numRow1st++;
			}
			else {
/*debug 0909*/
//				int direction=2*dir;
//				if(m_pParticleData->m_vNeighSize[4*index+direction+1]==-1)	direction++;
//				else m_pParticleData->m_vNeighSize[4*index+direction+1]=-1;

				{
//					m_pParticleData->m_vNeighSize[4*index+direction]=numRow;
					if(numRow>10)	printf("Notice: numRow=%zu for index = %zu\n",numRow,index);
//                	                size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
//					for(size_t i=0;i<numRow;i++)
//						m_pParticleData->m_vNeighList[40*index+10*direction+i]=neighbourList[index*maxNeiNum+i];
				}
/*
				if(fabs(m_pParticleData->m_vPositionY[index])<0.12&&fabs(m_pParticleData->m_vPositionX[index]-0.6)<0.6)
				{
					cout<<endl;
					cout<<"Number of neighbours for index "<<index<<" = "<<numRow<<endl;
					cout<<index<<" "<<m_pParticleData->m_vPositionX[index]<<" "<<m_pParticleData->m_vPositionY[index]<<" 0 0 "<<inPressure[index]<<" "<<inVelocity[index]<<" 0 0"<<endl;
					for(int i=0;i<numRow;i++)
					{
                                        	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
                                                int temp_index=neighbourList[index*maxNeiNum+i];
                                                cout<<temp_index<<" "<<m_pParticleData->m_vPositionX[temp_index]<<" "<<m_pParticleData->m_vPositionY[temp_index]<<" "<<m_pParticleData->m_vPositionX[temp_index]-m_pParticleData->m_vPositionX[index]<<" "<<m_pParticleData->m_vPositionY[temp_index]-m_pParticleData->m_vPositionY[index]<<" "<<inPressure[temp_index]<<" "<<inVelocity[temp_index]<<" "<<inPressure[temp_index]-inPressure[index]<<" "<<inVelocity[temp_index]-inVelocity[index]<<endl;
					}
                                        cout<<"1st derivative of p = "<<*p_d<<", 2nd derivative of p = "<<*p_dd<<endl;
                                        cout<<"1st derivative of vel = "<<*vel_d<<", 2nd derivative of vel = "<<*vel_dd<<endl;
				}*/
/*debug 0909*/
/*				if(index==3461){	//12/04 debug
					cout<<"Number of neighbours needed for 2nd order GFD: "<<numRow2nd<<endl;
					cout<<"Number of neighbours needed for 1st order GFD: "<<numRow1st<<endl;
					cout<<"Number of neighbours used::"<<numRow<<endl;
					cout<<"Neighbour list: index, x, y, pressure, velocity"<<endl;
					for(int i=0;i<numRow;i++)
						{
							size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
							int temp_index=neighbourList[index*maxNeiNum+i];
							cout<<temp_index<<" ("<<m_pParticleData->m_vPositionX[temp_index]<<", "<<m_pParticleData->m_vPositionY[temp_index]<<") "<<inPressure[temp_index]<<" "<<inVelocity[temp_index]<<endl;
						}
					cout<<"1st derivative of p = "<<*p_d<<", 2nd derivative of p = "<<*p_dd;
					cout<<"1st derivative of p = "<<*vel_d<<", 2nd derivative of p = "<<*vel_dd;
				}	//12/04 debug*/

				break;
			}
		}

	}
	
	//if(std::isnan(*p_d) || std::isnan(*p_dd) || std::isnan(*vel_d) || std::isnan(*vel_dd) ||
	//   std::isinf(*p_d) || std::isinf(*p_dd) || std::isinf(*vel_d) || std::isinf(*vel_dd))
	//   return false;
	//else
	//	return true;

	//cout<<"p_d="<<*p_d<<endl;
	//cout<<"p_dd="<<*p_dd<<endl;
	//cout<<"vel_d="<<*vel_d<<endl;
	//cout<<"vel_dd="<<*vel_dd<<endl;
	//cout<<"-----------------------------------------------------"<<endl;

}


void HyperbolicLPSolver::computeSpatialDer(int dir, size_t index, // input 
        int offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),
        const double* inPressure, const double* inVelocity,
        const int *neighbourList, const int *neighbourListSize,int additional,
        int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd, double* vel_e, double* p_e) { // output
	printf("Result is wrong! void HyperbolicLPSolver::computeSpatialDer(int dir, size_t index, int offset, double (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*), const double* inPressure, const double* inVelocity, const int *neighbourList, const int *neighbourListSize,int additional, int* LPFOrder, double* vel_d, double* vel_dd, double* p_d, double* p_dd, double* vel_e, double* p_e) needs to be updated!\n");
        //bool sufficientRank = false;

        // the initial value of the try-error process of finding a good A matrix with sufficient rank
        // these two variables will be increased if necessary in the while loop
        size_t numRow2nd = m_iNumRow2ndOrder+additional;
        size_t numRow1st = m_iNumRow1stOrder+additional;
	double distance;
        //cout<<"-------HyperbolicLPSolver::computeSpatialDer()-------"<<endl;
        //cout<<"numRow2nd="<<numRow2nd<<endl;
        //cout<<"numRow1st="<<numRow1st<<endl;

        while(true) {
        //while(!sufficientRank) {

                // decide row and column number based on LPFOrder (and numRow2nd, numRow1st) and the total number of neighbours
                size_t numRow, numCol;
                computeNumRowAndNumColAndLPFOrder(index, neighbourList, neighbourListSize, numRow2nd, numRow1st, // input
                                                                                  LPFOrder, &numRow, &numCol); // output        

                if(LPFOrder[index] == 0) {
                        *vel_d  = 0; *vel_dd = 0; *p_d  = 0; *p_dd = 0;
                        //return;
                        break;
                }

                // compute A
                double A[numRow*numCol];
                (this->*computeA)(index, neighbourList, LPFOrder, numRow, numCol, // input
                                                  A, &distance); // output
                double b[numRow];
                computeB(index, neighbourList, numRow, inPressure, // input: pressure first
                                 b); // output
                QRSolver qrSolver(numRow,numCol,A);

                double result[numCol];
                int info = qrSolver.solve(result,b);

                if(info!=0) { // then need to recompute A
                        //cout<<"rank="<<info<<endl;
                        if(LPFOrder[index]==2) numRow2nd++;
                        else if(LPFOrder[index]==1) numRow1st++;
                }
                else {
                        //sufficientRank = true;
			double result_p[numCol];
			for (size_t i=0;i<numCol; i++)
				result_p[i]=result[i];
                        *p_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
                        *p_dd = LPFOrder[index]==2? result[dir+offset]/distance/distance:0;

                        computeB(index, neighbourList, numRow, inVelocity, // input: velocity comes second 
                                         b); // output (rewrite array b)        

                        qrSolver.solve(result,b);

                        double result_vel[numCol];
                        for (size_t i=0;i<numCol; i++)
                                result_vel[i]=result[i];

                        *vel_d = result[dir]/distance; //dir=0 (x), dir=1(y), dir=2(z)
                        *vel_dd = LPFOrder[index]==2? result[dir+offset]/distance/distance:0;

                        if(std::isnan(*p_d) || std::isnan(*p_dd) || std::isnan(*vel_d) || std::isnan(*vel_dd) ||
                           std::isinf(*p_d) || std::isinf(*p_dd) || std::isinf(*vel_d) || std::isinf(*vel_dd)) {
                                if(LPFOrder[index]==2) numRow2nd++;
                                else if(LPFOrder[index]==1) numRow1st++;
                        }
                        else {
				int testRow=numRow+3;
				if(testRow>neighbourListSize[index])
					testRow=neighbourListSize[index];
		                double A1[testRow*numCol];
		                double b1[testRow];
		                (this->*computeA)(index, neighbourList, LPFOrder, testRow, numCol, // input
                                                  A1, &distance); // output
                		computeB(index, neighbourList, testRow, inPressure, // input: pressure first
                                		b1); // output
				compute_GFD_error(testRow,numCol,A1,b1,result_p,p_e);				
                                computeB(index, neighbourList, testRow, inVelocity, // input: pressure first
                                                b1); // output
                                compute_GFD_error(testRow,numCol,A1,b1,result_vel,vel_e);

                                break;
                        }
                }

        }


}

void HyperbolicLPSolver::computeSpatialDer(size_t index,  size_t offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),  const double* inPressure, const double* inVelocityU, const double* inVelocityV, const double* inVelocityW, const double* inVolume, const int *neighbourList, const int *neighbourListSize,
                                                          int *LPFOrder, double* Pd, double* Ud, double* Vd, double* Wd, double* Volumed, int number_of_derivative)
{
//	double threshold=0.1;
	double dis,first,second;
	double max_dis=0;
        for(int i=0;i<number_of_derivative;i++)
        {
        	Pd[i]=0.0;
                Ud[i]=0.0;
                Vd[i]=0.0;
                Wd[i]=0.0;
		Volumed[i]=0.0;
        }
	double distance;
        size_t numRow2nd = m_iNumRow2ndOrder;
        size_t numRow1st = m_iNumRow1stOrder;
        while(true) {
                size_t numRow, numCol;
                computeNumRowAndNumColAndLPFOrder(index, neighbourList, neighbourListSize, numRow2nd, numRow1st, // input

  LPFOrder, &numRow, &numCol); // output  
                if(LPFOrder[index] == 0) {
                        break;
                }

                double A[numRow*numCol];
                (this->*computeA)(index, neighbourList, LPFOrder, numRow, numCol, // input
                                                  A, &distance); // output

                double b[numRow];
                computeB(index, neighbourList, numRow, inPressure, // input: pressure first
                                 b); // output

                QRSolver qrSolver(numRow,numCol,A);

                double result[numCol];
                int info = qrSolver.solve(result,b);
                if(info!=0) { // then need to recompute A
                        //cout<<"index="<<index<<", numRow="<<numRow<<", rank="<<info<<endl;
                        if(LPFOrder[index]==2) numRow2nd++;
                        else if(LPFOrder[index]==1) numRow1st++;
                }
                else {

			first=0;
			second=0;
                        for(size_t i=0;i<numCol;i++)
			{
				if (i<offset)
				{
                                	Pd[i]=result[i]/distance;
					first=first+result[i];
				}
				else
				{
					Pd[i]=result[i]/distance/distance;
					second=second+result[i];
				}
			}
			dis=second/first;
			if (dis>max_dis)	max_dis=dis;

                        computeB(index, neighbourList, numRow, inVelocityU, // input: velocity comes second 
                                         b); // output (rewrite array b)        

                        qrSolver.solve(result,b);
                        first=0;
                        second=0;
                        for(size_t i=0;i<numCol;i++)
                        {
                                if (i<offset)
                                {
                                        Ud[i]=result[i]/distance;
                                        first=first+result[i];
                                }
                                else
                                {
                                        Ud[i]=result[i]/distance/distance;
                                        second=second+result[i];
                                }
                        }
                        dis=second/first;
                        if (dis>max_dis)        max_dis=dis;
                        computeB(index, neighbourList, numRow, inVelocityV, // input: velocity comes second 
                                         b); // output (rewrite array b)        

                        qrSolver.solve(result,b);
                        first=0;
                        second=0;
                        for(size_t i=0;i<numCol;i++)
                        {
                                if (i<offset)
                                {
                                        Vd[i]=result[i]/distance;
                                        first=first+result[i];
                                }
                                else
                                {
                                        Vd[i]=result[i]/distance/distance;
                                        second=second+result[i];
                                }
                        }
                        dis=second/first;
                        if (dis>max_dis)        max_dis=dis;

			if(m_iDimension==3)
			{
	                        computeB(index, neighbourList, numRow, inVelocityW, // input: velocity comes second 
        	                                 b); // output (rewrite array b)        
	
        	                qrSolver.solve(result,b);
            	        	first=0;
               	        	second=0;
                        	for(size_t i=0;i<numCol;i++)
                        	{
                                	if (i<offset)
                                	{
                                        	Wd[i]=result[i]/distance;
                                        	first=first+result[i];
                                	}
                                	else
                                	{
                                        	Wd[i]=result[i]/distance/distance;
                                        	second=second+result[i];
                                	}
                        	}
                        	dis=second/first;
                        	if (dis>max_dis)        max_dis=dis;
			}
                        computeB(index, neighbourList, numRow, inVolume, // input: volume comes last 
                                         b); // output (rewrite array b)        

                        qrSolver.solve(result,b);
                        first=0;
                        second=0;
                        for(size_t i=0;i<numCol;i++)
                        {
                                if (i<offset)
                                {
                                        Volumed[i]=result[i]/distance;
                                        first=first+result[i];
                                }
                                else
                                {
                                        Volumed[i]=result[i]/distance/distance;
                                        second=second+result[i];
                                }
                        }
                        dis=second/first;
                        if (dis>max_dis)        max_dis=dis;
//			if (max_dis>threshold)	LPFOrder[index]=1;
			int fail=0;
			for(int i=0; i<number_of_derivative;i++)
	                        if(std::isnan(Pd[i]) || std::isnan(Ud[i]) || std::isnan(Vd[i]) || std::isnan(Wd[i]) || std::isnan(Volumed[i]) ||
                           std::isinf(Pd[i]) || std::isinf(Ud[i]) || std::isinf(Vd[i]) || std::isinf(Wd[i])|| std::isnan(Volumed[i])) {
				fail=fail+1;
				}
			if(fail)
			{
                                if(LPFOrder[index]==2) numRow2nd++;
                                else if(LPFOrder[index]==1) numRow1st++;
			}
			else
			{
				break;
			}
		}
	}
}
//support function for density_derivative
void HyperbolicLPSolver::computeSpatialDer(size_t index, size_t offset, void (HyperbolicLPSolver::*computeA) (size_t, const int *, const int*, size_t, size_t,double*, double*),
                                                  const double* inVolume, const int *neighbourList, const int *neighbourListSize,
                                                  int *LPFOrder, double* volume_x, double* volume_y, double* volume_z)
{
	*volume_x=0;
	*volume_y=0;
	*volume_z=0;
	int order=LPFOrder[index];
	LPFOrder[index]=1;
        double distance;
        size_t numRow2nd = 4*m_iNumRow2ndOrder;
        size_t numRow1st = 4*m_iNumRow1stOrder;
	if (neighbourListSize[index]>(int)numRow1st)
		numRow1st=neighbourListSize[index];
        while(true) {
                size_t numRow, numCol;
                computeNumRowAndNumColAndLPFOrder(index, neighbourList, neighbourListSize, numRow2nd, numRow1st, // input

  LPFOrder, &numRow, &numCol); // output 

		if(LPFOrder[index]==0)
		{
			LPFOrder[index]=order;
			printf("Warning: cannot calculate derivative of volume for %zu th particle. SPH density estimator may be inaccurate.\n",index);
			break;
		}
                double A[numRow*numCol];
                (this->*computeA)(index, neighbourList, LPFOrder, numRow, numCol, // input
                                                  A, &distance); // output
                double b[numRow];
                computeB(index, neighbourList, numRow, inVolume, // input
                                 b); // output

                QRSolver qrSolver(numRow,numCol,A);

                double result[numCol];
                int info = qrSolver.solve(result,b);
                if(info!=0) { // then need to recompute A
                        //cout<<"index="<<index<<", numRow="<<numRow<<", rank="<<info<<endl;
                        if(LPFOrder[index]==2) numRow2nd++;
                        else if(LPFOrder[index]==1) numRow1st++;
                }
                else {
			*volume_x=result[0];
			*volume_y=result[1];
			if(offset==3)
				*volume_z=result[2];
			break;
		}
	}
	LPFOrder[index]=order;
/*	if(fabs(m_pParticleData->m_vPositionX[index])<0.05)
	{
//                cout<<m_pParticleData->m_iMaxNeighbourNum<<endl;
//		cout<<m_pParticleData->m_iMaxNeighbourNumInOneDir<<endl;
		cout<<endl;
		cout<<m_pParticleData->m_vPositionX[index]<<" "<<m_pParticleData->m_vPositionY[index]<<" "<<m_pParticleData->m_vPositionZ[index]<<" "<<inVolume[index]<<endl;;
		cout<<*volume_x<<" "<<*volume_y<<" "<<*volume_z<<endl;
		cout<<4*m_iNumRow1stOrder<<" "<<numRow1st<<endl;
		for(int j=0;j<numRow1st;j++)
		{
			int k=neighbourList[m_pParticleData->m_iMaxNeighbourNumInOneDir*index+j];
                cout<<k<<" "<<m_pParticleData->m_vPositionX[k]-m_pParticleData->m_vPositionX[index]<<" "<<m_pParticleData->m_vPositionY[k]-m_pParticleData->m_vPositionY[index]<<" "<<m_pParticleData->m_vPositionZ[k]-m_pParticleData->m_vPositionZ[index]<<" "<<inVolume[k]-inVolume[index]<<endl;

			}
	}*/
}

void HyperbolicLPSolver::compute_GFD_error(size_t numRow, size_t numCol, double *A, double *b, double *x, double *e)
{
	*e=0.0;
	double temp,e_temp;
	for(size_t i=0;i<numRow;i++)
	{
		e_temp=-b[i];
		for(size_t j=0;j<numCol;j++)
			e_temp=e_temp+x[j]*A[i+j*numRow];
		*e=*e+fabs(e_temp);
		temp=temp+fabs(b[i]);
	}
	if(temp<0.00001)
		*e=0.0;
	else
		*e=*e/temp;
//        if(std::isnan(*e) || std::isinf(*e)) {
//                printf("%d %d %f %f %f\n",numRow, numCol, e_temp,temp,(*e));
//
//                assert(false);
//        }
	
}

void HyperbolicLPSolver::computeNumRowAndNumColAndLPFOrder(size_t index, // input
    const int *neighbourList, const int *neighbourListSize, size_t numRow2nd, size_t numRow1st,
	int* LPFOrder, size_t *numRow, size_t *numCol) { // output
	
	// compute the numRow and numCol for matrix A
	size_t totalNeiNum = neighbourListSize[index];

	if(LPFOrder[index] == 2) {
		if(totalNeiNum >= numRow2nd) { // if the total number of neighbour >= current numRow2nd
			(*numRow) = numRow2nd;
			(*numCol) = m_iNumCol2ndOrder; 
		}
		else LPFOrder[index] = 1; // if no enough neighbour -> reduce LPFOrder to 1
	}

	if(LPFOrder[index] == 1) { 
		if(totalNeiNum >= numRow1st) { // if the total number of neighbour >= current numRow1st
			(*numRow) = numRow1st;
			(*numCol) = m_iNumCol1stOrder; 
		}
		else LPFOrder[index] = 0;
	}

	if(LPFOrder[index] == 0) {
		(*numRow) = 0;
		(*numCol) = 0;
	}

	//cout<<"-------HyperbolicLPSolver::computeNumRowAndNumColAndLPFOrder()-------"<<endl;
	//cout<<"index="<<index<<endl;
	//cout<<"totalNeiNum="<<totalNeiNum<<endl;
	//cout<<"LPFOrder="<<LPFOrder[index]<<endl;
	//cout<<"numRow="<<*numRow<<endl;
	//cout<<"numCol="<<*numCol<<endl;
	//cout<<"---------------------------------------------------------------------"<<endl;

}


void HyperbolicLPSolver::computeA2D(size_t index, const int *neighbourList, 
							  	    const int* LPFOrder, size_t numRow, size_t numCol,
								    double *A, double *dis) { // output 	
	//cout<<"--------------HyperbolicLPSolver::computeA2D()--------------"<<endl;
	//cout<<"index="<<index<<endl;
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
        double distance = sqrt((m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) * (m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) + (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]) * (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]));
	if(LPFOrder[index] == 1) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself 
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
			double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;

			A[i]            = h;
			A[i + 1*numRow] = k;	

			//cout<<A[i]<<"	"<<A[i + 1*numRow]<<endl;
		}
	}
	else if(LPFOrder[index] == 2) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
			double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;
			
			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow]	= 0.5*h*h;
			A[i + 3*numRow]	= 0.5*k*k;
			A[i + 4*numRow]	= h*k;
			

			//cout<<A[i]<<"	"<<A[i + 1*numRow]<<"	"<<A[i + 2*numRow]<<"	"
			//	<<A[i + 3*numRow]<<"	"<<A[i + 4*numRow]<<endl;
		}
	}	
	(*dis)=distance;	
	//cout<<"------------------------------------------------------------"<<endl;
}

/*void HyperbolicLPSolver::computeA2D(size_t index, const int *neighbourList,
                                                                    const int* LPFOrder, size_t numRow, size_t numCol,
                                                                    double *A, double *dis) { // output         
        //cout<<"--------------HyperbolicLPSolver::computeA2D()--------------"<<endl;
        //cout<<"index="<<index<<endl;

        size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
        double distance = sqrt((m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) * (m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) + (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]) * (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]));
        if(LPFOrder[index] == 1) {
                for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself 
                        int neiIndex = neighbourList[index*maxNeiNum+i];

                        double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
                        double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;

                        double r=sqrt(h*h+k*k)*distance;
			if(r<0.02)	r=0.02;
                        double weight=1/r;

                        A[i]            = h*weight;
                        A[i + 1*numRow] = k*weight;

                        //cout<<A[i]<<" "<<A[i + 1*numRow]<<endl;
                }
        }
        else if(LPFOrder[index] == 2) {
                for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
                        int neiIndex = neighbourList[index*maxNeiNum+i];
			
                        double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
                        double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;

			double r=sqrt(h*h+k*k)*distance;
                        if(r<0.02)      r=0.02;
			double weight=1/r;

                        A[i]            = h*weight;
                        A[i + 1*numRow] = k*weight;
                        A[i + 2*numRow] = 0.5*h*h*weight;
                        A[i + 3*numRow] = 0.5*k*k*weight;
                        A[i + 4*numRow] = h*k*weight;
                }
        }
        (*dis)=distance;
        //cout<<"------------------------------------------------------------"<<endl;
}*/


void HyperbolicLPSolver::computeA3D(size_t index, const int *neighbourList, 
								    const int* LPFOrder, size_t numRow, size_t numCol,
								    double *A, double *dis) { // output 	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
        double distance = sqrt((m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionX[index]) * (m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionX[index]) + (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionY[index]) * (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionY[index]) + (m_pParticleData->m_vPositionZ[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionZ[index]) * (m_pParticleData->m_vPositionZ[neighbourList[index*maxNeiNum]] - m_pParticleData->m_vPositionZ[index]));
	
	if(LPFOrder[index] == 1) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
			double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;
			double l = (m_pParticleData->m_vPositionZ[neiIndex] - m_pParticleData->m_vPositionZ[index])/distance;

			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow] = l;
		}
	}
	else if(LPFOrder[index] == 2) {
		for(size_t i=0; i<numRow; i++) { // Note that the neighbour list does not contain the particle itself
			int neiIndex = neighbourList[index*maxNeiNum+i];	
				
			double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
			double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;
			double l = (m_pParticleData->m_vPositionZ[neiIndex] - m_pParticleData->m_vPositionZ[index])/distance;

			A[i]            = h;
			A[i + 1*numRow] = k;
			A[i + 2*numRow] = l;
			A[i + 3*numRow] = 0.5*h*h;
			A[i + 4*numRow] = 0.5*k*k;
			A[i + 5*numRow] = 0.5*l*l;
			A[i + 6*numRow] = h*k;
			A[i + 7*numRow] = h*l;
			A[i + 8*numRow] = k*l;
		}
	}
	(*dis)=distance;
}


void HyperbolicLPSolver::computeB(size_t index, const int *neighbourList, size_t numRow, const double* inData, 
								  double *b) { // output 	
	
	//cout<<"--------------HyperbolicLPSolver::computeB()--------------"<<endl;
	//cout<<"index="<<index<<endl;

	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
	for(size_t i=0; i<numRow; i++) { 
		int neiIndex = neighbourList[index*maxNeiNum+i];	
		b[i] = inData[neiIndex] - inData[index];

		//cout<<b[i]<<endl;
	}	

	//cout<<"----------------------------------------------------------"<<endl;

}

/*void HyperbolicLPSolver::computeB(size_t index, const int *neighbourList, size_t numRow, const double* inData,
                                                                  double *b) { // output        

        //cout<<"--------------HyperbolicLPSolver::computeB()--------------"<<endl;
        //cout<<"index="<<index<<endl;

        size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNumInOneDir;
        double distance = sqrt((m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) * (m_pParticleData->m_vPositionX[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionX[index]) + (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]) * (m_pParticleData->m_vPositionY[neighbourList[index*maxNeiNum+numRow/4]] - m_pParticleData->m_vPositionY[index]));
        for(size_t i=0; i<numRow; i++) {
                int neiIndex = neighbourList[index*maxNeiNum+i];

                double h = (m_pParticleData->m_vPositionX[neiIndex] - m_pParticleData->m_vPositionX[index])/distance;
                double k = (m_pParticleData->m_vPositionY[neiIndex] - m_pParticleData->m_vPositionY[index])/distance;

                double r=sqrt(h*h+k*k)*distance;
                        if(r<0.02)      r=0.02;

		double weight=1/r;

                b[i] = (inData[neiIndex] - inData[index])*weight;

                //cout<<b[i]<<endl;
        }

        //cout<<"----------------------------------------------------------"<<endl;

}*/




void HyperbolicLPSolver::timeIntegration_van_leer(
	double realDt, double multiplier1st, double multiplier2nd, 
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double phi,
	double vel_d_0_first, double vel_d_1_first, double p_d_0_first, double p_d_1_first,
	double vel_d_0_second, double vel_d_1_second, double p_d_0_second, double p_d_1_second,
	double vel_dd_0_second, double vel_dd_1_second, double p_dd_0_second, double p_dd_1_second,
	double* outVolume, double* outVelocity, double* outPressure) { // output
	
	// TODO Note that this coeff K only works for Poly and Spoly EOS!!!!!!!
	double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume; 

	double flux_first = 0, flux_second = 0; 
	
	// P (first order flux)
	double Pt1st = -0.5*inVolume*K*(vel_d_0_first+vel_d_1_first) + 0.5*inVolume*sqrt(K)*(p_d_0_first-p_d_1_first);
	flux_first = multiplier1st*Pt1st;
	// P (second order flux)
	Pt1st = -0.5*inVolume*K*(vel_d_0_second+vel_d_1_second) + 0.5*inVolume*sqrt(K)*(p_d_0_second-p_d_1_second);
	double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0_second-vel_dd_1_second) 
	+ inVolume*inVolume*K*(p_dd_0_second+p_dd_1_second);
	flux_second = multiplier1st*Pt1st + multiplier2nd*Pt2nd; 
	// P flux (combined)
	double flux_p = flux_first + phi * (flux_second - flux_first);
	
	//Vt
	double flux_V = -flux_p/K;
	

	// VEL (first order flux) 
	double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0_first-vel_d_1_first) - 0.5*inVolume*(p_d_0_first+p_d_1_first);
	flux_first = multiplier1st*VELt1st;
	// VEL (second order flux)
	VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0_second-vel_d_1_second) - 0.5*inVolume*(p_d_0_second+p_d_1_second);
	double VELt2nd = inVolume*inVolume*K*(vel_dd_0_second+vel_dd_1_second) - 
	inVolume*inVolume*sqrt(K)*(p_dd_0_second-p_dd_1_second);
	flux_second = multiplier1st*VELt1st + multiplier2nd*VELt2nd;
	// VEL flux (combined)
	double flux_VEL = flux_first + phi * (flux_second - flux_first); 		

	// Note that the data pointers of in and out are different!!!!!!!
	(*outVolume)   = inVolume   + realDt*flux_V;
	(*outPressure) = inPressure + realDt*flux_p;
	(*outVelocity) = inVelocity + realDt*(flux_VEL+gravity);	
	
	if(std::isnan(*outVolume) || std::isinf(*outVolume) || 
	   std::isnan(*outPressure) || std::isinf(*outPressure) ||
	   std::isnan(*outVelocity) || std::isinf(*outVelocity)) {
		printf("%f %f %f %f\n",inVelocity, flux_first, flux_second, phi);
		printf("%f %f %f\n",(*outVolume),(*outPressure),(*outVelocity));

		assert(false);   
	}

}


void HyperbolicLPSolver::timeIntegration(
	double realDt, double multiplier1st, double multiplier2nd, 
	double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
	double* outVolume, double* outVelocity, double* outPressure) { // output
	
	// TODO Note that this coeff K only works for Poly and Spoly EOS!!!!!!!
	double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume; 
	//double K;
	//if(inSoundSpeed==0 && inVloume==0) K = 1;
	//else K = inSoundSpeed*inSoundSpeed/inVolume/inVolume;
	
	// Pt
	double Pt1st = -0.5*inVolume*K*(vel_d_0+vel_d_1) + 0.5*inVolume*sqrt(K)*(p_d_0-p_d_1);
	double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0-vel_dd_1) + inVolume*inVolume*K*(p_dd_0+p_dd_1);
	double Pt = multiplier1st*Pt1st + multiplier2nd*Pt2nd;
	
	// Vt
	double Vt = -Pt/K;
	

	// VELt
	double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0-vel_d_1) - 0.5*inVolume*(p_d_0+p_d_1);
	double VELt2nd = inVolume*inVolume*K*(vel_dd_0+vel_dd_1) - inVolume*inVolume*sqrt(K)*(p_dd_0-p_dd_1);
	double VELt = multiplier1st*VELt1st + multiplier2nd*VELt2nd;

	// Note that the data pointers of in and out are different!!!!!!!
	(*outVolume)   = inVolume   + realDt*Vt;
	(*outPressure) = inPressure + realDt*Pt;
	(*outVelocity) = inVelocity + realDt*(VELt+gravity);	
	
	if(std::isnan(*outVolume) || std::isinf(*outVolume) || 
	   std::isnan(*outPressure) || std::isinf(*outPressure) ||
	   std::isnan(*outVelocity) || std::isinf(*outVelocity)) {
		assert(false);   
	}

}

void HyperbolicLPSolver::timeIntegration(
        double realDt, double multiplier1st, double multiplier2nd, 
        double gravity, double inVolume, double inVelocity, double inPressure, double inSoundSpeed,
        double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0,
        double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1,
        double* outVolume, double* outVelocity, double* outPressure, double* volume_t, double* velocity_t, double* pressure_t) { // output

        double K = inSoundSpeed*inSoundSpeed/inVolume/inVolume;

        double Pt1st = -0.5*inVolume*K*(vel_d_0+vel_d_1) + 0.5*inVolume*sqrt(K)*(p_d_0-p_d_1);
        double Pt2nd = -inVolume*inVolume*pow(K,1.5)*(vel_dd_0-vel_dd_1) + inVolume*inVolume*K*(p_dd_0+p_dd_1);
        double Pt = multiplier1st*Pt1st + multiplier2nd*Pt2nd;

        double Vt = -Pt/K;


        double VELt1st = 0.5*inVolume*sqrt(K)*(vel_d_0-vel_d_1) - 0.5*inVolume*(p_d_0+p_d_1);
        double VELt2nd = inVolume*inVolume*K*(vel_dd_0+vel_dd_1) - inVolume*inVolume*sqrt(K)*(p_dd_0-p_dd_1);
        double VELt = multiplier1st*VELt1st + multiplier2nd*VELt2nd;

        (*volume_t)=Vt;
        (*velocity_t)=VELt;
        (*pressure_t)=Pt;

        (*outVolume)   = inVolume   + realDt*Vt;
        (*outPressure) = inPressure + realDt*Pt;
        (*outVelocity) = inVelocity + realDt*(VELt+gravity);

        if(std::isnan(*outVolume) || std::isinf(*outVolume) ||
           std::isnan(*outPressure) || std::isinf(*outPressure) ||
           std::isnan(*outVelocity) || std::isinf(*outVelocity)) {
                assert(false);
        }

}
void HyperbolicLPSolver::timeIntegration(double Dt, double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure){
//	double gamma=inSoundSpeed*inSoundSpeed/inVolume/inPressure;
	double gamma = m_pGamma;
	double Pinf=m_pPinf;

        // TODO Note that this functions only works for Poly EOS!!!!!!!
	if(m_iDimension==3)
	{
		double div=Ud[0]+Vd[1]+Wd[2];
		double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+Wd[2]*Wd[2]+2*Ud[1]*Vd[0]+2*Ud[2]*Wd[0]+2*Vd[2]*Wd[1];
		double Volumet=inVolume*div;
		double VelocityUt=-inVolume*Pd[0];
		double VelocityVt=-inVolume*Pd[1];
		double VelocityWt=-inVolume*Pd[2];
		double Pt=-gamma*(inPressure+Pinf)*div;
		double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-Volumed[2]*Pd[2]-inVolume*(Pd[3]+Pd[4]+Pd[5])-cross);
		double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[3]+Vd[6]+Wd[7])+Ud[0]*Pd[0]+Vd[0]*Pd[1]+Wd[0]*Pd[2]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[6]+Vd[4]+Wd[8])+Ud[1]*Pd[0]+Vd[1]*Pd[1]+Wd[1]*Pd[2]);
                double VelocityWtt=inVolume*((gamma-1)*Pd[2]*div+gamma*(inPressure+Pinf)*(Ud[7]+Vd[8]+Wd[5])+Ud[2]*Pd[0]+Vd[2]*Pd[1]+Wd[2]*Pd[2]);
		double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+Volumed[2]*Pd[2]+inVolume*(Pd[3]+Pd[4]+Pd[5])+cross);
		(*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;

                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;//TODO: MODIFIED GRAVITY DIRECTION
//                (*outVelocityW)=inVelocityW+Dt*VelocityWt+0.5*Dt*Dt*VelocityWtt;
                (*outVelocityW)=inVelocityW+Dt*VelocityWt+0.5*Dt*Dt*VelocityWtt;
                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
	        if(std::isnan(*outVolume) || std::isinf(*outVolume) ||
        	   std::isnan(*outPressure) || std::isinf(*outPressure) ||
          	 std::isnan(*outVelocityU) || std::isinf(*outVelocityU) ||std::isnan(*outVelocityV) || std::isinf(*outVelocityV) ||std::isnan(*outVelocityW) || std::isinf(*outVelocityW)) {
               	 	assert(false);
		}

	}
	if(m_iDimension==2)
	{
                double div=Ud[0]+Vd[1];
		double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+2*Ud[1]*Vd[0];
                double Volumet=inVolume*div;
                double VelocityUt=-inVolume*Pd[0];
                double VelocityVt=-inVolume*Pd[1];
                double Pt=-gamma*(inPressure+Pinf)*div;
                double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-inVolume*(Pd[2]+Pd[3])-cross);
                double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[2]+Vd[4])+Ud[0]*Pd[0]+Vd[0]*Pd[1]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[4]+Vd[3])+Ud[1]*Pd[0]+Vd[1]*Pd[1]);
                double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+inVolume*(Pd[2]+Pd[3])+cross);
                (*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;
                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
//                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;
                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
                if(std::isnan(*outVolume) || std::isinf(*outVolume) ||
                   std::isnan(*outPressure) || std::isinf(*outPressure) ||
                 std::isnan(*outVelocityU) || std::isinf(*outVelocityU) ||std::isnan(*outVelocityV) || std::isinf(*outVelocityV)) {
                        assert(false);
                }
		Pd[3]=gamma*(inPressure+Pinf)*(inVolume*(Pd[2]+Pd[3]));
		Pd[4]=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+cross);
		Pd[2]=-gamma*(inPressure+Pinf)*div;
//	      Pd[3]=VelocityUtt;
//              Pd[4]=Ptt;

	}
}

void HyperbolicLPSolver::timeIntegration(int index, double Dt, double gravity, double inVolume, double inVelocityU, double inVelocityV, double inVelocityW, double inPressure, double inSoundSpeed,
                                      double* Volumed, double* Ud, double* Vd, double *Wd, double *Pd,
                                                        double* outVolume, double* outVelocityU, double* outVelocityV, double* outVelocityW, double* outPressure){
//      double gamma=inSoundSpeed*inSoundSpeed/inVolume/inPressure;
        double gamma = m_pGamma;
        double Pinf=m_pPinf;

        // TODO Note that this functions only works for Poly EOS!!!!!!!
        if(m_iDimension==3)
        {

                double div=Ud[0]+Vd[1]+Wd[2];
                double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+Wd[2]*Wd[2]+2*Ud[1]*Vd[0]+2*Ud[2]*Wd[0]+2*Vd[2]*Wd[1];
                double Volumet=inVolume*div;
                double VelocityUt=-inVolume*Pd[0];
                double VelocityVt=-inVolume*Pd[1];
                double VelocityWt=-inVolume*Pd[2];
                double Pt=-gamma*(inPressure+Pinf)*div;
                double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-Volumed[2]*Pd[2]-inVolume*(Pd[3]+Pd[4]+Pd[5])-cross);
                double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[3]+Vd[6]+Wd[7])+Ud[0]*Pd[0]+Vd[0]*Pd[1]+Wd[0]*Pd[2]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[6]+Vd[4]+Wd[8])+Ud[1]*Pd[0]+Vd[1]*Pd[1]+Wd[1]*Pd[2]);
                double VelocityWtt=inVolume*((gamma-1)*Pd[2]*div+gamma*(inPressure+Pinf)*(Ud[7]+Vd[8]+Wd[5])+Ud[2]*Pd[0]+Vd[2]*Pd[1]+Wd[2]*Pd[2]);
                double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+Volumed[2]*Pd[2]+inVolume*(Pd[3]+Pd[4]+Pd[5])+cross);
                (*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;
                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;//TODO: MODIFIED GRAVITY DIRECTION
                (*outVelocityW)=inVelocityW+Dt*VelocityWt+0.5*Dt*Dt*VelocityWtt;
                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
                if(std::isnan(*outVolume) || std::isinf(*outVolume) ||
                   std::isnan(*outPressure) || std::isinf(*outPressure) ||
                 std::isnan(*outVelocityU) || std::isinf(*outVelocityU) ||std::isnan(*outVelocityV) || std::isinf(*outVelocityV) ||std::isnan(*outVelocityW) || std::isinf(*outVelocityW)) {
                        assert(false);
                }

        }
        if(m_iDimension==2)
        {

                double div=Ud[0]+Vd[1];
                double cross=Ud[0]*Ud[0]+Vd[1]*Vd[1]+2*Ud[1]*Vd[0];
                double Volumet=inVolume*div;
                double VelocityUt=-inVolume*Pd[0];
                double VelocityVt=-inVolume*Pd[1];
                double Pt=-gamma*(inPressure+Pinf)*div;
                double Volumett=inVolume*(div*div-Volumed[0]*Pd[0]-Volumed[1]*Pd[1]-inVolume*(Pd[2]+Pd[3])-cross);
                double VelocityUtt=inVolume*((gamma-1)*Pd[0]*div+gamma*(inPressure+Pinf)*(Ud[2]+Vd[4])+Ud[0]*Pd[0]+Vd[0]*Pd[1]);
                double VelocityVtt=inVolume*((gamma-1)*Pd[1]*div+gamma*(inPressure+Pinf)*(Ud[4]+Vd[3])+Ud[1]*Pd[0]+Vd[1]*Pd[1]);
                double Ptt=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+inVolume*(Pd[2]+Pd[3])+cross);
                (*outVolume)=inVolume+Dt*Volumet+0.5*Dt*Dt*Volumett;
                (*outVelocityU)=inVelocityU+Dt*VelocityUt+0.5*Dt*Dt*VelocityUtt;
                (*outVelocityV)=inVelocityV+Dt*VelocityVt+0.5*Dt*Dt*VelocityVtt+Dt*gravity;

                (*outPressure)=inPressure+Dt*Pt+0.5*Dt*Dt*Ptt;
                if(std::isnan(*outVolume) || std::isinf(*outVolume) ||
                   std::isnan(*outPressure) || std::isinf(*outPressure) ||
                 std::isnan(*outVelocityU) || std::isinf(*outVelocityU) ||std::isnan(*outVelocityV) || std::isinf(*outVelocityV)) {
                        assert(false);
                }
                Pd[3]=gamma*(inPressure+Pinf)*(inVolume*(Pd[2]+Pd[3]));
                Pd[4]=gamma*gamma*(inPressure+Pinf)*div*div+gamma*(inPressure+Pinf)*(Volumed[0]*Pd[0]+Volumed[1]*Pd[1]+cross);
                Pd[2]=-gamma*(inPressure+Pinf)*div;

        }
}

void HyperbolicLPSolver::printInvalidState(int phase, int dir, int index, 
	double positionX, double positionY, double positionZ, 
	double vel_d_0, double vel_dd_0, double p_d_0, double p_dd_0, 
	double vel_d_1, double vel_dd_1, double p_d_1, double p_dd_1) {
	
	cout<<"---------------------Invalid state-----------------------"<<endl;
	cout<<"phase="<<phase<<", direction="<<dir<<endl;
	cout<<"x["<<index<<"]="<<positionX<<", y["<<index<<"]="<<positionY
		<<", z["<<index<<"]="<<positionZ<<endl;	
	cout<<"---------------First Order Derivatives-------------------"<<endl;
	cout<<"vel_d_0="<<vel_d_0<<", vel_d_1="<<vel_d_1<<endl;
	cout<<"p_d_0="<<p_d_0<<", p_d_1="<<p_d_1<<endl;
	cout<<"---------------Second Order Derivatives------------------"<<endl;
	cout<<"vel_dd_0="<<vel_dd_0<<", vel_dd_1="<<vel_dd_1<<endl;
	cout<<"p_dd_0="<<p_dd_0<<", p_dd_1="<<p_dd_1<<endl;
	cout<<"---------------------------------------------------------"<<endl;
}


bool HyperbolicLPSolver::lowerLPFOrder(int index, const vector<int*>& LPFOrderOther, // input
	int* LPFOrder0, int* LPFOrder1) { // output
		
	//// lower the order in this direction (both sides)
	//if(LPFOrder0[index]>0) LPFOrder0[index]--;
	//if(LPFOrder1[index]>0) LPFOrder1[index]--;
	//
	//for(size_t i=0; i<LPFOrderOther.size(); i++) {
	//	if(LPFOrderOther[i][index]>0) return true;
	//}

	//if(LPFOrder0[index]==0 && LPFOrder1[index]==0) return false;
	//else return true;


	// lower the order in this direction (both sides)
	if(LPFOrder0[index]>0) LPFOrder0[index]--;
	if(LPFOrder1[index]>0) LPFOrder1[index]--;
	
	if(LPFOrder0[index]==0 && LPFOrder1[index]==0) {
		m_vDoNotUpdate[index]=true;
		for(size_t i=0; i<LPFOrderOther.size(); i++) 
			LPFOrderOther[i][index]=0;	
		return false;
	}

	return true;

}


void HyperbolicLPSolver::setBoundaryPressureAndVelocity(int phase) {
	
	// no need to update for the last phase
	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;

	// iteration start and end index
	size_t startIndex = m_pParticleData->getBoundaryStartIndex();
	size_t numParticle = m_pParticleData->getBoundaryNum();
	if(numParticle==0) return;
	
	//------- specifications of pointers and parameters-------
	
	// determine dir: x(0), y(1), or z(2)
	const int dir = (phase==-1)? -1:m_vDirSplitTable[m_iDirSplitOrder][phase];
	
	vector<const double*> position;
	position.push_back(m_pParticleData->m_vPositionX);
	position.push_back(m_pParticleData->m_vPositionY);
	if(m_iDimension==3) position.push_back(m_pParticleData->m_vPositionZ);		

	// input and also output data
	vector<double*> data;

	//-------------------------------------------------------
	
	if(phase==-1) {		
		// set input/output data
		data.push_back(m_pParticleData->m_vPressure);
		data.push_back(m_pParticleData->m_vVelocityU);
		data.push_back(m_pParticleData->m_vVelocityV);
		if(m_iDimension==3) data.push_back(m_pParticleData->m_vVelocityW);	
	}
	else {
		if(m_iDimension==2) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);
				if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityV); // next step is V
				else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityU); // next step is U
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
			}
			else assert(false);
		}
		else if(m_iDimension==3) {
			if(phase==0) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
			}
			else if(phase==2) {
				data.push_back(m_pParticleData->m_vTemp1Pressure);
				if(dir==0) { 
					if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp1VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp1VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp1VelocityW);
			}
			else if(phase==3) {
				data.push_back(m_pParticleData->m_vTemp2Pressure);	
				if(dir==0) { 
					if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp2VelocityW);
			}
			else assert(false);
		}	
	}
	
	compute0thOrderWeightedLPF(position, startIndex, numParticle, data);	
	
	for(size_t k=1; k<data.size(); k++) // k starts from 1: change velocity to the opposite of fluid particles
		for(size_t index=startIndex; index<startIndex+numParticle; index++)
			data[k][index] *= -1;
}


//void HyperbolicLPSolver::setGhostPressureAndVelocity(int phase) {
//	
//	// no need to update for the last phase
//	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;
//
//	// iteration start and end index
//	size_t startIndex = m_pParticleData->getGhostStartIndex();
//	size_t numParticle = m_pParticleData->getGhostNum();
//	if(numParticle==0) return;
////	cout<<"-------HyperbolicLPSolver::setGhostPressureAndVelocity()-------"<<endl;
////	cout<<"startIndex = "<<startIndex<<endl;
////	cout<<"numParticle = "<<numParticle<<endl;
////	cout<<"---------------------------------------------------------------"<<endl;
//	//------- specifications of pointers and parameters-------
//	
//	// determine dir: x(0), y(1), or z(2)
//	const int dir = (phase==-1)? -1:m_vDirSplitTable[m_iDirSplitOrder][phase];
//	
//	vector<const double*> position;
//	position.push_back(m_pParticleData->m_vPositionX);
//	position.push_back(m_pParticleData->m_vPositionY);
//	if(m_iDimension==3) position.push_back(m_pParticleData->m_vPositionZ);		
//
//	// input and also output data
//	vector<double*> data;
//
//	//-------------------------------------------------------
//	
//	if(phase==-1) {		
//		// set input/output data
//		data.push_back(m_pParticleData->m_vPressure);
//		data.push_back(m_pParticleData->m_vVelocityU);
//		data.push_back(m_pParticleData->m_vVelocityV);
//		if(m_iDimension==3) data.push_back(m_pParticleData->m_vVelocityW);	
//	}
//	else {
//		if(m_iDimension==2) {
//			if(phase==0) {
//				data.push_back(m_pParticleData->m_vTemp1Pressure);	
//			}
//			else if(phase==1) {
//				data.push_back(m_pParticleData->m_vTemp2Pressure);
//				if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityV); // next step is V
//				else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityU); // next step is U	
//			}	
//			else assert(false);
//		}
//		else if(m_iDimension==3) {
//			if(phase==0) {
//				data.push_back(m_pParticleData->m_vTemp1Pressure);	
//			}
//			else if(phase==1) {
//				data.push_back(m_pParticleData->m_vTemp2Pressure);	
//			}
//			else if(phase==2) {
//				data.push_back(m_pParticleData->m_vTemp1Pressure);
//				if(dir==0) { 
//					if(m_iDirSplitOrder==3)
//						data.push_back(m_pParticleData->m_vTemp1VelocityW);
//					else if(m_iDirSplitOrder==5)
//						data.push_back(m_pParticleData->m_vTemp1VelocityV);
//					else assert(false);
//				}
//				else if(dir==1) {
//					if(m_iDirSplitOrder==1)
//						data.push_back(m_pParticleData->m_vTemp1VelocityW);
//					else if(m_iDirSplitOrder==4)
//						data.push_back(m_pParticleData->m_vTemp1VelocityU);
//					else assert(false);
//				}
//				else if(dir==2) {
//					if(m_iDirSplitOrder==0)
//						data.push_back(m_pParticleData->m_vTemp1VelocityV);
//					else if(m_iDirSplitOrder==2)
//						data.push_back(m_pParticleData->m_vTemp1VelocityU);
//					else assert(false);
//				}
//				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp1VelocityU);
//				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp1VelocityV);
//				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp1VelocityW);
//			}
//			else if(phase==3) {
//				data.push_back(m_pParticleData->m_vTemp2Pressure);	
//				if(dir==0) { 
//					if(m_iDirSplitOrder==2)
//						data.push_back(m_pParticleData->m_vTemp1VelocityV);
//					else if(m_iDirSplitOrder==4)
//						data.push_back(m_pParticleData->m_vTemp1VelocityW);
//					else assert(false);
//				}
//				else if(dir==1) {
//					if(m_iDirSplitOrder==0)
//						data.push_back(m_pParticleData->m_vTemp1VelocityU);
//					else if(m_iDirSplitOrder==5)
//						data.push_back(m_pParticleData->m_vTemp1VelocityW);
//					else assert(false);
//				}
//				else if(dir==2) {
//					if(m_iDirSplitOrder==1)
//						data.push_back(m_pParticleData->m_vTemp1VelocityU);
//					else if(m_iDirSplitOrder==3)
//						data.push_back(m_pParticleData->m_vTemp1VelocityV);
//					else assert(false);
//				}
//				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
//				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
//				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp2VelocityW);
//			}
//			else assert(false);
//		}	
//	}
//	
//	compute0thOrderWeightedLPF(position, startIndex, numParticle, data);	
//		
//}


void HyperbolicLPSolver::setMirrorPressureAndVelocity(int phase) {

	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;

	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t boundaryStartIndex = m_pParticleData->m_iBoundaryStartIndex + m_pParticleData->m_iInflowNum;
	size_t boundaryEndIndex = m_pParticleData->m_iBoundaryStartIndex + m_pParticleData->m_iBoundaryNum;
	

	if(phase==0)
	{
                        for(size_t index=fluidEndIndex; index<boundaryStartIndex; index++)//inflow
                        {
                                m_pParticleData->m_vTemp1Volume[index]=m_pParticleData->m_vVolume[index];
                                m_pParticleData->m_vTemp1Pressure[index]=m_pParticleData->m_vPressure[index];
                                m_pParticleData->m_vTemp1VelocityU[index]=m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp1VelocityV[index]=m_pParticleData->m_vVelocityV[index];
                                if(m_iDimension==3) m_pParticleData->m_vTemp1VelocityW[index]=m_pParticleData->m_vVelocityW[index];
				m_pParticleData->m_vTemp1SoundSpeed[index]=m_pParticleData->m_vSoundSpeed[index];
                                m_pParticleData->m_vTemp2Volume[index]=m_pParticleData->m_vVolume[index];
                                m_pParticleData->m_vTemp2Pressure[index]=m_pParticleData->m_vPressure[index];
                                m_pParticleData->m_vTemp2VelocityU[index]=m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp2VelocityV[index]=m_pParticleData->m_vVelocityV[index];
                                if(m_iDimension==3)m_pParticleData->m_vTemp2VelocityW[index]=m_pParticleData->m_vVelocityW[index];
                                m_pParticleData->m_vTemp2SoundSpeed[index]=m_pParticleData->m_vSoundSpeed[index];

                        }

	}
	if(m_iDimension==3) {
		// pressure and velocity
		if(phase==0 || phase==2) {
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) {	
				size_t fIndex = m_vMirrorIndex[index-boundaryStartIndex];
				assert(fIndex>=fluidStartIndex && fIndex<boundaryStartIndex);
				m_pParticleData->m_vTemp1Pressure[index] = m_pParticleData->m_vTemp1Pressure[fIndex];
//                                m_pParticleData->m_vTemp1Pressure[index] = m_pParticleData->m_vPressure[index];
                                m_pParticleData->m_vTemp1VelocityU[index] = m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp1VelocityV[index] = m_pParticleData->m_vVelocityV[index];
                                m_pParticleData->m_vTemp1VelocityW[index] = m_pParticleData->m_vVelocityW[index];
//				if(phase==2) {
//					m_pParticleData->m_vTemp1VelocityU[index] = m_pParticleData->m_vTemp1VelocityU[fIndex];
//					m_pParticleData->m_vTemp1VelocityV[index] = m_pParticleData->m_vTemp1VelocityV[fIndex];
//					m_pParticleData->m_vTemp1VelocityW[index] = m_pParticleData->m_vTemp1VelocityW[fIndex];	
//				}
			}
		}
		else if(phase==1 || phase==3) {
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) {	
				size_t fIndex = m_vMirrorIndex[index-boundaryStartIndex];
				assert(fIndex>=fluidStartIndex && fIndex<boundaryStartIndex);	
//                                m_pParticleData->m_vTemp2Pressure[index] = m_pParticleData->m_vPressure[index];
				m_pParticleData->m_vTemp2Pressure[index] = m_pParticleData->m_vTemp2Pressure[fIndex];
                                m_pParticleData->m_vTemp2VelocityU[index] = m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp2VelocityV[index] = m_pParticleData->m_vVelocityV[index];
                                m_pParticleData->m_vTemp2VelocityW[index] = m_pParticleData->m_vVelocityW[index];
//				if(phase==3) {
//					m_pParticleData->m_vTemp2VelocityU[index] = m_pParticleData->m_vTemp2VelocityU[fIndex];
//					m_pParticleData->m_vTemp2VelocityV[index] = m_pParticleData->m_vTemp2VelocityV[fIndex];
//					m_pParticleData->m_vTemp2VelocityW[index] = m_pParticleData->m_vTemp2VelocityW[fIndex];	
//				}	
			}	
		}
 
	}
	else if(m_iDimension==2) {	
		// pressure and velocity
		if(phase==0) {
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) {	
				size_t fIndex = m_vMirrorIndex[index-boundaryStartIndex];
				assert(fIndex>=fluidStartIndex && fIndex<boundaryStartIndex);
				m_pParticleData->m_vTemp1Pressure[index] = m_pParticleData->m_vTemp1Pressure[fIndex];
                                m_pParticleData->m_vTemp1VelocityU[index] = m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp1VelocityV[index] = m_pParticleData->m_vVelocityV[index];
			}
		}
		else if(phase==1) {
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) {	
				size_t fIndex = m_vMirrorIndex[index-boundaryStartIndex];
				assert(fIndex>=fluidStartIndex && fIndex<boundaryStartIndex);	
				m_pParticleData->m_vTemp2Pressure[index] = m_pParticleData->m_vTemp2Pressure[fIndex];
                                m_pParticleData->m_vTemp2VelocityU[index] = m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vTemp2VelocityV[index] = m_pParticleData->m_vVelocityV[index];
//				m_pParticleData->m_vTemp2VelocityU[index] = m_pParticleData->m_vTemp2VelocityU[fIndex];
//				m_pParticleData->m_vTemp2VelocityV[index] = m_pParticleData->m_vTemp2VelocityV[fIndex];
			}	
		
		}
		
	}


}


void HyperbolicLPSolver::setGhostVelocity(int phase) {
		
//	if((m_iDimension==2 && (phase==0 || phase==2)) || 
//	   (m_iDimension==3 && (phase==0 || phase==1 || phase==4))) return;
	size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
	size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	
	if(m_iDimension==3) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			if(m_vFillGhost[index]) {
				size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
				
				for(int i=0; i<m_pParticleData->m_vNeighbourListSize[index]; i++) {	
							
					size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];

					if(neiIndex>=ghostStartIndex) {// ghost
                                                m_pParticleData->m_vTemp1Pressure[neiIndex] = 0;//Need to set pressure as well!
                                                m_pParticleData->m_vTemp2Pressure[neiIndex] = 0;
						m_pParticleData->m_vTemp1VelocityU[neiIndex] = m_pParticleData->m_vTemp1VelocityU[index];
						m_pParticleData->m_vTemp2VelocityU[neiIndex] = m_pParticleData->m_vTemp2VelocityU[index];
						m_pParticleData->m_vTemp1VelocityV[neiIndex] = m_pParticleData->m_vTemp1VelocityV[index];
						m_pParticleData->m_vTemp2VelocityV[neiIndex] = m_pParticleData->m_vTemp2VelocityV[index];
						m_pParticleData->m_vTemp1VelocityW[neiIndex] = m_pParticleData->m_vTemp1VelocityW[index];
						m_pParticleData->m_vTemp2VelocityW[neiIndex] = m_pParticleData->m_vTemp2VelocityW[index];
					}
				}
			}
		}
	}
	else if(m_iDimension==2) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			if(m_vFillGhost[index]) {
				size_t neiListStartIndex = index*m_pParticleData->m_iMaxNeighbourNum;
				for(int i=0; i<m_pParticleData->m_vNeighbourListSize[index]; i++) {
						
					size_t neiIndex = m_pParticleData->m_vNeighbourList[neiListStartIndex+i];

					if(neiIndex>=ghostStartIndex) {// ghost
						m_pParticleData->m_vTemp1Pressure[neiIndex] = 0;//Need to set pressure as well!
                                                m_pParticleData->m_vTemp2Pressure[neiIndex] = 0;
						m_pParticleData->m_vTemp1VelocityU[neiIndex] = m_pParticleData->m_vTemp1VelocityU[index];
						m_pParticleData->m_vTemp2VelocityU[neiIndex] = m_pParticleData->m_vTemp2VelocityU[index];
						m_pParticleData->m_vTemp1VelocityV[neiIndex] = m_pParticleData->m_vTemp1VelocityV[index];
						m_pParticleData->m_vTemp2VelocityV[neiIndex] = m_pParticleData->m_vTemp2VelocityV[index];	
					}
				}
			}
		}
	}


}


void HyperbolicLPSolver::setGhostPressureAndVelocity(int phase) {
	
	// no need to update for the last phase
	if((m_iDimension==2 && phase==2) || (m_iDimension==3 && phase==4)) return;

	// iteration start and end index
	size_t startIndex = m_pParticleData->getGhostStartIndex();
	size_t numParticle = m_pParticleData->getGhostNum();
	if(numParticle==0) return;
//	cout<<"-------HyperbolicLPSolver::setGhostPressureAndVelocity()-------"<<endl;
//	cout<<"startIndex = "<<startIndex<<endl;
//	cout<<"numParticle = "<<numParticle<<endl;
//	cout<<"---------------------------------------------------------------"<<endl;
	//------- specifications of pointers and parameters-------
	
	// determine dir: x(0), y(1), or z(2)
	const int dir = (phase==-1)? -1:m_vDirSplitTable[m_iDirSplitOrder][phase];
	
	vector<const double*> position;
	position.push_back(m_pParticleData->m_vPositionX);
	position.push_back(m_pParticleData->m_vPositionY);
	if(m_iDimension==3) position.push_back(m_pParticleData->m_vPositionZ);		

	// input and also output data
	vector<double*> data;

	//-------------------------------------------------------
	
	if(phase==-1) {		
		// set input/output data
		//data.push_back(m_pParticleData->m_vPressure);
		fill_n(m_pParticleData->m_vPressure+startIndex,numParticle,0);
		fill_n(m_pParticleData->m_vTemp1Pressure+startIndex,numParticle,0);
		fill_n(m_pParticleData->m_vTemp2Pressure+startIndex,numParticle,0);
		data.push_back(m_pParticleData->m_vVelocityU);
		data.push_back(m_pParticleData->m_vVelocityV);
		if(m_iDimension==3) data.push_back(m_pParticleData->m_vVelocityW);	
	}
	else {
		if(m_iDimension==2) {
			if(phase==0) {
				return;
				//data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				//data.push_back(m_pParticleData->m_vTemp2Pressure);
				if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityV); // next step is V
				else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityU); // next step is U	
			}	
			else assert(false);
		}
		else if(m_iDimension==3) {
			if(phase==0) {
				return;
				//data.push_back(m_pParticleData->m_vTemp1Pressure);	
			}
			else if(phase==1) {
				return;
				//data.push_back(m_pParticleData->m_vTemp2Pressure);	
			}
			else if(phase==2) {
				//data.push_back(m_pParticleData->m_vTemp1Pressure);
				if(m_iDirSplitOrder==2 || m_iDirSplitOrder==4)
					data.push_back(m_pParticleData->m_vTemp1VelocityU);
				else if(m_iDirSplitOrder==0 || m_iDirSplitOrder==5)
					data.push_back(m_pParticleData->m_vTemp1VelocityV);
				else if(m_iDirSplitOrder==1 || m_iDirSplitOrder==3)
					data.push_back(m_pParticleData->m_vTemp1VelocityW);
				else assert(false);
	
				/*	
				if(dir==0) { 
					if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else assert(false);
				}
				*/
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp1VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp1VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp1VelocityW);
			}
			else if(phase==3) {
				//data.push_back(m_pParticleData->m_vTemp2Pressure);	
				if(m_iDirSplitOrder==0 || m_iDirSplitOrder==1)
					data.push_back(m_pParticleData->m_vTemp2VelocityU);
				else if(m_iDirSplitOrder==2 || m_iDirSplitOrder==3)
					data.push_back(m_pParticleData->m_vTemp2VelocityV);
				else if(m_iDirSplitOrder==4 || m_iDirSplitOrder==5)
					data.push_back(m_pParticleData->m_vTemp2VelocityW);
				else assert(false);	
				/*
				if(dir==0) { 
					if(m_iDirSplitOrder==2)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else if(m_iDirSplitOrder==4)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==1) {
					if(m_iDirSplitOrder==0)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==5)
						data.push_back(m_pParticleData->m_vTemp1VelocityW);
					else assert(false);
				}
				else if(dir==2) {
					if(m_iDirSplitOrder==1)
						data.push_back(m_pParticleData->m_vTemp1VelocityU);
					else if(m_iDirSplitOrder==3)
						data.push_back(m_pParticleData->m_vTemp1VelocityV);
					else assert(false);
				}
				*/
				//if(dir==0)		data.push_back(m_pParticleData->m_vTemp2VelocityU);
				//else if(dir==1) data.push_back(m_pParticleData->m_vTemp2VelocityV);
				//else if(dir==2) data.push_back(m_pParticleData->m_vTemp2VelocityW);
			}
			else assert(false);
		}	
	}
	
	compute0thOrderWeightedLPF(position, startIndex, numParticle, data);	
		
}

void HyperbolicLPSolver::compute0thOrderWeightedLPF(vector<const double*>& position, 
													size_t startIndex, size_t numParticle,
													vector<double*>& data) {
	
	//cout<<"-------HyperbolicLPSolver::computeOthOrderWeightedLPF()-------"<<endl;	
	
	//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
	//double startTime = omp_get_wtime();
	
	//cout<<"startIndex = "<<startIndex<<endl;
	//cout<<"numParticle = "<<numParticle<<endl;
	
	//int n = omp_get_max_threads();
	//vector<ofstream> ofs(n);
	//for(int i=0; i<n; i++) {
	//	string s = "./omp_debug/n_" + to_string(i);
	//	ofs[i].open(s.c_str());
	//}

	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for(size_t index=startIndex; index<startIndex+numParticle; index++) {	
		
		//int tid = omp_get_thread_num(); 
		//ofs[tid]<<"tid="<<tid<<" index="<<index<<endl;	

		
		vector<double> data_w2_total(data.size(),0);	
		double w2_total = 0;
		for(int i=0; i<m_pParticleData->m_vNeighbourListSize[index]; i++) { 
			int neiIndex = m_pParticleData->m_vNeighbourList[index*m_pParticleData->m_iMaxNeighbourNum+i];
			
			// make sure neighbour is a fluid particle
			assert(neiIndex >= (int)m_pParticleData->m_iFluidStartIndex && 
			       neiIndex < (int)m_pParticleData->m_iFluidStartIndex + (int)m_pParticleData->m_iFluidNum);	
			
			double dis2 = 0;
			for(size_t j=0; j<position.size(); j++) {
				double h = ( position[j][neiIndex]-position[j][index] );
				dis2 += (h*h);
			}	

			double weight = exp(-dis2);
			w2_total += (weight*weight);
			
			for(size_t k=0; k<data.size(); k++) 
				data_w2_total[k] += (data[k][neiIndex]*weight*weight);	
			
		}
		
		if(w2_total != 0) {
			for(size_t k=0; k<data.size(); k++)
				data[k][index] = data_w2_total[k]/w2_total;
		}
		else {
			for(size_t k=0; k<data.size(); k++) 
				data[k][index] = 0;
		}
			
	}

	//double elipsedTime = omp_get_wtime() - startTime;
	//printf("Compute zeroth order LPF takes %.16g seconds\n", elipsedTime);	
	
	//cout<<"--------------------------------------------------------------"<<endl;

}


//void HyperbolicLPSolver::compute0thOrderWeightedLPF(vector<const double*>& position, 
//													size_t startIndex, size_t numParticle,
//													vector<double*>& data) {
//	
//	//cout<<"-------HyperbolicLPSolver::computeOthOrderWeightedLPF()-------"<<endl;
//	
//	// alias
//	const int *neighbourList = m_pParticleData->m_vNeighbourList;
//	const int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;
//	const size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;
//	const size_t endIndex = startIndex + numParticle;	
//	
//	//cout<<"omp_get_max_threads() = "<<omp_get_max_threads()<<endl;
//	//double startTime = omp_get_wtime();
//	
//	//cout<<"startIndex = "<<startIndex<<endl;
//	//cout<<"numParticle = "<<numParticle<<endl;
//	
//	int n = omp_get_max_threads();
//	vector<ofstream> ofs(n);
//	for(int i=0; i<n; i++) {
//		string s = "./omp_debug/n_" + to_string(i);
//		ofs[i].open(s.c_str());
//	}
//
//	#ifdef _OPENMP
//	#pragma omp parallel for
//	#endif
//	for(size_t index=startIndex; index<endIndex; index++) {	
//		
//		//int tid = omp_get_thread_num(); 
//		//ofs[tid]<<"tid="<<tid<<" index="<<index<<endl;
//		
//
//		const size_t totalNumNei = neighbourListSize[index];
//		if(totalNumNei == 0) { // does not have neighbour
//			for(size_t k=0; k<data.size(); k++) {data[k][index] = 0;}
//			continue; 
//		}
//
//		// do 0th order weighted LPF
//		double data_w2_total[data.size()];
//		fill_n(data_w2_total,data.size(),0);
//		double w2_total = 0;
//		for(size_t i=0; i<totalNumNei; i++) { 
//			size_t neiIndex = neighbourList[index*maxNeiNum+i];
//			
//			// make sure neighbour is a fluid particle
//			assert(neiIndex >= m_pParticleData->m_iFluidStartIndex && 
//			       neiIndex < m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum);	
//			
//			double dis2 = 0;
//			for(size_t j=0; j<position.size(); j++) {
//				double h = ( position[j][neiIndex]-position[j][index] );
//				dis2 += (h*h);
//			}	
//
//			double weight = exp(-dis2);
//			w2_total += (weight*weight);
//			
//			for(size_t k=0; k<data.size(); k++) 
//				data_w2_total[k] += (data[k][neiIndex]*weight*weight);	
//			
//		}
//		
//		if(w2_total != 0) {
//			for(size_t k=0; k<data.size(); k++)
//				data[k][index] = data_w2_total[k]/w2_total;
//		}
//		else assert(false);
//			
//	}
//
//	//double elipsedTime = omp_get_wtime() - startTime;
//	//printf("Compute zeroth order LPF takes %.16g seconds\n", elipsedTime);	
//	
//	//cout<<"--------------------------------------------------------------"<<endl;
//
//}



void HyperbolicLPSolver::updateFluidState() {
		
	swap(m_pParticleData->m_vTemp1Volume, m_pParticleData->m_vVolume);
	swap(m_pParticleData->m_vTemp1Pressure, m_pParticleData->m_vPressure);
	swap(m_pParticleData->m_vTemp1SoundSpeed, m_pParticleData->m_vSoundSpeed);	
	
}


void HyperbolicLPSolver::moveFluidParticle() {
		
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	if(m_iDimension==2) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vPositionX[index] += 0.5 * m_fDt * 
			(1*m_pParticleData->m_vVelocityU[index] + 1*m_pParticleData->m_vTemp1VelocityU[index]); // 0.5 (old + new)	
			
			m_pParticleData->m_vPositionY[index] += 0.5 * m_fDt * 
			(1*m_pParticleData->m_vVelocityV[index] + 1*m_pParticleData->m_vTemp1VelocityV[index]);
		}	
	}
	else if(m_iDimension==3) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif	
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			m_pParticleData->m_vPositionX[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index]); // 0.5 (old + new)	
			
			m_pParticleData->m_vPositionY[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityV[index] + m_pParticleData->m_vTemp1VelocityV[index]);
			
			m_pParticleData->m_vPositionZ[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityW[index] + m_pParticleData->m_vTemp1VelocityW[index]);
		}
	}

}


void HyperbolicLPSolver::moveFluidParticleAdjusted() {
		
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	double eps = 0.5;

	if(m_iDimension==2) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			vector<double> result(2,0);
			computeLocalWeightedAverageVelocity(index,result);
			
			m_pParticleData->m_vPositionX[index] += m_fDt * (0.5 * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index])
			+ eps*result[0]); 	
			
			m_pParticleData->m_vPositionY[index] += m_fDt * (0.5 *
			(m_pParticleData->m_vVelocityV[index] + m_pParticleData->m_vTemp1VelocityV[index])
			+ eps*result[1]);
		}	
	}
	else if(m_iDimension==3) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
			vector<double> result(3,0);
			computeLocalWeightedAverageVelocity(index,result);

			m_pParticleData->m_vPositionX[index] += m_fDt * (0.5 * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index])
			+ eps*result[0]); 	
			
			m_pParticleData->m_vPositionY[index] += m_fDt * (0.5 *
			(m_pParticleData->m_vVelocityV[index] + m_pParticleData->m_vTemp1VelocityV[index])
			+ eps*result[1]);
			
			m_pParticleData->m_vPositionZ[index] += m_fDt * (0.5 *
			(m_pParticleData->m_vVelocityW[index] + m_pParticleData->m_vTemp1VelocityW[index])
			+ eps*result[2]);
		}
	}

}


void HyperbolicLPSolver::updateFluidVelocity() {
	
	swap(m_pParticleData->m_vTemp1VelocityU, m_pParticleData->m_vVelocityU);
	swap(m_pParticleData->m_vTemp1VelocityV, m_pParticleData->m_vVelocityV);	
	if(m_iDimension==3)	swap(m_pParticleData->m_vTemp1VelocityW, m_pParticleData->m_vVelocityW);

}

/*****************************************************************************
*
* The following functions are for testing purposes only
*
*****************************************************************************/

bool compByDist(const pair<double,int>& l, const pair<double,int>& r) {
	return l.first < r.first;
}

void HyperbolicLPSolver::testNeighbourSearch() {
	
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();
	
	// build the search structure
	m_pNeighbourSearcher->buildSearchStructure(positionX, positionY, positionZ, totalNumParticle);	
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;

	// the entire neighbour list (to be updated in the following loops)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	
	//double neiListDistTemp[maxNeiNum]; // a temp array for dist between a particle and its neighbours
//cout<<"1"<<endl;
	double* neighbourListDist = new double[m_pParticleData->m_iCapacity*maxNeiNum];	
	//double neighbourListDist[m_pParticleData->m_iCapacity*maxNeiNum]; 
//cout<<"2"<<endl;
	//
	//int* bruteForceNeighbourList = new int[m_pParticleData->m_iCapacity*maxNeiNum];
	//int* bruteForceNeighbourListSize = new int[m_pParticleData->m_iCapacity];
	int bruteForceNeighbourList[m_pParticleData->m_iCapacity*maxNeiNum];
//cout<<"3"<<endl;
	int bruteForceNeighbourListSize[m_pParticleData->m_iCapacity];
//cout<<"4"<<endl;
    //double bruteForceNeighbourListDist[m_pParticleData->m_iCapacity*maxNeiNum];
	double* bruteForceNeighbourListDist = new double[m_pParticleData->m_iCapacity*maxNeiNum];
//cout<<"5"<<endl;
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();		
	
	// fluid
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		
		size_t count = 0;
		vector<pair<double,int>> neiDistIndex;
		
		//-------brute-force neighbour search-------
		for(size_t j=fluidStartIndex; j<fluidStartIndex+totalNumParticle; j++) {
			double xd = positionX[index]-positionX[j];
			double yd = positionY[index]-positionY[j];
			double zd = positionZ[index]-positionZ[j];
			double dist = sqrt(xd*xd+yd*yd+zd*zd);
			if(dist <= m_fNeiSearchRadius && j!=index) {	
				neiDistIndex.push_back({dist,j});
				count++; // found a neighbour
			}
			if(count == maxNeiNum) break;
		}
		sort(neiDistIndex.begin(),neiDistIndex.end(),compByDist);
		for(size_t k=0; k<count; k++) {
			bruteForceNeighbourList[neiListStartIndex+k] = neiDistIndex[k].second;
			bruteForceNeighbourListDist[neiListStartIndex+k] = neiDistIndex[k].first;
		}
		bruteForceNeighbourListSize[index] = count;
		//------------------------------------------	
	
		
		//-------begin octree search-------
		size_t numNeiFound;
		m_pNeighbourSearcher->searchNeighbour(positionX[index], positionY[index], positionZ[index], m_fNeiSearchRadius, 
											  neighbourList+neiListStartIndex,neighbourListDist+neiListStartIndex,
											  numNeiFound,index); // output
		
		neighbourListSize[index] = numNeiFound;
		//---------------------------------
	}	
	
	// print debug info
	ofstream ofs1("bruteForceSearch");
	ofstream ofs2("octreeSearch");
	//cout<<"-------HyperbolicLPSolver::testNeighbourSearch()-------"<<endl;
	//cout<<"totalNumParticle="<<totalNumParticle<<endl;
	//cout<<"maxNeiNum="<<maxNeiNum<<endl;
	//cout<<"The neighbour list for Fluid particles:"<<endl;
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 	
		size_t neiListStartIndex = index*maxNeiNum;
		
		//cout<<"-----Brute Force SEARCH------"<<endl;
		ofs1<<"neighbourListSize["<<index<<"]="<<bruteForceNeighbourListSize[index]<<endl;
		ofs1<<"-----------------------------"<<endl;
		for(size_t k=neiListStartIndex; k<neiListStartIndex+bruteForceNeighbourListSize[index]; k++) 
			ofs1<<"neiIndex="<<bruteForceNeighbourList[k]<<"   "<<bruteForceNeighbourListDist[k]<<endl;
		ofs1<<"-----------------------------"<<endl;

		//cout<<"---------OCTREE SEARCH-------"<<endl;
		ofs2<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
		ofs2<<"-----------------------------"<<endl;
		for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++) 
			ofs2<<"neiIndex="<<neighbourList[k]<<"   "<<neighbourListDist[k]<<endl;
		ofs2<<"-----------------------------"<<endl;

	}	
	//cout<<"-----------------------------------------------------------------"<<endl;
	
	delete[] neighbourListDist;
	delete[] bruteForceNeighbourListDist;
}



void HyperbolicLPSolver::searchNeighbourBruteForce(int index, int* neighbourList, int* neighbourListSize) {

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();			
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
		
	size_t count = 0;
	vector<pair<double,int>> neiDistIndex(maxNeiNum);
	
	// WILL SEARCH FOR "ALL" PARTICLES FOR NEIGHBOURS!!!
	//-------brute-force neighbour search-------
	for(size_t j=0; j<totalNumParticle; j++) {
		if(j == (size_t)index) continue; // DO NOT INCLUDE THE PARTICLE ITSELF AS A NEIGHBOUR!!!

		double xd = positionX[index]-positionX[j];
		double yd = positionY[index]-positionY[j];
		double zd = positionZ[index]-positionZ[j];
		double dist = sqrt(xd*xd+yd*yd+zd*zd);
		if(dist <= m_fNeiSearchRadius) {	
			neiDistIndex[count] = {dist,j};	
			count++;
		}
		if(count == maxNeiNum) { //THIS IS NOT AN ERROR, JUST LIMIT THE SIZE OF NEIGHBOURS LIST
			cout<<"-------HyperbolicSolver::searchNeighbourBruteForce()-------"<<endl;
			cout<<"particle "<<index<<" has more than "<<maxNeiNum<<" neighbours!!!"<<endl;
			cout<<"-----------------------------------------------------------"<<endl;
			break;
		}
	}

	// SORT THE NEIGHBOURS IN ASCENDING ORDER OF THE DISTANCE FROM PARTICLE "INDEX"
	sort(neiDistIndex.begin(),neiDistIndex.begin()+count,compByDist);
	
	for(size_t k=0; k<count; k++) {
		neighbourList[k] = neiDistIndex[k].second;
		//cout<<"neiDistIndex["<<k<<"].first="<<neiDistIndex[k].first<<endl;
		//cout<<"neiDistIndex["<<k<<"].second="<<neiDistIndex[k].second<<endl;
		//cout<<"(size_t)neighbourList+k="<<(size_t)neighbourList+k<<endl;
		//cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	}
	neighbourListSize[0] = count;
	//------------------------------------------			
	
	//cout<<"-------HyperbolicLPSolver::searchNeighbourBruteForce()-------"<<endl;
	//cout<<"neighbourList="<<neighbourList<<endl;
	//cout<<"neighbourListSize[0]="<<neighbourListSize[0]<<endl;
	//cout<<"-------------------------------------------------------------"<<endl<<endl;

}


void HyperbolicLPSolver::searchNeighbourBruteForce(int index, int* neighbourList, 
												   double* neighbourListDist, int* neighbourListSize) {

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case
	const size_t totalNumParticle = m_pParticleData->getTotalNum();			
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
		
	size_t count = 0;
	vector<pair<double,int>> neiDistIndex(maxNeiNum);
	
	// WILL SEARCH FOR "ALL" PARTICLES FOR NEIGHBOURS!!!
	//-------brute-force neighbour search-------
	for(size_t j=0; j<totalNumParticle; j++) {
		if(j == (size_t)index) continue; // DO NOT INCLUDE THE PARTICLE ITSELF AS A NEIGHBOUR!!!

		double xd = positionX[index]-positionX[j];
		double yd = positionY[index]-positionY[j];
		double zd = positionZ[index]-positionZ[j];
		double dist = sqrt(xd*xd+yd*yd+zd*zd);
		if(dist <= m_fNeiSearchRadius) {	
			neiDistIndex[count] = {dist,j};	
			count++;
		}
		if(count == maxNeiNum) { //THIS IS NOT AN ERROR, JUST LIMIT THE SIZE OF NEIGHBOURS LIST
			cout<<"-------HyperbolicSolver::searchNeighbourBruteForce()-------"<<endl;
			cout<<"particle "<<index<<" has more than "<<maxNeiNum<<" neighbours!!!"<<endl;
			cout<<"-----------------------------------------------------------"<<endl;
			break;
		}
	}

	// SORT THE NEIGHBOURS IN ASCENDING ORDER OF THE DISTANCE FROM PARTICLE "INDEX"
	sort(neiDistIndex.begin(),neiDistIndex.begin()+count,compByDist);
	
	for(size_t k=0; k<count; k++) {
		neighbourList[k]     = neiDistIndex[k].second;
		neighbourListDist[k] = neiDistIndex[k].first;
		//cout<<"neiDistIndex["<<k<<"].first="<<neiDistIndex[k].first<<endl;
		//cout<<"neiDistIndex["<<k<<"].second="<<neiDistIndex[k].second<<endl;
		//cout<<"(size_t)neighbourList+k="<<(size_t)neighbourList+k<<endl;
		//cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	}
	neighbourListSize[0] = count;
	//------------------------------------------			
	
	//cout<<"-------HyperbolicLPSolver::searchNeighbourBruteForce()-------"<<endl;
	//cout<<"neighbourList="<<neighbourList<<endl;
	//cout<<"neighbourListSize[0]="<<neighbourListSize[0]<<endl;
	//cout<<"-------------------------------------------------------------"<<endl<<endl;

}


void HyperbolicLPSolver::searchNeighbourBruteForce(double x0, double y0, double z0, 
	int* neighbourList, size_t& numNeiFound) {

	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;
	const double *positionZ = m_pParticleData->m_vPositionZ; // is all zero for the 2D case			
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
		
	size_t count = 0;
	vector<pair<double,int>> neiDistIndex(maxNeiNum);
	
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	// WILL ONLY SEARCH FLUID PARTICLES AS NEIGHOBURS!!!
	//-------brute-force neighbour search-------
	for(size_t j=fluidStartIndex; j<fluidEndIndex; j++) {
		double xd = x0-positionX[j];
		double yd = y0-positionY[j];
		double zd = z0-positionZ[j];
		double dist = sqrt(xd*xd+yd*yd+zd*zd);
		if(dist <= m_fNeiSearchRadius) {	
			neiDistIndex[count] = {dist,j};
			count++;
		}
		if(count == maxNeiNum) { 
			cout<<"-------HyperbolicSolver::searchNeighbourBruteForce()-------"<<endl;
			cout<<"A ghost particle has more than "<<maxNeiNum<<" neighbours!!!"<<endl;
			cout<<"-----------------------------------------------------------"<<endl;
			break;
		}
	}
	sort(neiDistIndex.begin(),neiDistIndex.begin()+count,compByDist);
	for(size_t k=0; k<count; k++) 
		neighbourList[k] = neiDistIndex[k].second;
	numNeiFound = count;
	//------------------------------------------			

}


void HyperbolicLPSolver::generateGhostParticleByBruteForceNeighbourSearch() {		
	
	// use each fluid bounding box to generate ghost particles
	const vector<BoundingBox*>& fluidBoxes = m_pParticleData->m_vFluidBoundingBox;
	
	const double h_r = 0.5*m_fAvgParticleSpacing;
	const size_t capacity = m_pParticleData->m_iCapacity;	
	
	// temp space for valid ghost particles
	const size_t n = capacity - (m_pParticleData->m_iFluidNum + m_pParticleData->m_iBoundaryNum);
	vector<double> xGhost(n,0), yGhost(n,0), zGhost(n,0);
		
	int neiListTemp[m_pParticleData->m_iMaxNeighbourNum]; // the index of the neighbours of the ghost particle	
	
	size_t ghostCount = 0;
	if(m_iDimension==2) {

		for(size_t p=0; p<fluidBoxes.size(); p++) {
			
			int tag = fluidBoxes[p]->getObjectTag();

			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
				
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	

			// compute the location of particles	
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
						
						size_t numNeiFound;
						searchNeighbourBruteForce(x,y,0,neiListTemp,numNeiFound);
						
						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,tag)) continue;		
						if(ghostCount >= n) assert(false); // exceed array size capacity
						xGhost[ghostCount] = x;
						yGhost[ghostCount] = y;
						ghostCount++;
					
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						
						size_t numNeiFound;
						searchNeighbourBruteForce(x,y,0,neiListTemp,numNeiFound);

						//cout<<numNeiFound<<endl;
						if(!isValidGhostParticle(x,y,0,neiListTemp,numNeiFound,tag)) continue;		
						if(ghostCount >= n) assert(false); // exceed array size capacity	
						xGhost[ghostCount] = x;
						yGhost[ghostCount] = y;
						ghostCount++;
					}
				}
			}	
		}	
	}
	else if(m_iDimension==3) {
	
		for(size_t p=0; p<fluidBoxes.size(); p++) {
			
			int tag = fluidBoxes[p]->getObjectTag();

			double xmin = fluidBoxes[p]->getXmin();	
			double xmax = fluidBoxes[p]->getXmax();
			double ymin = fluidBoxes[p]->getYmin();	
			double ymax = fluidBoxes[p]->getYmax();
			double zmin = fluidBoxes[p]->getZmin();	
			double zmax = fluidBoxes[p]->getZmax();

			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			//get parameters of hexagonal packing
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
			
			// compute the location of particles
			for(size_t i=l0; i<=l1; i++) { 
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);	
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,tag)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;		

							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,tag)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
							}
						}
					}
						
				} 
				else { //even-numbered layers
					for(size_t j=m0_even; j<=m1_even; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows
							for(size_t k=nn0_odd; k<=nn1_odd; k++) { 
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,tag)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								
								size_t numNeiFound;
								searchNeighbourBruteForce(x,y,z,neiListTemp,numNeiFound);	

								if(!isValidGhostParticle(x,y,z,neiListTemp,numNeiFound,tag)) continue;		
								if(ghostCount >= n) assert(false); // exceed array size capacity
								xGhost[ghostCount] = x;
								yGhost[ghostCount] = y;
								zGhost[ghostCount] = z;
								ghostCount++;	
							}
						}
					}	
				}    
			}	
		}	
	
	}
	
	size_t ghostIndex = m_pParticleData->getGhostStartIndex();
	for(size_t count=0; count<ghostCount; count++) {
		m_pParticleData->m_vPositionX[ghostIndex] = xGhost[count];
		m_pParticleData->m_vPositionY[ghostIndex] = yGhost[count];
		m_pParticleData->m_vPositionZ[ghostIndex] = zGhost[count];
		ghostIndex++;
	}

	m_pParticleData->m_iGhostNum = ghostIndex - m_pParticleData->m_iGhostStartIndex;
	m_pParticleData->m_iTotalNum = m_pParticleData->m_iFluidNum + 
								   m_pParticleData->m_iBoundaryNum + 
								   m_pParticleData->m_iGhostNum;
	
	cout<<"-------HyperbolicLPSolver::generateGhostParticleByBruteForceNeighbourSearch()-------"<<endl;
	cout<<"m_pParticleData->m_iGhostStartIndex="<<m_pParticleData->m_iGhostStartIndex<<endl;
	cout<<"m_pParticleData->m_iGhostNum="<<m_pParticleData->m_iGhostNum<<endl;
	cout<<"m_pParticleData->m_iTotalNum="<<m_pParticleData->m_iTotalNum<<endl;
	cout<<"------------------------------------------------------------------------------------"<<endl;
}

void HyperbolicLPSolver::verifyGhostParticleBySumOfDist(ofstream& ofs) {

	size_t start = m_pParticleData->getGhostStartIndex();
	size_t end   = m_pParticleData->getGhostStartIndex() + m_pParticleData->m_iGhostNum;
	ofs<<"# of ghost particles = "<<m_pParticleData->m_iGhostNum<<endl;
	double distSum = 0;
	for(size_t i=start; i<end; i++) {
		double x = m_pParticleData->m_vPositionX[i];
		double y = m_pParticleData->m_vPositionY[i];
		double z = m_pParticleData->m_vPositionZ[i];
		double dist = sqrt(x*x+y*y+z*z);
		distSum += dist;
	}
	ofs<<"Sum of distance to origin = "<<fixed<<setprecision(16)<<setw(24)<<distSum<<endl;

}


void HyperbolicLPSolver::searchNeighbourForAllParticleByBruteForceNeighbourSearch() {	
	
	// the entire neighbour list (to be updated in the following loops)
	int *neighbourList = m_pParticleData->m_vNeighbourList;
	int *neighbourListSize = m_pParticleData->m_vNeighbourListSize;	

	size_t boundaryStartIndex = m_pParticleData->getBoundaryStartIndex();
	size_t boundaryEndIndex = boundaryStartIndex + m_pParticleData->getBoundaryNum();
	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	size_t ghostStartIndex = m_pParticleData->getGhostStartIndex();
	size_t ghostEndIndex = ghostStartIndex + m_pParticleData->getGhostNum();
	
	size_t maxNeiNum = m_pParticleData->m_iMaxNeighbourNum;	
	
	//cout<<"---HyperbolicLPSolver::searchNeighbourForAllParticleByBruteForceNeighbourSearch()---"<<endl;

	// fluid
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);
			

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;
		
	}
	
	// boundary
	for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;			
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);

		// boundary particle take only fluid particles as neighbours
		int incr = 0;
		for(int k=0; k<neighbourListSize[index]; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;
	}


	// ghost
	for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
		
		size_t neiListStartIndex = index*maxNeiNum;
		searchNeighbourBruteForce(index,neighbourList+neiListStartIndex,neighbourListSize+index);

		// boundary particle take only fluid particles as neighbours
		int incr = 0;
		for(int k=0; k<neighbourListSize[index]; k++) {
			size_t neiI = (size_t)neighbourList[neiListStartIndex+k];
			if(neiI>=fluidStartIndex && neiI<fluidEndIndex) { // only fluid particle
				neighbourList[neiListStartIndex+incr] = neiI;
				incr++;
			} 
		}
		neighbourListSize[index] = incr;

		// PRINT DEBUG INFO
		//cout<<"index="<<index<<endl;
		//cout<<"neiListStartIndex="<<neiListStartIndex<<endl;
		//cout<<"neighbourList+neiListStartIndex="<<neighbourList+neiListStartIndex<<endl;
		//cout<<"neighbourListSize="<<neighbourListSize[index]<<endl;
		//for(int k=0; k<neighbourListSize[index]; k++) 
		//	cout<<"neighbourList["<<neiListStartIndex+k<<"]="<<neighbourList[neiListStartIndex+k]<<endl;		    
	}


	// PRINT DEBUG INFO
	//cout<<"-------HyperbolicLPSolver::searchNeighbourForAllParticle()-------"<<endl;
	//cout<<"m_pParticleData->getTotalNum()="<<m_pParticleData->getTotalNum()<<endl;
	//cout<<"maxNeiNum="<<maxNeiNum<<endl;
	//cout<<"Fluid particles:"<<endl;
	//for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	//cout<<"Boundary particles:"<<endl;
	//for(size_t index=boundaryStartIndex; index<boundaryEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	//cout<<"Ghost particles:"<<endl;
	//for(size_t index=ghostStartIndex; index<ghostEndIndex; index++) { 
	//	
	//	size_t neiListStartIndex = index*maxNeiNum;			
	//	cout<<"neighbourListSize["<<index<<"]="<<neighbourListSize[index]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//	for(size_t k=neiListStartIndex; k<neiListStartIndex+neighbourListSize[index]; k++)
	//		cout<<"neighbourList["<<k<<"]="<<neighbourList[k]<<endl;
	//	cout<<"-----------------------------"<<endl;
	//}
	cout<<"-----------------------------------------------------------------"<<endl;
	
}



bool HyperbolicLPSolver::checkUpwindNeighbourList() {
	
	const double *positionX = m_pParticleData->m_vPositionX;
	const double *positionY = m_pParticleData->m_vPositionY;		
	const double *positionZ = m_pParticleData->m_vPositionZ;
	if(m_iDimension==2) positionZ = nullptr;

	// get the upwind neighbour lists (to be updated in the following)
	int *neighbourListRight = m_pParticleData->m_vNeighbourListRight;
	int *neighbourListLeft = m_pParticleData->m_vNeighbourListLeft;
	int *neighbourListNorth = m_pParticleData->m_vNeighbourListNorth;
	int *neighbourListSouth = m_pParticleData->m_vNeighbourListSouth;
	int *neighbourListUp, *neighbourListDown;
	if(m_iDimension==3) {
		neighbourListUp = m_pParticleData->m_vNeighbourListUp;
		neighbourListDown = m_pParticleData->m_vNeighbourListDown;
	}

	// get the size of neighbour lists (to be updated in the following)
	int *neighbourListRightSize = m_pParticleData->m_vNeighbourListRightSize;
	int *neighbourListLeftSize = m_pParticleData->m_vNeighbourListLeftSize;
	int *neighbourListNorthSize = m_pParticleData->m_vNeighbourListNorthSize;
	int *neighbourListSouthSize = m_pParticleData->m_vNeighbourListSouthSize;
	int *neighbourListUpSize, *neighbourListDownSize;
	if(m_iDimension==3) {
		neighbourListUpSize = m_pParticleData->m_vNeighbourListUpSize;
		neighbourListDownSize = m_pParticleData->m_vNeighbourListDownSize;
	}
	
	size_t maxNeiNumInOneDir = m_pParticleData->m_iMaxNeighbourNumInOneDir;

	size_t fluidStartIndex = m_pParticleData->getFluidStartIndex();
	size_t fluidEndIndex = fluidStartIndex + m_pParticleData->getFluidNum();
	
	for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		
		size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
			
		// print the upwind neighbour lists
		//if(index%100==0 && index<=1000) cout<<"neighbourListRightSize["<<index<<"]="<<neighbourListRightSize[index]<<endl;
		for(int k=0; k<neighbourListRightSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListRight[j];
			assert(positionX[neiIndex] > positionX[index]);
			//cout<<"neighbourListRight["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionX[neiIndex]=%.16g\n",positionX[neiIndex]);
			//printf("positionX[index]=%.16g\n",positionX[index]);
		}
		//if(index%100==0 && index<=1000)cout<<"neighbourListLeftSize["<<index<<"]="<<neighbourListLeftSize[index]<<endl;
		for(int k=0; k<neighbourListLeftSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListLeft[j];
			assert(positionX[neiIndex] < positionX[index]);
			//cout<<"neighbourListLeft["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionX[neiIndex]=%.16g\n",positionX[neiIndex]);
			//printf("positionX[index]=%.16g\n",positionX[index]);
		}
		//if(index%100==0 && index<=1000) cout<<"neighbourListNorthSize["<<index<<"]="<<neighbourListNorthSize[index]<<endl;
		for(int k=0; k<neighbourListNorthSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListNorth[j];
			assert(positionY[neiIndex] > positionY[index]);
			//cout<<"neighbourListNorth["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionY[neiIndex]=%.16g\n",positionY[neiIndex]);
			//printf("positionY[index]=%.16g\n",positionY[index]);
		}
		//if(index%100==0 && index<=1000) cout<<"neighbourListSouthSize["<<index<<"]="<<neighbourListSouthSize[index]<<endl;
		for(int k=0; k<neighbourListSouthSize[index]; k++) {
			int j = neiListInOneDirStartIndex+k;
			int neiIndex = neighbourListSouth[j];
			assert(positionY[neiIndex] < positionY[index]);
			//cout<<"neighbourListSOuth["<<j<<"]="<<neiIndex<<endl;	
			//printf("positionY[neiIndex]=%.16g\n",positionY[neiIndex]);
			//printf("positionY[index]=%.16g\n",positionY[index]);
		}
	}	
	if(m_iDimension==3) {
		for(size_t index=fluidStartIndex; index<fluidEndIndex; index++) {
		
			size_t neiListInOneDirStartIndex = index*maxNeiNumInOneDir;
			
			// print the upwind neighbour lists
			//if(index%100==0 && index<=1000) cout<<"neighbourListUpSize["<<index<<"]="<<neighbourListUpSize[index]<<endl;
			for(int k=0; k<neighbourListUpSize[index]; k++) {
				int j = neiListInOneDirStartIndex+k;
				int neiIndex = neighbourListUp[j];
				assert(positionZ[neiIndex] > positionZ[index]);
				//cout<<"neighbourListUp["<<j<<"]="<<neiIndex<<endl;	
				//printf("positionZ[neiIndex]=%.16g\n",positionZ[neiIndex]);
				//printf("positionZ[index]=%.16g\n",positionZ[index]);
			}
			//if(index%100==0 && index<=1000) cout<<"neighbourListDownSize["<<index<<"]="<<neighbourListDownSize[index]<<endl;
			for(int k=0; k<neighbourListDownSize[index]; k++) {
				int j = neiListInOneDirStartIndex+k;
				int neiIndex = neighbourListDown[j];
				assert(positionZ[neiIndex] < positionZ[index]);
				//cout<<"neighbourListDown["<<j<<"]="<<neiIndex<<endl;	
				//printf("positionZ[neiIndex]=%.16g\n",positionZ[neiIndex]);
				//printf("positionZ[index]=%.16g\n",positionZ[index]);
			}	
		}	
	}	

	return true;
}



string HyperbolicLPSolver::rightFlush(size_t writeStep, size_t numDigits) {
	
	assert(pow(10,numDigits) >= writeStep);

	string result;

	if(writeStep == 0) numDigits--;
	for(size_t i=writeStep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writeStep); 
	
	return result;

}

int HyperbolicLPSolver::writeResult(double time, size_t writeStep, size_t startIndex, size_t numParticle) {
	
	// Create an output file the name "filename"
	string filename = "vtkDebug" + rightFlush(writeStep, 7) + ".vtk";
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Error opening file: %s\n",filename.c_str()); 
		return 1;
	}
	size_t endIndex = startIndex + numParticle;


	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);
	if(m_pParticleData->getDimension()==2) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vPositionX[i], m_pParticleData->m_vPositionY[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i = startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vPositionX[i], m_pParticleData->m_vPositionY[i], m_pParticleData->m_vPositionZ[i]);
	}
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	
	fprintf(outfile,"VECTORS Velocity double\n");
	if(m_pParticleData->getDimension()==2) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vVelocityU[i], m_pParticleData->m_vVelocityV[i], 0.);
	}
	else if(m_pParticleData->getDimension()==3) {
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g %.16g %.16g\n",m_pParticleData->m_vVelocityU[i], m_pParticleData->m_vVelocityV[i], m_pParticleData->m_vVelocityW[i]);
	}
	
	fprintf(outfile,"SCALARS LPFOrder_right int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderRight[i]);
	
	fprintf(outfile,"SCALARS LPFOrder_left int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderLeft[i]);

	fprintf(outfile,"SCALARS LPFOrder_north int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderNorth[i]);

	fprintf(outfile,"SCALARS LPFOrder_south int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderSouth[i]);
		

	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS LPFOrder_up int\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderUp[i]);

		fprintf(outfile,"SCALARS LPFOrder_down int\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%d\n",m_pParticleData->m_vLPFOrderDown[i]);	
	}


	fprintf(outfile,"SCALARS velocity_u double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityU[i]);
	
	fprintf(outfile,"SCALARS velocity_v double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityV[i]);
	
	if(m_pParticleData->getDimension()==3) {
		fprintf(outfile,"SCALARS velocity_w double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(size_t i=startIndex; i<endIndex; i++)
			fprintf(outfile,"%.16g\n",m_pParticleData->m_vVelocityW[i]);
	}
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vSoundSpeed[i]);
	
	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",m_pParticleData->m_vPressure[i]);
		
	fprintf(outfile,"SCALARS density double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",1./m_pParticleData->m_vVolume[i]);
		
        fprintf(outfile,"SCALARS phi double\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(size_t i=startIndex; i<endIndex; i++)
                fprintf(outfile,"%.16g\n",m_pParticleData->m_vPhi[i]);
	fclose(outfile);
	
	return 0;
	
	
}

////////////////////////////////////////////////////////////////////////////////////////
// End of HyperbolicLPSolver
////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////
// Start of HyperbolicLPSolver1D
////////////////////////////////////////////////////////////////////////////////////////


HyperbolicLPSolver1D::HyperbolicLPSolver1D(Initializer& init, ParticleData* pData) {
	
	srand(time(0));
	
	// from arg list
	m_pParticleData = pData; 
	m_pEOS = init.getEOS();
	
	// get parameters from init
	m_iDimension = init.getDimension();
	m_iLPFOrder = init.getLPFOrder();  
	m_iNumRow1stOrder = init.getNumRow1stOrder(); // NOT USED
	m_iNumRow2ndOrder = init.getNumRow2ndOrder(); // NOT USED
	m_iNumCol1stOrder = init.getNumCol1stOrder(); // NOT USED
	m_iNumCol2ndOrder = init.getNumCol2ndOrder(); // NOT USED
	m_fAvgParticleSpacing = init.getInitParticleSpacing();
	m_fInitParticleSpacing = m_fAvgParticleSpacing; 
	m_fInvalidPressure = init.getInvalidPressure(); 
	m_fInvalidDensity = init.getInvalidDensity();
	m_iUseCriticalPressure = init.getUseCriticalPressure();
	m_fCriticalPressure = init.getCriticalPressure();
        m_iIfSPH=init.getIfSPH();
//	cout<<"m_iUseCriticalPressure="<<m_iUseCriticalPressure<<endl;
//	cout<<"m_fCriticalPressure="<<m_fCriticalPressure<<endl;
        if(m_iIfSPH==1)
                printf("SPH density estimator\n");
        else if(m_iIfSPH==3)
		printf("SPH / PDE density estimator switch\n");
	else
                printf("PDE density updator\n");

	m_iIfDebug = init.getIfDebug();
	debug.open(init.getDebugfileName(), std::ofstream::out | std::ofstream::app);

	m_sBoundaryType = init.getBoundaryObjTypes()[0]; // 1D ONLY
	m_iUseLimiter = init.getUseLimiter(); // 1D ONLY 
	m_fThresholdP = init.getThresholdP(); // 1D ONLY

	if(m_iLPFOrder == 1) m_iUseLimiter = false;
	if(m_sBoundaryType == "periodic") { // 1D ONLY
		m_fPeriodicLeftBoundary =  m_pParticleData->m_vPositionX[0] 
		                           - 0.5*m_fAvgParticleSpacing;
		m_fPeriodicRightBoundary = m_pParticleData->m_vPositionX[m_pParticleData->m_iTotalNum-1] 
		                           + 0.5*m_fAvgParticleSpacing;
		m_fDisBetweenPeriodicBoundaries = m_fPeriodicRightBoundary - m_fPeriodicLeftBoundary;
	}
	else {
		m_fPeriodicLeftBoundary = 0; 
		m_fPeriodicRightBoundary = 0;
		m_fDisBetweenPeriodicBoundaries = 0;	
	}

	if(m_iIfDebug) {
		debug<<"m_fThresholdP="<<m_fThresholdP<<endl;	
		debug<<"m_fPeriodicLeftBoundary="<<m_fPeriodicLeftBoundary<<endl;
		debug<<"m_fPeriodicRightBoundary="<<m_fPeriodicRightBoundary<<endl;
		debug<<"m_fDisBetweenPeriodicBoundaries="<<m_fDisBetweenPeriodicBoundaries<<endl;
	}

	// for completeness initialize to zero	
	m_fDt = 0;
	
	computeSetupsForNextIteration();
		

}

void HyperbolicLPSolver1D::computeSetupsForNextIteration() {	

        if(m_iIfSPH) SPHDensityEstimatorForFluidParticle();
	
	cavitation();

	// update boundary particles
	if(m_sBoundaryType == "free") {
		updateFreeBoundaryLocation();
		updateFreeBoundaryPressureAndVelocity();
	}
	else if(m_sBoundaryType == "solid") {
		updateSolidBoundaryPressureAndVelocity();	
	}
	
	// to determine the dt for next step
	computeMinParticleSpacing();
	computeMaxSoundSpeed();
	computeMaxFluidVelocity();
	// update values of divided difference for limiter
	if(m_iUseLimiter) 
		updateLimiter();
}

void HyperbolicLPSolver1D::SPHDensityEstimatorForFluidParticle() {
	double SmoothingLengthMultiplier=1.5;
	double ErrorTolerance=0.2;
	size_t  np = m_pParticleData->m_iTotalNum;

        double* x     = m_pParticleData->m_vPositionX;
        double* volume     = m_pParticleData->m_vVolume;
        double* mass     = m_pParticleData->m_vMass;

		
        size_t begin, end;
        if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") {
                begin = 1;
                end = np-2;
        }
        else if(m_sBoundaryType == "periodic") {
                begin = 0;
                end = np-1;
        }
        else {
                cout<<"ERROR: Invalid boundary type!"<<endl;
                exit(1);
        }

	for(size_t i = begin; i <=end; i++) {
		double SmoothingLength = SmoothingLengthMultiplier*volume[i]*mass[i];
		double c=0.666667/SmoothingLength;
		double c1=2.0*c;
		double rho_estimate=0.0;
		double rho_estimatel=c*mass[i];
                double rho_estimater=c*mass[i];
		bool breakl=0;
		bool breakr=0;

		for (size_t j=i+1; (j<np)&&(!breakr); j++){
			double dis = (x[j]-x[i])/SmoothingLength;
			if (dis<1)
				rho_estimater=rho_estimater+c1*(1.0-1.5*dis*dis+0.75*dis*dis*dis)*mass[j];
			else if(dis<2)
				rho_estimater=rho_estimater+c1/4.0*(2-dis)*(2-dis)*(2-dis)*mass[j];
			else
				breakr=1;
		}
                for (int j=i-1; (j>=0)&&(!breakl); j--){
                        double dis = (x[i]-x[j])/SmoothingLength;
                        if (dis<1)
                                rho_estimatel=rho_estimatel+c1*(1.0-1.5*dis*dis+0.75*dis*dis*dis)*mass[j];
                        else if(dis<2)
                                rho_estimatel=rho_estimatel+c1/4.0*(2-dis)*(2-dis)*(2-dis)*mass[j];
                        else
                                breakl=1;
//			if(i==1)
//				cout<<j<<" "<<rho_estimatel<<endl;
                }
		double errorl=fabs((1.0/rho_estimatel-volume[i])/volume[i]);
                double errorr=fabs((1.0/rho_estimater-volume[i])/volume[i]);
		if ((errorl<ErrorTolerance&&errorr<ErrorTolerance))
			rho_estimate=(rho_estimatel+rho_estimater)/2.0;
		else if(errorl<errorr)
			rho_estimate=rho_estimatel;
		else
			rho_estimate=rho_estimater;

//		if(i==1)
//		{
//			cout<<mass[0]<<" "<<x[0]<<" "<<(x[1]-x[0])/SmoothingLength<<endl;
//			cout<<volume[i]<<" "<<rho_estimatel<<" "<<rho_estimater<<" "<<rho_estimate<<endl;
//		}

		volume[i]=1.0/rho_estimate;
	}
}


void HyperbolicLPSolver1D::cavitation() {
	if(m_iUseCriticalPressure==0) return; // no cavitation model
	
	size_t np      = m_pParticleData->m_iTotalNum;
	//cout<<"np="<<np<<endl;

	size_t begin, end;
	if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") { 
		begin = 1;
		end = np-2;
	}	
	else if(m_sBoundaryType == "periodic") {
		begin = 0;
		end = np-1;
	}
	else {
		cout<<"ERROR: Invalid boundary type!"<<endl;	
		exit(1);
	}
	
	size_t count = 0;
	
	for(size_t index=begin; index<=end; index++) { // update inner particles
		if(m_pParticleData->m_vPressure[index] < m_fCriticalPressure) {
			m_pParticleData->m_vPressure[index]=0;
			m_pParticleData->m_vSoundSpeed[index] = 
			m_pEOS->getSoundSpeed(m_pParticleData->m_vPressure[index],1./m_pParticleData->m_vVolume[index]);
			count++;
		}	
	}	
	
	cout<<"Cavitation: number of particles: "<<count<<endl;
}

int HyperbolicLPSolver1D::solve(double dt) {
	return solveUpwind(dt);
//	return solveNewUpwind(dt);
//	return solveUpwindPredictorCorrector(dt);
}
int HyperbolicLPSolver1D::solveDefault(double dt) {
	
//	return solve_van_leer(dt);

	
	if(m_iIfDebug) debug<<"--------------HyperbolicLPSolver1D::solve()--------------"<<endl;
	
	// dt for this time step 
	m_fDt = dt;

	//alias	
	double* up_old = m_pParticleData->m_vVelocityU;
	double* V_old  = m_pParticleData->m_vVolume;
	double* p_old  = m_pParticleData->m_vPressure;
	double* cs_old = m_pParticleData->m_vSoundSpeed;
	double* up     = m_pParticleData->m_vTemp1VelocityU;
	double* V      = m_pParticleData->m_vTemp1Volume;
	double* p      = m_pParticleData->m_vTemp1Pressure;
	double* cs     = m_pParticleData->m_vTemp1SoundSpeed;
	double* xp     = m_pParticleData->m_vPositionX;
	
	size_t np      = m_pParticleData->m_iTotalNum;
	//cout<<"np="<<np<<endl;

	size_t begin, end;
	if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") { 
		begin = 1;
		end = np-2;
	}	
	else if(m_sBoundaryType == "periodic") {
		begin = 0;
		end = np-1;
	}
	else {
		cout<<"ERROR: Invalid boundary type!"<<endl;	
		exit(1);
	} 
	
	// update inner particles
	for(size_t i=begin; i<=end; i++) {
		
		int left_order=0, right_order=0;
		computeLPFOrder(i, left_order, right_order);
		if(m_iIfDebug) {assert(left_order!=0); assert(right_order!=0);}
		
		//cout<<"i="<<i<<"	left_order="<<left_order<<"	right_order="<<right_order<<endl;

		double ux_right=0, ux_left=0, px_right=0, px_left=0;
		double uxx_right=0, uxx_left=0, pxx_right=0, pxx_left=0;
		computeSpatialDer(left_order,  0,  i,  up_old, xp, ux_left,  uxx_left); // 0:left nei
		computeSpatialDer(left_order,  0,  i,  p_old,  xp, px_left,  pxx_left); // 0:left nei
		computeSpatialDer(right_order, 1,  i,  up_old, xp, ux_right, uxx_right); // 1: right nei
		computeSpatialDer(right_order, 1,  i,  p_old,  xp, px_right, pxx_right); // 1: right nei
		
		//cout<<"i="<<i<<"	ux_left="<<px_left<<"	uxx_left="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	px_left="<<px_left<<"	pxx_left="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	ux_right="<<px_left<<"	uxx_right="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	px_right="<<px_left<<"	pxx_right="<<pxx_left<<endl;
		
		// update V, up, p
		timeIntegration(i, V_old, up_old, p_old, cs_old, 
		                ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right,
						V, up, p); // output
		
		// check if there is invalid state
		bool isInvalid = (p[i] < m_fInvalidPressure || V[i] < m_fInvalidDensity);
		if(isInvalid) {	
			printInvalidState(i, ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right);
			cout<<"ERROR: Invalid state occurs!"<<endl;
			exit(1);
		}
		
		// update sound speed
		cs[i] = m_pEOS->getSoundSpeed(p[i],1./V[i]);

	}	
	
	updateFluidState();
	moveFluidParticle();

	computeSetupsForNextIteration();
	
	if(m_iIfDebug) debug<<"-------------------------------------------------------"<<endl;
	return 0;

	
}


int HyperbolicLPSolver1D::solve_van_leer(double dt) {
	
	if(m_iIfDebug) debug<<"--------------HyperbolicLPSolver1D::solve_van_leer()--------------"<<endl;
	
	// dt for this time step 
	m_fDt = dt;

	//alias	
	double* up_old = m_pParticleData->m_vVelocityU;
	double* V_old  = m_pParticleData->m_vVolume;
	double* p_old  = m_pParticleData->m_vPressure;
	double* cs_old = m_pParticleData->m_vSoundSpeed;
	double* up     = m_pParticleData->m_vTemp1VelocityU;
	double* V      = m_pParticleData->m_vTemp1Volume;
	double* p      = m_pParticleData->m_vTemp1Pressure;
	double* cs     = m_pParticleData->m_vTemp1SoundSpeed;
	double* xp     = m_pParticleData->m_vPositionX;
	
	size_t np      = m_pParticleData->m_iTotalNum;
	//cout<<"np="<<np<<endl;

	size_t begin, end;
	if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") { 
		begin = 1;
		end = np-2;
	}	
	else if(m_sBoundaryType == "periodic") {
		begin = 0;
		end = np-1;
	}
	else {
		cout<<"ERROR: Invalid boundary type!"<<endl;	
		exit(1);
	} 
	
	// update inner particles
	for(size_t i=begin; i<=end; i++) {	
	
		double ux_right1=0, ux_left1=0, px_right1=0, px_left1=0; // 1st order first derivatives
		double ux_right2=0, ux_left2=0, px_right2=0, px_left2=0; // 2nd order first derivatives
		double uxx_right=0, uxx_left=0, pxx_right=0, pxx_left=0; // 2nd order second derivatives
		
		computeSpatialDer(1, 0,  i,  up_old, xp, ux_left1,  uxx_left); // 0:left nei
		computeSpatialDer(1, 1,  i,  up_old, xp, ux_right1, uxx_right); // 1: right nei
		
		double theta_u;
		//method1
		//if(ux_right1 == 0) theta_u = -1; 
		//else theta_u = ux_left1/ux_right1;
		
		//method2
		if(ux_right1==0 || ux_left1==0) theta_u = -1;	
		else theta_u = 0.5*(ux_right1/ux_left1 + ux_left1/ux_right1);

		double phi_u;
		if(theta_u <= 0) phi_u = 0;
		else phi_u = 2*theta_u/(1+theta_u);
		if((m_sBoundaryType == "free" || m_sBoundaryType == "solid") && (i==begin || i==end)) {
			ux_left2 = 0; uxx_left = 0;
			ux_right2 = 0; uxx_right = 0;
		}
		else {
			computeSpatialDer(2, 0,  i,  up_old, xp, ux_left2,  uxx_left); // 0:left nei
			computeSpatialDer(2, 1,  i,  up_old, xp, ux_right2, uxx_right); // 1: right nei
		}

		computeSpatialDer(1, 0,  i,  p_old, xp, px_left1,  pxx_left); // 0:left nei
		computeSpatialDer(1, 1,  i,  p_old, xp, px_right1, pxx_right); // 1: right nei
		
		double theta_p;
		
		//method1
		//if(px_right1 == 0) theta_p = -1; 
		//else theta_p = px_left1/px_right1;
		
		//method2
		if(px_right1==0 || px_left1==0) theta_p = -1;	
		else theta_p = 0.5*(px_right1/px_left1 + px_left1/px_right1);


		double phi_p;
		if(theta_p <= 0) phi_p = 0;
		else phi_p = 2*theta_p/(1+theta_p);
		
		if((m_sBoundaryType == "free" || m_sBoundaryType == "solid") && (i==begin || i==end)) {
			px_left2 = 0; pxx_left = 0;
			px_right2 = 0; pxx_right = 0;
		}
		else {
			computeSpatialDer(2, 0,  i,  p_old, xp, px_left2,  pxx_left); // 0:left nei
			computeSpatialDer(2, 1,  i,  p_old, xp, px_right2, pxx_right); // 1: right nei	
		}

		//double phi = max(phi_u,phi_p);
		double phi = min(phi_u,phi_p);	
		//double phi = 0.5*(phi_u+phi_p);

		//cout<<"i="<<i<<"	ux_left="<<px_left<<"	uxx_left="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	px_left="<<px_left<<"	pxx_left="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	ux_right="<<px_left<<"	uxx_right="<<pxx_left<<endl;
		//cout<<"i="<<i<<"	px_right="<<px_left<<"	pxx_right="<<pxx_left<<endl;
		
		// update V, up, p
		timeIntegration_van_leer(i, V_old, up_old, p_old, cs_old, phi, 
		                         ux_left1, px_left1, ux_right1, px_right1, 
								 ux_left2, uxx_left, px_left2, pxx_left, ux_right2, uxx_right, px_right2, pxx_right,
								 V, up, p); // output
		
		// check if there is invalid state
		bool isInvalid = (p[i] < m_fInvalidPressure || V[i] < m_fInvalidDensity);
		if(isInvalid) {	
			//printInvalidState(i, ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right);
			cout<<"ERROR: Invalid state occurs!"<<endl;
			exit(1);
		}
		
		// update sound speed
		cs[i] = m_pEOS->getSoundSpeed(p[i],1./V[i]);

	}	
	
	updateFluidState();
	moveFluidParticle();

	computeSetupsForNextIteration();
	
	if(m_iIfDebug) debug<<"-------------------------------------------------------"<<endl;
	return 0;
}

int HyperbolicLPSolver1D::solveUpwind(double dt) {

//      return solve_van_leer(dt);

        if(m_iIfDebug) debug<<"--------------HyperbolicLPSolver1D::solve()--------------"<<endl;

        // dt for this time step 
        m_fDt = dt;

        //alias 
        double* up_old = m_pParticleData->m_vVelocityU;
        double* V_old  = m_pParticleData->m_vVolume;
        double* p_old  = m_pParticleData->m_vPressure;
        double* cs_old = m_pParticleData->m_vSoundSpeed;
        double* up     = m_pParticleData->m_vTemp1VelocityU;
        double* V      = m_pParticleData->m_vTemp1Volume;
        double* p      = m_pParticleData->m_vTemp1Pressure;
        double* cs     = m_pParticleData->m_vTemp1SoundSpeed;
        double* xp     = m_pParticleData->m_vPositionX;

        size_t np      = m_pParticleData->m_iTotalNum;
        //cout<<"np="<<np<<endl;

        size_t begin, end;
        if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") {
                begin = 1;
                end = np-2;
        }      
        else if(m_sBoundaryType == "periodic") {
                begin = 0;
                end = np-1;
        }
        else {
                cout<<"ERROR: Invalid boundary type!"<<endl;
                exit(1);
        }

        // update inner particles
        for(size_t i=begin; i<=end; i++) {

                int left_order=0, right_order=0;
                computeLPFOrder(i, left_order, right_order);
                if(m_iIfDebug) {assert(left_order!=0); assert(right_order!=0);}

                //cout<<"i="<<i<<"      left_order="<<left_order<<"     right_order="<<right_order<<endl;

                double ux_right=0, ux_left=0, px_right=0, px_left=0;
                double uxx_right=0, uxx_left=0, pxx_right=0, pxx_left=0;
		computeSpatialDerUpwind(left_order, 0, i, up_old, p_old, xp, ux_left, px_left);//1: left nei
                computeSpatialDerUpwind(right_order, 1, i, up_old, p_old, xp, ux_right, px_right);//1: right nei

                //cout<<"i="<<i<<"      ux_left="<<px_left<<"   uxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_left="<<px_left<<"   pxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      ux_right="<<px_left<<"  uxx_right="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_right="<<px_left<<"  pxx_right="<<pxx_left<<endl;

                // update V, up, p
                timeIntegration(i, V_old, up_old, p_old, cs_old,
                                ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right,
                                                V, up, p); // output

                // check if there is invalid state
                bool isInvalid = (p[i] < m_fInvalidPressure || V[i] < m_fInvalidDensity || std::isnan(p[i]) || std::isnan(V[i]));
                if(isInvalid) {
                        printInvalidState(i, ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right);
                        cout<<"ERROR: Invalid state occurs!"<<endl;
                        exit(1);
                }

                // update sound speed
                cs[i] = m_pEOS->getSoundSpeed(p[i],1./V[i]);


		//debug 0329
/*		if(i==79)
		{
			printf("x_78 x_79 x_80 Deltax_left Deltax_right: \n%10f %10f %10f %10f %10f\n",xp[i-1],xp[i],xp[i+1],xp[i]-xp[i-1],xp[i+1]-xp[i]);
			printf("V_78 V_79 V_80 DeltaV_left DeltaV_right dV_left dV_right: \n%10f %10f %10f %10f %10f %10f %10f\n",V_old[i-1],V_old[i],V_old[i+1],V_old[i]-V_old[i-1],V_old[i+1]-V_old[i],(V_old[i]-V_old[i-1])/(xp[i]-xp[i-1]),(V_old[i+1]-V_old[i])/(xp[i+1]-xp[i]));
			printf("u_78 u_79 u_80 Deltau_left Deltau_right du_left du_right: \n%10f %10f %10f %10f %10f %10f %10f\n",up_old[i-1],up_old[i],up_old[i+1],up_old[i]-up_old[i-1],up_old[i+1]-up_old[i],(up_old[i]-up_old[i-1])/(xp[i]-xp[i-1]),(up_old[i+1]-up_old[i])/(xp[i+1]-xp[i]));
			printf("P_78 P_79 P_80 DeltaP_left DeltaP_right dP_left dP_right: \n%10f %10f %10f %10f %10f %10f %10f\n",p_old[i-1],p_old[i],p_old[i+1],p_old[i]-p_old[i-1],p_old[i+1]-p_old[i],(p_old[i]-p_old[i-1])/(xp[i]-xp[i-1]),(p_old[i+1]-p_old[i])/(xp[i+1]-xp[i]));
			double du_left=(up_old[i]-up_old[i-1])/(xp[i]-xp[i-1]);
			double du_right=(up_old[i+1]-up_old[i])/(xp[i+1]-xp[i]);
			double dP_left=(p_old[i]-p_old[i-1])/(xp[i]-xp[i-1]);
			double dP_right=(p_old[i+1]-p_old[i])/(xp[i+1]-xp[i]);
			double K=1.667*p_old[i]/V_old[i];
			printf("v/2, du_left+du_right, v/2*(du_left+du_right, -v/2sqrt(k), dp_right-dp_left), -v/2sqrt(k)*(dp_right-dp_left) vt: \n%10f %10f %10f %10f %10f %10f %10f\n",V_old[i]/2.0, du_left+du_right, V_old[i]/2.0*(du_left+du_right), -V_old[i]/2.0/sqrt(K), dP_right-dP_left, -V_old[i]/2.0/sqrt(K)*(dP_right-dP_left),V_old[i]/2.0*(du_left+du_right)-V_old[i]/2.0/sqrt(K)*(dP_right-dP_left));
			printf("Vsqrt(K)/2 du_right-du_left, Vsqrt(K)/2(du_right-du_left), -V/2, dP_right+dP_left, -V/2(dP_right+dP_left), Ut: \n%10f %10f %10f %10f %10f %10f %10f\n",V_old[i]*sqrt(K)/2.0, du_right-du_left, V_old[i]*sqrt(K)/2.0*(du_right-du_left),-V_old[i]/2.0, dP_right+dP_left, -V_old[i]/2.0*(dP_right+dP_left), V_old[i]*sqrt(K)/2.0*(du_right-du_left)-V_old[i]/2.0*(dP_right+dP_left));
			printf("-VK/2 du_right+du_left, -VK/2(du_right+du_left), Vsqrt(K)/2, dP_right-dP_left, Vsqrt(K)/2(dP_right-dP_left), Pt: \n%10f %10f %10f %10f %10f %10f %10f\n", -V_old[i]*K/2.0, du_right+du_left, -V_old[i]*K/2.0*(du_right+du_left), V_old[i]*sqrt(K)/2.0, dP_right-dP_left, V_old[i]*sqrt(K)/2.0*(dP_right-dP_left), -V_old[i]*K/2.0*(du_right+du_left)+V_old[i]*sqrt(K)/2.0*(dP_right-dP_left));
		}*/
        }

        updateFluidState();
        moveFluidParticle();
//	reorderFluidParticle();
        computeSetupsForNextIteration();
        if(m_iIfDebug) debug<<"-------------------------------------------------------"<<endl;
        return 0;


}



int HyperbolicLPSolver1D::solveUpwindPredictorCorrector(double dt) {

//      return solve_van_leer(dt);


        if(m_iIfDebug) debug<<"--------------HyperbolicLPSolver1D::solve()--------------"<<endl;

        // dt for this time step 
        m_fDt = dt;

        //alias 
        double* up_old = m_pParticleData->m_vVelocityU;
        double* V_old  = m_pParticleData->m_vVolume;
        double* p_old  = m_pParticleData->m_vPressure;
        double* cs_old = m_pParticleData->m_vSoundSpeed;

        double* up_pc   = m_pParticleData->m_vTemp2VelocityU;
        double* V_pc    = m_pParticleData->m_vTemp2Volume;
        double* p_pc    = m_pParticleData->m_vTemp2Pressure;
        double* cs_pc   = m_pParticleData->m_vTemp2SoundSpeed;
//        double* xp_old   = m_pParticleData->m_vTemp1PositionX;


        double* up     = m_pParticleData->m_vTemp1VelocityU;
        double* V      = m_pParticleData->m_vTemp1Volume;
        double* p      = m_pParticleData->m_vTemp1Pressure;
        double* cs     = m_pParticleData->m_vTemp1SoundSpeed;

        double* xp     = m_pParticleData->m_vPositionX;

        size_t np      = m_pParticleData->m_iTotalNum;
        //cout<<"np="<<np<<endl;

        size_t begin, end;
        if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") {
                begin = 1;
                end = np-2;
        }
        else if(m_sBoundaryType == "periodic") {
                begin = 0;
                end = np-1;
        }
        else {
                cout<<"ERROR: Invalid boundary type!"<<endl;
                exit(1);
        }

        // update inner particles
        for(size_t i=begin; i<=end; i++) {

                int left_order=0, right_order=0;
                computeLPFOrder(i, left_order, right_order);
                if(m_iIfDebug) {assert(left_order!=0); assert(right_order!=0);}

                //cout<<"i="<<i<<"      left_order="<<left_order<<"     right_order="<<right_order<<endl;

                double ux_right=0, ux_left=0, px_right=0, px_left=0;
                double uxx_right=0, uxx_left=0, pxx_right=0, pxx_left=0;
                computeSpatialDerUpwind(left_order, 0, i, up_old, p_old, xp, ux_left, px_left);//1: left nei
                computeSpatialDerUpwind(right_order, 1, i, up_old, p_old, xp, ux_right, px_right);//1: right nei

                //cout<<"i="<<i<<"      ux_left="<<px_left<<"   uxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_left="<<px_left<<"   pxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      ux_right="<<px_left<<"  uxx_right="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_right="<<px_left<<"  pxx_right="<<pxx_left<<endl;

                // update V, up, p
                timeIntegration(0.5*dt, i, V_old, up_old, p_old, cs_old,
                                ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right,
                                                V_pc, up_pc, p_pc); // output

                // check if there is invalid state
                bool isInvalid = (p_pc[i] < m_fInvalidPressure || V_pc[i] < m_fInvalidDensity);
                if(isInvalid) {
                        printInvalidState(i, ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right);
                        cout<<"ERROR: Invalid state occurs!"<<endl;
                        exit(1);
                }

                // update sound speed
                cs_pc[i] = m_pEOS->getSoundSpeed(p_pc[i],1./V_pc[i]);

        }
        updateFluidState_pc1();
        moveFluidParticle_pc1();
        computeSetupsForNextIteration();

        for(size_t i=begin; i<=end; i++) {

                int left_order=0, right_order=0;
                computeLPFOrder(i, left_order, right_order);
                if(m_iIfDebug) {assert(left_order!=0); assert(right_order!=0);}

                //cout<<"i="<<i<<"      left_order="<<left_order<<"     right_order="<<right_order<<endl;

                double ux_right=0, ux_left=0, px_right=0, px_left=0;
                double uxx_right=0, uxx_left=0, pxx_right=0, pxx_left=0;
                computeSpatialDerUpwind(left_order, 0, i, up_pc, p_pc, xp, ux_left, px_left);//1: left nei
                computeSpatialDerUpwind(right_order, 1, i, up_pc, p_pc, xp, ux_right, px_right);//1: right nei

                //cout<<"i="<<i<<"      ux_left="<<px_left<<"   uxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_left="<<px_left<<"   pxx_left="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      ux_right="<<px_left<<"  uxx_right="<<pxx_left<<endl;
                //cout<<"i="<<i<<"      px_right="<<px_left<<"  pxx_right="<<pxx_left<<endl;

                // update V, up, p
                timeIntegration(i, V_old, up_old, p_old, cs_old,
                                ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right,
                                                V, up, p); // output

                // check if there is invalid state
                bool isInvalid = (p[i] < m_fInvalidPressure || V[i] < m_fInvalidDensity);
                if(isInvalid) {
                        printInvalidState(i, ux_left, uxx_left, px_left, pxx_left, ux_right, uxx_right, px_right, pxx_right);
                        cout<<"ERROR: Invalid state occurs!"<<endl;
                        exit(1);
                }

                // update sound speed
                cs[i] = m_pEOS->getSoundSpeed(p[i],1./V[i]);

        }
        updateFluidState_pc2();
        moveFluidParticle_pc2();
//        reorderFluidParticle();
        computeSetupsForNextIteration();

        if(m_iIfDebug) debug<<"-------------------------------------------------------"<<endl;
        return 0;


}

int HyperbolicLPSolver1D::solveNewUpwind(double dt) {

//      return solve_van_leer(dt);

        if(m_iIfDebug) debug<<"--------------HyperbolicLPSolver1D::solve()--------------"<<endl;

        // dt for this time step 
        m_fDt = dt;

        //alias 
        double* up_old = m_pParticleData->m_vVelocityU;
        double* V_old  = m_pParticleData->m_vVolume;
        double* p_old  = m_pParticleData->m_vPressure;
        double* cs_old = m_pParticleData->m_vSoundSpeed;
        double* up     = m_pParticleData->m_vTemp1VelocityU;
        double* V      = m_pParticleData->m_vTemp1Volume;
        double* p      = m_pParticleData->m_vTemp1Pressure;
        double* cs     = m_pParticleData->m_vTemp1SoundSpeed;
        double* xp     = m_pParticleData->m_vPositionX;
	double* mass   = m_pParticleData->m_vMass;

        size_t np      = m_pParticleData->m_iTotalNum;
        //cout<<"np="<<np<<endl;

        size_t begin, end;
        if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") {
                begin = 1;
                end = np-2;
        }
        else if(m_sBoundaryType == "periodic") {
                begin = 0;
                end = np-1;
        }
        else {
                cout<<"ERROR: Invalid boundary type!"<<endl;
                exit(1);
        }

        // update inner particles
        for(size_t i=begin; i<=end; i++) {

		int order = 1;
		double ux=0, uxx=0, px=0, pxx=0;
		computeSpatialDerCenter(order, i, up_old, p_old, xp, ux, uxx, px, pxx);
                timeIntegrationCenter(i, V_old, up_old, p_old, cs_old,mass[i], 
                                ux, uxx, px, pxx,
                                                V, up, p); // output
//		if(i==799)
//                        cout<<i<<" x: "<<xp[i]<<" V_old: "<<V_old[i]<<" up_old: "<<up_old[i]<<" p_old: "<<p_old[i]<<" ux: "<<ux<<" uxx: "<<uxx<<" px: "<<px<<" pxx: "<<pxx<<" V: "<<V[i]<<" up: "<<up[i]<<" p: "<<p[i]<<endl;
                // check if there is invalid state
                bool isInvalid = (p[i] < m_fInvalidPressure || V[i] < m_fInvalidDensity || std::isnan(p[i]) || std::isnan(V[i]));
                if(isInvalid) {
			cout<<i<<" x: "<<xp[i]<<" V_old: "<<V_old[i]<<" up_old: "<<up_old[i]<<" p_old: "<<p_old[i]<<" ux: "<<ux<<" uxx: "<<uxx<<" px: "<<px<<" pxx: "<<pxx<<" V: "<<V[i]<<" up: "<<up[i]<<" p: "<<p[i]<<endl;
                        cout<<"ERROR: Invalid state occurs!"<<endl;
                        exit(1);
                }

                // update sound speed
                cs[i] = m_pEOS->getSoundSpeed(p[i],1./V[i]);
        }

        updateFluidState();
        moveFluidParticle();
//      reorderFluidParticle();
        computeSetupsForNextIteration();
        if(m_iIfDebug) debug<<"-------------------------------------------------------"<<endl;
        return 0;


}




void HyperbolicLPSolver1D::computeLPFOrder(size_t index, int& left_order, int& right_order) {
	
	if(m_iLPFOrder == 1) {
		left_order = 1, right_order = 1; 
	}
	else if(m_iLPFOrder == 2) {
		// alias
		const double* dd3_left  = m_pParticleData->m_vDD3Left;
		const double* dd3_right = m_pParticleData->m_vDD3Right;
		size_t  np              = m_pParticleData->m_iTotalNum;  

		left_order = 2, right_order = 2; // default

		if(m_iUseLimiter) {
			
			//double thres = 10;
			if(fabs(dd3_left[index])  > m_fThresholdP) // If the left 2 cells contain discontinuity, use 1st order GFD
				left_order = 1;
			if(fabs(dd3_right[index]) > m_fThresholdP) // If the right 2 cells contain discontinuity, use 1st order GFD
				right_order = 1;
		
			/*
			//double thres1 = 1, thres2 = 10, thres;
			//thres = xp[i]>0? thres2:thres1;
			double dd3left = xp[i]>0? dd3_left[i]:fabs(dd3_left[i]);
			double dd3right = xp[i]>0? dd3_right[i]:fabs(dd3_right[i]);
			if(dd3left > thres) // If the left 2 cells contain discontinuity, use 1st order GFD
				left_order = 1;
			if(dd3right > thres) // If the right 2 cells contain discontinuity, use 1st order GFD
				right_order = 1;
			*/
			//if(fabs(dd3_left[i]) > thres) // If the left 2 cells contain discontinuity, use 1st order GFD
			//	left_order = 1;
			//if(fabs(dd3_right[i]) > thres) // If the right 2 cells contain discontinuity, use 1st order GFD
			//	right_order = 1;
		}
		
		if(m_sBoundaryType == "free" || m_sBoundaryType == "solid") {
			if(index == 1) left_order = 1;
			else if(index == np-2) right_order = 1;
		}

	}
	else {
		cout<<"ERROR: m_iLPFOrder!=1 && m_iLPFOrder!=2!"<<endl;
		exit(1);
	}

}


void HyperbolicLPSolver1D::computeSpatialDer(int order, int direction, int i, 
											 const double *local_u, const double *local_x,
											 double& dudx, double& dudxx) { //output
	
	//alias
	int np = (int)m_pParticleData->m_iTotalNum;

	int num_nei = order;
	
	//double *u_nei, *xp_nei;
	//u_nei = new double[num_nei+1];
	//xp_nei = new double[num_nei+1];
	double u_nei[num_nei+1];
	double xp_nei[num_nei+1];
	
	// First element must be itself
	u_nei[0]  = local_u[i];
	xp_nei[0] = local_x[i];

	// Find neighbors - assuming periodic boundary
	if(direction == 0) //Wave propagating to the right, find neighbors on the left 
	{ 
		for(int j=0; j<num_nei; j++) 
		{
			if(i-j-1 >= 0) // if there are num_nei neighbors on the left
			{
				xp_nei[j+1] = local_x[i-j-1];
				u_nei[j+1] = local_u[i-j-1];
			}
			else // else take those rightmost particles as neighbors
			{
				int index_temp = np + (i-j-1); //Corresponding particle on the rightmost
				double xp_temp = local_x[index_temp] - m_fDisBetweenPeriodicBoundaries;
				if(xp_temp >= local_x[i]) assert(false);
				xp_nei[j+1] = xp_temp;
				u_nei[j+1] = local_u[index_temp];
			}
		}
	}
	else if(direction == 1) //Wave propagating to the left, find neighbors on the right
	{
		for(int j=0; j<num_nei; j++) 
		{
			if(i+j+1 <= np-1) // if there are num_nei neighbors on the right
			{
				xp_nei[j+1] = local_x[i+j+1];
				u_nei[j+1] = local_u[i+j+1];
			}
			else // else take those leftmost particles as neighbors
			{
				int index_temp = (i+j+1) - np; //Corresponding particle on the leftmost
				double xp_temp = local_x[index_temp] + m_fDisBetweenPeriodicBoundaries;
				if(xp_temp <= local_x[i]) assert(false);
				xp_nei[j+1] = xp_temp;
				u_nei[j+1] = local_u[index_temp];
			}
		}
	}

	//for(int k=0; k<=num_nei; k++) {
	//	cout<<"xp_nei["<<k<<"]="<<xp_nei[k]<<endl;
	//	cout<<"u_nei["<<k<<"]="<<u_nei[k]<<endl;
	//}

	if(order == 1) {
		if(direction == 0)
			dudx = (u_nei[0] - u_nei[1]) / (xp_nei[0] - xp_nei[1]);
		else if(direction == 1)
			dudx = (u_nei[1] - u_nei[0]) / (xp_nei[1] - xp_nei[0]);
		dudxx = 0;
	}
	else if(order == 2 || order == 3){
		//solveByCramer(order, num_nei, u_nei, xp_nei, dudx, dudxx);
		solveByQR(order, num_nei, u_nei, xp_nei, dudx, dudxx); 
	}
		

	//delete[] u_nei;
	//delete[] xp_nei;
}

void HyperbolicLPSolver1D::computeSpatialDerCenter(int order, int i, const double * local_u, const double *local_p, const double *local_x, double& dudx, double& dudxx, double& dpdx, double& dpdxx) { //output

	if(order>1)
	{
		cout<<"Error: try to use second order upwind method."<<endl;
		exit(1);
	}
	bool GhostParticlesForVacancy=1;
        double VacancyMultiplier = 4.0;
	int num_nei = 4;
	bool Vacancy_left=0;
	bool Vacancy_right=0;
        //alias
        int np = (int)m_pParticleData->m_iTotalNum;

	double mindx=local_x[np-1]-local_x[0];
	if (i<np&&(local_x[i+1]-local_x[i])<mindx)	mindx=local_x[i+1]-local_x[i];
        if (i>0&&(local_x[i]-local_x[i-1])<mindx)      mindx=local_x[i]-local_x[i-1];


	double averagedx=0;
	int i1=i+5;
	int i2=i-5;
	if (i1>=np) i1=np-1;
	if (i2<0) i2=0;
	averagedx=(local_x[i1]-local_x[i2])/(i1-i2);

        double u_nei[num_nei+1];
	double p_nei[num_nei+1];
        double xp_nei[num_nei+1];

        // First element must be itself
        u_nei[0]  = local_u[i];
	p_nei[0]  = local_p[i];
        xp_nei[0] = local_x[i];

	//Vacancy=..
		Vacancy_left=0;
		Vacancy_right=0;
		if(i>0)
		{
                                xp_nei[1] = local_x[i-1];
                                u_nei[1] = local_u[i-1];
				p_nei[1] = local_p[i-1];
//				if((local_x[i]-local_x[i-1])>mindx*VacancyMultiplier)
				if((local_x[i]-local_x[i-1])>averagedx*VacancyMultiplier)
					Vacancy_left=GhostParticlesForVacancy;
		}
		if(i==0 ||Vacancy_left)
		{
			xp_nei[1]=2*local_x[i]-local_x[i+1];
			u_nei[1]=local_u[i];
			p_nei[1]=0.0;
			cout<<"Use ghost particle as left neighbour for "<<i<<"th particle at x = "<<local_x[i]<<"."<<endl;
		}

                if(i<(np-1))
                {
                                xp_nei[2] = local_x[i+1];
                                u_nei[2] = local_u[i+1];
                                p_nei[2] = local_p[i+1];
//                                if((local_x[i+1]-local_x[i])>mindx*VacancyMultiplier)
                                if((local_x[i+1]-local_x[i])>averagedx*VacancyMultiplier)
                                        Vacancy_right=GhostParticlesForVacancy;
		}
                if(i==(np-1) ||Vacancy_right)
                {
                        xp_nei[2]=2*local_x[i]-local_x[i-1];
                        u_nei[2]=local_u[i];
                        p_nei[2]=0.0;
                        cout<<"Use ghost particle as right neighbour for "<<i<<"th particle at x = "<<local_x[i]<<"."<<endl;
                }
		if(i>1 && i<(np-2) &&Vacancy_left==0 &&Vacancy_right==0)
		{
			xp_nei[3] = local_x[i-2];
                        u_nei[3] = local_u[i-2];
                        p_nei[3] = local_p[i-2];
                        xp_nei[4] = local_x[i+2];
                        u_nei[4] = local_u[i+2];
                        p_nei[4] = local_p[i+2];
		}
		else
			num_nei=2;

                solveByQR(2, num_nei, u_nei, xp_nei, dudx, dudxx);
                solveByQR(2, num_nei, p_nei, xp_nei, dpdx, dpdxx);

}

void HyperbolicLPSolver1D::computeSpatialDerUpwind(int order, int direction, int i, const double * local_u, const double *local_p, const double *local_x, double& dudx, double& dpdx) { //output

        if(order>1)
        {
                cout<<"Error: try to use second order upwind method."<<endl;
                exit(1);
        }
        bool GhostParticlesForVacancy=1;
        double VacancyMultiplier = 4.0;
        bool Vacancy=0;
        //alias
        int np = (int)m_pParticleData->m_iTotalNum;

        double mindx=local_x[np-1]-local_x[0];
        if (i<np&&(local_x[i+1]-local_x[i])<mindx)      mindx=local_x[i+1]-local_x[i];
        if (i>0&&(local_x[i]-local_x[i-1])<mindx)      mindx=local_x[i]-local_x[i-1];


        double averagedx=0;
        int i1=i+5;
        int i2=i-5;
        if (i1>=np) i1=np-1;
        if (i2<0) i2=0;
        averagedx=(local_x[i1]-local_x[i2])/(i1-i2);

        double u_nei[2];
        double p_nei[2];
        double xp_nei[2];

        // First element must be itself
        u_nei[0]  = local_u[i];
        p_nei[0]  = local_p[i];
        xp_nei[0] = local_x[i];

        //Vacancy=..

        // Find neighbors - assuming periodic boundary
        if(direction == 0) //Wave propagating to the right, find neighbors on the left 
        {

                if(i>0)
                {
                                xp_nei[1] = local_x[i-1];
                                u_nei[1] = local_u[i-1];
                                p_nei[1] = local_p[i-1];
//                              if((local_x[i]-local_x[i-1])>mindx*VacancyMultiplier)
                                if((local_x[i]-local_x[i-1])>averagedx*VacancyMultiplier)
                                        Vacancy=GhostParticlesForVacancy;
                }
                if(i==0 ||Vacancy)
                {
                        xp_nei[1]=2*local_x[i]-local_x[i+1];
                        u_nei[1]=local_u[i];
                        p_nei[1]=0.0;
                        cout<<"Use ghost particle as left neighbour for "<<i<<"th particle at x = "<<local_x[i]<<"."<<endl;
                }
        }
        else if(direction == 1) //Wave propagating to the left, find neighbors on the right
        {

                if(i<(np-1))
                {
                                xp_nei[1] = local_x[i+1];
                                u_nei[1] = local_u[i+1];
                                p_nei[1] = local_p[i+1];
//                                if((local_x[i+1]-local_x[i])>mindx*VacancyMultiplier)
                                if((local_x[i+1]-local_x[i])>averagedx*VacancyMultiplier)
                                        Vacancy=GhostParticlesForVacancy;
                }
                if(i==(np-1) ||Vacancy)
                {
                        xp_nei[1]=2*local_x[i]-local_x[i-1];
                        u_nei[1]=local_u[i];
                        p_nei[1]=0.0;
                        cout<<"Use ghost particle as right neighbour for "<<i<<"th particle at x = "<<local_x[i]<<"."<<endl;
                }
        }

        dudx = (u_nei[0] - u_nei[1]) / (xp_nei[0] - xp_nei[1]);
        dpdx = (p_nei[0] - p_nei[1]) / (xp_nei[0] - xp_nei[1]);
}



void HyperbolicLPSolver1D::solveByQR(int order, int num_nei, const double *local_u, const double *local_x, 
									double &dudx, double &dudxx) { //output
	using std::isnan; // avoid ambiguity
	
	int num_row = num_nei;
    int num_col = order;
    
    if(order == 2) { // 2nd order GFD
	int info=1;
	num_row=order;
	while(info)
	{
//		cout<<local_x[0]<<" "<<num_row<<endl;
                double A[num_row*num_col];
                double b[num_row];
                double normal = 0.5*(fabs(local_x[1]-local_x[0])+fabs(local_x[2]-local_x[0]));
                for(int k=1; k<num_row+1 ; k++) {
                        double h = (local_x[k] - local_x[0])/normal;
                    A[(k-1)] = h;
                    A[(k-1)+1*(num_row)] = 0.5 * h * h;
                    b[k-1] = local_u[k] - local_u[0];
                }
                unique_ptr<LSSolver> solver(new QRSolver(num_row,num_col,A));
                double result[num_col];
                info = solver->solve(result,b);
                if(info == 0) { // QR succeeds
                        dudx  = result[0]/normal;
                        dudxx = result[1]/normal/normal;
                        if(isnan(dudx) || isnan(dudxx)) {
                                cout<<"ERROR: nan derivatives!"<<endl;
                                exit(1);
                        }
			if(num_row>order)
				cout<<num_row <<" neighbours are used to calculate the derivative of particle x = "<<local_x[0]<<endl;
                }
		else{ //QR fails, increase the number of neighbors
			num_row++;
//                        cout<<A[0]<<" "<<A[1]<<" "<<A[2]<<" "<<A[3]<<endl;
			if(num_row>num_nei){
				cout<<local_x[0]<<" "<<num_nei<<""<<endl;
				cout<<local_x[1]-local_x[0]<<" "<<local_x[2]-local_x[0]<<" "<<local_x[3]-local_x[0]<<" "<<local_x[4]-local_x[0]<<endl;
				cout<<A[0]<<" "<<A[1+4]<<" "<<A[2+8]<<" "<<A[3+12]<<endl;
	                        cout<<"ERROR: QR decomposition has insufficient rank for order = "<<order<<"!"<<endl;
				exit(1);
			}
		}
	}		
		//double* A = new double[num_row*num_col];
		//double* b = new double[num_row];
/*		double A[num_row*num_col];
		double b[num_row];
		double normal = local_x[1]-local_x[0];
		for(int k=1; k<num_nei+1 ; k++) {
			double h = (local_x[k] - local_x[0])/normal;
		    A[(k-1)] = h;
		    A[(k-1)+1*(num_row)] = 0.5 * h * h;
		    b[k-1] = local_u[k] - local_u[0];
		}
		unique_ptr<LSSolver> solver(new QRSolver(num_row,num_col,A));
		double result[num_col];
		int info = solver->solve(result,b);
		if(info == 0) { // QR succeeds
			dudx  = result[0]/normal;
			dudxx = result[1]/normal/normal;	
			if(isnan(dudx) || isnan(dudxx)) {
				cout<<"ERROR: nan derivatives!"<<endl;
				exit(1);
			}
		}
		else { // QR fails
			cout<<"ERROR: QR decomposition has insufficient rank for order = "<<order<<"!"<<endl;
			cout<<local_u[0]<<", "<<local_u[1]<<", "<<local_u[2]<<endl;
			cout<<local_x[0]<<", "<<local_x[1]<<", "<<local_x[2]<<endl;
			cout<<A[0]<<" "<<A[1]<<" "<<A[2]<<" "<<A[3]<<endl;
			exit(1);
		}*/

		//Lapack_qr solver(num_row, num_col, A);
		//int status = solver.QR_decomposition();
		//if(status == 0) { // QR succeeds
		//	double* der_results = new double[num_col];
		//	solver.solve_ls(b, der_results);
		//	dudx = der_results[0];
		//	dudxx = der_results[1];
		//	delete[] der_results;
		//	// Check point
		//	if(isnan(dudx) || isnan(dudxx)) {
		//		cout<<"CHECK POINT: Lapack_qr::solve_ls() gets nan derivatives!!!"<<endl;
		//		assert(false);
		//	}
		//}
		//else {
		//	// Check point
		//	cout<<"CHECK POINT: QR fails (status != 0)!!!"<<endl;
		//	assert(false);
		//}
		//
		//delete[] A;
		//delete[] b;
		
	}
	else if (order == 3) { // 3rd order GFD
		
		//double* A = new double[num_row*num_col];
		//double* b = new double[num_row];
		double A[num_row*num_col];
		double b[num_row];
		for(int k=1; k<num_nei+1 ; k++) {
			double h = (local_x[k] - local_x[0]);
		  A[(k-1)] = h;
		  A[(k-1)+1*(num_row)] = 0.5 * h * h;
		  A[(k-1)+2*(num_row)] = 1. / 6. * h * h * h;
		  b[k-1] = local_u[k] - local_u[0];
		}
		unique_ptr<LSSolver> solver(new QRSolver(num_row,num_col,A));
		double result[num_col];
		int info = solver->solve(result,b);
		if(info == 0) { // QR succeeds
			dudx  = result[0];
			dudxx = result[1];	
			if(isnan(dudx) || isnan(dudxx)) {
				cout<<"ERROR: nan derivatives!"<<endl;
				exit(1);
			}
		}
		else { // QR fails
			cout<<"ERROR: QR decomposition has insufficient rank for order = "<<order<<"!"<<endl;
			exit(1);
		}
		
		
		//Lapack_qr solver(num_row, num_col, A);
		//int status = solver.QR_decomposition();
		//if(status == 0) { // QR succeeds
		//	double* der_results = new double[num_col];
		//	solver.solve_ls(b, der_results);
		//	dudx = der_results[0];
		//	dudxx = der_results[1];
		//	delete[] der_results;
		//	// Check point
		//	if(isnan(dudx) || isnan(dudxx)) {
		//		cout<<"CHECK POINT: Lapack_qr::solve_ls() gets nan derivatives!!!"<<endl;
		//		assert(false);
		//	}
		//}
		//else {
		//	// Check point
		//	cout<<"CHECK POINT: QR fails (status != 0)!!!"<<endl;
		//	assert(false);
		//}
		//
		//delete[] A;
		//delete[] b;
			
	}													

	
}


double HyperbolicLPSolver1D::constantWeight(double h) {
	return 1;
}

void HyperbolicLPSolver1D::solveByCramer(int order, int num_nei, const double *local_u, const double *local_x, 
									     double &dudx, double &dudxx) { //output
	
	// set up a function pointer to the weight function
	double (HyperbolicLPSolver1D::*weight)(double) = &HyperbolicLPSolver1D::constantWeight;

	if(order == 2) { // 2nd order GFD
		
		//-------------  | a11  a12 | |   du/dx   | = b1 -------------------
		//-------------  | a21  a22 | |d^2 u/d x^2| = b2 -------------------
		double a11=0; double a12=0;
		double a21=0; double a22=0;
		double b1=0;  double b2=0;

		for(int k=1; k<num_nei+1 ; k++) {
			double qq = (local_x[k] - local_x[0]);
			//double ww = weight(qq); // weight
		    double ww = (this->*weight)(qq); // weight	
			ww *= ww;
		  
			a11 += qq*qq*ww;
			a21 += qq*qq*qq*ww/2.;
			a22 += qq*qq*qq*qq*ww/4.;
		
			double uu = local_u[k] - local_u[0];
			b1 += uu*qq*ww;	
			b2 += uu*qq*qq*ww/2.;

			//cout<<"ww="<<ww<<endl;
			//cout<<"qq="<<qq<<endl;
			//cout<<"uu="<<uu<<endl;
		}
		
		a12 = a21;
		double delta = a11*a22 - a12*a21;
		dudx  = (a22*b1 - a12 *b2)/delta;
		dudxx = (a11*b2 - a21 *b1)/delta;
		
//		if (time-dt == 0) {
//			printf("%.16g %.16g %.16g %.16g\n",local_x[0],delta,dudx,dudxx);
//		}
		
	}
	else if (order == 3) { // 3rd order GFD
		
		//-------------  | a11  a12 a13 | |   du/dx   | = b1 -------------------
		//-------------  | a21  a22 a23 | |d^2 u/d x^2| = b2 -------------------
		//-------------  | a31  a32 a33 | |d^3 u/d x^3| = b3 -------------------
		double a11=0; double a12=0; double a13=0;
		double a21=0; double a22=0; double a23=0;
		double a31=0; double a32=0; double a33=0;
		double b1 =0; double b2 =0; double b3 =0;

		for(int k=1; k<num_nei+1 ; k++) {
			double qq = (local_x[k] - local_x[0]); 
		    //double ww = weight(qq); // weight
		    double ww = (this->*weight)(qq); // weight
			ww *= ww;
		  
			a11 += qq*qq*ww;
			a12 += qq*qq*qq*ww/2.;
			a13 += qq*qq*qq*qq*ww/6.;
			a22 += qq*qq*qq*qq*ww/4.;
			a23 += qq*qq*qq*qq*qq*ww/12.;
			a33 += qq*qq*qq*qq*qq*qq*ww/36.;
			
			double uu = local_u[k] - local_u[0];
			b1 += uu*qq*ww;	
			b2 += uu*qq*qq*ww/2.;
			b3 += uu*qq*qq*qq*ww/6.;
		}
		
		a21 = a12;
		a31 = a13;
		a32 = a23;
		double delta = a11*a22*a33 + a21*a32*a13 + a12*a23*a31 - a13*a22*a31 - a23*a32*a11 - a12*a21*a33;
		dudx  = (b1*a22*a33 + b2*a32*a13 + a12*a23*b3 - a13*a22*b3 - a23*a32*b1 - a12*b2*a33)/delta;
		dudxx = (a11*b2*a33 + a21*b3*a13 + b1*a23*a31 - a13*b2*a31 - a23*b3*a11 - b1*a21*a33)/delta;
	}
	
	using std::isnan; // avoid ambiguity
	if(isnan(dudx) || isnan(dudxx)) {
		cout<<"ERROR: nan derivatives!"<<endl;
		exit(1);
	}

}

void HyperbolicLPSolver1D::timeIntegration(int i, const double* V_old, const double* up_old,
										   const double* p_old, const double* cs_old,	
										   double ux_left,  double uxx_left,
										   double px_left,  double pxx_left,
										   double ux_right, double uxx_right,
										   double px_right, double pxx_right,
										   double* V, double* up, double* p) { // output
	//alias
	double dt = m_fDt;

	double eosCoeff;
	//if(eos_choice == poly || eos_choice == spoly) // poly || spoly
		eosCoeff = cs_old[i]*cs_old[i] / (V_old[i]*V_old[i]);
	//else assert(false);
		
	double tmp1 = 0.5*V_old[i];
	double tmp2 = -0.5*V_old[i]/sqrt(eosCoeff);
	double tmp3 = 0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
	double tmp4 = -0.25*V_old[i]*V_old[i];
	V[i]  = V_old[i]  + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
			dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	
	tmp1 = 0.5*V_old[i]*sqrt(eosCoeff);
	tmp2 = -0.5*V_old[i];
	tmp3 = 0.25*V_old[i]*V_old[i]*eosCoeff;
	tmp4 = -0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
	up[i] = up_old[i] + dt * (tmp1*(ux_right-ux_left)+tmp2*(px_right+px_left)) +
			dt * dt * (tmp3*(uxx_right+uxx_left)+tmp4*(pxx_right-pxx_left));
	
	tmp1 = -0.5*V_old[i]*eosCoeff;
	tmp2 = 0.5*V_old[i]*sqrt(eosCoeff);
	tmp3 = -0.25*V_old[i]*V_old[i]*pow(eosCoeff,1.5);
	tmp4 = 0.25*V_old[i]*V_old[i]*eosCoeff;
	p[i]  = p_old[i] + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
			dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	

	//cout<<"dt="<<dt<<"	eosCoeff="<<eosCoeff<<endl;
	//cout<<"V_old="<<V_old[i]<<"	V="<<V[i]<<endl;
	//cout<<"up_old="<<up_old[i]<<" up="<<up[i]<<endl;
	//cout<<"p_old="<<p_old[i]<<"	p="<<p[i]<<endl;

}

void HyperbolicLPSolver1D::timeIntegration(double dt, int i, const double* V_old, const double* up_old,
                                                                                   const double* p_old, const double* cs_old,
                                                                                   double ux_left,  double uxx_left,
                                                                                   double px_left,  double pxx_left,
                                                                                   double ux_right, double uxx_right,
                                                                                   double px_right, double pxx_right,
                                                                                   double* V, double* up, double* p) { // output
        //alias
//        double dt = m_fDt;

        double eosCoeff;
        //if(eos_choice == poly || eos_choice == spoly) // poly || spoly
                eosCoeff = cs_old[i]*cs_old[i] / (V_old[i]*V_old[i]);
        //else assert(false);

        double tmp1 = 0.5*V_old[i];
        double tmp2 = -0.5*V_old[i]/sqrt(eosCoeff);
        double tmp3 = 0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
        double tmp4 = -0.25*V_old[i]*V_old[i];
        V[i]  = V_old[i]  + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
                        dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));

        tmp1 = 0.5*V_old[i]*sqrt(eosCoeff);
        tmp2 = -0.5*V_old[i];
        tmp3 = 0.25*V_old[i]*V_old[i]*eosCoeff;
        tmp4 = -0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
        up[i] = up_old[i] + dt * (tmp1*(ux_right-ux_left)+tmp2*(px_right+px_left)) +
                        dt * dt * (tmp3*(uxx_right+uxx_left)+tmp4*(pxx_right-pxx_left));

        tmp1 = -0.5*V_old[i]*eosCoeff;
        tmp2 = 0.5*V_old[i]*sqrt(eosCoeff);
        tmp3 = -0.25*V_old[i]*V_old[i]*pow(eosCoeff,1.5);
        tmp4 = 0.25*V_old[i]*V_old[i]*eosCoeff;
        p[i]  = p_old[i] + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
                        dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
}

void HyperbolicLPSolver1D::timeIntegrationCenter(int i, const double* V_old, const double* up_old,
                                                                                   const double* p_old, const double* cs_old, double mass,  
                                                                                   double ux,  double uxx,
                                                                                   double px,  double pxx,
                                                                                   double* V, double* up, double* p) { // output
        //alias
	double h=mass*V_old[i];
        double dt = m_fDt;

        double eosCoeff;
        //if(eos_choice == poly || eos_choice == spoly) // poly || spoly
                eosCoeff = cs_old[i]*cs_old[i] / (V_old[i]*V_old[i]);
        //else assert(false);

        double tmp1 = V_old[i];
        double tmp2 = -0.5*V_old[i]/sqrt(eosCoeff)*h;
        V[i]  = V_old[i]  + dt * (tmp1*ux+tmp2*pxx);

        tmp1 = 0.5*V_old[i]*sqrt(eosCoeff)*h;
        tmp2 = -V_old[i];
        up[i] = up_old[i] + dt * (tmp1*uxx+tmp2*px);

        tmp1 = -V_old[i]*eosCoeff;
        tmp2 = 0.5*V_old[i]*sqrt(eosCoeff)*h;
        p[i]  = p_old[i] + dt * (tmp1*ux+tmp2*pxx);
}

void HyperbolicLPSolver1D::timeIntegration_van_leer(
	int i, const double* V_old, const double* up_old, const double* p_old, const double* cs_old, double phi,	
	double ux_left1, double px_left1, double ux_right1, double px_right1, 
	double ux_left2,  double uxx_left, double px_left2,  double pxx_left,
	double ux_right2, double uxx_right, double px_right2, double pxx_right,
	double* V, double* up, double* p) { // output	
	
	//alias
	double dt = m_fDt;

	double eosCoeff;
	//if(eos_choice == poly || eos_choice == spoly) // poly || spoly
		eosCoeff = cs_old[i]*cs_old[i] / (V_old[i]*V_old[i]);
	//else assert(false);
	
	double flux = 0;

	double tmp1 = 0.5*V_old[i];
	double tmp2 = -0.5*V_old[i]/sqrt(eosCoeff);
	double tmp3 = 0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
	double tmp4 = -0.25*V_old[i]*V_old[i];
	//double Vt1  = dt * (tmp1*(ux_left1+ux_right1)+tmp2*(px_right1-px_left1));
	//double Vt2  = dt * (tmp1*(ux_left2+ux_right2)+tmp2*(px_right2-px_left2)) +
	//		      dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	//flux = Vt1 + phi * (Vt2-Vt1);
	//V[i] = V_old[i] + flux;
	double Vt1  = (tmp1*(ux_left1+ux_right1)+tmp2*(px_right1-px_left1));
	double Vt2  = (tmp1*(ux_left2+ux_right2)+tmp2*(px_right2-px_left2)) +
			      dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	flux = Vt1 + phi * (Vt2-Vt1);
	V[i] = V_old[i] + dt * flux;
	//V[i]  = V_old[i]  + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
	//		dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	
	tmp1 = 0.5*V_old[i]*sqrt(eosCoeff);
	tmp2 = -0.5*V_old[i];
	tmp3 = 0.25*V_old[i]*V_old[i]*eosCoeff;
	tmp4 = -0.25*V_old[i]*V_old[i]*sqrt(eosCoeff);
	//double ut1 = dt * (tmp1*(ux_right1-ux_left1)+tmp2*(px_right1+px_left1));
	//double ut2 = dt * (tmp1*(ux_right2-ux_left2)+tmp2*(px_right2+px_left2)) +
	//		     dt * dt * (tmp3*(uxx_right+uxx_left)+tmp4*(pxx_right-pxx_left));
	//flux = ut1 + phi * (ut2-ut1);
	//up[i] = up_old[i] + flux;
	double ut1 = (tmp1*(ux_right1-ux_left1)+tmp2*(px_right1+px_left1));
	double ut2 = (tmp1*(ux_right2-ux_left2)+tmp2*(px_right2+px_left2)) +
			     dt * (tmp3*(uxx_right+uxx_left)+tmp4*(pxx_right-pxx_left));
	flux = ut1 + phi * (ut2-ut1);
	up[i] = up_old[i] + dt * flux;
	//up[i] = up_old[i] + dt * (tmp1*(ux_right-ux_left)+tmp2*(px_right+px_left)) +
	//		dt * dt * (tmp3*(uxx_right+uxx_left)+tmp4*(pxx_right-pxx_left));
	
	tmp1 = -0.5*V_old[i]*eosCoeff;
	tmp2 = 0.5*V_old[i]*sqrt(eosCoeff);
	tmp3 = -0.25*V_old[i]*V_old[i]*pow(eosCoeff,1.5);
	tmp4 = 0.25*V_old[i]*V_old[i]*eosCoeff;
	//double pt1 = dt * (tmp1*(ux_left1+ux_right1)+tmp2*(px_right1-px_left1));
	//double pt2 = dt * (tmp1*(ux_left2+ux_right2)+tmp2*(px_right2-px_left2)) +
	//		     dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	//flux = pt1 + phi * (pt2-pt1);
	//p[i] = p_old[i] + flux;
	double pt1 = (tmp1*(ux_left1+ux_right1)+tmp2*(px_right1-px_left1));
	double pt2 = (tmp1*(ux_left2+ux_right2)+tmp2*(px_right2-px_left2)) +
			     dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	flux = pt1 + phi * (pt2-pt1);
	p[i] = p_old[i] + dt * flux;
	//p[i]  = p_old[i] + dt * (tmp1*(ux_left+ux_right)+tmp2*(px_right-px_left)) +
	//		dt * dt * (tmp3*(uxx_right-uxx_left)+tmp4*(pxx_right+pxx_left));
	
	//cout<<"dt="<<dt<<"	eosCoeff="<<eosCoeff<<endl;
	//cout<<"V_old="<<V_old[i]<<"	V="<<V[i]<<endl;
	//cout<<"up_old="<<up_old[i]<<" up="<<up[i]<<endl;
	//cout<<"p_old="<<p_old[i]<<"	p="<<p[i]<<endl;

}



void HyperbolicLPSolver1D::printInvalidState(int i, 
					   double ux_left, double uxx_left, double px_left, double pxx_left, 
					   double ux_right, double uxx_right, double px_right, double pxx_right) {
	
/*	debug<<"-----------------HyperbolicLPSolver1D::printInvalidState()------------------"<<endl;
	debug<<"index="<<i<<"  x="<<m_pParticleData->m_vPositionX[i]<<endl;
	debug<<"p_old = "<<m_pParticleData->m_vPressure[i]
	     <<"  p = "<<m_pParticleData->m_vTemp1Pressure[i]<<endl;
	debug<<"V_old = "<<m_pParticleData->m_vVolume[i]
	     <<"  V = "<<m_pParticleData->m_vTemp1Volume[i]<<endl;
	debug<<"up_old = "<<m_pParticleData->m_vVelocityU[i]
	     <<"  up = "<<m_pParticleData->m_vTemp1VelocityU[i]<<endl;
	debug<<"ux_left = "<<ux_left<<setw(15)<<"ux_right = "<<ux_right<<setw(15)<<endl;
	debug<<"px_left = "<<px_left<<setw(15)<<"px_right = "<<px_right<<setw(15)<<endl;
	debug<<"uxx_left = "<<uxx_left<<setw(15)<<"uxx_right = "<<uxx_right<<setw(15)<<endl;
	debug<<"pxx_left = "<<pxx_left<<setw(15)<<"pxx_right = "<<pxx_right<<setw(15)<<endl;
	debug<<"----------------------------------------------------------------------------"<<endl;
*/
        cout<<"-----------------HyperbolicLPSolver1D::printInvalidState()------------------"<<endl;
        cout<<"index="<<i<<"  x="<<m_pParticleData->m_vPositionX[i]<<endl;
        cout<<"p_old = "<<m_pParticleData->m_vPressure[i]
             <<"  p = "<<m_pParticleData->m_vTemp1Pressure[i]<<endl;
        cout<<"V_old = "<<m_pParticleData->m_vVolume[i]
             <<"  V = "<<m_pParticleData->m_vTemp1Volume[i]<<endl;
        cout<<"up_old = "<<m_pParticleData->m_vVelocityU[i]
             <<"  up = "<<m_pParticleData->m_vTemp1VelocityU[i]<<endl;
        cout<<"ux_left = "<<ux_left<<setw(15)<<"ux_right = "<<ux_right<<setw(15)<<endl;
        cout<<"px_left = "<<px_left<<setw(15)<<"px_right = "<<px_right<<setw(15)<<endl;
        cout<<"uxx_left = "<<uxx_left<<setw(15)<<"uxx_right = "<<uxx_right<<setw(15)<<endl;
        cout<<"pxx_left = "<<pxx_left<<setw(15)<<"pxx_right = "<<pxx_right<<setw(15)<<endl;
        cout<<"----------------------------------------------------------------------------"<<endl;


}


void HyperbolicLPSolver1D::updateFluidState() {
		
	swap(m_pParticleData->m_vTemp1Volume,     m_pParticleData->m_vVolume);
	swap(m_pParticleData->m_vTemp1Pressure,   m_pParticleData->m_vPressure);
	swap(m_pParticleData->m_vTemp1SoundSpeed, m_pParticleData->m_vSoundSpeed);	
	swap(m_pParticleData->m_vTemp1VelocityU,  m_pParticleData->m_vVelocityU);
}

void HyperbolicLPSolver1D::updateFluidState_pc1() {

        swap(m_pParticleData->m_vTemp2Volume,     m_pParticleData->m_vVolume);
        swap(m_pParticleData->m_vTemp2Pressure,   m_pParticleData->m_vPressure);
        swap(m_pParticleData->m_vTemp2SoundSpeed, m_pParticleData->m_vSoundSpeed);
        swap(m_pParticleData->m_vTemp2VelocityU,  m_pParticleData->m_vVelocityU);
}

void HyperbolicLPSolver1D::updateFluidState_pc2() {

        swap(m_pParticleData->m_vTemp1Volume,     m_pParticleData->m_vVolume);
        swap(m_pParticleData->m_vTemp1Pressure,   m_pParticleData->m_vPressure);
        swap(m_pParticleData->m_vTemp1SoundSpeed, m_pParticleData->m_vSoundSpeed);
        swap(m_pParticleData->m_vTemp1VelocityU,  m_pParticleData->m_vVelocityU);
}

void HyperbolicLPSolver1D::moveFluidParticle() {
	
	for(size_t index=0; index<m_pParticleData->m_iTotalNum; index++)
		m_pParticleData->m_vPositionX[index] += 0.5 * m_fDt * 
			(m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp1VelocityU[index]); // 0.5 (old + new)
}

void HyperbolicLPSolver1D::moveFluidParticle_pc1() {

        for(size_t index=0; index<m_pParticleData->m_iTotalNum; index++)
	{
		m_pParticleData->m_vTemp1PositionX[index] = m_pParticleData->m_vPositionX[index];
                m_pParticleData->m_vPositionX[index] += 0.25 * m_fDt *
                        (m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp2VelocityU[index]); // 0.5 (old + new)
	}
}

void HyperbolicLPSolver1D::moveFluidParticle_pc2() {

        for(size_t index=0; index<m_pParticleData->m_iTotalNum; index++)
                m_pParticleData->m_vPositionX[index] = m_pParticleData->m_vTemp1PositionX[index] + 0.5 * m_fDt *
                        (m_pParticleData->m_vVelocityU[index] + m_pParticleData->m_vTemp2VelocityU[index]); // 0.5 (old + new)
}




void HyperbolicLPSolver1D::reorderFluidParticle() {

	bool swap=1;
	double temp;

	while(swap)
	{
	swap=0;
		for(size_t index=1; index<m_pParticleData->m_iTotalNum; index++)
		{
			if(m_pParticleData->m_vPositionX[index]<m_pParticleData->m_vPositionX[index-1])
				{
				cout<<"Warning: particle penetrate each other. Reorder particle "<<index-1<<" and particle "<<index<<endl;
				swap=1;
				temp=m_pParticleData->m_vPositionX[index];
				m_pParticleData->m_vPositionX[index]=m_pParticleData->m_vPositionX[index-1];
				m_pParticleData->m_vPositionX[index-1]=temp;
                                temp=m_pParticleData->m_vVolume[index];
                                m_pParticleData->m_vVolume[index]=m_pParticleData->m_vVolume[index-1];
                                m_pParticleData->m_vVolume[index-1]=temp;
                                temp=m_pParticleData->m_vPressure[index];
                                m_pParticleData->m_vPressure[index]=m_pParticleData->m_vPressure[index-1];
                                m_pParticleData->m_vPressure[index-1]=temp;
                                temp=m_pParticleData->m_vSoundSpeed[index];
                                m_pParticleData->m_vSoundSpeed[index]=m_pParticleData->m_vSoundSpeed[index-1];
                                m_pParticleData->m_vSoundSpeed[index-1]=temp;
                                temp=m_pParticleData->m_vVelocityU[index];
                                m_pParticleData->m_vVelocityU[index]=m_pParticleData->m_vVelocityU[index-1];
                                m_pParticleData->m_vVelocityU[index-1]=temp;
                                temp=m_pParticleData->m_vMass[index];
                                m_pParticleData->m_vMass[index]=m_pParticleData->m_vMass[index-1];
                                m_pParticleData->m_vMass[index-1]=temp;
				}
		}
	}

	bool merge=1;
	double merge_tolerance=10;

	while(merge)
	{
		merge=0;
		for(size_t index=2; index<(m_pParticleData->m_iTotalNum-1); index++)
		{
			if((merge_tolerance*(m_pParticleData->m_vPositionX[index]-m_pParticleData->m_vPositionX[index-1])<(m_pParticleData->m_vPositionX[index-1]-m_pParticleData->m_vPositionX[index-2]))&&(merge_tolerance*(m_pParticleData->m_vPositionX[index]-m_pParticleData->m_vPositionX[index-1])<(m_pParticleData->m_vPositionX[index+1]-m_pParticleData->m_vPositionX[index])))
			{
				cout<<"Warning: particles are very close. Merge particle "<<index-1<<" and particle "<<index<<endl;
				merge=1;
				m_pParticleData->m_vPositionX[index-1]=(m_pParticleData->m_vMass[index-1]*m_pParticleData->m_vPositionX[index-1]+m_pParticleData->m_vMass[index]*m_pParticleData->m_vPositionX[index])/(m_pParticleData->m_vMass[index-1]+m_pParticleData->m_vMass[index]);
                                m_pParticleData->m_vVolume[index-1]=(m_pParticleData->m_vMass[index-1]*m_pParticleData->m_vVolume[index-1]+m_pParticleData->m_vMass[index]*m_pParticleData->m_vVolume[index])/(m_pParticleData->m_vMass[index-1]+m_pParticleData->m_vMass[index]);
                                m_pParticleData->m_vPressure[index-1]=(m_pParticleData->m_vMass[index-1]*m_pParticleData->m_vPressure[index-1]+m_pParticleData->m_vMass[index]*m_pParticleData->m_vPressure[index])/(m_pParticleData->m_vMass[index-1]+m_pParticleData->m_vMass[index]);
                                m_pParticleData->m_vVelocityU[index-1]=(m_pParticleData->m_vMass[index-1]*m_pParticleData->m_vVelocityU[index-1]+m_pParticleData->m_vMass[index]*m_pParticleData->m_vVelocityU[index])/(m_pParticleData->m_vMass[index-1]+m_pParticleData->m_vMass[index]);
				m_pParticleData->m_vSoundSpeed[index-1]=m_pEOS->getSoundSpeed(m_pParticleData->m_vPressure[index-1],1./m_pParticleData->m_vVolume[index-1]);
				m_pParticleData->m_vMass[index-1]=m_pParticleData->m_vMass[index-1]+m_pParticleData->m_vMass[index];
				m_pParticleData->m_iTotalNum=m_pParticleData->m_iTotalNum-1;
				for(size_t index2=index; index2<m_pParticleData->m_iTotalNum; index2++)
				{
					m_pParticleData->m_vPositionX[index2]=m_pParticleData->m_vPositionX[index2+1];
                                        m_pParticleData->m_vVolume[index2]=m_pParticleData->m_vVolume[index2+1];
                                        m_pParticleData->m_vVelocityU[index2]=m_pParticleData->m_vVelocityU[index2+1];
                                        m_pParticleData->m_vPressure[index2]=m_pParticleData->m_vPressure[index2+1];
                                        m_pParticleData->m_vSoundSpeed[index2]=m_pParticleData->m_vSoundSpeed[index2+1];
                                        m_pParticleData->m_vMass[index2]=m_pParticleData->m_vMass[index2+1];
				}
			}
		}
	}
}



void HyperbolicLPSolver1D::updateFreeBoundaryLocation() {
	
	// alias
	double* positionX = m_pParticleData->m_vPositionX;
	size_t totalNum = m_pParticleData->m_iTotalNum;
	
	// update location for free boundary by reconstruction using density
	double dis = (positionX[2] - positionX[1]) * 1.0;
	positionX[0] = positionX[1] - dis;
	dis = (positionX[totalNum-2] - positionX[totalNum-3]) * 1.0;
	positionX[totalNum-1] = positionX[totalNum-2] + dis;	

}


void HyperbolicLPSolver1D::updateFreeBoundaryPressureAndVelocity() {
	
	// alias
	double* pressure = m_pParticleData->m_vPressure;
	double* velocityU = m_pParticleData->m_vVelocityU;
	double* volume = m_pParticleData->m_vVolume;
	double* soundSpeed = m_pParticleData->m_vSoundSpeed;
	
	size_t totalNum = m_pParticleData->m_iTotalNum;
	
	// Assign pressure and velocity for computation of spatial derivatives
	pressure[0] = pressure[totalNum-1] = 1e-9; // 1e-9 is the pressure of vacuum
	velocityU[0] = velocityU[1];
	velocityU[totalNum-1] = velocityU[totalNum-2];	
	
	// Assign the following not for computation; just for consistent values for plotting
	volume[0] = volume[1];
	volume[totalNum-1] = volume[totalNum-2];
	soundSpeed[0] = soundSpeed[1];
	soundSpeed[totalNum-1] = soundSpeed[totalNum-2];	

}


void HyperbolicLPSolver1D::updateSolidBoundaryPressureAndVelocity() {

	// alias
	double* pressure = m_pParticleData->m_vPressure;
	double* velocityU = m_pParticleData->m_vVelocityU;
	double* volume = m_pParticleData->m_vVolume;
	double* soundSpeed = m_pParticleData->m_vSoundSpeed;
	
	size_t totalNum = m_pParticleData->m_iTotalNum;
	
	// Assign pressure and velocity for computation of spatial derivatives
	pressure[0] = pressure[1];
	pressure[totalNum-1] = pressure[totalNum-2];
	velocityU[0] = -velocityU[1];
	velocityU[totalNum-1] = -velocityU[totalNum-2];	
	
	// Assign the following not for computation; just for consistent values for plotting
	volume[0] = volume[1];
	volume[totalNum-1] = volume[totalNum-2];
	soundSpeed[0] = soundSpeed[1];
	soundSpeed[totalNum-1] = soundSpeed[totalNum-2];

}


void HyperbolicLPSolver1D::computeMinParticleSpacing() {
	
	// alias
	const double *positionX = m_pParticleData->m_vPositionX;	
	size_t totalNum = m_pParticleData->m_iTotalNum;

	// initial value
	m_fMinParticleSpacing = numeric_limits<double>::max();	
	
	for(size_t i=1; i<totalNum; i++) {
		double dis = fabs(positionX[i]-positionX[i-1]);
		m_fMinParticleSpacing = min(m_fMinParticleSpacing, dis);
	}			
	
	if(m_fMinParticleSpacing == numeric_limits<double>::max()) {
		cout<<"ERROR: Cannot find valid minimum inter-particle distance!"<<endl;
		exit(1);
	}
	
	if(m_iIfDebug) {
		debug.precision(16);
		//debug<<"-------HyperbolicLPSolver1D::computeMinParticleSpacing()-------"<<endl;
		debug<<"m_fMinParticleSpacing="<<m_fMinParticleSpacing<<endl;
		//debug<<"---------------------------------------------------------------"<<endl;
	}

}


void HyperbolicLPSolver1D::computeMaxSoundSpeed() {
	
	// alias
	const double* soundSpeed = m_pParticleData->m_vSoundSpeed;
	size_t totalNum = m_pParticleData->m_iTotalNum;

	//initial value
	m_fMaxSoundSpeed = -1;
	
	size_t startIndex = m_sBoundaryType=="periodic"? 0:1;
	size_t endIndex = m_sBoundaryType=="periodic"? totalNum:totalNum-1;
	for(size_t i=startIndex; i<endIndex; i++) 
		m_fMaxSoundSpeed = max(m_fMaxSoundSpeed, soundSpeed[i]);	
		
	if(m_fMaxSoundSpeed == -1) {
		cout<<"ERROR: Cannot find valid maximum sound speed!"<<endl;
		exit(1);	
	}
	
	if(m_iIfDebug) {
		debug.precision(16);
		//debug<<"-------HyperbolicLPSolver1D::computeMaxSoundSpeed()-------"<<endl;
		debug<<"m_fMaxSoundSpeed="<<m_fMaxSoundSpeed<<endl;
		//debug<<"----------------------------------------------------------"<<endl;
	}
}



void HyperbolicLPSolver1D::computeMaxFluidVelocity() {
	
	// alias
	const double* velocityU = m_pParticleData->m_vVelocityU;
	size_t totalNum = m_pParticleData->m_iTotalNum;

	// initial value
	m_fMaxFluidVelocity = -1;

	size_t startIndex = m_sBoundaryType=="periodic"? 0:1;
	size_t endIndex = m_sBoundaryType=="periodic"? totalNum:totalNum-1;
	for(size_t i=startIndex; i<endIndex; i++) 
		m_fMaxFluidVelocity = max(m_fMaxFluidVelocity, fabs(velocityU[i]));	
		
	if(m_fMaxFluidVelocity == -1) {
		cout<<"ERROR: Cannot find valid maximum fluid velocity!"<<endl;
		exit(1);	
	}
	
	if(m_iIfDebug) {
		debug.precision(16);
		//debug<<"-------HyperbolicLPSolver1D::computeMaxFluidVelocity()-------"<<endl;
		debug<<"m_fMaxFluidVelocity="<<m_fMaxFluidVelocity<<endl;
		//debug<<"-------------------------------------------------------------"<<endl;
	}

}


void HyperbolicLPSolver1D::updateLimiter() {
	
	// alias so that do not need to modify original 1D code
	double* xp        = m_pParticleData->m_vPositionX;
	double* p         = m_pParticleData->m_vPressure;
	double* dd1       = m_pParticleData->m_vDD1;
	double* dd2_left  = m_pParticleData->m_vDD2Left;
	double* dd2_right = m_pParticleData->m_vDD2Right;
	double* dd3_left  = m_pParticleData->m_vDD3Left;
	double* dd3_right = m_pParticleData->m_vDD3Right;
	double* cumP      = m_pParticleData->m_vCumP;
	double* xp_m      = m_pParticleData->m_vPositionXm;	
	int     np        = m_pParticleData->m_iTotalNum;

	// Calculate middle values for xp's 
	// and their corresponding cumulative pressure cumP
	double dis   = 0.5*(xp[1]-xp[0]);
	double xp05  = xp[0] + dis;
	double xp_05 = xp[0] - dis;
	//double xp05  = 0.5*(xp[0]+xp[1]);
	//double xp_05 = 0.5*(xp[0]+(xp[np-1]-dis_between_periodic_boundaries)); 
	xp_m[0] = xp05;
	cumP[0] = p[0]*(xp05-xp_05);

	for(int i=1; i<np-1; i++) {
		xp05  = 0.5*(xp[i]+xp[i+1]);
		xp_05 = 0.5*(xp[i]+xp[i-1]);
		xp_m[i] = xp05;
		cumP[i] = cumP[i-1] + p[i]*(xp05-xp_05);
	}
	dis = 0.5*(xp[np-1]-xp[np-2]);
	xp05 = xp[np-1] + dis;
	xp_05 = xp[np-1] - dis;
	//xp05  = 0.5*(xp[np-1]+(xp[0]+dis_between_periodic_boundaries)); 
	//xp_05 = 0.5*(xp[np-1]+xp[np-2]);
	xp_m[np-1] = xp05;
	cumP[np-1] = cumP[np-2] + p[np-1]*(xp05-xp_05);
	
	// compute dd1
	for(int i=1; i<np; i++) 
		dd1[i] = (cumP[i] - cumP[i-1]) / (xp_m[i] - xp_m[i-1]);  
	dd1[0] = dd1[1];
	
	// compute dd2_left
	for(int i=2; i<np; i++) 
		dd2_left[i] = (dd1[i]-dd1[i-1]) / (xp_m[i]-xp_m[i-2]);		
	dd2_left[0] = dd2_left[2];
	dd2_left[1] = dd2_left[2];
	
	// compute dd2_right
	for(int i=1; i<np-1; i++) 
		dd2_right[i] = (dd1[i+1]-dd1[i]) / (xp_m[i+1]-xp_m[i-1]);
	dd2_right[0] = dd2_right[1];
	dd2_right[np-1] = dd2_right[np-2];
	
	// compute dd3_left
	for(int i=3; i<np; i++) 
		dd3_left[i] = (dd2_left[i]-dd2_left[i-1]) / (xp_m[i]-xp_m[i-3]);
	dd3_left[0] = dd3_left[3];
	dd3_left[1] = dd3_left[3];
	dd3_left[2] = dd3_left[3];
	
	// compute dd3_right
	for(int i=1; i<np-2; i++) 
		dd3_right[i] = (dd2_right[i+1]-dd2_right[i]) / (xp_m[i+2]-xp_m[i-1]);
	dd3_right[0] = dd3_right[1];
	dd3_right[np-1] = dd3_right[np-3];
	dd3_right[np-2] = dd3_right[np-3];

}

////////////////////////////////////////////////////////////////////////////////////////
// End of HyperbolicLPSolver1D
////////////////////////////////////////////////////////////////////////////////////////
