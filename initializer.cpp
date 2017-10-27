#include "initializer.h"
#include "geometry.h"
#include "state.h"
#include "boundary.h"
#include "eos.h"
#include "hexagonal_packing.h"
#include "neighbour_searcher.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <omp.h>
//#include <unordered_map>

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of Initializer
////////////////////////////////////////////////////////////////////////////////////////

Initializer::Initializer(const string& inputfileName, bool ifDebug, const string& debugfileName)
: m_iIfDebug(ifDebug), m_sDebugfileName(debugfileName) {	
//	cout<<"readInputfile"<<endl;	
	readInputfile(inputfileName); // read inputfile
//	cout<<"setparams"<<endl; 
	setParams();
//        cout<<"seteos"<<endl; 	
	setEOS();	
//        cout<<"setobjs"<<endl; 	
	setObjs();
//        cout<<"setboundingbox"<<endl; 
	setBoundingBox(); // set up fluid bounding box for geometry init and assigning object tag	
//        cout<<"initgeometryandstate"<<endl; 	
	initGeometryAndState();	

	setObjectTag();

	setLocalParSpacing();
//	setLocalParSpacingTemp();

	m_iIfRestart = false;
}


Initializer::Initializer(const string& inputfileName, const string& datafileName, 
bool ifDebug, const string& debugfileName)
: m_iIfDebug(ifDebug), m_sDebugfileName(debugfileName) {
	
	readInputfile(inputfileName);
	
	setParams();
	
	setEOS();
		
	readDatafile(datafileName); // read datafile
	
	//modifyLocalParSpacing();
	//modifyInitParticleSpacing();
	//modifyInitNeighbourSearchRadius();
	//modifyInitContactLength();
	
	//modifyNumParticleWithinSearchRadius();

	setBoundingBox(m_vObjectTag,m_iFluidNum);

	if(m_iDimension != 1) setBoundingBoxStartIndex();
	
	m_iIfRestart = true;
}


void Initializer::setLocalParSpacing() {

	for(size_t i=0; i<m_iFluidNum; i++) 
		m_vLocalParSpacing[i] = m_fInitParticleSpacing;
		
}

void Initializer::setLocalParSpacingTemp() {
	
	double factor;
	double denDiff = 2.0;
	if(m_iDimension==2) factor = sqrt(denDiff);
	else if (m_iDimension==3) factor=std::cbrt(denDiff);
	for(size_t i=0; i<m_iFluidNum; i++) {
		if(m_vPositionX[i]<=0) m_vLocalParSpacing[i] = m_fInitParticleSpacing;
		else m_vLocalParSpacing[i] = m_fInitParticleSpacing * factor;
	}

}


Initializer::~Initializer() {
	for(auto obj:m_vFluidObj) delete obj;
//	for(auto obj:m_vBoundaryObj) delete obj;	
	for(auto state:m_vFluidObjState) delete state;
//	for(auto box:m_vFluidBoundingBox) delete box;
//	for(auto box:m_vBoundaryBoundingBox) delete box;
}


void Initializer::readDatafile(const string& datafileName) {
	
	ofstream save(m_sFilenameSaveInit, ofstream::app); // save init info

	ifstream ifs(datafileName);
	
	double v1, v2, v3;
	string s;
	size_t num;
	
	getline(ifs,s); // Skip 1 line
	ifs >> s >> s >> s >> s >> v1;
	save<<"(RESTART CHANGE): m_fStartTime has been changed from "<<m_fStartTime;
	m_fStartTime = v1;
	save<<" to "<<m_fStartTime<<endl;	
	
	save<<"(RESTART CHANGE): m_iWriteStep has been changed from "<<m_iWriteStep;
	m_iWriteStep = m_fStartTime/m_fWriteTimeInterval;
	save<<" to "<<m_iWriteStep<<endl;

	for(int line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	ifs >> s >> num >> s; // Read the number of particles (fluid + boundary)
	save<<"(RESTART) number of particles (boundary + fluid) "<<num<<endl;
	m_iCapacity = (size_t)(m_fTimesCapacity*num + m_fAdditionalCapacity); // set capacity of data arays
	initParticleDataMemory(); // Allocate memory

	for(size_t i=0; i<num; i++) { // Read location
		ifs >> v1 >> v2 >> v3;
		m_vPositionX[i] = v1;
		m_vPositionY[i] = v2;
		m_vPositionZ[i] = v3;
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<num; i++) { // Read velocity
		ifs >> v1 >> v2 >> v3;
		m_vVelocityU[i] = v1;
		m_vVelocityV[i] = v2;
		if(m_iDimension == 3) m_vVelocityW[i] = v3;
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines	
	for(size_t i=0; i<num; i++) { // Read pressure
		ifs >> v1;
		m_vPressure[i] = v1;
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<num; i++) { // Read volume
		ifs >> v1;
		m_vVolume[i] = v1;
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<num; i++) { // Read sound speed
		ifs >> v1;
		m_vSoundSpeed[i] = v1;
	}
	
	//map<int,vector<double>> counter; // counter for each tag
	//counter {tag:[num,xmin,xmax,ymin,ymax,zmin,zmax]}
//	m_iBoundaryNum = 0;
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<num; i++) { // Read object tag
		ifs >> v1;
		m_vObjectTag[i] = (int)v1;
//		if(m_vObjectTag[i]==0) m_iBoundaryNum++;
	//	int v = int(v1);
	//	m_vObjectTag[i] = v;
	//	if(counter.find(v) != counter.end()) {
	//		counter[v][0] += 1;
	//		counter[v][1] = min(counter[v][1],m_vPositionX[i]);
	//		counter[v][2] = max(counter[v][2],m_vPositionX[i]);
	//		counter[v][3] = min(counter[v][3],m_vPositionY[i]);
	//		counter[v][4] = max(counter[v][4],m_vPositionY[i]);
	//		counter[v][5] = min(counter[v][5],m_vPositionZ[i]);
	//		counter[v][6] = max(counter[v][6],m_vPositionZ[i]);
	//	}
	//	else {
	//		counter.insert({v, vector<double>(7,0)});
	//		counter[v][0] = 1;
	//		counter[v][1] = m_vPositionX[i];
	//		counter[v][2] = m_vPositionX[i];
	//		counter[v][3] = m_vPositionY[i];
	//		counter[v][4] = m_vPositionY[i];
	//		counter[v][5] = m_vPositionZ[i];
	//		counter[v][6] = m_vPositionZ[i];
	//	}	
	}
	//if(counter.find(0) != counter.end()) // no boundary particle
	//	m_iBoundaryNum = counter[0][0];
	//else
	//	m_iBoundaryNum = 0;
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<num; i++) { // Read sound speed
		ifs >> v1;
		m_vLocalParSpacing[i] = v1;
	}
	
//	m_iFluidNum = num - m_iBoundaryNum;
	m_iFluidNum = num;
	m_iFluidStartIndex = 0;
//	m_iBoundaryStartIndex = m_iFluidStartIndex + m_iFluidNum;
//	m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum + m_iBoundaryNum;
	
	//setBoundingBox(counter); // set bounding boxes for fluid & boundary objects

	save<<"m_iFluidNum "<<m_iFluidNum<<endl;
//	save<<"m_iBoundaryNum "<<m_iBoundaryNum<<endl;
	save<<"m_iCapacity "<<m_iCapacity<<endl;
	save<<"m_iFluidStartIndex "<<m_iFluidStartIndex<<endl;
//	save<<"m_iBoundaryStartIndex "<<m_iBoundaryStartIndex<<endl;
//	save<<"m_iGhostStartIndex "<<m_iGhostStartIndex<<endl;


}




void Initializer::readInputfile(const string& inputfileName) {
	
	ifstream ifs(inputfileName); // read the input file
	m_sFilenameSaveInit = "save_init_param"; // filename for saving info for restart

	vector<string> lines;
	string s;
	ofstream save(m_sFilenameSaveInit); // save info for restart
	
	while(getline(ifs, s)) {
		//ofs<<s<<endl; 
		lines.push_back(s);
	}

	istringstream iss;
	size_t i = 0; // input file line number
	
	iss.str(lines[i++]);
	iss>>m_iNumThreads;
	save<<"m_iNumThreads "<<m_iNumThreads<<endl;
/*	#ifdef _OPENMP
	if(m_iNumThreads<=1) {
		cerr<<"(ERROR): Multithread mode: m_iNumThreads="<<m_iNumThreads<<endl;
		cerr<<"Please modify the inputfile."<<endl;
		assert(false);
	}
	#else
	if(m_iNumThreads>1) {
		cerr<<"(ERROR): Single thread mode: m_iNumThreads="<<m_iNumThreads<<endl;
		cerr<<"Please modify the inputfile."<<endl;
		assert(false);
	}
	#endif
*/
	iss.str(lines[i++]);
	iss>>m_fStartTime;
	save<<"m_fStartTime "<<m_fStartTime<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fEndTime;
	save<<"m_fEndTime "<<m_fEndTime<<endl;

	iss.str(lines[i++]);
	iss>>m_fWriteTimeInterval;
	save<<"m_fWriteTimeInterval "<<m_fWriteTimeInterval<<endl;

	iss.str(lines[i++]);
	iss>>m_fCFLCoeff;
	save<<"m_fCFLCoeff "<<m_fCFLCoeff<<endl;

	iss.str(lines[i++]);
	iss>>m_iDimension;
	save<<"m_iDimension "<<m_iDimension<<endl;

	iss.str(lines[i++]);
	iss>>m_iFluidObjNum;
	save<<"m_iFluidObjNum "<<m_iFluidObjNum<<endl;
 
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vFluidObjNames.push_back(tmpS);
		save<<"m_vFluidObjectNames"<<j<<" "<<m_vFluidObjNames[j]<<endl;
	}
	
	for(int j=0; j<m_iFluidObjNum; j++) {
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vFluidObjStateNames.push_back(tmpS);
		save<<"m_vFluidObjectStateNames"<<j<<" "<<m_vFluidObjStateNames[j]<<endl;
	}

	
	iss.str(lines[i++]);
	iss>>m_iBoundaryObjNum;
	save<<"m_iBoundaryObjNum "<<m_iBoundaryObjNum<<endl;
/*	
	if(m_iDimension == 1) { // In the 1D case only boundary type is relevant (only one type is allowed)
		string tmpS;
		iss.str(lines[i++]);
		iss>>tmpS;
		m_vBoundaryObjTypes.push_back(tmpS);
		save<<"m_vBoundaryObjectTypes[0] "<<m_vBoundaryObjTypes[0]<<endl;
	}
	else { // In 2D & 3D both boundary geometry name and boundary type should be specified
*/
		for(int j=0; j<m_iBoundaryObjNum; j++) {
			string tmpS;
		
			iss.str(lines[i++]);
			iss>>tmpS;
			m_vBoundaryObjNames.push_back(tmpS);
			save<<"m_vBoundaryObjectNames"<<j<<" "<<m_vBoundaryObjNames[j]<<endl;
			
			iss.str(lines[i++]);
			iss>>tmpS;
			m_vBoundaryObjTypes.push_back(tmpS);
			save<<"m_vBoundaryObjectTypes"<<j<<" "<<m_vBoundaryObjTypes[j]<<endl;
		}
//	}

	iss.str(lines[i++]);
	iss>>m_iRandomDirSplitOrder;
	save<<"m_iRandomDirSplitOrder "<<m_iRandomDirSplitOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iLPFOrder;
	save<<"m_iLPFOrder "<<m_iLPFOrder<<endl;

	iss.str(lines[i++]);
	iss>>m_iEOSChoice;
	save<<"m_iEOSChoice "<<m_iEOSChoice<<endl;

	iss.str(lines[i++]);
	iss>>m_fGamma;
	save<<"m_fGamma "<<m_fGamma<<endl;

	iss.str(lines[i++]);
	iss>>m_fPinf;
	save<<"m_fPinf "<<m_fPinf<<endl;

	iss.str(lines[i++]);
	iss>>m_fEinf;
	save<<"m_fEinf "<<m_fEinf<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fInitParticleSpacing;
	save<<"m_fInitParticleSpacing "<<m_fInitParticleSpacing<<endl;	

	iss.str(lines[i++]);
	iss>>m_fGravity;
	save<<"m_fGravity "<<m_fGravity<<endl;

	iss.str(lines[i++]);
	iss>>m_iMovingBoxForGhostParticle;
	save<<"m_iMovingBoxForGhostParticle "<<m_iMovingBoxForGhostParticle<<endl;	
	
	iss.str(lines[i++]);
	iss>>m_iUseLimiter;
	save<<"m_iUseLimiter "<<m_iUseLimiter<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fThresholdP;
	save<<"m_fThresholdP "<<m_fThresholdP<<endl;
	
	iss.str(lines[i++]);
	iss>>m_iUseCriticalPressure;
	save<<"m_iUseCriticalPressure "<<m_iUseCriticalPressure<<endl;
	
	iss.str(lines[i++]);
	iss>>m_fCriticalPressure;
	save<<"m_fCriticalPressure "<<m_fCriticalPressure<<endl;	
	

	iss.str(lines[i++]);
	iss>>m_iVariableNeiSearchRadius;
	save<<"m_iVariableNeiSearchRadius "<<m_iVariableNeiSearchRadius<<endl;

        iss.str(lines[i++]);
        iss>>m_iIfLaxWendroff;
        save<<"m_iIfLaxWendroff "<<m_iIfLaxWendroff<<endl;

        iss.str(lines[i++]);
        iss>>m_iIfNoSplit;
        save<<"m_iIfNoSplit "<<m_iIfNoSplit<<endl;

        iss.str(lines[i++]);
        iss>>m_iIfSPH;
        save<<"m_iIfSPH "<<m_iIfSPH<<endl;

}


void Initializer::setParams() {
			
	ofstream save(m_sFilenameSaveInit, ofstream::app); 

	m_iWriteStep = 0; // The default value
	if(m_iDimension==2)
		m_fTimesNeiSearchRadius = 4;//4
        if(m_iDimension==3)
                m_fTimesNeiSearchRadius = 2.5;//4

	m_fNeiSearchRadius = m_fTimesNeiSearchRadius*m_fInitParticleSpacing;

	m_fTimesContactLength = 1.1;
	m_fContactLength = 1.1*m_fInitParticleSpacing;
	
	m_fTimesBoundingBox = m_fTimesNeiSearchRadius+1;
	
	if(m_iDimension==2)
		m_fTimesCapacity = 3;
	if(m_iDimension==3)
		m_fTimesCapacity = 2;//1.5

	m_fAdditionalCapacity = 500000;
	if(m_iDimension==3)
		m_fAdditionalCapacity = 4000000;	
	// compute the number of particles within m_fNeiSearchRadius based on the initial packing of particles
	computeNumParticleWithinSearchRadius();
	
	if(m_iDimension==3) {
		m_iMaxNeighbourNum = 3*m_iNumParticleWithinSearchRadius;
		cout<<"Maxinum number of neighbour for each particles: "<<m_iMaxNeighbourNum<<endl;
//		m_iMaxNeighbourNumInOneDir = m_iMaxNeighbourNum/2;
		m_iMaxNeighbourNumInOneDir = m_iMaxNeighbourNum;				
		m_iNumRow2ndOrder = 54; 
		m_iNumRow1stOrder = 8; 
		m_iNumCol2ndOrder = 9; 
		m_iNumCol1stOrder = 3;

//		if(m_iMaxNeighbourNum>(2*m_iNumRow2ndOrder))
//			m_iMaxNeighbourNum=2*m_iNumRow2ndOrder;
//		if(m_iMaxNeighbourNumInOneDir>(8*m_iNumRow1stOrder))
//			m_iMaxNeighbourNumInOneDir=8*m_iNumRow1stOrder; 
	}
	else if(m_iDimension==2) {
		m_iMaxNeighbourNum = 3*m_iNumParticleWithinSearchRadius; 
		m_iMaxNeighbourNumInOneDir = m_iMaxNeighbourNum/2; 
				
		m_iNumRow2ndOrder = 36;//36 
		m_iNumRow1stOrder = 3; 
		m_iNumCol2ndOrder = 5; 
		m_iNumCol1stOrder = 2;	
	}
	else { // 1D
		m_iMaxNeighbourNum = m_iLPFOrder==1? 2:4; 
		m_iMaxNeighbourNumInOneDir = m_iLPFOrder==1? 1:2;
				
		m_iNumRow2ndOrder = 2; 
		m_iNumRow1stOrder = 1; 
		m_iNumCol2ndOrder = 2; 
		m_iNumCol1stOrder = 1;	
	
	}

	// (eos chocie=1:poly 2:spoly)	
	if(m_iEOSChoice==1) 
		m_fInvalidPressure = 0;
	else if(m_iEOSChoice==2) 	
		m_fInvalidPressure = -m_fPinf;
	m_fInvalidDensity = 0;
	
	m_iTreeDepth = 7;
        if(m_iDimension==3) m_iTreeDepth = 5;
	
	// save info for restart
	save<<"m_iWriteStep "<<m_iWriteStep<<endl;
	save<<"m_fTimesNeiSearchRadius "<<m_fTimesNeiSearchRadius<<endl;
	save<<"m_fTimesCapacity "<<m_fTimesCapacity<<endl;
	save<<"m_fAdditionalCapacity "<<m_fAdditionalCapacity<<endl;
	save<<"m_fNeiSearchRadius "<<m_fNeiSearchRadius<<endl;
	//save<<"m_iNumParticleWithinSearchRadius "<<m_iNumParticleWithinSearchRadius<<endl;
	save<<"m_iMaxNeighbourNum "<<m_iMaxNeighbourNum<<endl;
	save<<"m_iMaxNeighbourNumInOneDir "<<m_iMaxNeighbourNumInOneDir<<endl;
	save<<"m_iNumRow2ndOrder "<<m_iNumRow2ndOrder<<endl;
	save<<"m_iNumRow1stOrder "<<m_iNumRow1stOrder<<endl;
	save<<"m_iNumCol2ndOrder "<<m_iNumCol2ndOrder<<endl;
	save<<"m_iNumCol1stOrder "<<m_iNumCol1stOrder<<endl;	
	save<<"m_fInvalidPressure "<<m_fInvalidPressure<<endl;
	save<<"m_fInvalidDensity "<<m_fInvalidDensity<<endl;
	save<<"m_iTreeDepth "<<m_iTreeDepth<<endl;		
	save<<"m_fContactLength "<<m_fContactLength<<endl;
	
}


void Initializer::setEOS() {
	 
	if(m_iEOSChoice == 1) // Polytropic gas EOS
		m_pEOS = new PolytropicGasEOS(m_fGamma);
	else if(m_iEOSChoice==2) // Stiffened Polytropic gas EOS
		m_pEOS = new StiffPolytropicGasEOS(m_fGamma,m_fPinf,m_fEinf);
	else {
		cout<<"The choice of EOS does not exist!!! Please correct the input file."<<endl;
		assert(false);
	}

}


void Initializer::setObjs() {
	
	// create fluid geometry objects 
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObj.push_back(GeometryFactory::instance().createGeometry(m_vFluidObjNames[j]));
	
	// create fluid state objects
	for(int j=0; j<m_iFluidObjNum; j++) 
		m_vFluidObjState.push_back(StateFactory::instance().createState(m_vFluidObjStateNames[j]));

	if(m_iDimension != 1) { // // create boundary objects (for 2D & 3D)	
		for(int j=0; j<m_iBoundaryObjNum; j++)
			if(m_vBoundaryObjTypes[j]!="free")
				m_vBoundaryObj.push_back(BoundaryFactory::instance().createBoundary(m_vBoundaryObjNames[j]));
			//m_vBoundaryObj.push_back(GeometryFactory::instance().createGeometry(m_vBoundaryObjNames[j]));	
	}

}


void Initializer::setBoundingBox() {
	
	if(m_iDimension == 1) return; // bounding box for geometry initialization is only for 2D & 3D
//	std::cout<<1<<std::endl;
	ofstream save(m_sFilenameSaveInit, ofstream::app);
//        std::cout<<2<<std::endl;

	double sp  = m_fTimesBoundingBox*m_fInitParticleSpacing;
//        std::cout<<2<<std::endl;

	for(size_t i=0; i<m_vFluidObj.size(); i++) {
		double xmin, xmax, ymin, ymax, zmin, zmax;
//		std::cout<<i<<std::endl;
		m_vFluidObj[i]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		BoundingBox* box;
		if(m_iDimension==3)
//			box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, zmin-sp, zmax+sp);
			box = new BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		else
//			box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, 0, 0);	
			box = new BoundingBox(xmin, xmax, ymin, ymax, 0, 0); 
		m_vFluidBoundingBox.push_back(box);
//                std::cout<<i<<std::endl;
		
		save<<"m_vFluidBoundingBox("<<i<<"):"<<endl;
		save<<"m_fXmin="<<m_vFluidBoundingBox[i]->getXmin()<<endl;
		save<<"m_fXmax="<<m_vFluidBoundingBox[i]->getXmax()<<endl;
		save<<"m_fYmin="<<m_vFluidBoundingBox[i]->getYmin()<<endl;
		save<<"m_fYmax="<<m_vFluidBoundingBox[i]->getYmax()<<endl;
		save<<"m_fZmin="<<m_vFluidBoundingBox[i]->getZmin()<<endl;
		save<<"m_fZmax="<<m_vFluidBoundingBox[i]->getZmax()<<endl;
		save<<"-----------------"<<endl;
	}
/*
	for(size_t i=0; i<m_vBoundaryObj.size(); i++) {
		double xmin, xmax, ymin, ymax, zmin, zmax;
		m_vBoundaryObj[i]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
		BoundingBox* box;
		if(m_iDimension==3)
			box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, zmin-sp, zmax+sp);
		else
			box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, 0, 0);	
		m_vBoundaryBoundingBox.push_back(box);
		
		save<<"m_vBoundaryBoundingBox("<<i<<"):"<<endl;
		save<<"m_fXmin="<<m_vBoundaryBoundingBox[i]->getXmin()<<endl;
		save<<"m_fXmax="<<m_vBoundaryBoundingBox[i]->getXmax()<<endl;
		save<<"m_fYmin="<<m_vBoundaryBoundingBox[i]->getYmin()<<endl;
		save<<"m_fYmax="<<m_vBoundaryBoundingBox[i]->getYmax()<<endl;
		save<<"m_fZmin="<<m_vBoundaryBoundingBox[i]->getZmin()<<endl;
		save<<"m_fZmax="<<m_vBoundaryBoundingBox[i]->getZmax()<<endl;
		save<<"-----------------"<<endl;
		
	}
*/

}


void Initializer::setBoundingBox(const int* tag, size_t num) {
/*	
	double sp = m_fTimesBoundingBox*m_fInitParticleSpacing;

	double xmin=m_vPositionX[0], xmax=m_vPositionX[0], 
	ymin=m_vPositionY[0], ymax=m_vPositionY[0], 
	zmin=m_vPositionZ[0], zmax=m_vPositionZ[0];
	int tprev = tag[0];
	size_t n = 1;
	for(size_t i=1; i<num; i++) {
		if(tag[i]==tprev) {
			xmin=min(xmin,m_vPositionX[i]);
			xmax=max(xmax,m_vPositionX[i]);
			ymin=min(ymin,m_vPositionY[i]);
			ymax=max(ymax,m_vPositionY[i]);
			zmin=min(zmin,m_vPositionZ[i]);
			zmax=max(zmax,m_vPositionZ[i]);
			n++;
		}
		else {
			BoundingBox* box;
			if(m_iDimension==3) 
				box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, zmin-sp, zmax+sp);
			else
				box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, 0, 0);
			box->setNumber(n);
			box->setObjectTag(tprev);
			if(tprev==0) m_vBoundaryBoundingBox.push_back(box);
			else m_vFluidBoundingBox.push_back(box);
			xmin=m_vPositionX[i];
			xmax=m_vPositionX[i];
			ymin=m_vPositionY[i];
			ymax=m_vPositionY[i];
			zmin=m_vPositionZ[i];
			zmax=m_vPositionZ[i];
			n=1;
			tprev=tag[i];
		}

	}
	BoundingBox* box;
	if(m_iDimension==3) 
		box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, zmin-sp, zmax+sp);
	else
		box = new BoundingBox(xmin-sp, xmax+sp, ymin-sp, ymax+sp, 0, 0);
	box->setNumber(n);
	box->setObjectTag(tprev);
	if(tprev==0) m_vBoundaryBoundingBox.push_back(box);
	else m_vFluidBoundingBox.push_back(box);	

	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) {
		cout<<"m_vFluidBoundingBox["<<p<<"]:"<<endl;
		cout<<"m_iObjectTag="<<m_vFluidBoundingBox[p]->getObjectTag()<<endl;
		cout<<"m_iNumber="<<m_vFluidBoundingBox[p]->getNumber()<<endl;
		cout<<"m_fXmin="<<m_vFluidBoundingBox[p]->getXmin()<<endl;
		cout<<"m_fXmax="<<m_vFluidBoundingBox[p]->getXmax()<<endl;
		cout<<"m_fYmin="<<m_vFluidBoundingBox[p]->getYmin()<<endl;
		cout<<"m_fYmax="<<m_vFluidBoundingBox[p]->getYmax()<<endl;
		cout<<"m_fZmin="<<m_vFluidBoundingBox[p]->getZmin()<<endl;
		cout<<"m_fZmax="<<m_vFluidBoundingBox[p]->getZmax()<<endl;
	}

*/	
}


void Initializer::initGeometryAndState() {
	
	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart
	
	if(m_iDimension == 1) {
		
		bool saveData = false;
		// Compute the number of particles	
//		m_iBoundaryNum = 0; 
		m_iFluidNum = initGeometryAndState1D(saveData);
		m_iCapacity = m_iFluidNum;
		
		// use m_iCapacity to create memory
		initParticleDataMemory();
		
		saveData = true;
		// Actually initialize particle location and state this time	
		initGeometryAndState1D(saveData);	
		
		// the following are NOT relevant in 1D; In 1D only m_iFluidNum is used
		m_iFluidStartIndex = 0;
//		m_iBoundaryStartIndex = 0;
//		m_iGhostStartIndex = 0;
//		m_iBoundaryNum = 0;	

	}
	else {
		bool saveData = false;
		 
		m_iFluidNum = initGeometryAndStateOnHexPacking(saveData);// Compute the number of particles
//		m_iFluidNum = initGeometryAndStateOnHexPackingTemp(saveData);//for 2d shock tube

		m_iCapacity = (size_t)(m_fTimesCapacity*m_iFluidNum + m_fAdditionalCapacity);

		
		initParticleDataMemory();// use m_iCapacity to create memory
		
		saveData = true;
		
		initGeometryAndStateOnHexPacking(saveData);// Actually initialize particle location and state this time
//		initGeometryAndStateOnHexPackingTemp(saveData);//for 2d shock tube

		// assign the start index of fluid, boundary, and ghost particles
		// this variables are only relevant in 2D & 3D; In 1D only m_iFluidNum is used
		m_iFluidStartIndex = 0;
//		m_iBoundaryStartIndex = m_iFluidStartIndex + m_iFluidNum;
//		m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum;
//		m_iGhostStartIndex = m_iFluidStartIndex + m_iFluidNum + m_iBoundaryNum;
		cout<<"Number of fluid particles: "<<m_iFluidNum<<endl;
		cout<<"Allocated storage for all particles: "<<m_iCapacity<<endl;			
			
	}	

	save<<"m_iFluidNum "<<m_iFluidNum<<endl;
//	save<<"m_iBoundaryNum "<<m_iBoundaryNum<<endl;
	save<<"m_iCapacity "<<m_iCapacity<<endl;
	save<<"m_iFluidStartIndex "<<m_iFluidStartIndex<<endl;
//	save<<"m_iBoundaryStartIndex "<<m_iBoundaryStartIndex<<endl;
//	save<<"m_iGhostStartIndex "<<m_iGhostStartIndex<<endl;
	
}


void Initializer::initParticleDataMemory() {
	
	// location
	m_vPositionX = new double[m_iCapacity];
	m_vPositionY = new double[m_iCapacity];
	m_vPositionZ = new double[m_iCapacity];
	fill_n(m_vPositionX,m_iCapacity,0);
	fill_n(m_vPositionY,m_iCapacity,0);
	fill_n(m_vPositionZ,m_iCapacity,0);	
	
	// velocity
	m_vVelocityU = new double[m_iCapacity];
	fill_n(m_vVelocityU,m_iCapacity,0);
	if(m_iDimension==2 || m_iDimension==3) {
		m_vVelocityV = new double[m_iCapacity];
		fill_n(m_vVelocityV,m_iCapacity,0);
	}
	else m_vVelocityV = nullptr;
	if(m_iDimension==3) {
		m_vVelocityW = new double[m_iCapacity];
		fill_n(m_vVelocityW,m_iCapacity,0);
	}
	else m_vVelocityW = nullptr; 
	
	// states
	m_vVolume = new double[m_iCapacity];
	fill_n(m_vVolume,m_iCapacity,0);
	
	m_vPressure = new double[m_iCapacity];
	fill_n(m_vPressure,m_iCapacity,0);	

        m_vMass = new double[m_iCapacity];
        fill_n(m_vMass,m_iCapacity,0);
	
	m_vSoundSpeed = new double[m_iCapacity];
	fill_n(m_vSoundSpeed,m_iCapacity,0);
	
	// object tags
	m_vObjectTag = new int[m_iCapacity];
	fill_n(m_vObjectTag,m_iCapacity,0);

	// local inter-particle spacing
	m_vLocalParSpacing = new double[m_iCapacity];
	fill_n(m_vLocalParSpacing,m_iCapacity,0);
	
}


//size_t Initializer::initGeometryAndState1D(bool saveData) {
//	
//	// alias
//	const vector<Geometry*>& objs = m_vFluidObj;
//	const vector<State*>& states = m_vFluidObjState;	
//
//	double xmin, xmax, ymin, ymax, zmin, zmax;
//	objs[0]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
//	size_t numParticle = (size_t)((xmax - xmin)/m_fInitParticleSpacing) + 1;
//		
//	if(saveData) { 	
//		for(size_t i=0; i<numParticle; i++) { 
//			double x = xmin + i*m_fInitParticleSpacing;
//			m_vPositionX[i] = x; 
//			m_vVolume[i]    = 1./states[0]->density(x,0,0);
//			m_vPressure[i]  = states[0]->pressure(x,0,0);
//			double tmpY, tmpZ;
//			states[0]->velocity(x,0,0,m_vVelocityU[i],tmpY,tmpZ);
//			m_vSoundSpeed[i] = m_pEOS->getSoundSpeed(m_vPressure[i],1./m_vVolume[i]);
//		}
//	}
//
//	//cout<<"xmin="<<xmin<<"	xmax="<<xmax<<"	numParticle="<<numParticle<<endl;
//	return numParticle;
//}

size_t Initializer::initGeometryAndState1D(bool saveData) {
	
	// alias
	const vector<Geometry*>& objs = m_vFluidObj;
	const vector<State*>& states = m_vFluidObjState;	
	
	size_t numParticle = 0;
	size_t i = 0;
	for(size_t p=0; p<objs.size(); p++) {
		
		double xmin, xmax, ymin, ymax, zmin, zmax;
		objs[p]->getBoundingBox(xmin, xmax, ymin, ymax, zmin, zmax);
//		cout<<xmin<<" "<<xmax<<endl;
		size_t numParticleI = (size_t)((xmax - xmin)/m_fInitParticleSpacing) + 1;
		size_t end = i+numParticleI;

		if(saveData) { 	
			for(; i<end; i++) { 
				double x = xmin + (i-numParticle) *m_fInitParticleSpacing;
				m_vPositionX[i] = x; 
				m_vVolume[i]    = 1./states[p]->density(x,0,0);
				m_vPressure[i]  = states[p]->pressure(x,0,0);
				double tmpY, tmpZ;
				states[p]->velocity(x,0,0,m_vVelocityU[i],tmpY,tmpZ);
				m_vSoundSpeed[i] = m_pEOS->getSoundSpeed(m_vPressure[i],1./m_vVolume[i]);
				m_vMass[i]=m_fInitParticleSpacing/m_vVolume[i];
			}
		}
		numParticle += numParticleI;
	
	}
	//cout<<"xmin="<<xmin<<"	xmax="<<xmax<<"	numParticle="<<numParticle<<endl;
	return numParticle;
}


size_t Initializer::initGeometryAndStateOnHexPackingTemp(bool saveData) {

	double h_r = 0.5*m_fInitParticleSpacing;
	
	vector<BoundingBox*>& boxes = m_vFluidBoundingBox;
	const vector<Geometry*>& objs = m_vFluidObj;
	const vector<State*>& states = m_vFluidObjState;

	size_t numParticle = 0;

	if(m_iDimension==2) {

		for(size_t p=0; p<objs.size(); p++) {
		
			size_t numParticleI = 0;
			
			//left
			double xmin = boxes[p]->getXmin();	
			//double xmax = boxes[p]->getXmax();
			double xmax = 0;
			double ymin = boxes[p]->getYmin();	
			double ymax = boxes[p]->getYmax();
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	

			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; // level function of fluid object	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;	
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; 	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;	
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				}
			}

			
			
			//right
			//double xmin = boxes[p]->getXmin();
			xmin = h_r/2.;
			xmax = boxes[p]->getXmax();
			ymin = boxes[p]->getYmin();	
			ymax = boxes[p]->getYmax();
			
			h_r *= std::sqrt(2.0);

			HexagonalPacking2D hex2DD(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			//size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2DD.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	

			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2DD.computeX(0,k);
						double y = hex2DD.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; // level function of fluid object	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;	
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2DD.computeX(1,k);
						double y = hex2DD.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; 	
						if(saveData) {
							m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;	
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				}
			}
			boxes[p]->setNumber(numParticleI);
		}	
	}
	else assert(false);	
	
	return numParticle;
}



size_t Initializer::initGeometryAndStateOnHexPacking(bool saveData) {

	double h_r = 0.5*m_fInitParticleSpacing;
	
	vector<BoundingBox*>& boxes = m_vFluidBoundingBox;
	const vector<Geometry*>& objs = m_vFluidObj;
	const vector<State*>& states = m_vFluidObjState;

	size_t numParticle = 0;
	srand(1);
	if(m_iDimension==2) {

		for(size_t p=0; p<objs.size(); p++) {
			h_r = 0.5*m_fInitParticleSpacing;		
			size_t numParticleI = 0;

			double xmin = boxes[p]->getXmin();	
			double xmax = boxes[p]->getXmax();
			double ymin = boxes[p]->getYmin();	
			double ymax = boxes[p]->getYmax();
//                        cout<<ymin<<endl;

/*			double xt=(xmin+xmax)/2;
			double yt=(ymin+ymax)/2;
			int counter=0;
			while((!objs[p]->operator()(xt,yt,0))&&counter<10000){
				xt=xmin+((double)rand()/(double)RAND_MAX)*(xmax-xmin);
                                yt=ymin+((double)rand()/(double)RAND_MAX)*(ymax-ymin);
				counter++;
			}
			if (counter==10000)
			{
				cout<<"Error in initialization: cannot find particles in domain"<<endl;
				return(-1);
			}
			
			double dt=states[p]->density(xt,yt,0);
*/
//			h_r=h_r/std::sqrt(dt);	
						
			HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
			// get parameters of hexagonal packing
			size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
			hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);
//			n1_odd=n1_odd*2;
//			n1_even=n1_even*2;
			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t j=m0; j<=m1; j++) { 
				if((j+1)%2 != 0) { // odd-numbered rows 
					for(size_t k=n0_odd; k<=n1_odd; k++) { 
						double x = hex2D.computeX(0,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; // level function of fluid object	
//						if((!objs[p]->operator()(8*x,y,0))&&x>0)	continue;
						if(saveData) {
//							if(x>0)
//								m_vPositionX[numParticle] = 8*x;
//							else
								m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;
							
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
							m_vMass[numParticle] = sqrt(3)*2*h_r*h_r/m_vVolume[numParticle];
							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				} 
				else{ // even-numbered rows
					for(size_t k=n0_even; k<=n1_even; k++) {
						double x = hex2D.computeX(1,k);
						double y = hex2D.computeY(j);
						if(!objs[p]->operator()(x,y,0)) continue; 	
//                                                if((!objs[p]->operator()(8*x,y,0))&&x>0)        continue;
						if(saveData) {
//                                                        if(x>0)
//                                                                m_vPositionX[numParticle] = 8*x;
//                                                        else
								m_vPositionX[numParticle] = x;
							m_vPositionY[numParticle] = y;
							
							m_vVolume[numParticle] = 1./states[p]->density(x,y,0);
							m_vPressure[numParticle] = states[p]->pressure(x,y,0);
                                                        m_vMass[numParticle] = sqrt(3)*2*h_r*h_r/m_vVolume[numParticle];

							double tmpZ;
							states[p]->velocity(x,y,0,
								m_vVelocityU[numParticle],m_vVelocityV[numParticle],tmpZ);
							m_vSoundSpeed[numParticle] = 
								m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
							
						}	
						numParticle++;
						numParticleI++;
					}
				}
			}
			boxes[p]->setNumber(numParticleI);
		}	
	}
	else if(m_iDimension==3) {
	
		for(size_t p=0; p<objs.size(); p++) {
			size_t numParticleI = 0;	

			double xmin = boxes[p]->getXmin();	
			double xmax = boxes[p]->getXmax();
			double ymin = boxes[p]->getYmin();	
			double ymax = boxes[p]->getYmax();
			double zmin = boxes[p]->getZmin();	
			double zmax = boxes[p]->getZmax();
			
			HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
			//get parameters of hexagonal packing
			size_t l0,l1;
			size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
			size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
			hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
								n0_odd, n1_odd, n0_even, n1_even, 
								nn0_odd, nn1_odd, nn0_even, nn1_even);	
			
			// compute the location of particles
			// and compute the number of fluid particles
			for(size_t i=l0; i<=l1; i++) { 
				if((i+1)%2 != 0) { //odd-numbered layers
					for(size_t j=m0_odd; j<=m1_odd; j++) { 
						if((j+1)%2 != 0) { //odd-numbered rows 
							for(size_t k=n0_odd; k<=n1_odd; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue;	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									
									m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
                                                                        m_vMass[numParticle] = sqrt(2)*4*h_r*h_r*h_r/m_vVolume[numParticle];
									m_vPressure[numParticle] = states[p]->pressure(x,y,z);
									states[p]->velocity(x,y,z,
										m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
									m_vSoundSpeed[numParticle] = 
										m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									
								}	
								numParticle++;
								numParticleI++;	
								
							}
						} 
						else{ //even-numbered rows
							for(size_t k=n0_even; k<=n1_even; k++) {
								double x = hex3D.computeX(1,k);
								double y = hex3D.computeY(0,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue;	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									
									m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
									m_vPressure[numParticle] = states[p]->pressure(x,y,z);
                                                        		m_vMass[numParticle] = sqrt(2)*4*h_r*h_r*h_r/m_vVolume[numParticle];
									states[p]->velocity(x,y,z,
										m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
									m_vSoundSpeed[numParticle] = 
										m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									
								}
								numParticle++;
								numParticleI++;
								
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
								if(!objs[p]->operator()(x,y,z)) continue;
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									
									m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
									m_vPressure[numParticle] = states[p]->pressure(x,y,z);
                                                                        m_vMass[numParticle] = sqrt(2)*4*h_r*h_r*h_r/m_vVolume[numParticle];
									states[p]->velocity(x,y,z,
										m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
									m_vSoundSpeed[numParticle] = 
										m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									
								}	
								numParticle++;
								numParticleI++;
								
							}
						} 
						else { //even-numbered rows
							for(size_t k=nn0_even; k<=nn1_even; k++) {
								double x = hex3D.computeX(0,k);
								double y = hex3D.computeY(1,j);
								double z = hex3D.computeZ(i);
								if(!objs[p]->operator()(x,y,z)) continue; 	
								if(saveData) {
									m_vPositionX[numParticle] = x;
									m_vPositionY[numParticle] = y;
									m_vPositionZ[numParticle] = z;
									
									m_vVolume[numParticle] = 1./states[p]->density(x,y,z);
									m_vPressure[numParticle] = states[p]->pressure(x,y,z);
                                                                        m_vMass[numParticle] = sqrt(2)*4*h_r*h_r*h_r/m_vVolume[numParticle];
									states[p]->velocity(x,y,z,
										m_vVelocityU[numParticle],m_vVelocityV[numParticle],m_vVelocityW[numParticle]);
									m_vSoundSpeed[numParticle] = 
										m_pEOS->getSoundSpeed(m_vPressure[numParticle],1./m_vVolume[numParticle]);
									
								}
								numParticle++;
								numParticleI++;
								
							}
						}
					}	
				}    
			}
			
			boxes[p]->setNumber(numParticleI);
			
		}	
	
	}
	
	return numParticle;
}


void Initializer::computeNumParticleWithinSearchRadius() {
	
	size_t result = 0;
        srand(1);

	double h_r = 0.5*m_fInitParticleSpacing;

	if(m_iDimension==2) {
		
		// a small rectangular space 
		double xmin = -1.5*m_fNeiSearchRadius;
		double xmax =  1.5*m_fNeiSearchRadius;
		double ymin = -1.5*m_fNeiSearchRadius;
		double ymax =  1.5*m_fNeiSearchRadius;		
		
		HexagonalPacking2D hex2D(xmin,xmax,ymin,ymax,h_r);
		// get parameters of hexagonal packing
		size_t m0, m1, n0_odd, n1_odd, n0_even, n1_even;
		hex2D.getParameters(m0, m1, n0_odd, n1_odd, n0_even, n1_even);	
			
		//compute the location of particles 
		for(size_t j=m0; j<=m1; j++) { 
			if((j+1)%2 != 0) { //odd-numbered rows 
				for(size_t k=n0_odd; k<=n1_odd; k++) { 
					double x = hex2D.computeX(0,k);
					double y = hex2D.computeY(j);	
					if(sqrt(x*x+y*y)<=m_fNeiSearchRadius) result++;
				}
			} 
			else{ //even-numbered rows
				for(size_t k=n0_even; k<=n1_even; k++) {
					double x = hex2D.computeX(1,k);
					double y = hex2D.computeY(j);	
					if(sqrt(x*x+y*y)<=m_fNeiSearchRadius) result++;
				}
			}
		}
	}
	else if(m_iDimension==3) {
		
		// a small cube 
		double xmin = -1.5*m_fNeiSearchRadius;
		double xmax =  1.5*m_fNeiSearchRadius;
		double ymin = -1.5*m_fNeiSearchRadius;
		double ymax =  1.5*m_fNeiSearchRadius;
		double zmin = -1.5*m_fNeiSearchRadius;
		double zmax =  1.5*m_fNeiSearchRadius;
		
		HexagonalPacking3D hex3D(xmin, xmax, ymin, ymax, zmin, zmax, h_r);
		//get parameters of hexagonal packing
		size_t l0,l1;
		size_t m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
		size_t nn0_odd, nn1_odd, nn0_even, nn1_even; 
		hex3D.getParameters(l0, l1, m0_odd, m1_odd, m0_even, m1_even, 
						    n0_odd, n1_odd, n0_even, n1_even, 
							nn0_odd, nn1_odd, nn0_even, nn1_even);	
		
		// compute the location of particles
		// and compute the number of fluid particles
		for(size_t i=l0; i<=l1; i++) { 
			if((i+1)%2 != 0) { //odd-numbered layers
				for(size_t j=m0_odd; j<=m1_odd; j++) { 
					if((j+1)%2 != 0) { //odd-numbered rows 
						for(size_t k=n0_odd; k<=n1_odd; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					} 
					else{ //even-numbered rows
						for(size_t k=n0_even; k<=n1_even; k++) {
							double x = hex3D.computeX(1,k);
							double y = hex3D.computeY(0,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
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
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					} 
					else { //even-numbered rows
						for(size_t k=nn0_even; k<=nn1_even; k++) {
							double x = hex3D.computeX(0,k);
							double y = hex3D.computeY(1,j);
							double z = hex3D.computeZ(i);	
							if(sqrt(x*x+y*y+z*z)<=m_fNeiSearchRadius) result++;
						}
					}
				}	
			}    
		}	
	}

	m_iNumParticleWithinSearchRadius = result;
	
}


void Initializer::setBoundingBoxStartIndex() {
	
	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart
	
	size_t tmp=0;	
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) {
		// set the start index of fluid bound box
		m_vFluidBoundingBox[p]->setStartIndex(m_iFluidStartIndex+tmp); 
		size_t num = m_vFluidBoundingBox[p]->getNumber();	
		tmp += num;
		save<<"m_vFluidBoundingBox["<<p<<"]->m_iStartIndex "<<m_vFluidBoundingBox[p]->getStartIndex()<<endl;
		save<<"m_vFluidBoundingBox["<<p<<"]->m_iNumber "<<num<<endl;
	}
	/*
	tmp=0;	
	for(size_t p=0; p<m_vBoundaryBoundingBox.size(); p++) {
		// set the start index of boundary bound box
		m_vBoundaryBoundingBox[p]->setStartIndex(m_iBoundaryStartIndex+tmp); 
		size_t num = m_vBoundaryBoundingBox[p]->getNumber();	
		tmp += num;
		save<<"m_vBoundaryBoundingBox["<<p<<"]->m_iStartIndex "<<m_vBoundaryBoundingBox[p]->getStartIndex()<<endl;
		save<<"m_vBoundaryBoundingBox["<<p<<"]->m_iNumber "<<num<<endl;
	}
	*/
	
}


 
void Initializer::setObjectTag() {
	
	if(m_iDimension==1) return;

	// set the fluid object tag (fluid tag = the object num: 1, 2, 3,...)
	size_t tmp=0;
	for(size_t p=0; p<m_vFluidBoundingBox.size(); p++) { 
		m_vFluidBoundingBox[p]->setObjectTag(p+1);
		size_t num = m_vFluidBoundingBox[p]->getNumber();	
		for(size_t index=m_iFluidStartIndex+tmp; index<m_iFluidStartIndex+tmp+num; index++) {
			m_vObjectTag[index] = p+1;
		}
		tmp += num;
	}
	
	/*
	tmp=0;
	//set the boundary object tag, which is 0
	for(size_t p=0; p<m_vBoundaryBoundingBox.size(); p++) { 
		m_vBoundaryBoundingBox[p]->setObjectTag(0);
		size_t num = m_vBoundaryBoundingBox[p]->getNumber();	
		for(size_t index=m_iBoundaryStartIndex+tmp; index<m_iBoundaryStartIndex+tmp+num; index++) {
			m_vObjectTag[index] = 0;
		}
		tmp += num;
	}
	*/
	
}

void Initializer::readParamfile(const string& paramfileName) {
	
	/*
	m_iNumThreads 1
m_fStartTime 0
m_fEndTime 0.05
m_fWriteTimeInterval 0.0005
m_fCFLCoeff 0.999
m_iDimension 2
m_iFluidObjNum 1
m_vFluidObjectNames0 jet2dexp
m_vFluidObjectStateNames0 jet2dexpstate
m_iBoundaryObjNum 0
m_iRandomDirSplitOrder 1
m_iLPFOrder 2
m_iEOSChoice 1
m_fGamma 5
m_fPinf 0
m_fEinf 0
m_fInitParticleSpacing 0.05
m_fGravity 0
m_iMovingBoxForGhostParticle 1
m_iUseLimiter 0
m_fThresholdP 10
m_iUseCriticalPressure 0
m_fCriticalPressure 0
m_fNeiSearchRadius 0.3
m_iNumParticleWithinSearchRadius 127
m_iMaxNeighbourNum 889
m_iMaxNeighbourNumInOneDir 444
m_iNumRow2ndOrder 7
m_iNumRow1stOrder 3
m_iNumCol2ndOrder 5
m_iNumCol1stOrder 2
m_fInvalidPressure 0
m_fInvalidVolume 0
m_iTreeDepth 5
m_fContactLength 0.055
m_iFluidNum 55252
m_iBoundaryNum 0
m_iCapacity 193382
m_iFluidStartIndex 0
m_iBoundaryStartIndex 55252
m_iGhostStartIndex 55252
m_vFluidBoundingBox[0]->m_iStartIndex 0
m_vFluidBoundingBox[0]->m_iNumber 55252	
	*/
	ifstream ifs(paramfileName);
	
	unordered_map<string,string> table;
	
	string s1, s2;
	while((ifs >> s1 >> s2)) table[s1] = s2;
	
	cout<<"table.size()="<<table.size()<<endl;
	for(auto i:table) {
		cout<<i.first<<" "<<i.second<<endl;
	}

	if(table.find("m_iNumThreads") != table.end())
		m_iNumThreads = stoi(table["m_iNumThreads"]);
	else {
		cout<<"Cannot find m_iNumThreads"<<endl;
		assert(false);
	}

}






////////////////////////////////////////////////////////////////////////////////////////
// deprecated functions
////////////////////////////////////////////////////////////////////////////////////////

void Initializer::modifyLocalParSpacing() {
	
	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart
	save<<"(RESTART CHANGE): m_fInitParticleSpacing has been changed from "<<m_fInitParticleSpacing;	

	NeighbourSearcher* neiSearcher = new OctreeSearcher(m_iCapacity,m_iMaxNeighbourNum,m_iTreeDepth);

	neiSearcher->buildSearchStructure(m_vPositionX, m_vPositionY, m_vPositionZ, m_iFluidStartIndex, m_iFluidNum);
	
	#ifdef _OPENMP
	#pragma omp parallel  
	{
	int tid = omp_get_thread_num();	
	#endif

	int* neiList = new int[m_iCapacity*m_iMaxNeighbourNum]; // the index of the neighbours
	double neiDist[m_iMaxNeighbourNum]; // the distance from neighbours
	size_t numNei; // the number of neighbours found
	size_t neiListSize[m_iCapacity];
	double neiDistNearest[m_iCapacity];
	double initSearchRadius = 5*m_fInitParticleSpacing;
	
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=m_iFluidStartIndex; index<m_iFluidStartIndex+m_iFluidNum; index++) { 
		size_t neiListStartIndex = index*m_iMaxNeighbourNum;		
		
		#ifdef _OPENMP
		neiSearcher->searchNeighbour(
			m_vPositionX[index],m_vPositionY[index],m_vPositionZ[index],initSearchRadius, 
			neiList+neiListStartIndex,neiDist,numNei,tid,index); // output
			
		#else
		neiSearcher->searchNeighbour(
			m_vPositionX[index],m_vPositionY[index],m_vPositionZ[index],initSearchRadius, 
			neiList+neiListStartIndex,neiDist,numNei,index); // output
		#endif	
		
		assert(numNei>0);
		
		double localRadius = neiDist[0] * m_fTimesNeiSearchRadius;// define the local neighoburhood

		#ifdef _OPENMP
		neiSearcher->searchNeighbour(
			m_vPositionX[index],m_vPositionY[index],m_vPositionZ[index],localRadius, 
			neiList+neiListStartIndex,neiDist,numNei,tid,index); // output
			
		#else
		neiSearcher->searchNeighbour(
			m_vPositionX[index],m_vPositionY[index],m_vPositionZ[index],localRadius, 
			neiList+neiListStartIndex,neiDist,numNei,index); // output
		#endif
		
		assert(numNei>0);

		neiListSize[index] = numNei;
		neiDistNearest[index] = neiDist[0];
		
	}
	
	

	// use the defined neighbourhood/radius to search neighbours 
	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=m_iFluidStartIndex; index<m_iFluidStartIndex+m_iFluidNum; index++) { 
		size_t neiListStartIndex = index*m_iMaxNeighbourNum;		
		
		double sumDist = 0;
		for(size_t j=0; j<neiListSize[index]; j++) {
			size_t neiIndex = neiList[neiListStartIndex+j];
			sumDist += neiDistNearest[neiIndex];
		}
		
		sumDist += neiDistNearest[index];

		m_vLocalParSpacing[index] = sumDist/(double)(neiListSize[index]+1);
		
	}	
	
//	double sum = 0;
//	#ifdef _OPENMP
//	#pragma omp for default(shared) reduction(+:sum)
//	#endif
//	for(size_t index=m_iFluidStartIndex; index<m_iFluidStartIndex+m_iFluidNum; index++) { 	
//		sum += neiDistNearest[index];
//	}
//	m_fInitParticleSpacing = sum/(double)(m_iFluidNum);

	delete[] neiList;

	#ifdef _OPENMP  
	} // end of parallel region
	#endif
	
	save<<" to "<<m_fInitParticleSpacing<<endl;
	
	delete neiSearcher;

}

void Initializer::modifyNumParticleWithinSearchRadius() {

	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart
	
	#ifdef _OPENMP
	omp_set_num_threads(min(omp_get_max_threads(), m_iNumThreads));
	size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
	vector<double> numNei(numThreads,0);
	#else
	double numNei = 0;
	#endif
	
	#ifdef _OPENMP
	#pragma omp parallel 
	{ // begin openmp parallel region
	int tid = omp_get_thread_num();
	#endif

	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t i=m_iFluidStartIndex; i<m_iFluidStartIndex+m_iFluidNum; i++) {
		size_t num = 0;	
		double x0 = m_vPositionX[i];
		double y0 = m_vPositionY[i];
		double z0 = m_vPositionZ[i];
		for(size_t j=m_iFluidStartIndex; j<m_iFluidStartIndex+m_iFluidNum; j++) {
			if(i==j) continue;
			double x = m_vPositionX[j];
			double y = m_vPositionY[j];
			double z = m_vPositionZ[j];
			double dist = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
			if(dist <= m_fNeiSearchRadius) num++;
		}
		#ifdef _OPENMP
		numNei[tid] += num;
		#else
		numNei += num;
		#endif
	}
	
	#ifdef _OPENMP
	} // end openmp parallel region
	#endif


	save<<"(RESTART CHANGE): m_iNumParticleWithinSearchRadius has been changed from "
	<<m_iNumParticleWithinSearchRadius;
	
	#ifdef _OPENMP
	double sum = 0;
	for(size_t i=0; i<numThreads; i++) sum += numNei[i];
	m_iNumParticleWithinSearchRadius = sum/(double)(m_iFluidNum);
	#else
	m_iNumParticleWithinSearchRadius = numNei/(double)(m_iFluidNum);
	#endif
	
	save<<" to "<<m_iNumParticleWithinSearchRadius<<endl;

}


void Initializer::modifyInitParticleSpacing() {
	/*	
	NeighbourSearcher* neiSearcher = new OctreeSearcher(m_iCapacity,m_iMaxNeighbourNum,m_iTreeDepth);

	neiSearcher->buildSearchStructure(m_vPositionX, m_vPositionY, m_vPositionZ, m_iFluidStartIndex, m_iFluidNum);
	
	#ifdef _OPENMP
	#pragma omp parallel  
	{
	int tid = omp_get_thread_num();	
	#endif

	int neiListTemp[m_iMaxNeighbourNum]; // the index of the neighbours
	double neiListDistTemp[m_iMaxNeighbourNum]; // the distance from neighbours
	size_t numNeiFound; // the number of neighbours found

	#ifdef _OPENMP
	#pragma omp for
	#endif
	for(size_t index=m_iFluidStartIndex; index<m_iFluidStartIndex+m_iFluidNum; index++) { 
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
		
		neighbourListSize[index] = numNeiFound;		    	
	}
	#ifdef _OPENMP  
	} // end of parallel region
	#endif
	*/


	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart
	
//	#ifdef _OPENMP
//	omp_set_num_threads(min(omp_get_max_threads(), m_iNumThreads));
//	size_t numThreads = min(omp_get_max_threads(), m_iNumThreads);
//	vector<double> sumDist(numThreads,0);
//	#else
	double sumDist = 0;
//	#endif
	
//	#ifdef _OPENMP
//	#pragma omp parallel 
//	{ // begin openmp parallel region
//	int tid = omp_get_thread_num();
//	#endif

//	#ifdef _OPENMP
//	#pragma omp for
//	#endif
	for(size_t i=m_iFluidStartIndex; i<m_iFluidStartIndex+m_iFluidNum; i++) {
		
		double minDist = numeric_limits<double>::max();
		double x0 = m_vPositionX[i];
		double y0 = m_vPositionY[i];
		double z0 = m_vPositionZ[i];
		for(size_t j=m_iFluidStartIndex; j<m_iFluidStartIndex+m_iFluidNum; j++) {
			if(i==j) continue;
			double x = m_vPositionX[j];
			double y = m_vPositionY[j];
			double z = m_vPositionZ[j];
			double dist = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
			minDist = min(minDist,dist);
		}
//		#ifdef _OPENMP
//		sumDist[tid] += minDist;
//		#else
		sumDist += minDist;
//		#endif
	}
	
//	#ifdef _OPENMP
//	} // end openmp parallel region
//	#endif


	save<<"(RESTART CHANGE): m_fInitParticleSpacing has been changed from "<<m_fInitParticleSpacing;
	
//	#ifdef _OPENMP
//	double sum = 0;
//	for(size_t i=0; i<numThreads; i++) sum += sumDist[i];
//	m_fInitParticleSpacing = sum/(double)(m_iFluidNum);
//	#else
	m_fInitParticleSpacing = sumDist/(double)(m_iFluidNum);
//	#endif
	
	save<<" to "<<m_fInitParticleSpacing<<endl;

}

void Initializer::modifyInitNeighbourSearchRadius() {

	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart	

	save<<"(RESTART CHANGE): m_fNeiSearchRadius has been changed from "<<m_fNeiSearchRadius;
	m_fNeiSearchRadius = m_fTimesNeiSearchRadius*m_fInitParticleSpacing;
	save<<" to "<<m_fNeiSearchRadius<<endl;

}


void Initializer::modifyInitContactLength() {

	ofstream save(m_sFilenameSaveInit, ofstream::app); // save info for restart	

	save<<"(RESTART CHANGE): m_fContactLength has been changed from "<<m_fContactLength;
	m_fContactLength = m_fTimesContactLength*m_fInitParticleSpacing;
	save<<" to "<<m_fContactLength<<endl;

}








////////////////////////////////////////////////////////////////////////////////////////
// End of Initializer
////////////////////////////////////////////////////////////////////////////////////////
