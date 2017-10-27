#include "initializer.h"
#include "neighbour_searcher.h"
#include "particle_data.h"
#include "particle_viewer.h"
#include "lp_solver.h"
#include "time_controller.h"

#include <iostream>
#include <cassert>
#include <sys/stat.h> // mkdir()
#include <cstring> // strerror()
#include <vector>
using namespace std;

//#define LP_DEBUG
// execute as ./lp -i <inputfile name> -o <output directory name> 

int main(int argc, const char* argv[]) {	
	
	string inputfileName;
	string datafileName; // restart datafile
	string outputfileNameAll, outputfileNameFluid;		
	string debugfileName = "debug"; 
	bool ifDebug = false;
#ifdef LP_DEBUG
        cout<<"Begin"<<endl;
#endif
	
	// Assign command line arguments to variables
	for(int i=1; i<argc; i++) { // argv[0] is name of the executable 
		if(!strcmp(argv[i], "-i")) { // input; strcmp returns 0 if two c_str are the same
			if(i+1 >= argc || argv[i+1][0]=='-') { // no inputfile name following -i option
				cerr<<"ERROR: Needs to provide inputfile name!"<<endl;
				exit(1);
			}
			inputfileName = argv[i+1];
			cout<<"input file name: "<<inputfileName<<endl;
			if(i+2 < argc && argv[i+2][0]!='-') { // there is second input file (restart data file)
				datafileName = argv[i+2];
				cout<<"resatrt data file name: "<<datafileName<<endl;
			}
		}
		else if(!strcmp(argv[i], "-o")) { // output
			if(i+1 >= argc || argv[i+1][0]=='-') { // no outputfile name following -o option
				cerr<<"ERROR: Needs to provide outputfile name!"<<endl;
				exit(1);
			}
			if(mkdir(argv[i+1], 0777) == -1) { //mkdir takes const char*
				cerr<<"ERROR:"<<strerror(errno)<<endl;
				exit(1);
			}
			outputfileNameAll = string(argv[i+1]) + "/out";
			outputfileNameFluid = string(argv[i+1]) + "/out_fluid";
			cout<<"output file name: "<<outputfileNameAll<<endl;
			cout<<"output file name: "<<outputfileNameFluid<<endl;
		}
		else if(!strcmp(argv[i], "-d")) { // debug
			ifDebug = true;
			if(i+1 >= argc || argv[i+1][0]=='-') // no debugfile name following -d option is fine, use default name
				continue;
			debugfileName = argv[i+1];
			cout<<"debug file name: "<<debugfileName<<endl;
		}	
	}
	if(inputfileName.empty()) {
		cerr<<"ERROR: Needs to provide -i option!"<<endl;
		exit(1);
	}
	if(outputfileNameAll.empty() || outputfileNameFluid.empty()) {
		cerr<<"ERROR: Needs to provide -o option!"<<endl;
		exit(1);
	}
#ifdef LP_DEBUG
	cout<<"Initialization"<<endl;
#endif
	Initializer* init;
	if(datafileName.empty())
		init = new Initializer(inputfileName, ifDebug, debugfileName);
	else // restart
		init = new Initializer(inputfileName, datafileName, ifDebug, debugfileName);
#ifdef LP_DEBUG
	cout<<"ParticleData"<<endl;;
#endif
	// load boundary and fluid geoemtry/state into particleData
	ParticleData* pData = new ParticleData(*init);
#ifdef LP_DEBUG
	cout<<"ParticleSolver"<<endl;
#endif	
	//TODO: do this in command line using eg. -all -fluid or in input file as options entered by user
	//const string outputfileNameAll = string(argv[4]) + "/out";		
	//const string outputfileNameFluid = string(argv[4]) + "/out_fluid";
	ParticleViewer* pViewerAll, *pViewerFluid;
#ifdef LP_DEBUG
        cout<<"ParticleSolver1"<<endl;
#endif	
	NeighbourSearcher* neiSearcher;
#ifdef LP_DEBUG
        cout<<"ParticleSolver2"<<endl;
#endif	
	LPSolver* lpSolver;
#ifdef LP_DEBUG
        cout<<"ParticleSolver3"<<endl;
#endif	
	if(init->getDimension() == 1) {
		pViewerAll = new TXTParticleViewer1D(pData, "all", outputfileNameAll);	
		pViewerFluid = new TXTParticleViewer1D(pData, "fluid", outputfileNameFluid);
		lpSolver = new HyperbolicLPSolver1D(*init, pData);	
	}
	else {	
#ifdef LP_DEBUG
        cout<<"VTKParticleSolver"<<endl;
#endif
		pViewerAll = new VTKParticleViewer(pData, "all", outputfileNameAll);	
		pViewerFluid = new VTKParticleViewer(pData, "fluid", outputfileNameFluid);
#ifdef LP_DEBUG
        cout<<"OctreeSearcher"<<endl;
#endif		
	neiSearcher = new OctreeSearcher(*init); // neighbour is used only in 2D and 3D
#ifdef LP_DEBUG
        cout<<"LPSolver"<<endl;
#endif
	lpSolver = new HyperbolicLPSolver(*init, pData, neiSearcher);
	}
#ifdef LP_DEBUG
        cout<<"ParticleViewer"<<endl;
#endif
	vector<ParticleViewer*> viewers({pViewerAll,pViewerFluid});	

	TimeController* timeControl = new DefaultTimeController(*init, lpSolver, viewers);

	delete init; // initialization finished
#ifdef LP_DEBUG
	cout<<"Begin to solve"<<endl;
#endif	
	return timeControl->solve(); 

}
