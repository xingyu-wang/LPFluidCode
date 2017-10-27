#include "vtkfile_processor.h"
#include "neighbour_searcher.h"
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include "omp.h"
#include <iomanip> //setw

using namespace std;


VTKFileProcessor::VTKFileProcessor(): numParticle(0), time(0), ns(nullptr) {}

VTKFileProcessor::~VTKFileProcessor() {delete ns;}


string VTKFileProcessor::rightFlush(size_t writeStep, size_t numDigits) {
	
	assert(pow(10,numDigits) > writeStep);

	string result;

	if(writeStep == 0) numDigits--;
	for(size_t i=writeStep; i!=0; i /= 10) numDigits--;

	for( ; numDigits>0; numDigits--) result.push_back('0');
	
	result += to_string(writeStep); 
	
	return result;

}


void VTKFileProcessor::read(const std::string& filename) {
	
	dealloc();

	ifstream ifs(filename);
	ifs.precision(9);	
	string s;
	
	getline(ifs,s); // Skip 1 line
	ifs >> s >> s >> s >> s >> time;	

	for(int line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	ifs >> s >> numParticle >> s; // Read the number of particles (fluid + boundary)
	
	alloc(); // allocate space 

	for(size_t i=0; i<numParticle; i++) { // Read location
		ifs >> x[i] >> y[i] >> z[i];
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<numParticle; i++) { // Read velocity
		ifs >> u[i] >> v[i] >> w[i];
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines	
	for(size_t i=0; i<numParticle; i++) { // Read pressure
		ifs >> pressure[i];
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<numParticle; i++) { // Read volume
		ifs >> volume[i];
	}

        for(size_t line=1; line<=2; line++) getline(ifs,s); // Skip 1 lines
        for(size_t i=0; i<numParticle; i++) getline(ifs,s); // Skpi volume_voronoi

        for(size_t line=1; line<=2; line++) getline(ifs,s); // Skip 1 lines
        for(size_t i=0; i<numParticle; i++) getline(ifs,s); // Skpi density

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<numParticle; i++) { // Read sound speed
		ifs >> soundSpeed[i];
	}

        for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
        for(size_t i=0; i<numParticle; i++) { // Read mass
                ifs >> mass[i];
        }	

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<numParticle; i++) { // Read object tag
		ifs >> objectTag[i];	
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=0; i<numParticle; i++) { // Read particle spacing
		ifs >> localParSpacing[i];	
	}

	cout<<"Number of particles: "<<numParticle<<endl;
//	cout<< x[0] <<" "<<y[0]<<" "<<z[0]<<endl;
	std::cout.precision(9);
//        cout<<volume[0]<<endl;
//        cout<<volume[numParticle-1]<<endl;
//        cout<<soundSpeed[0]<<endl;
//        cout<<soundSpeed[numParticle-1]<<endl;
//	cout<<mass[0]<<endl;
//	cout<<mass[numParticle-1]<<endl;
//        cout<<localParSpacing[0]<<endl;
//        cout<<localParSpacing[numParticle-1]<<endl;

}

void VTKFileProcessor::interpolation_to_mesh(double xLeft, double xRight, double yLeft, double yRight, double zLeft, double zRight, int xN, int yN, int zN, const std::string& filename){
        if(numParticle==0) {
                cout<<"Have not read a vtkfile yet, there is no data in this class!!!"<<endl;
                assert(false);
        }

        size_t maxNeiNum = 10000;
        int treeDepth = 5;
        double total_volume = (xRight-xLeft)*(yRight-yLeft)*(zRight-zLeft);
        int numMesh = xN*yN*zN;
        double spacing=cbrt(total_volume/numParticle);
        double radius = 3.*spacing;
        if(ns==nullptr) {
                this->ns = new OctreeSearcher(numParticle,maxNeiNum,treeDepth);
                ns->buildSearchStructure(&x[0], &y[0], &z[0],numParticle);
        }

        assert(spacing!=0);
	vector<double> massMesh(numMesh);

        int neiList[maxNeiNum];
        double neiListDist[maxNeiNum];
        for(int i=0; i<numMesh; i++) {
		if((10*i)%numMesh==0) cout<<"Interpolating... "<<(100*i)/numMesh<<"% completed."<<endl;
		int xi=i%xN;
		int yi=(i/xN)%yN;
		int zi=i/xN/yN;
                double xM = xLeft + (xRight-xLeft)*xi/xN;
                double yM = yLeft + (yRight-yLeft)*yi/yN;
                double zM = zLeft + (zRight-zLeft)*zi/zN;
                size_t numNeiFound = 0;
//		cout<<xM<<" "<<yM<<" "<<zM<<endl;
                ns->searchNeighbour(xM,yM,zM,radius, neiList, neiListDist, numNeiFound);
                //cout<<"xCen="<<xCen<<"        numNeiFound="<<numNeiFound<<endl;
//		cout<<numNeiFound<<endl;
		if(numNeiFound==0){ //cout<<"Cannot find neighbour for grid point "<< xM<<" "<<yM<<" "<<zM<<endl;
			massMesh[i]=0;
//		assert(numNeiFound);
//
//		massMesh[i]=mass[neiList[0]];
		}
		else{
//		cout<<numNeiFound<<endl;
		int numAve=10;
		if(numNeiFound<numAve)	numAve=numNeiFound;
		double tempSum=0;
		for (size_t j=0; j<numAve; j++)
//			tempSum+=pressure[neiList[j]];
			tempSum+=1.0/volume[neiList[j]];
		massMesh[i]=tempSum/numAve;
//		cout<<massMesh[i]<<endl;
		}
/*                int indexL = -1, indexR = -1, indexExact = -1;
                for(size_t j=0; j<numNeiFound; j++) {
                        int neiIndex = neiList[j];
                        if(x[neiIndex]==xCen) {
                                indexExact = neiIndex;
                                break;
                        }
                        else if(x[neiIndex]>xCen) {
                                if(indexR==-1) indexR = neiIndex;
                        }
                        else {
                                if(indexL==-1) indexL = neiIndex;
                        }

                        if(indexR!=-1 && indexL!=-1) break;
                }


                if(indexExact!=-1) {
                        result_x[i] = xCen;
                        result_volume[i]     = volume[indexExact];
                        result_pressure[i]   = pressure[indexExact];
                        result_soundSpeed[i] = soundSpeed[indexExact];
                        //cout<<"xCen="<<xCen<<", p="<<result_pressure[i]<<endl;

                }
                else if(indexL!=-1 && indexR!=-1) {
                        result_x[i] = xCen;
                        double xDiffL = xCen-x[indexL], xDiffR = x[indexR]-xCen;
                        result_volume[i] = linear_interpolate(xDiffL,volume[indexL],xDiffR,volume[indexR]);
                        result_pressure[i] = linear_interpolate(xDiffL,pressure[indexL],xDiffR,pressure[indexR]);
                        result_soundSpeed[i] = linear_interpolate(xDiffL,soundSpeed[indexL],xDiffR,soundSpeed[indexR]);
                        //cout<<"xCen="<<xCen<<", xleft="<<x[indexL]<<", xright="<<x[indexR]<<endl;
                        //cout<<"pL="<<pressure[indexL]<<", pR="<<pressure[indexR]<<", p="<<result_pressure[i]<<endl;
                }
                else assert(false);*/
        }

//      delete ns;

//        vector<string> titles = {"x", "volume", "pressure", "sound_speed"};
//        vector<double*> inputs = {&massMesh[0]};

//        writeMesh(filename,titles,xLeft,xRight,yLeft,yRight,zLeft,zRight,xN,yN,zN,&inputs);
	cout<<"Writing..."<<endl;
	writeMesh(&(massMesh[0]),xLeft,xRight,yLeft,yRight,zLeft,zRight,xN,yN,zN,filename);
}

void VTKFileProcessor::writeMesh(double* massMesh, double xLeft, double xRight, double yLeft, double yRight, double zLeft, double zRight, int xN, int yN, int zN, const std::string& filename){
        FILE *outfile;
        outfile = fopen(filename.c_str(), "w");
        if(outfile==nullptr) {
                printf("Unable to open output file: %s\n",filename.c_str());
                assert(false);
        }

        size_t startIndex = 0;
	size_t numMesh = xN*yN*zN;
        size_t endIndex = startIndex + numMesh;

        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"vtk output\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
	fprintf(outfile,"DIMENSIONS %d %d %d\n",xN,yN,zN);
	fprintf(outfile,"X_COORDINATES %d float\n",xN);
	for(size_t xi=0;xi<xN;xi++)
		fprintf(outfile,"%f ", xLeft+(xRight-xLeft)*xi/xN);
        fprintf(outfile,"\nY_COORDINATES %d float\n",yN);
        for(size_t yi=0;yi<yN;yi++)
                fprintf(outfile,"%f ", yLeft+(yRight-yLeft)*yi/yN);
        fprintf(outfile,"\nZ_COORDINATES %d float\n",zN);
        for(size_t zi=0;zi<zN;zi++)
                fprintf(outfile,"%f ", zLeft+(zRight-zLeft)*zi/zN);
//	fprintf(outfile,"\nCELL_DATA 1\n");
        fprintf(outfile,"\nPOINT_DATA %d\n",numMesh);
//	fprintf(outfile,"FIELD FieldData 1\n");
//	fprintf(outfile,"mass 1 %d float\n",numMesh);
	fprintf(outfile,"SCALARS scalars float\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t index=startIndex;index<endIndex;index++)
	{
//                if((10*index)%numMesh==0) cout<<"Writing... "<<(100*index)/numMesh<<"% completed."<<endl;
		fprintf(outfile,"%.16g ",massMesh[index]);
		if(((index-startIndex)%xN)==(xN-1)) fprintf(outfile,"\n");
	}
        fclose(outfile);
	cout<<"Writing completed"<<endl;

}


void VTKFileProcessor::combine(const std::string& filename) {
	
	ifstream ifs(filename);
	
	double time_;
	size_t num;
	string s;
	
	getline(ifs,s); // Skip 1 line
	ifs >> s >> s >> s >> s >> time_;	
	
	assert(time==time_); // combined file should be at the same time

	for(int line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	ifs >> s >> num >> s; // Read the number of particles (fluid + boundary)
	
	size_t oldNumParticle = numParticle;

	augment(num);

	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read location
		ifs >> x[i] >> y[i] >> z[i];
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read velocity
		ifs >> u[i] >> v[i] >> w[i];
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines	
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read pressure
		ifs >> pressure[i];
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read volume
		ifs >> volume[i];
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read sound speed
		ifs >> soundSpeed[i];
	}
	
	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read object tag
		ifs >> objectTag[i];	
	}

	for(size_t line=1; line<=3; line++) getline(ifs,s); // Skip 2 lines
	for(size_t i=oldNumParticle; i<numParticle; i++) { // Read object tag
		ifs >> localParSpacing[i];	
	}
		
}

void VTKFileProcessor::write(const std::string& filename) {	
	
	FILE *outfile;
	outfile = fopen(filename.c_str(), "w");
	if(outfile==nullptr) {
		printf("Unable to open output file: %s\n",filename.c_str()); 
		assert(false);
	}
	
	size_t startIndex = 0;
	size_t endIndex = startIndex + numParticle;

	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.16g\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	
	fprintf(outfile,"POINTS %ld double\n",numParticle);	
	for(size_t i = startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g %.16g %.16g\n",x[i], y[i], z[i]);
	
	fprintf(outfile,"POINT_DATA %ld\n",numParticle);
	fprintf(outfile,"VECTORS Velocity double\n");	
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g %.16g %.16g\n",u[i], v[i], w[i]);

	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",pressure[i]);
		
	fprintf(outfile,"SCALARS volume double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",volume[i]);
	
	fprintf(outfile,"SCALARS sound_speed double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",soundSpeed[i]);

	fprintf(outfile,"SCALARS object_tag int\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%d\n",objectTag[i]);	
	
	fprintf(outfile,"SCALARS local_par_spacing double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for(size_t i=startIndex; i<endIndex; i++)
		fprintf(outfile,"%.16g\n",localParSpacing[i]);

	fclose(outfile);
	
}


void VTKFileProcessor::rotate(double angle) {

	for(size_t i=0; i<numParticle; i++) {
		double xtmp, ytmp;
		xtmp = x[i]*cos(angle) - y[i]*sin(angle);
		ytmp = x[i]*sin(angle) + y[i]*cos(angle);
		x[i]=xtmp;
		y[i]=ytmp;

		double utmp, vtmp;
		utmp = u[i]*cos(angle) - v[i]*sin(angle);
		vtmp = u[i]*sin(angle) + v[i]*cos(angle);
		u[i]=utmp;
		v[i]=vtmp;
	}

}

void VTKFileProcessor::translate(double xshift, double yshift, double zshift) {

	for(size_t i=0; i<numParticle; i++) {
		x[i] += xshift;
		y[i] += yshift;
		z[i] += zshift;
	}


}

void VTKFileProcessor::changeTag(int tag) {

	for(size_t i=0; i<numParticle; i++) {
		objectTag[i] = tag;
	}

}

void VTKFileProcessor::changeVelocity(double u_, double v_, double w_) {

	for(size_t i=0; i<numParticle; i++) {
		u[i] += u_;
		v[i] += v_;
		w[i] += w_;
	}

}



void VTKFileProcessor::reduce2dY(double yCen, double xLeft, double xRight, double spacing, const string& fname) {
	
	if(numParticle==0) {
		cout<<"Have not read a vtkfile yet, there is no data in this class!!!"<<endl;
		assert(false);
	}
	
	size_t maxNeiNum = 500;
	int treeDepth = 5;
	double radius = 3.*spacing;
	if(ns==nullptr) {
		this->ns = new OctreeSearcher(numParticle,maxNeiNum,treeDepth);
		ns->buildSearchStructure(&x[0], &y[0], &z[0],numParticle);
	}
	
	assert(spacing!=0);
	int K = (xRight-xLeft)/spacing;
	vector<double> result_x(K+1); 
	vector<double> result_volume(K+1);
	vector<double> result_pressure(K+1);
	vector<double> result_soundSpeed(K+1);
	
	int neiList[maxNeiNum];
	double neiListDist[maxNeiNum];
	for(int i=0; i<=K; i++) {
		double xCen = xLeft + i*spacing;	
		size_t numNeiFound = 0;
		ns->searchNeighbour(xCen,yCen,0,radius, neiList, neiListDist, numNeiFound);
		//cout<<"xCen="<<xCen<<"	numNeiFound="<<numNeiFound<<endl;

		int indexL = -1, indexR = -1, indexExact = -1;
		for(size_t j=0; j<numNeiFound; j++) {
			int neiIndex = neiList[j];
			if(x[neiIndex]==xCen) {
				indexExact = neiIndex;
				break;
			}	
			else if(x[neiIndex]>xCen) {
				if(indexR==-1) indexR = neiIndex;
			}
			else {
				if(indexL==-1) indexL = neiIndex;
			}

			if(indexR!=-1 && indexL!=-1) break;
		}
		
		
		if(indexExact!=-1) {
			result_x[i] = xCen;
			result_volume[i]     = volume[indexExact];
			result_pressure[i]   = pressure[indexExact];
			result_soundSpeed[i] = soundSpeed[indexExact];
			//cout<<"xCen="<<xCen<<", p="<<result_pressure[i]<<endl;

		}
		else if(indexL!=-1 && indexR!=-1) {
			result_x[i] = xCen;
			double xDiffL = xCen-x[indexL], xDiffR = x[indexR]-xCen;
			result_volume[i] = linear_interpolate(xDiffL,volume[indexL],xDiffR,volume[indexR]);
			result_pressure[i] = linear_interpolate(xDiffL,pressure[indexL],xDiffR,pressure[indexR]);	
			result_soundSpeed[i] = linear_interpolate(xDiffL,soundSpeed[indexL],xDiffR,soundSpeed[indexR]);
			//cout<<"xCen="<<xCen<<", xleft="<<x[indexL]<<", xright="<<x[indexR]<<endl;
			//cout<<"pL="<<pressure[indexL]<<", pR="<<pressure[indexR]<<", p="<<result_pressure[i]<<endl;
		}
		else assert(false);
	}	
	
//	delete ns;

	vector<string> titles = {"x", "volume", "pressure", "sound_speed"};
	vector<double*> inputs = {&result_x[0], &result_volume[0], &result_pressure[0], &result_soundSpeed[0]};
	
	write1d(fname,titles,inputs,K+1);

	
}


void VTKFileProcessor::reduce3dYZ(double yCen, double zCen, double xLeft, double xRight, double spacing, const string& fname) {
	
	if(numParticle==0) {
		cout<<"Have not read a vtkfile yet, there is no data in this class!!!"<<endl;
		assert(false);
	}
	
	size_t maxNeiNum = 1500;
	int treeDepth = 5;
	double radius = 3.*spacing;
	if(ns==nullptr) {
		this->ns = new OctreeSearcher(numParticle,maxNeiNum,treeDepth);
		ns->buildSearchStructure(&x[0], &y[0], &z[0],numParticle);
	}
	
	assert(spacing!=0);
	int K = (xRight-xLeft)/spacing;
	vector<double> result_x(K+1); 
	vector<double> result_volume(K+1);
	vector<double> result_pressure(K+1);
	vector<double> result_soundSpeed(K+1);
	
	int neiList[maxNeiNum];
	double neiListDist[maxNeiNum];
	for(int i=0; i<=K; i++) {
		double xCen = xLeft + i*spacing;	
		size_t numNeiFound = 0;
		ns->searchNeighbour(xCen,yCen,zCen,radius, neiList, neiListDist, numNeiFound);
		//cout<<"xCen="<<xCen<<"	numNeiFound="<<numNeiFound<<endl;

		int indexL = -1, indexR = -1, indexExact = -1;
		for(size_t j=0; j<numNeiFound; j++) {
			int neiIndex = neiList[j];
			if(x[neiIndex]==xCen) {
				indexExact = neiIndex;
				break;
			}	
			else if(x[neiIndex]>xCen) {
				if(indexR==-1) indexR = neiIndex;
			}
			else {
				if(indexL==-1) indexL = neiIndex;
			}

			if(indexR!=-1 && indexL!=-1) break;
		}
		
		
		if(indexExact!=-1) {
			result_x[i] = xCen;
			result_volume[i]     = volume[indexExact];
			result_pressure[i]   = pressure[indexExact];
			result_soundSpeed[i] = soundSpeed[indexExact];
			//cout<<"xCen="<<xCen<<", p="<<result_pressure[i]<<endl;

		}
		else if(indexL!=-1 && indexR!=-1) {
			result_x[i] = xCen;
			double xDiffL = xCen-x[indexL], xDiffR = x[indexR]-xCen;
			result_volume[i] = linear_interpolate(xDiffL,volume[indexL],xDiffR,volume[indexR]);
			result_pressure[i] = linear_interpolate(xDiffL,pressure[indexL],xDiffR,pressure[indexR]);	
			result_soundSpeed[i] = linear_interpolate(xDiffL,soundSpeed[indexL],xDiffR,soundSpeed[indexR]);
			//cout<<"xCen="<<xCen<<", xleft="<<x[indexL]<<", xright="<<x[indexR]<<endl;
			//cout<<"pL="<<pressure[indexL]<<", pR="<<pressure[indexR]<<", p="<<result_pressure[i]<<endl;
		}
		else {
			result_x[i] = xCen;
			result_volume[i]     = 0;
			result_pressure[i]   = 0;
			result_soundSpeed[i] = 0;	
		}
	}	
	


	vector<string> titles = {"x", "volume", "pressure", "sound_speed"};
	vector<double*> inputs = {&result_x[0], &result_volume[0], &result_pressure[0], &result_soundSpeed[0]};
	
	write1d(fname,titles,inputs,K+1);

	
}



void VTKFileProcessor::write1d(const string& fname, const vector<string>& titles, 
const vector<double*>& inputs, size_t num) {

	ofstream ofs(fname);	
	if(!ofs.is_open()) {
		cout<<"Unable to open file: "<<fname<<endl;	
		assert(false);
	}
	
	// print the titles
	for(const auto& ti: titles) 
		ofs<<ti<<setw(24);
	ofs<<endl;
	
	for(size_t i=0; i<num; i++) {
		ofs.precision(16);
		ofs<<left<<setw(24);
		for(size_t j=0; j<inputs.size(); j++) {
			ofs<<inputs[j][i]<<setw(24);
		}
		ofs<<endl;
	}

}


double VTKFileProcessor::linear_interpolate(double diff0, double value0, double diff1, double value1) {
	double sum = diff0+diff1;
	assert(sum!=0);
	return value0*(diff1/sum) + value1*(diff0/sum);
}



void VTKFileProcessor::alloc() {
	x = vector<double> (numParticle,0);
	y = vector<double> (numParticle,0);
	z = vector<double> (numParticle,0);
	u = vector<double> (numParticle,0);
	v = vector<double> (numParticle,0);
	w = vector<double> (numParticle,0);
	volume = vector<double> (numParticle,0);
	pressure = vector<double> (numParticle,0);
	soundSpeed = vector<double> (numParticle,0);
	objectTag = vector<int> (numParticle,0);
	localParSpacing = vector<double> (numParticle,0);
	mass = vector<double> (numParticle,0);
}

void VTKFileProcessor::dealloc() {
	x.clear();
	y.clear();
	z.clear();
	u.clear();
	v.clear();
	w.clear();
	volume.clear();
	pressure.clear();
	soundSpeed.clear();
	objectTag.clear();
	localParSpacing.clear();
	mass.clear();
}


void VTKFileProcessor::augment(size_t added_size) {
	assert(added_size > 0);
	numParticle += added_size;
	x.resize(numParticle);
	y.resize(numParticle);
	z.resize(numParticle);
	u.resize(numParticle);
	v.resize(numParticle);
	w.resize(numParticle);
	volume.resize(numParticle);
	pressure.resize(numParticle);
	soundSpeed.resize(numParticle);
	objectTag.resize(numParticle);
	localParSpacing.resize(numParticle);
	mass.resize(numParticle);
}
