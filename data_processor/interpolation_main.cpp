#include <iostream>
#include "vtkfile_processor.h"
#include <cmath>

const double PI = atan(1) * 4.;

using namespace std;

int main(int argc, const char* argv[]) {

	double xLeft= -0, xRight = 0.20;
        double yLeft= -0.06, yRight = 0.06;
        double zLeft= -0.06, zRight = 0.06;
        double xN = 100, yN = 100, zN = 100;

        string ifname0 = "/home/xingyu/lp_output/lp_nozzle_rothe_delete_3d_symmetric/out_fluid0000050.vtk";
        string ofname0 = "/home/xingyu/vtktest/rothe_delete_3d_symmetric.vtk";
        VTKFileProcessor p0;
        p0.read(ifname0);
        p0.interpolation_to_mesh(xLeft, xRight, yLeft, yRight, zLeft, zRight, xN, yN, zN, ofname0);
        cout<<"completed"<<endl;

/*	double xLeft= -0.5, xRight = 0.5;
        double yLeft= -0.5, yRight = 0.5;
        double zLeft= -0.5, zRight = 0.5;
	double xN = 50, yN = 50, zN = 50;

      	string ifname0 = "/home/xingyu/lp_output/lp_rayleightaylor3d_lw/out_fluid0000000.vtk";
        string ofname0 = "/home/xingyu/vtktest/rayleightaylor3d_mesh_00.vtk";
        VTKFileProcessor p0;
        p0.read(ifname0);
        p0.interpolation_to_mesh(xLeft, xRight, yLeft, yRight, zLeft, zRight, xN, yN, zN, ofname0);
        cout<<"t=0 completed"<<endl;

	string ifname1 = "/home/xingyu/lp_output/lp_rayleightaylor3d_lw/out_fluid0000030.vtk";
	string ofname1 = "/home/xingyu/vtktest/rayleightaylor3d_mesh_30.vtk";
	VTKFileProcessor p1;
	p1.read(ifname1);
	p1.interpolation_to_mesh(xLeft, xRight, yLeft, yRight, zLeft, zRight, xN, yN, zN, ofname1);
	cout<<"t=30 completed"<<endl;

        string ifname2 = "/home/xingyu/lp_output/lp_rayleightaylor3d_lw/out_fluid0000060.vtk";
        string ofname2 = "/home/xingyu/vtktest/rayleightaylor3d_mesh_60.vtk";
        VTKFileProcessor p2;
        p2.read(ifname2);
        p2.interpolation_to_mesh(xLeft, xRight, yLeft, yRight, zLeft, zRight, xN, yN, zN, ofname2);
        cout<<"t=60 completed"<<endl;	

        string ifname3 = "/home/xingyu/lp_output/lp_rayleightaylor3d_lw/out_fluid0000085.vtk";
        string ofname3 = "/home/xingyu/vtktest/rayleightaylor3d_mesh_85.vtk";
        VTKFileProcessor p3;
        p3.read(ifname3);
        p3.interpolation_to_mesh(xLeft, xRight, yLeft, yRight, zLeft, zRight, xN, yN, zN, ofname3);
        cout<<"t=85 completed"<<endl;*/

	return 0;
}
