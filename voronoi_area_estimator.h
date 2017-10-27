#ifndef __VORONOI_AREA_ESTIMATOR__
#define __VORONOI_AREA_ESTIMATOR__

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullSet.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullPoint.h"

#include <cstdio>  
#include <ostream>
#include <stdexcept>
#include <vector>


using std::cerr;
using std::cin;
using std::cout;
using std::endl;

using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullQh;
using orgQhull::RboxPoints;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexList;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;


class VoronoiAreaEstimator{

private:
	int dimension;
	int number_of_points;
	const double* x_coord;
	const double* y_coord;
	const double* z_coord;
	const double* mass;
	double* area;
	

public:
	VoronoiAreaEstimator(int d, int N, const double* x, const double* y, const double* z, const double* m, double* a): dimension(d),number_of_points(N), x_coord(x), y_coord(y), z_coord(z), area(a) {};
	virtual ~VoronoiAreaEstimator() {};

	int ComputeVoronoiArea();

};
#endif
