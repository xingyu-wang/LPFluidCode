#include "voronoi_area_estimator.h"
#include <ctime>
#include "omp.h"
int VoronoiAreaEstimator::ComputeVoronoiArea() {
    QHULL_LIB_CHECK

    if(dimension == 3){
      std::cout<<"3D density estimator is not implemented yet."<<std::endl;
      return -1;
    }

    try{
//      double startTime;
//        startTime = omp_get_wtime();

      RboxPoints rbox;
      Qhull qhull;


std::ostringstream os;
os.precision(10);
os<<dimension<<" ";

os<<number_of_points<<" ";

for (int i=0;i<number_of_points;i++)
{
  os<<x_coord[i]<<" ";
  os<<y_coord[i]<<" ";
  area[i]=0.0;
}

std::istringstream is;

//std::cout << os.str() ;

is.str(os.str());

      rbox.appendPoints(is);
//        printf("Prepare data for Voronoi mesh takes %.16g seconds\n", omp_get_wtime() - startTime);
//        startTime = omp_get_wtime();

      qhull.runQhull(rbox, "v");
//        printf("Generate Voronoi mesh takes %.16g seconds\n", omp_get_wtime() - startTime);
//        startTime = omp_get_wtime();
      QhullFacetList facets=qhull.facetList();
  double radius2,length2,area_t;
  int count=0;
  for (QhullFacetList::iterator it = facets.begin(); it !=facets.end();++it)
  {
//    cout<<"Facet "<<count<<endl;
    count++;
    if (!(*it).isGood()) continue;      
    QhullFacet f = *it;

    QhullPoint center=f.getCenter();
      double *c_coord = center.coordinates();
//      std::cout<<"Center coordinates"<<std::endl;
//      std::cout<<c_coord[0]<<" "<<c_coord[1]<<" "<<c_coord[2]<<std::endl;
    QhullVertexSet vSet = f.vertices();
    int id[3];
    int count_vertex=0;
    for (QhullVertexSet::iterator vIt = vSet.begin(); (vIt != vSet.end()&&(count_vertex<3)); ++vIt)
    {
//      cout<<count_vertex<<endl;
      QhullVertex v = *vIt;
      QhullPoint p = v.point();
      id[count_vertex]=p.id();
//      cout<<count_vertex<<" "<<id[count_vertex]<<" "<<x_coord[id[count_vertex]]<<" "<<y_coord[id[count_vertex]]<<std::endl;
      count_vertex++;
    }

//    std::cout<<id[0]<<" "<<x_coord[id[0]]<<" "<<y_coord[id[0]]<<std::endl;
//    std::cout<<id[1]<<" "<<x_coord[id[1]]<<" "<<y_coord[id[1]]<<std::endl;
//    std::cout<<id[2]<<" "<<x_coord[id[2]]<<" "<<y_coord[id[2]]<<std::endl;

    radius2=(c_coord[0]-x_coord[id[0]])*(c_coord[0]-x_coord[id[0]])+(c_coord[1]-y_coord[id[0]])*(c_coord[1]-y_coord[id[0]]);

    length2=(x_coord[id[1]]-x_coord[id[0]])*(x_coord[id[1]]-x_coord[id[0]])+(y_coord[id[1]]-y_coord[id[0]])*(y_coord[id[1]]-y_coord[id[0]]);
    if (radius2>length2/4)
      area_t=sqrt((radius2-length2/4)*length2/4)/2;
    else
      area_t=0.0;

    area[id[0]]=area[id[0]]+area_t;
    area[id[1]]=area[id[1]]+area_t;

    length2=(x_coord[id[2]]-x_coord[id[0]])*(x_coord[id[2]]-x_coord[id[0]])+(y_coord[id[2]]-y_coord[id[0]])*(y_coord[id[2]]-y_coord[id[0]]);
    if (radius2>length2/4)
      area_t=sqrt((radius2-length2/4)*length2/4)/2;
    else
      area_t=0.0;

    area[id[0]]=area[id[0]]+area_t;
    area[id[2]]=area[id[2]]+area_t;

    length2=(x_coord[id[2]]-x_coord[id[1]])*(x_coord[id[2]]-x_coord[id[1]])+(y_coord[id[2]]-y_coord[id[1]])*(y_coord[id[2]]-y_coord[id[1]]);
    if (radius2>length2/4)
      area_t=sqrt((radius2-length2/4)*length2/4)/2;
    else
      area_t=0.0;

    area[id[1]]=area[id[1]]+area_t;
    area[id[2]]=area[id[2]]+area_t;
  }


//   for (int i=0;i<number_of_points;i++)
//    {
//      std::cout<<area[i]<<std::endl;
//    }
//        printf("Calculate Voronoi mesh area takes %.16g seconds\n", omp_get_wtime() - startTime);

      return 0;
    }catch(QhullError &e){
        cerr << e.what() << std::endl;
        return e.errorCode();
    }
  return 0;
}//voronoi_area_estimator

