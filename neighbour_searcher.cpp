/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Oct. 15 2014
 * 
 * \brief
 *      
 */


#include <algorithm>
#include <iostream>
#include "neighbour_searcher.h"
#include "initializer.h"
#include <omp.h>
#include <cassert>


////////////////////////////////////////////////////////////////////////////////////////
// Start of OctreeSearcher
////////////////////////////////////////////////////////////////////////////////////////

OctreeSearcher::OctreeSearcher(const Initializer& init) {
	m_iMaxParticleNum = init.getCapacity();
	m_iMaxNeighborNum = init.getMaxNeighbourNum()-init.getNumParticleWithinSearchRadius();
//	std::cout<<init.getMaxNeighbourNum()<<std::endl;
//	std::cout<<init.getNumParticleWithinSearchRadius()<<std::endl;

//	m_iMinNeighborNum = init.getNumParticleWithinSearchRadius()/2;
	m_iMinNeighborNum = init.getNumRow2ndOrder();
	int treeDepth = init.getTreeDepth(); 

	m_pOctree = new Octree(treeDepth, m_iMaxParticleNum);
	//m_pSearchResult = new SearchResult[m_iMaxNeighborNum];
	
#ifdef _OPENMP
	m_iTheadNum = omp_get_max_threads();
	
	m_pSearchResult = new SearchResult*[m_iTheadNum];
	for(int i=0 ; i<m_iTheadNum ; i++) {
		m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
	}

        m_pSearchResultTemp = new SearchResult*[m_iTheadNum];
        for(int i=0 ; i<m_iTheadNum ; i++) {
                m_pSearchResultTemp[i] = new SearchResult[m_iMaxNeighborNum];
        }
#else
	m_pSearchResulttemp = new SearchResult[m_iMaxNeighborNum];
#endif


	// print debug info
	std::cout<<"-------OctreeSearcher::OctreeSearcher()-------"<<std::endl;
	std::cout<<"m_iMaxParticleNum = "<<m_iMaxParticleNum<<std::endl;
	std::cout<<"m_iMaxNeighborNum = "<<m_iMaxNeighborNum<<std::endl;
	std::cout<<"m_iMinNeighborNum = "<<m_iMinNeighborNum<<std::endl;
	std::cout<<"treeDepth = "<<treeDepth<<std::endl;
	std::cout<<"----------------------------------------------"<<std::endl;
}



OctreeSearcher::OctreeSearcher(size_t maxParticleNum, size_t maxNeighbourNum, int treedepth)
: m_iMaxParticleNum(maxParticleNum), m_iMaxNeighborNum(maxNeighbourNum)
{
  m_pOctree = new Octree(treedepth, m_iMaxParticleNum);

#ifdef _OPENMP
  m_iTheadNum = omp_get_max_threads();
  /*
  m_pSearchResult.resize(m_iTheadNum);
  for(int i=0 ; i<m_iTheadNum ; i++) {
    //m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
    m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
  }
  */
  m_pSearchResult = new SearchResult*[m_iTheadNum];
  for(int i=0 ; i<m_iTheadNum ; i++) {
    m_pSearchResult[i] = new SearchResult[m_iMaxNeighborNum];
  }
#else
  m_pSearchResult = new SearchResult[m_iMaxNeighborNum];
#endif
}



int OctreeSearcher::buildSearchStructure(const double* x, const double* y, const double* z, size_t begin, size_t numParticles) {
  return m_pOctree->buildOctree(x+begin, y+begin, z+begin, numParticles);
}

int OctreeSearcher::buildSearchStructure(const double* x, const double* y, const double* z, const double* m, size_t begin, size_t numParticles) {
  return m_pOctree->buildOctree(x+begin, y+begin, z+begin, m+begin, numParticles);
}

int OctreeSearcher::buildSearchStructure(const double* x, const double* y, const double* z, const double* m, const double* vv, size_t begin, size_t numParticles) {
  return m_pOctree->buildOctree(x+begin, y+begin, z+begin, m+begin, vv+begin, numParticles);
}


/*
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double radius, int* result, size_t& result_length) {
  return m_pOctree->searchNeighbor(x, y, z, radius, result, result_length);
}
*/


#ifdef _OPENMP
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int tid, int index) {

	result_length=0;
	double radius=initial_radius;
	double min_radius=-1.0, max_radius=-1.0;
	int lower_bound=0,upper_bound=0;

	while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
	{
		if(upper_bound*lower_bound==1)
			radius=(min_radius+max_radius)/2.0;
		if(upper_bound>lower_bound)
			radius=max_radius/2.0;
		if(upper_bound<lower_bound)
			radius=min_radius*2.0;
		result_length=m_iMaxNeighborNum;
	  m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult[tid], result_length);
		if(result_length>m_iMaxNeighborNum)
		{
			max_radius=radius;
			upper_bound=1;
		}
		if(result_length<m_iMinNeighborNum)
		{
			min_radius=radius;
			lower_bound=1;
		}
  }

	sort(m_pSearchResult[tid], m_pSearchResult[tid]+result_length);

  if(index > -1) {
    size_t i=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[tid][i].index) {
        i++;
        break;
      } else {
        result[i] = m_pSearchResult[tid][i].index;
        distance[i] = m_pSearchResult[tid][i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[tid][i].index;
      distance[i-1] = m_pSearchResult[tid][i].distance;
    }

    result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[tid][i].index;
      distance[i] = m_pSearchResult[tid][i].distance;
    }
  }
//	if(radius>initial_radius)
//		std::cout<<"Search radius enlarged for "<<index<<"th particle. Initial search radius = "<<initial_radius<<", search radius = "<<radius<<"."<<std::endl;
//	if(radius<initial_radius)
//		std::cout<<"Search radius shrunk for "<<index<<"th particle. Initial search radius = "<<initial_radius<<", search radius = "<<radius<<"."<<std::endl;

  return 0;
}

#else
int OctreeSearcher::searchNeighbour(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int index) {

	result_length=0;
	double radius=initial_radius;
	double min_radius=-1.0, max_radius=-1.0;
	int lower_bound=0,upper_bound=0;

	while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
	{
		if(upper_bound*lower_bound==1)
			radius=(min_radius+max_radius)/2.0;
		if(upper_bound>lower_bound)
			radius=max_radius/2.0;
		if(upper_bound<lower_bound)
			radius=min_radius*2.0;
		result_length=m_iMaxNeighborNum;
  	m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult, result_length);
		if(result_length>m_iMaxNeighborNum)
		{
			max_radius=radius;
			upper_bound=1;
		}
		if(result_length<m_iMinNeighborNum)
		{
			min_radius=radius;
			lower_bound=1;
		}
  }

  sort(m_pSearchResult, m_pSearchResult+result_length);

  if(index > -1) {
    size_t i=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[i].index) {
        i++;
        break;
      } else {
        result[i] = m_pSearchResult[i].index;
        distance[i] = m_pSearchResult[i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[i].index;
      distance[i-1] = m_pSearchResult[i].distance;
    }

    result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[i].index;
      distance[i] = m_pSearchResult[i].distance;
    }
  }
//	if(radius>initial_radius)
//		std::cout<<"Search radius enlarged for "<<index<<"th particle. Initial search radius = "<<initial_radius<<", search radius = "<<radius<<"."<<std::endl;
//	if(radius<initial_radius)
//		std::cout<<"Search radius shrunk for "<<index<<"th particle. Initial search radius = "<<initial_radius<<", search radius = "<<radius<<"."<<std::endl;


  return 0;
}

#endif


  void OctreeSearcher::setMaxParticleNum(size_t maxParticleNum) {
	m_iMaxParticleNum = maxParticleNum;
    m_pOctree->setMaxParticleNum(m_iMaxParticleNum);
  }


#ifdef _OPENMP
int OctreeSearcher::searchNeighbourQuadrant(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int tid, int index) {
//Search neighbour in a sphere, should give good neighbourhood for uniform particle distribution
	double fixed_radius=initial_radius;
	{
		result_length=0;
		double min_radius=-1.0, max_radius=-1.0;
		int lower_bound=0,upper_bound=0;
//Search radius is adaptively adjusted if too many / too few neighbours are found
		while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
		{
//			if(index%100000==0)
//				cout<<index<<" "<<result_length<<endl;
			if(upper_bound*lower_bound==1)
				fixed_radius=(min_radius+max_radius)/2.0;
			if(upper_bound>lower_bound)
				fixed_radius=max_radius/2.0;
			if(upper_bound<lower_bound)
				fixed_radius=min_radius*2.0;
			result_length=m_iMaxNeighborNum;
			m_pOctree->searchNeighbor(x, y, z, fixed_radius, m_pSearchResultTemp[tid], result_length);
			if(result_length>m_iMaxNeighborNum)
			{
				max_radius=fixed_radius;
				upper_bound=1;
			}
			if(result_length<m_iMinNeighborNum)
			{
				min_radius=fixed_radius;
				lower_bound=1;
			}
		}
	}
	
	size_t result_lengthtemp=result_length;
	int count=0;
	size_t temp_length=0;
	int number_of_quadrant=4;
	double d_angle=2.0*M_PI/number_of_quadrant;
	double angle=45.0/180.0*M_PI;
	int radius_multiplier=16;
	for(int quadrant_index=0;quadrant_index<number_of_quadrant;quadrant_index++)
	{
		double radius=fixed_radius;
		double min_radius=-1.0, max_radius=-1.0;
		int lower_bound=0,upper_bound=0;
		temp_length=0;
//Select the neighbours in the quadrant from the spherical neighbourhood
		for(size_t indextemp=0;indextemp<result_lengthtemp;indextemp++){
			if(m_pSearchResultTemp[tid][indextemp].angle>=angle && m_pSearchResultTemp[tid][indextemp].angle<(angle+d_angle))
				m_pSearchResult[tid][count+temp_length++]=m_pSearchResultTemp[tid][indextemp];
			if(temp_length>m_iMaxNeighborNum/number_of_quadrant)
				break;
		}
//If the particle distribution is non-uniform and the spherical neighbourhood is not good enough, search for neighbours in the quadrant individually using octree
//Search radius is adaptively adjusted if too many / too few neighbours are found
		while((temp_length>(m_iMaxNeighborNum/number_of_quadrant))||((temp_length<(m_iMinNeighborNum/number_of_quadrant))&&(radius<(radius_multiplier*fixed_radius))))
		{
//                        if(index%100000==0)
//                                cout<<index<<" "<<temp_length<<" "<<angle<<endl;
			if(upper_bound*lower_bound==1)
				radius=(min_radius+max_radius)/2.0;
			if(upper_bound>lower_bound)
				radius=max_radius/2.0;
			if(upper_bound<lower_bound)
				radius=min_radius*2.0;
			temp_length=m_iMaxNeighborNum/number_of_quadrant;
			m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult[tid], temp_length, 1*count, 1.0*angle, angle+d_angle);
			if(temp_length>(m_iMaxNeighborNum/number_of_quadrant))
			{
				max_radius=radius;
				upper_bound=1;
			}
			if(temp_length<(m_iMinNeighborNum/number_of_quadrant))
			{
				min_radius=radius;
				lower_bound=1;
			}
		}
		angle=angle+d_angle;
		count=count+temp_length;
	}
	result_length=count;
//The output neighbour list is sorted by the distance

	sort(m_pSearchResult[tid], m_pSearchResult[tid]+result_length);
//The particle is not a neighbour of itself, so remove it from the neighbour list
  if(index > -1) {
    size_t i=0;
		int flag=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[tid][i].index) {
        i++;
	flag=1;
        break;
      } else {
        result[i] = m_pSearchResult[tid][i].index;
        distance[i] = m_pSearchResult[tid][i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[tid][i].index;
      distance[i-1] = m_pSearchResult[tid][i].distance;
    }
	if(flag==1)
    		result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[tid][i].index;
      distance[i] = m_pSearchResult[tid][i].distance;
    }
  }

  return 0;
}

#else
int OctreeSearcher::searchNeighbourQuadrant(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int index) {

	double fixed_radius=initial_radius;
	{
		result_length=0;
		double min_radius=-1.0, max_radius=-1.0;
		int lower_bound=0,upper_bound=0;

		while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
		{
			if(upper_bound*lower_bound==1)
				fixed_radius=(min_radius+max_radius)/2.0;
			if(upper_bound>lower_bound)
				fixed_radius=max_radius/2.0;
			if(upper_bound<lower_bound)
				fixed_radius=min_radius*2.0;
			result_length=m_iMaxNeighborNum;
			m_pOctree->searchNeighbor(x, y, z, fixed_radius, m_pSearchResult, result_length);
			if(result_length>m_iMaxNeighborNum)
			{
				max_radius=fixed_radius;
				upper_bound=1;
			}
			if(result_length<m_iMinNeighborNum)
			{
				min_radius=fixed_radius;
				lower_bound=1;
			}
		}
	}


	int count=0;
	size_t temp_length=0;
	int number_of_quadrant=4;
	double d_angle=2.0*M_PI/number_of_quadrant;
	double angle=45.0/180.0*M_PI;
	int radius_multiplier=16;
	for(int quadrant_index=0;quadrant_index<number_of_quadrant;quadrant_index++)
	{
		double radius=fixed_radius;
		double min_radius=-1.0, max_radius=-1.0;
		int lower_bound=0,upper_bound=0;
		temp_length=0;

		while((temp_length>(m_iMaxNeighborNum/number_of_quadrant))||((temp_length<(m_iMinNeighborNum/number_of_quadrant))&&(radius<(radius_multiplier*fixed_radius))))
		{
			if(upper_bound*lower_bound==1)
				radius=(min_radius+max_radius)/2.0;
			if(upper_bound>lower_bound)
				radius=max_radius/2.0;
			if(upper_bound<lower_bound)
				radius=min_radius*2.0;
			temp_length=m_iMaxNeighborNum/number_of_quadrant;
			m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult, temp_length, 1*count, 1.0*angle, angle+d_angle);
			if(temp_length>(m_iMaxNeighborNum/number_of_quadrant))
			{
				max_radius=radius;
				upper_bound=1;
			}
			if(temp_length<(m_iMinNeighborNum/number_of_quadrant))
			{
				min_radius=radius;
				lower_bound=1;
			}
		}
		angle=angle+d_angle;
		count=count+temp_length;
	}
	result_length=count;

  sort(m_pSearchResult, m_pSearchResult+result_length);

  if(index > -1) {
    size_t i=0;
		int flag=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[i].index) {
        i++;
				flag=1;
        break;
      } else {
        result[i] = m_pSearchResult[i].index;
        distance[i] = m_pSearchResult[i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[i].index;
      distance[i-1] = m_pSearchResult[i].distance;
    }
		if(flag==1)
    	result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[i].index;
      distance[i] = m_pSearchResult[i].distance;
    }
  }


  return 0;
}

#endif

#ifdef _OPENMP
int OctreeSearcher::searchNeighbourDirection(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int tid, int index) {

        double fixed_radius=initial_radius;
        {
                result_length=0;
                double min_radius=-1.0, max_radius=-1.0;
                int lower_bound=0,upper_bound=0;

                while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
                {
                        if(upper_bound*lower_bound==1)
                                fixed_radius=(min_radius+max_radius)/2.0;
                        if(upper_bound>lower_bound)
                                fixed_radius=max_radius/2.0;
                        if(upper_bound<lower_bound)
                                fixed_radius=min_radius*2.0;
                        result_length=m_iMaxNeighborNum;
                        m_pOctree->searchNeighbor(x, y, z, fixed_radius, m_pSearchResult[tid], result_length);
                        if(result_length>m_iMaxNeighborNum)
                        {
                                max_radius=fixed_radius;
                                upper_bound=1;
                        }
                        if(result_length<m_iMinNeighborNum)
                        {
                                min_radius=fixed_radius;
                                lower_bound=1;
                        }
                }
        }

        int count=0;
        size_t temp_length=0;
	int number_of_dir=6;
        int radius_multiplier=16;
	int dir;
        for(int dir_index=0;dir_index<number_of_dir;dir_index++)
        {
		dir=(2*(dir_index%2)-1)*(dir_index/2+1);
                double radius=fixed_radius;
                double min_radius=-1.0, max_radius=-1.0;
                int lower_bound=0,upper_bound=0;
                temp_length=0;

                while((temp_length>(m_iMaxNeighborNum/number_of_dir))||((temp_length<(m_iMinNeighborNum/number_of_dir))&&(radius<(radius_multiplier*fixed_radius))))
                {
                        if(upper_bound*lower_bound==1)
                                radius=(min_radius+max_radius)/2.0;
                        if(upper_bound>lower_bound)
                                radius=max_radius/2.0;
                        if(upper_bound<lower_bound)
                                radius=min_radius*2.0;
                        temp_length=m_iMaxNeighborNum/number_of_dir;
                        m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult[tid], temp_length, count, dir);
                        if(temp_length>(m_iMaxNeighborNum/number_of_dir))
                        {
                                max_radius=radius;
                                upper_bound=1;
                        }
                        if(temp_length<(m_iMinNeighborNum/number_of_dir))
                        {
                                min_radius=radius;
                                lower_bound=1;
                        }
                }
                count=count+temp_length;
        }
        result_length=count;

        sort(m_pSearchResult[tid], m_pSearchResult[tid]+result_length);

/*  if(index > -1) {
    size_t i=0;
                int flag=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[tid][i].index) {
        i++;
                                flag=1;
        break;
      } else {
        result[i] = m_pSearchResult[tid][i].index;
        distance[i] = m_pSearchResult[tid][i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[tid][i].index;
      distance[i-1] = m_pSearchResult[tid][i].distance;
    }
                if(flag==1)
        result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[tid][i].index;
      distance[i] = m_pSearchResult[tid][i].distance;
    }
  }*/

  size_t i=0;
  size_t j=0;
  for(i=0; i<result_length; i++) {
      if(index == m_pSearchResult[tid][i].index) {
//	result_length--;
        continue;
      } else {
        result[j] = m_pSearchResult[tid][i].index;
        distance[j] = m_pSearchResult[tid][i].distance;
        j++;
      }
    }
  result_length=j;

  return 0;
}

#else
int OctreeSearcher::searchNeighbourDirection(const double x, const double y, const double z, const double initial_radius, int* result, double* distance, size_t& result_length, int index) {

        double fixed_radius=initial_radius;
        {
                result_length=0;
                double min_radius=-1.0, max_radius=-1.0;
                int lower_bound=0,upper_bound=0;

                while((result_length>m_iMaxNeighborNum)||(result_length<m_iMinNeighborNum))
                {
                        if(upper_bound*lower_bound==1)
                                fixed_radius=(min_radius+max_radius)/2.0;
                        if(upper_bound>lower_bound)
                                fixed_radius=max_radius/2.0;
                        if(upper_bound<lower_bound)
                                fixed_radius=min_radius*2.0;
                        result_length=m_iMaxNeighborNum;
                        m_pOctree->searchNeighbor(x, y, z, fixed_radius, m_pSearchResult, result_length);
                        if(result_length>m_iMaxNeighborNum)
                        {
                                max_radius=fixed_radius;
                                upper_bound=1;
                        }
                        if(result_length<m_iMinNeighborNum)
                        {
                                min_radius=fixed_radius;
                                lower_bound=1;
                        }
                }
        }

        int count=0;
        size_t temp_length=0;
        int number_of_dir=6;
        int radius_multiplier=16;
	int dir;
        for(int dir_index=0;dir_index<number_of_dir;dir_index++)
        {
                dir=(2*(dir_index%2)-1)*(dir_index/2+1);
                double radius=fixed_radius;
                double min_radius=-1.0, max_radius=-1.0;
                int lower_bound=0,upper_bound=0;
                temp_length=0;

                while((temp_length>(m_iMaxNeighborNum/number_of_dir))||((temp_length<(m_iMinNeighborNum/number_of_dir))&&(radius<(radius_multiplier*fixed_radius))))
                {
                        if(upper_bound*lower_bound==1)
                                radius=(min_radius+max_radius)/2.0;
                        if(upper_bound>lower_bound)
                                radius=max_radius/2.0;
                        if(upper_bound<lower_bound)
                                radius=min_radius*2.0;
                        temp_length=m_iMaxNeighborNum/number_of_dir;
                        m_pOctree->searchNeighbor(x, y, z, radius, m_pSearchResult, temp_length, count, dir);
                        if(temp_length>(m_iMaxNeighborNum/number_of_dir))
                        {
                                max_radius=radius;
                                upper_bound=1;
                        }
                        if(temp_length<(m_iMinNeighborNum/number_of_dir))
                        {
                                min_radius=radius;
                                lower_bound=1;
                        }
                }
                count=count+temp_length;
        }
        result_length=count;

  sort(m_pSearchResult, m_pSearchResult+result_length);

/*  if(index > -1) {
    size_t i=0;
    int flag=0;
    for(i=0 ; i<result_length ; i++) {
      if(index == m_pSearchResult[i].index) {
        i++;
                                flag=1;
        break;
      } else {
        result[i] = m_pSearchResult[i].index;
        distance[i] = m_pSearchResult[i].distance;
      }
    }
    for(; i<result_length ; i++) {
      result[i-1] = m_pSearchResult[i].index;
      distance[i-1] = m_pSearchResult[i].distance;
    }
                if(flag==1)
        result_length--;

  } else {
    for(size_t i=0 ; i<result_length ; i++) {
      result[i] = m_pSearchResult[i].index;
      distance[i] = m_pSearchResult[i].distance;
    }
  }*/

  size_t i=0;
  size_t j=0;
  for(i=0; i<result_length; i++) {
      if(index == m_pSearchResult[i].index) {
//	result_length--;
	continue;
      } else {
        result[j] = m_pSearchResult[i].index;
        distance[j] = m_pSearchResult[i].distance;
	j++;
      }
    }
  result_length=j;

  return 0;
}

#endif























int OctreeSearcher::densityEstimator(const int index, const double x, const double y, const double z, const double radius, double* count_density, double dir_x, double dir_y, double dir_z) {

        return m_pOctree->densityEstimator(index, x, y, z, radius, count_density, dir_x, dir_y, dir_z);
}

int OctreeSearcher::densityEstimator(const double x, const double y, const double z, const double radius, double* count_density, double dir_x, double dir_y, double dir_z) {

        return m_pOctree->densityEstimator(x, y, z, radius, count_density, dir_x, dir_y, dir_z);
}

int OctreeSearcher::VoronoiDensityEstimator(const int index, const double x, const double y, const double z, const double radius, double* density) {

	return m_pOctree->VoronoiDensityEstimator(index, x, y, z, radius, density);
}

int OctreeSearcher::densityEstimator(const int index, const double x, const double y, const double z, const double radius, double* count_density, const double* volume, double* vmin, double* vmax) {

        return m_pOctree->densityEstimator(index, x, y, z, radius, count_density, volume, vmin, vmax);
}

////////////////////////////////////////////////////////////////////////////////////////
// End of OctreeSearcher
////////////////////////////////////////////////////////////////////////////////////////

  



  
