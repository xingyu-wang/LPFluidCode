/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \author Gaurish Telang <gaurish108@gaurish108> 
 * \date Mon Jan. 30 2012
 * 
 * \brief Implementations of Octree class.
 *  
 */

#include "octree.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <numeric> 
#include <fstream>
#include <cmath>
#include <omp.h>
#include <new>
#include <cassert>

//#define OUTPUT_TIME



using namespace std;

bool compareSearchResultAngle(const SearchResult& sr1, const SearchResult& sr2) {
  return (sr1.angle<sr2.angle);
}

inline bool is_node_intersect_search_region(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const double& search_x, const double& search_y, const double& search_z, const double& radius) {
  //The following calculation has been fashioned such that:
  //if final value of squared_dmin == 0 then the point lies inside the node. 
  //if final value of squared_dmin !=0  then the point lies outside the node, AND 
  //             tells the SQUARE of the minimum distance of the point to the points on the node boundary(surface).


  double squared_dmin = 0.0;
  double temp; //Used as a temporary variable to store some intermediate results.
 
  //Process the x cooridinates
  if( search_x < min_x ) { 
    temp = search_x - min_x; 
    squared_dmin += temp*temp;
  } else if( search_x > max_x )	{ 
    temp         = search_x - max_x ; 
    squared_dmin += temp*temp                  ;
  }   


  //Process the Y-coorindtaes
  if( search_y < min_y ) { 
    temp = search_y - min_y ; 
    squared_dmin += temp*temp;
  } else if( search_y > max_y ) { 
    temp = search_y - max_y ; 
    squared_dmin += temp*temp;
  }   

  //Process the Z-coorindtaes
  if( search_z < min_z ) 
  { 
    temp          = search_z - min_z; 
    squared_dmin += temp*temp;
  } else if( search_z > max_z ) { 
    temp          = search_z - max_z; 
    squared_dmin += temp*temp;
  }   

  if (squared_dmin <= radius*radius) return true;
  else return false;

}

/*
inline bool is_node_intersect_search_region(const node& treenode, const double& search_x, const double& search_y, const double& search_z, const double& radius) {
  return is_node_intersect_search_region(treenode.vmin[0], treenode.vmax[0], treenode.vmin[1], treenode.vmax[1], treenode.vmin[2], treenode.vmax[2], search_x, search_y, search_z, radius);
} 
*/





Octree::Octree(int treedepth, size_t maxParticleNum)
: m_iMaxParticleNum(maxParticleNum), m_iMaxDepth(treedepth)
  , m_vNodeKey(0), m_vDepth(0), m_vFirstParticleIndex(0), m_vNumberOfContainedParticles(0), m_vFirstChildIndex(0), m_vNumberOfChildren(0)
  , m_vLowerLimitOfX(0), m_vLowerLimitOfY(0), m_vLowerLimitOfZ(0), m_vUpperLimitOfX(0),  m_vUpperLimitOfY(0),  m_vUpperLimitOfZ(0)
{

#ifdef _OPENMP
  //m_vParticleKeyIndex.resize(m_iMaxParticleNum);
  m_vParticleKeyIndex = new KeyIndex [ m_iMaxParticleNum ];
#else
  m_vParticleKeyIndex = new KeyIndex [ m_iMaxParticleNum ];
#endif

  
  //initialize(numParticles);
}


void Octree::setMaxParticleNum(size_t maxParticleNum) {
	m_iMaxParticleNum = maxParticleNum;
	KeyIndex* tmp = m_vParticleKeyIndex;
	try {
		m_vParticleKeyIndex = new KeyIndex [m_iMaxParticleNum];
	}
	catch(std::bad_alloc& ba) {
		std::cerr<<"std::bad_alloc caught during data array augmentation: "<<ba.what()<<std::endl;
		assert(false);
	}
	delete[] tmp;
}

Octree::~Octree() {
  //delete[] m_vCoordX;
  //delete[] m_vCoordY;
  //delete[] m_vCoordZ;

  //delete m_vParticleNumList;

  delete[] m_vParticleKeyIndex;


  if(m_vNodeKey != 0) {
    clearOctreeStructure();
  }
  



}

/*
int Octree::initialize(size_t numParticles) {

  if(numParticles > m_iMaxParticleNum) {
    return 1; // error // TODO:
  }

  m_iTotalNumberOfParticles = numParticles;

  for(size_t i=0 ; i<m_iTotalNumberOfParticles ; i++) {
    m_vParticleKeyIndex[i].index = i;
  }

  return 0;
}
*/


//This might need to be uncommented later. One of a few possible constructors.
// Octree::Octree(const double *xp, const double *yp, const double *zp, int treedepth, int numOfParticles)
// : m_iTotalNumberOfParticles(numOfParticles), m_iMaxDepth(treedepth)
// {
//   buildOctree(xp, yp, zp);
// }

int Octree::buildOctree(const double* x, const double* y, const double* z, const double* m, size_t numParticles) {
	buildOctree(x,y,z,numParticles);
	mass = m;
	return 0;
}

int Octree::buildOctree(const double* x, const double* y, const double* z, const double* m, const double* vv, size_t numParticles) {
        buildOctree(x,y,z,numParticles);
        mass = m;
	VolumeVoronoi=vv;
	return 0;
}

int Octree::buildOctree(const double* x, const double* y, const double* z, size_t numParticles) {

  m_vCoordX = x;
  m_vCoordY = y;
  m_vCoordZ = z;
  m_iTotalNumberOfParticles = numParticles;

  m_dBoundingBox_min_x = *std::min_element(m_vCoordX, m_vCoordX + m_iTotalNumberOfParticles);
  m_dBoundingBox_max_x = *std::max_element(m_vCoordX, m_vCoordX + m_iTotalNumberOfParticles);
  m_dBoundingBox_min_y = *std::min_element(m_vCoordY, m_vCoordY + m_iTotalNumberOfParticles);
  m_dBoundingBox_max_y = *std::max_element(m_vCoordY, m_vCoordY + m_iTotalNumberOfParticles);
  m_dBoundingBox_min_z = *std::min_element(m_vCoordZ, m_vCoordZ + m_iTotalNumberOfParticles);
  m_dBoundingBox_max_z = *std::max_element(m_vCoordZ, m_vCoordZ + m_iTotalNumberOfParticles);

  //Create an array which will hold the morton key and index of each particle.
  //After we fill this array, it will be sorted by the "key" field in the "KeyIndex" struct.
  //The operator < has been overloaded for this purpose
  


  //int index = 0;

#ifdef OUTPUT_TIME
        double time1;
  time1 = omp_get_wtime();
#endif
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=0 ; i<m_iTotalNumberOfParticles ; i++) {
    m_vParticleKeyIndex[i].index = i;
    m_vParticleKeyIndex[i].key = computeKey(m_vCoordX[i], m_vCoordY[i], m_vCoordZ[i]);
    //(m_vParticleKeyIndex.at(i)).key = computeKey(m_vCoordX[i], m_vCoordY[i], m_vCoordZ[i]);
  }
#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "Key Computation Time : " << time1 << std::endl;
#endif
  //Now sort the key-index array
#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
#ifdef _OPENMP
  //thrust::sort(m_vParticleKeyIndex.begin(), m_vParticleKeyIndex.begin() + m_iTotalNumberOfParticles);
  std::sort(m_vParticleKeyIndex, m_vParticleKeyIndex + m_iTotalNumberOfParticles);
#else
  std::sort(m_vParticleKeyIndex, m_vParticleKeyIndex + m_iTotalNumberOfParticles);
#endif
#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "Key Sorting Time : " << time1 << std::endl;
#endif  

  //Now that we have the sorted key-index array we can begin the tree construction level by level
  //For that we need to calculate three helper vectors. bitmasks, baseaddresses and numnodes. 
  //numnodes[i]     = number of nodes at level i of the octree. zero level correposnds to root nodes. hence numnodes[0] = 1
  //basaddresses[i] = position in the octree arrays of first node at level i
  //bitmasks[i]     = Bit-Mask for level i. used in the computation of the numnodes vector 
  std::vector<int> bitmasks(m_iMaxDepth+1);
  std::vector<int> baseaddresses(m_iMaxDepth+1);
  std::vector<int> numnodes(m_iMaxDepth+1);

  //CALCULATE THE BITMASKS VECTOR
  bitmasks[m_iMaxDepth] = (1 << (3 * m_iMaxDepth)) - 1 ;//Deepest level bitmask 
  bitmasks[    0      ] = 0                       ;//Root    level bitmask 
  //Compute remaining level bitmasks
  for (int k = m_iMaxDepth -1 ; k >= 1 ; --k) {
    int shiftbits = 3 * ( m_iMaxDepth - k ) ;
    bitmasks[k] = ( bitmasks[ m_iMaxDepth ] >> shiftbits ) << shiftbits;  
  } 
#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
  //COMPUTE THE NUMNODES VECTOR
  numnodes[ 0 ] =  1 ;//Root level has only 1 node i.e. the root node itself.
  for (int k = 1; k < m_iMaxDepth + 1; ++k) {
    int count=1;
    
    #ifdef _OPENMP
	#pragma omp parallel for default(shared) reduction(+:count)
    #endif
	for ( int i = 1; i < m_iTotalNumberOfParticles ; ++i) {
      if ( (m_vParticleKeyIndex[i].key & bitmasks[k]) != (m_vParticleKeyIndex[i-1].key & bitmasks[k] )) {
	      ++count;
	    }
    }

    numnodes[k]=count;   
  }
#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "NUMNODES Computation Time : " << time1 << std::endl;
#endif
  //COMPUTE THE BASEADDRESSES VECTOR.
  baseaddresses[0]=0; 
  baseaddresses[1]=1;
  for (unsigned int k = 2; k < baseaddresses.size(); ++k) {
   baseaddresses[k]=baseaddresses[k-1]+numnodes[k-1];
  }   

  //With the bitmasks, numnodes and baseaddresses vectors we can now begin the tree computation.
  
  //This is the length of the arrays used in descirbing octrees. 
  m_iTreeLength = std::accumulate(numnodes.begin(),numnodes.end(),0);

  if(m_vNodeKey != 0) {
    clearOctreeStructure();
  }
  //Now we allocate the arrays using the new operator and m_iTreeLength
  m_vNodeKey			= new int[m_iTreeLength];
  m_vDepth			= new int[m_iTreeLength];
  m_vFirstParticleIndex		= new int[m_iTreeLength];
  m_vNumberOfContainedParticles = new int[m_iTreeLength];
  m_vFirstChildIndex		= new int[m_iTreeLength];
  m_vNumberOfChildren		= new int[m_iTreeLength];  

  m_vLowerLimitOfX		= new double[m_iTreeLength];     
  m_vUpperLimitOfX		= new double[m_iTreeLength]; 

  m_vLowerLimitOfY		= new double[m_iTreeLength]; 
  m_vUpperLimitOfY		= new double[m_iTreeLength]; 

  m_vLowerLimitOfZ		= new double[m_iTreeLength]; 
  m_vUpperLimitOfZ		= new double[m_iTreeLength]; 
#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
  //Start building the last level
  
  //Construct the deepest level of the tree first.
  m_vNodeKey[baseaddresses[m_iMaxDepth]]   = m_vParticleKeyIndex[0].key;
  m_vDepth[baseaddresses[m_iMaxDepth]]     = m_iMaxDepth;

  m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]]	= 0;
 //tree[baseaddressed[m_iMaxDepth]].pnum calculated in the for loop below.

  m_vNumberOfChildren[baseaddresses[m_iMaxDepth]]	= 0;
  m_vFirstChildIndex[baseaddresses[m_iMaxDepth]]	= (-1);
  

  int j   = 1; //This counter changes on encountering a new key.  
  int fix = 0; //Used in the pnum calculation. Records the index where the key change last took place.  
               //pnum = i-fix where i and fix are positions where successive key changes take place.  

  for(int i = 1; i < m_iTotalNumberOfParticles ; ++i) {
    if(m_vParticleKeyIndex[i].key!=m_vParticleKeyIndex[i-1].key)	{
      m_vNodeKey[baseaddresses[m_iMaxDepth]+j] = m_vParticleKeyIndex[i].key;  
      m_vDepth[baseaddresses[m_iMaxDepth]+j] = m_iMaxDepth; 
      m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]+j]	= i; 
      m_vNumberOfChildren[baseaddresses[m_iMaxDepth]+j] = 0; //No child nodes. Deepest level.
      m_vFirstChildIndex[baseaddresses[m_iMaxDepth]+j]   		=(-1); //Special sentinel value. No children since leaf node.    
      m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]	= i - fix; //Justified above!
      ++j;  
      fix	= i; //change origin of the measure tape. 
    }//END computation on encountering a new label
  }
  m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]  = m_iTotalNumberOfParticles - fix;//This cannot be tackled within above for loop.
  //END DEEPEST LEVEL NODE CONSTRUCTION.
#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "Leaf Construction Time : " << time1 << std::endl;
  time1 = omp_get_wtime();
#endif
  //Now that the deepest level has been constructed lets us start filling in the rest of the octree nodes.
  //Fill the other levels.
  for (int level=m_iMaxDepth-1; level>=0; --level) {
    j=1;
    fix=0;    //used in cnum calculation. Records the index where key change LAST took place. 
              //cnum=i-fix. 'i' and 'fix' are positions where successive MASKED key changes take place
          
    //LEVEL CONSTRUCTION BEGINS. WE MUST SCAN LEVEL+1 SECTION TO INSERT APPRoPRIATE VALUES IN THE LEVEL SECTION.
    m_vNodeKey[baseaddresses[level]]			= ( m_vNodeKey[baseaddresses[level+1]] & bitmasks[level])     ;
    m_vDepth[baseaddresses[level]]			= level;
    m_vFirstParticleIndex[baseaddresses[level]]	= 0                         ;
    m_vFirstChildIndex[baseaddresses[level]]		= baseaddresses[level+1]     ;

    int sum=m_vNumberOfContainedParticles[baseaddresses[level+1]];//This variable used in pnum calculations.
    m_vNumberOfContainedParticles[baseaddresses[level]]	= sum  ; 

    for (int i = 1; i < numnodes[level+1] ; ++i) {
      if( ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] ) != ( m_vNodeKey[baseaddresses[level+1] + i-1] & bitmasks[level] )  )//Relate i and i-1
      {
         m_vNodeKey[baseaddresses[level]+j]			= ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] )  ; //Give it the bitmasked key.  
         m_vDepth[baseaddresses[level]+j]				= level                       ;//Assign the depth
         m_vFirstParticleIndex[baseaddresses[level]+j]		= m_vFirstParticleIndex[baseaddresses[level+1]+i]                              ;
         m_vNumberOfContainedParticles[baseaddresses[level]+j-1]	= sum;//Calculate the difference where succesive changes take place in keys masked with a level mask.
         m_vFirstChildIndex[baseaddresses[level]+j]		= baseaddresses[level+1]+i      ;//Since we encountered a new masked level key. we record the place where this is encountered/
         m_vNumberOfChildren[baseaddresses[level]+j-1]		= i-fix                           ;//Will have to keep some kind of running sum
       
         ++j  ;
         fix=i;
         sum = m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//initializing sum to pnum of the new key encountered.              
       }

       else sum += m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//This is executed only if baseaddresses[level+1]+i and baseaddresses[level+1]+i-1 share the same parent . 
                                                  //Hence we increment the sum counter which is keeping track of the number of particles within the parent node.
    }//end inner for 

    m_vNumberOfChildren[baseaddresses[level]+j-1]                = numnodes[level+1]-fix   ;
    m_vNumberOfContainedParticles[baseaddresses[level]+j-1]      = sum                     ;
  }//end outer for

//END LEVEL CONSTRUCTION.
#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "Non-leaf Construction Time : " << time1 << std::endl;
#endif
  //Finally we fill in the vertex information required by all the nodes.
  //This is a separate computatin because of its length.
  //Loop over all the tree vector from left to right
  //for center & vertex computation. The only information necessary 
  //for the entire computation here is the DEPTH and the KEY of all the octree nodes.
  //These have already been computed above.

        //In the previous code, the vertex location of each node is computed as the limit of bounding box of the domain plus a series of displacement. Each level higher than the current level offers a displacement for each node, so the time complexity is O(m_iMaxDepth*m_iTreeLength). For some problem, the time to compute bounding box can take ~50% of the total time to build octree. In this code, the vertex location of each node is computed as the vertex location of its parent plus only one displacement given by the current level, so the new time complexity is O(m_iTreeLength). This modification reduces the time to compute vertex location by ~80%, compared to the previous code.
#ifdef OUTPUT_TIME
        time1= omp_get_wtime();
#endif
    double xwidth =  m_dBoundingBox_max_x - m_dBoundingBox_min_x   ;
    double ywidth =  m_dBoundingBox_max_y - m_dBoundingBox_min_y   ;
    double zwidth =  m_dBoundingBox_max_z - m_dBoundingBox_min_z   ;

    m_vLowerLimitOfX[0] = m_dBoundingBox_min_x;
    m_vLowerLimitOfY[0] = m_dBoundingBox_min_y;
    m_vLowerLimitOfZ[0] = m_dBoundingBox_min_z;

    m_vUpperLimitOfX[0] = m_dBoundingBox_max_x;
    m_vUpperLimitOfY[0] = m_dBoundingBox_max_y;
    m_vUpperLimitOfZ[0] = m_dBoundingBox_max_z;

        for(int level=0; level<m_iMaxDepth; ++level)
        {
                xwidth=0.5*xwidth;
                ywidth=0.5*ywidth;
                zwidth=0.5*zwidth;
    #pragma omp parallel for
                for(int i=0; i < numnodes[level] ; ++i)
                {
                        for(int j=m_vFirstChildIndex[baseaddresses[level]+i]; j<m_vFirstChildIndex[baseaddresses[level]+i]+m_vNumberOfChildren[baseaddresses[level]+i]; ++j)
                        {
                                unsigned int bit_triplet_extracted = ( m_vNodeKey[j] >> (3*(m_iMaxDepth-level-1)) ) & 0x7;
        //Case 000
        if (bit_triplet_extracted ==       0)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
          }
        //Case 001
        else if (bit_triplet_extracted  == 1)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
          }
         //Case 010
        else if (bit_triplet_extracted  == 2)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
          }
         //Case 011
        else if (bit_triplet_extracted  == 3)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
          }
          //Case 100
        else if (bit_triplet_extracted ==  4)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
          }
        //Case 101
        else if (bit_triplet_extracted  == 5)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
          }
          //Case 110 
        else if (bit_triplet_extracted  == 6)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
          }
         //Case 111
        else if (bit_triplet_extracted  == 7)
          {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
          }

                        }//j
                }//i
        }//level

#ifdef OUTPUT_TIME
        time1= omp_get_wtime()-time1;
        std::cout << "Bounding box Construction Time : " << time1 << std::endl;
#endif
  //Print the octree information to a file. 
  /*
  std::ofstream octree_file;
  octree_file.open("printed_octree_vector.txt")			;
  octree_file << "BOUNDING_BOX INFORMATION :" << std::endl	;
  octree_file << "Min_x:  "  << m_dBoundingBox_min_x   << std::endl			;
  octree_file << "Max_x:  "  << m_dBoundingBox_max_x<< std::endl			;

  octree_file << "\nMin_y:  " << m_dBoundingBox_min_y<< std::endl			;
  octree_file << "Max_y:  "   << m_dBoundingBox_max_y<< std::endl			;

  octree_file << "\nMin_z:  "  << m_dBoundingBox_min_z<< std::endl			;
  octree_file << "Max_z:  "   << m_dBoundingBox_max_z<< std::endl			;
  
  octree_file << "\nLEVEL\t"<<"BITMASK\t"<<"BASEADDRESSES\t"<<"NUMNODES"<<std::endl;
  for (int i = 0; i < m_iMaxDepth+1; ++i)
    {
      octree_file << i << "\t" << bitmasks[i] << "\t" << baseaddresses[i] << "\t\t" << numnodes[i] << std::endl;
    }
   octree_file << "Length of the octree vector is m_iTreeLength = " << m_iTreeLength << std::endl;  
  for (int level = 0; level <= m_iMaxDepth; ++level)
  {
    octree_file<<"============================================================================================================================================================"<<std::endl;
    octree_file<<"RANK\t"<<"DEPTH\t"<< "KEY\t"<<"PIDX\t"<<"PNUM\t # "<<"CIDX\t"<<"CNUM\t"<<"GLOBAL POSN\t\t"<<"VMIN[0]\tVMIN[1]\tVMIN[2]\t\t\tVMAX[0]\tVMAX[1]\tVMAX[2]"<<std::endl;
    octree_file<<"============================================================================================================================================================"<<std::endl; 
     for (int i = 0; i < numnodes[level]; ++i)
	{
	  octree_file<<std::setprecision(5)<< i      <<"\t"
              << m_vDepth[baseaddresses[level]+i]    <<"\t"  
	      << m_vNodeKey[baseaddresses[level]+i]      <<"\t"
	      << m_vFirstParticleIndex[baseaddresses[level]+i] <<"\t"
	      << m_vNumberOfContainedParticles[baseaddresses[level]+i] <<"\t # "
	      << m_vFirstChildIndex[baseaddresses[level]+i] <<"\t"
	      << m_vNumberOfChildren[baseaddresses[level]+i] <<"\t\t"
		     <<        baseaddresses[level]+i     <<"\t\t" 

		   << m_vLowerLimitOfX[baseaddresses[level]+i] << "\t"
		   << m_vLowerLimitOfY[baseaddresses[level]+i] << "\t"
		   << m_vLowerLimitOfZ[baseaddresses[level]+i] << "\t\t\t"

		   << m_vUpperLimitOfX[baseaddresses[level]+i] << "\t"
		   << m_vUpperLimitOfY[baseaddresses[level]+i] << "\t"
		   << m_vUpperLimitOfZ[baseaddresses[level]+i] << "\t" <<std::endl;

       	}
      octree_file<<'\n';
  }
  */


  return 0;

}//End of the buildOctree method
 


void Octree::clearOctreeStructure() {
  delete[] m_vNodeKey;
  delete[] m_vDepth;
  delete[] m_vFirstParticleIndex;
  delete[] m_vNumberOfContainedParticles;
  delete[] m_vFirstChildIndex;
  delete[] m_vNumberOfChildren;

  delete[] m_vLowerLimitOfX;
  delete[] m_vUpperLimitOfX;

  delete[] m_vLowerLimitOfY;
  delete[] m_vUpperLimitOfY;

  delete[] m_vLowerLimitOfZ;
  delete[] m_vUpperLimitOfZ; 

  m_vNodeKey = 0;
  m_vDepth = 0;
  m_vFirstParticleIndex = 0;
  m_vNumberOfContainedParticles = 0;
  m_vFirstChildIndex = 0;
  m_vNumberOfChildren = 0;

  m_vLowerLimitOfX = 0;
  m_vUpperLimitOfX = 0;

  m_vLowerLimitOfY = 0;
  m_vUpperLimitOfY = 0;

  m_vLowerLimitOfZ = 0;
  m_vUpperLimitOfZ = 0;
}


inline uint32_t Octree::computeKey(const double& x, const double& y, const double& z)
{
  uint32_t    	key=0;//Key which will be returned. 
                      //Excluding this bit, the key can containg upto 30 bits 
                      //corresponding to a max-tree depth of 10. 



  

  




  double	left_x=m_dBoundingBox_min_x ,right_x=m_dBoundingBox_max_x ;
  double	left_y=m_dBoundingBox_min_y ,right_y=m_dBoundingBox_max_y ;
  double	left_z=m_dBoundingBox_min_z ,right_z=m_dBoundingBox_max_z ;

  //Midpoint of the various intervals in the following for loop;
  double	midpt_x;
  double        midpt_y;
  double        midpt_z;

  for (int i = 1; i <= m_iMaxDepth; ++i) {
    //Compute midpoints.
    midpt_x=(left_x+right_x)/2.0;
    midpt_y=(left_y+right_y)/2.0;
    midpt_z=(left_z+right_z)/2.0;

    //Now we consider 8 cases. EXACTLY ONE of the following will be true within this for loop. 
    //Hence we place it within an if-elsef-if condtruct which guarantees mutual exclusion
    
    //left_x, left_y et al. will get modified within these if-else constructs 
    //since the boundary of the containing interval is continously shrinking.

      //Case 1.
     if (x<midpt_x && y<midpt_y && z<midpt_z) {
         // -----------------------------abc becomes --------------------------abc000 
 	    key=(key << 3);//insert three zeros at the end
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=left_z;
         right_z=midpt_z;
       }
   

      //Case 2.
     else if (x>=midpt_x &&   \
              y<midpt_y &&        \
              z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc100 
 	 key=(key << 3) + (1<<2);
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     //Case 3
        else if (x<midpt_x &&   \
                 y>=midpt_y &&        \
                 z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc010 
 	 key=(key << 3) + (1<<1);
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     
     //Case 4
      else if ( x<midpt_x &&   \
                y<midpt_y &&        \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc001 
 	 key=(key << 3) + 1;
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=midpt_z;
         right_z=right_z;
       }

     //Case 5
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc110 
 	 key=(key << 3) + (1<<2) + (1<<1) ;
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y=midpt_y;
         right_y=right_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     //Case 6
      else if ( x>=midpt_x &&   \
                y<midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc101 
 	 key=(key << 3) + (1<<2) + 1 ;
 
         left_x =midpt_x;
         right_x=right_x;
       
         left_y =left_y;
         right_y=midpt_y;

         left_z=midpt_z;
         right_z=right_z;
       }

     //Case 7
      else if ( x<midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc011 
 	 key=(key << 3) + (1<<1) + 1 ;
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z =midpt_z;
         right_z=right_z;
       }

     //Case 8
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc111 
 	 key=(key << 3) + (1<<2) + (1<<1) +1 ;
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z=midpt_z;
         right_z=right_z;
       }

    }//End for loop
  return key ;
}



int Octree::searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, int* result, size_t& result_length) {

  std::deque<int> nodelist;
  nodelist.push_back(0);

  result_length = 0;
  int result_index = 0;

  double radius_sq = radius*radius;


  // before leaf node
  //for(int i=0 ; i<treedepth ; i++) {
  while(1) {
    
    // If nodelist is empty, then there is no leaf node. So there is no leaf node which intersects with the sphere.
    if(nodelist.empty()) return 0;

    // popping from front
    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    // search child nodes which are intersect with sphere haveing radius "radius".
    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }

  // When above while-loop end normally(by the break in the loop), the list is not empty and it has only leaf node.
  // For leaf nodes
  while(!nodelist.empty()) {

    //node cnode = tree[nodelist[0]];
    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];
    
    int n;
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = search_x - m_vCoordX[n];
      double dy = search_y - m_vCoordY[n];
      double dz = search_z - m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius_sq) {
        //result->push_back(m_vParticleKeyIndex[i].index);//push back the index of the particle's coorindates in the x,y,z arrays onto the neghbourlist . 
        result[result_index] = m_vParticleKeyIndex[i].index;
        result_index++;
      }
    }
    nodelist.pop_front();
  }

  result_length = result_index;
  return 0;
   
}





int Octree::searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, SearchResult* result, size_t& result_length) {

  
  std::deque<int> nodelist;
  nodelist.push_back(0);

//  result_length = 0;
  size_t result_index = 0;

  double radius_sq = radius*radius;


  // before leaf node
  //for(int i=0 ; i<treedepth ; i++) {
  while(1) {
    
    // If nodelist is empty, then there is no leaf node. So there is no leaf node which intersects with the sphere.
    if(nodelist.empty()) return 0;

    // popping from front
    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    // search child nodes which are intersect with sphere haveing radius "radius".
    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }

  // When above while-loop end normally(by the break in the loop), the list is not empty and it has only leaf node.
  // For leaf nodes
  while(!nodelist.empty()) {

    //node cnode = tree[nodelist[0]];
    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];
    
    int n;
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = - search_x + m_vCoordX[n];
      double dy = - search_y + m_vCoordY[n];
      double dz = - search_z + m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius_sq) {
        //result->push_back(m_vParticleKeyIndex[i].index);//push back the index of the particle's coorindates in the x,y,z arrays onto the neghbourlist . 
				if(result_index==result_length){
					result_length = result_length + 1;
					return 0;
				}
        result[result_index].index = m_vParticleKeyIndex[i].index;
        result[result_index].distance = sqrt(squared_distance);
	if(dx==0)
		result[result_index].angle=(dy>0)?0.5*M_PI: 1.5*M_PI;
	else
		result[result_index].angle = atan(dy/dx) + ((dx<0)?M_PI: ((dy<0)?2*M_PI: 0));
	if(result[result_index].angle<(45.0/180.0*M_PI))
		result[result_index].angle+=2.0*M_PI;
//	if(rand()%10000==0)
//		cout<<dx<<" "<<dy<<" "<<result[result_index].angle<<endl;
        result_index++;
      }
    }
    nodelist.pop_front();
  }

  result_length = result_index;
  return 0;



}

//Following are 3 geometric support function for quadrant search algorithm in 2D.
inline bool left_test(const double a_x, const double a_y, const double b_x, const double b_y, const double c_x, const double c_y) {
  //if left_test=true then c lies in the left of ab
  //if left_test=false then c lies in the right of ab or on ab
	return (((b_x-a_x)*(c_y-a_y)-(c_x-a_x)*(b_y-a_y))>0.0);
}

inline bool incone_test(const double a_x, const double a_y, const double b_x, const double b_y, const double alpha, const double beta){
	//if incone_test=true then b lies in the cone which center is a, right arm is alpha and left arm is beta.
	//Only works when 0<beta-alpha<M_PI
	double c_x=a_x+cos(alpha);
	double c_y=a_y+sin(alpha);
	double d_x=a_x+cos(beta);
	double d_y=a_y+sin(beta);
	return (left_test(a_x, a_y, c_x, c_y, b_x, b_y)&& left_test(a_x, a_y, b_x, b_y, d_x, d_y));
}

inline bool incone_test(const double a_x, const double a_y, const double a_z, const double b_x, const double b_y, const double b_z, const int dir){
        //if incone_test=true then b lies in the dir direction from a
        //dir=1: x+; -1: x-; 2: y+; -2: y-; 3:z+; -3: z-
	if (dir==1)
		return ((fabs(b_x-a_x)>=fabs(b_y-a_y))&&(fabs(b_x-a_x)>=fabs(b_z-a_z))&&((b_x-a_x)>=0));
        else if (dir==-1)
                return ((fabs(b_x-a_x)>=fabs(b_y-a_y))&&(fabs(b_x-a_x)>=fabs(b_z-a_z))&&((b_x-a_x)<=0));
        else if (dir==2)
                return ((fabs(b_y-a_y)>=fabs(b_x-a_x))&&(fabs(b_y-a_y)>=fabs(b_z-a_z))&&((b_y-a_y)>=0));
        else if (dir==-2)
                return ((fabs(b_y-a_y)>=fabs(b_x-a_x))&&(fabs(b_y-a_y)>=fabs(b_z-a_z))&&((b_y-a_y)<=0));
        else if (dir==3)
                return ((fabs(b_z-a_z)>=fabs(b_x-a_x))&&(fabs(b_z-a_z)>=fabs(b_y-a_y))&&((b_z-a_z)>=0));
        else if (dir==-3)
                return ((fabs(b_z-a_z)>=fabs(b_x-a_x))&&(fabs(b_z-a_z)>=fabs(b_y-a_y))&&((b_z-a_z)<=0));
	else
		return 0;
}

inline bool is_node_intersect_search_region(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const double& search_x, const double& search_y, const double& search_z, const double& radius, const double alpha, const double beta) {
	//if output=true then rectangle [min_x, max_x]*[min_y,max_y] intersects the fan which center is (search_x, search_y), radius is radius, right arm is alpha and left arm is beta.

if(!is_node_intersect_search_region(min_x, max_x, min_y, max_y, 0.0, 0.0, search_x, search_y, 0.0, radius)) return false;

if(incone_test(search_x,search_y,min_x,min_y,alpha,beta)) return true;
if(incone_test(search_x,search_y,min_x,max_y,alpha,beta)) return true;
if(incone_test(search_x,search_y,max_x,min_y,alpha,beta)) return true;
if(incone_test(search_x,search_y,max_x,max_y,alpha,beta)) return true;

double sinalpha=sin(alpha),cosalpha=cos(alpha);
double b_x=search_x+cosalpha;
double b_y=search_y+sinalpha;
if(((left_test(search_x,search_y,b_x,b_y,min_x,min_y)+left_test(search_x,search_y,b_x,b_y,max_x,min_y))==1)&&((min_y-search_y)*sinalpha>0)) return true;
if(((left_test(search_x,search_y,b_x,b_y,min_x,max_y)+left_test(search_x,search_y,b_x,b_y,max_x,max_y))==1)&&((max_y-search_y)*sinalpha>0)) return true;
if(((left_test(search_x,search_y,b_x,b_y,min_x,min_y)+left_test(search_x,search_y,b_x,b_y,min_x,max_y))==1)&&((min_x-search_x)*cosalpha>0)) return true;
if(((left_test(search_x,search_y,b_x,b_y,max_x,min_y)+left_test(search_x,search_y,b_x,b_y,max_x,max_y))==1)&&((max_x-search_x)*cosalpha>0)) return true;

return false;
}

inline bool is_node_intersect_search_region(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const double& search_x, const double& search_y, const double& search_z, const double& radius, int dir) {

if(!is_node_intersect_search_region(min_x, max_x, min_y, max_y, min_z, max_z, search_x, search_y, search_z, radius)) return false;

double min_dx;
double min_dy;
double min_dz;

if (min_x>search_x)
	min_dx=min_x-search_x;
else if (max_x<search_x)
	min_dx=search_x-max_x;
else
	min_dx=0.;

if (min_y>search_y)
        min_dy=min_y-search_y;
else if (max_y<search_y)
        min_dy=search_y-max_y;
else
        min_dy=0.;

if (min_z>search_z)
        min_dz=min_z-search_z;
else if (max_z<search_z)
        min_dz=search_z-max_z;
else
        min_dz=0.;

if (dir==1)
	return ((max_x>=search_x)&&((max_x-search_x)>=min_dy)&&((max_x-search_x)>=min_dz));
else if (dir==-1)
        return ((min_x<=search_x)&&((search_x-min_x)>=min_dy)&&((search_x-min_x)>=min_dz));
else if (dir==2)
        return ((max_y>=search_y)&&((max_y-search_y)>=min_dx)&&((max_y-search_y)>=min_dz));
else if (dir==-2)
        return ((min_y<=search_y)&&((search_y-min_y)>=min_dx)&&((search_y-min_y)>=min_dz));
else if (dir==3)
        return ((max_z>=search_z)&&((max_z-search_z)>=min_dx)&&((max_z-search_z)>=min_dy));
else if (dir==-3)
        return ((min_z<=search_z)&&((search_z-min_z)>=min_dx)&&((search_z-min_z)>=min_dy));
else
	return 0;
}

int Octree::searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, SearchResult* result, size_t& result_length, int count, double alpha, double beta) {
//Function for 2d neighbour search in a fan with center = (search_x, search_y, search_z), radius = radius and angle = [alpha, beta]
  std::deque<int> nodelist;
  nodelist.push_back(0);

  size_t result_index = 0;

  double radius_sq = radius*radius;

//Locate octree leaf nodes that may contain the neighbours, using a top-to-botoom recursive algorithm testing the intersection of the subtree and the fan
  // before leaf node
  while(1) {
    
    // If nodelist is empty, then there is no leaf node. So there is no leaf node which intersects with the sphere.
    if(nodelist.empty()) 
		{
			result_length = result_index;
			return 0;
		}

    // popping from front
    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    // search child nodes which are intersect with the fan
    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius, alpha, beta)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }

  // When above while-loop end normally(by the break in the loop), the list is not empty and it has only leaf node.
  // For leaf nodes given by the previous loop, test if each particle inside are in the fine or not
  while(!nodelist.empty()) {

    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];
    
    int n;
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = search_x - m_vCoordX[n];
      double dy = search_y - m_vCoordY[n];
      double dz = search_z - m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( (squared_distance < radius_sq)&&incone_test(dx,dy,0.0,0.0,alpha,beta)) {
//Test if the number of neighbours are greater than the capacity of the result list. If so, return an error
				if(result_index==result_length){
					result_length = result_length + 1;
					return 0;
				}
        result[result_index+count].index = m_vParticleKeyIndex[i].index;//put the index of the neighbour in the result
        result[result_index+count].distance = sqrt(squared_distance);
        result_index++;
      }
    }
    nodelist.pop_front();
  }

  result_length = result_index;
  return 0;



}


int Octree::searchNeighbor(const double search_x, const double search_y, const double search_z, const double radius, SearchResult* result, size_t& result_length, int count, int dir) {
//dir=1: x+; -1: x-; 2: y+; -2: y-; 3:z+; -3: z-
  std::deque<int> nodelist;
  nodelist.push_back(0);

//  result_length = 0;
  size_t result_index = 0;

  double radius_sq = radius*radius;


  // before leaf node
  //for(int i=0 ; i<treedepth ; i++) {
  while(1) {

    // If nodelist is empty, then there is no leaf node. So there is no leaf node which intersects with the sphere.
    if(nodelist.empty())
                {
                        result_length = result_index;
                        return 0;
                }

    // popping from front
    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    // search child nodes which are intersect with sphere haveing radius "radius".
    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius, dir)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }

  // When above while-loop end normally(by the break in the loop), the list is not empty and it has only leaf node.
  // For leaf nodes
  while(!nodelist.empty()) {

    //node cnode = tree[nodelist[0]];
    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];

    int n;
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = search_x - m_vCoordX[n];
      double dy = search_y - m_vCoordY[n];
      double dz = search_z - m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( (squared_distance < radius_sq)&&incone_test(dx,dy,dz,0.0,0.0,0.0,dir)) {
        //result->push_back(m_vParticleKeyIndex[i].index);//push back the index of the particle's coorindates in the x,y,z arrays onto the neghbourlist . 

                                if(result_index==result_length){
                                        result_length = result_length + 1;
                                        return 0;
                                }
        result[result_index+count].index = m_vParticleKeyIndex[i].index;
        result[result_index+count].distance = sqrt(squared_distance);
        result_index++;
      }
    }
    nodelist.pop_front();
  }

  result_length = result_index;
  return 0;



}






int Octree::densityEstimator(const double search_x, const double search_y, const double search_z, const double radius, double* density_count, double dir_x, double dir_y, double dir_z)
{
	return densityEstimator(0, search_x, search_y, search_z, radius, density_count, dir_x, dir_y, dir_z);
}



int Octree::densityEstimator(const int search_index, const double search_x, const double search_y, const double search_z, const double radius, double* density_count, double dir_x, double dir_y, double dir_z){

  std::deque<int> nodelist;
  nodelist.push_back(0);
  double old_density=(*density_count);
  (*density_count)=0;
//  int result_index = 0;

  double radius_sq = radius*radius;

  double mag=sqrt(dir_x*dir_x+dir_y*dir_y+dir_z*dir_z);
  double dirx, diry, dirz;
  if(mag<1e-08)
  {
    dirx=1.0;
    diry=0.0;
    dirz=0.0; 
  }
  else
  {
    dirx=dir_x/mag;
    diry=dir_y/mag;
    dirz=dir_z/mag;
  }
//  dirx=1.0;
//  diry=0.0;
//  dirz=0.0;

  while(1) {

    if(nodelist.empty()) return 0;

    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }
    double h=radius/2;
    double c=10.0/7.0/M_PI/h/h;
//    double cl=3/M_PI/radius/radius/radius;
    int n;
    double count_left=0;
    double count_right=0;
    double delta=0.2*h;
    int index0=0;
    int small_right=0;
    int small_left=0;
    int big_right=0;
    int big_left=0;
    int bad_right=0;
    int bad_left=0;


  while(!nodelist.empty()) {

    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];
    
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = - search_x + m_vCoordX[n];
      double dy = - search_y + m_vCoordY[n];
      double dz = - search_z + m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
	if (squared_distance<0.000001)
		index0=n;
      if( squared_distance < radius_sq) {
	double distance=sqrt(squared_distance);
	double u=distance/h;
	double dir=dirx*dx+diry*dy+dirz*dz;
	if (u<1){
//		(*density_count)=(*density_count)+c*(1-1.5*u*u+0.75*u*u*u);
		if(fabs(dir)<delta){
			count_left=count_left+c*(1-1.5*u*u+0.75*u*u*u)*(1-dir/delta)*mass[n];
        	        count_right=count_right+c*(1-1.5*u*u+0.75*u*u*u)*(1+dir/delta)*mass[n];
//                        count_left=count_left+c*(1-1.5*u*u+0.75*u*u*u)*(1-sin(0.5*M_PI*dir/delta))*mass[n];
//                        count_right=count_right+c*(1-1.5*u*u+0.75*u*u*u)*(1+sin(0.5*M_PI*dir/delta))*mass[n];
			if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
				big_right++;
				big_left++;
			}
			else
			{
				small_right++;
				small_left++;
			}
		}
		else if (dir<0){
                        count_left=count_left+c*(1-1.5*u*u+0.75*u*u*u)*2*mass[n];
                        if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
                                big_left++;
                        }
                        else
                        {
                                small_left++;
                        }
		}
		else{
                        count_right=count_right+c*(1-1.5*u*u+0.75*u*u*u)*2*mass[n];
                        if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
                                big_right++;
                        }
                        else
                        {
                                small_right++;
                        }

		}
	}
	else{
//		(*density_count)=(*density_count)+c/4*(2-u)*(2-u)*(2-u);
                if(fabs(dir)<delta){
                        count_left=count_left+c/4*(2-u)*(2-u)*(2-u)*(1-dir/delta)*mass[n];
                        count_right=count_right+c/4*(2-u)*(2-u)*(2-u)*(1+dir/delta)*mass[n];
//                        count_left=count_left+c/4*(2-u)*(2-u)*(2-u)*(1-sin(0.5*M_PI*dir/delta))*mass[n];
//                        count_right=count_right+c/4*(2-u)*(2-u)*(2-u)*(1+sin(0.5*M_PI*dir/delta))*mass[n];
                        if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
                                big_right++;
                                big_left++;
                        }
                        else
                        {
                                small_right++;
                                small_left++;
                        }
                }
                else if (dir<0){
                        count_left=count_left+2*c/4*(2-u)*(2-u)*(2-u)*mass[n];
                        if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
                                big_left++;
                        }
                        else
                        {
                                small_left++;
                        }
		}
                else{
                        count_right=count_right+2*c/4*(2-u)*(2-u)*(2-u)*mass[n];
                        if (fabs(mass[n]/mass[search_index]-1.0)>0.3){
                                big_right++;
                        }
                        else
                        {
                                small_right++;
                        }
		}
	}
//	(*density_count)=(*density_count)+cl*(radius-distance);
	
      }
	
    }
    nodelist.pop_front();
  }


	if(fabs(mass[index0]/mass[search_index]-1.0)<0.3)
	{
		bad_left=big_left;
		bad_right=big_right;
	}
	else
	{
		bad_left=small_left;
		bad_right=small_right;
	}

	double sigma=0.20;
	double error_left=fabs(1/count_left-1/old_density)*old_density;
	double error_right=fabs(1/count_right-1/old_density)*old_density;
//	if (((error_left<sigma)&&(error_right<sigma)))
//	if (1)
/*	if(bad_left+bad_right)
	{
		cout<<search_x<<" "<<search_y<<endl;
		cout<<radius<<endl;
		cout<<small_left<<" "<<big_left<<" "<<bad_left<<endl;
                cout<<small_right<<" "<<big_right<<" "<<bad_right<<endl;
		cout<<count_left<<" "<<count_right<<endl;
	}*/
	if ((bad_left==bad_right)&&(error_left<sigma)&&(error_right<sigma))
//	if(bad_left==bad_right)
		(*density_count)=(count_left+count_right)/2;
	else if((bad_left<bad_right)||((bad_left==bad_right)&&(error_left<error_right))){
//	else if(bad_left<bad_right) {
//	else if(error_left<error_right){
//	else if (mass[index0]>0.0001){
		(*density_count)=count_left;
/*		cout<<"left "<<search_x<<" "<<search_y<<endl;
		cout<<radius<<endl;
		cout<<old_density<<endl;
		cout<<count_left<<" "<<count_right<<endl;
                cout<<error_left<<" "<<error_right<<endl;
		cout<<endl;*/
	}
	else{
		(*density_count)=count_right;
/*                cout<<"right "<<search_x<<" "<<search_y<<endl;
                cout<<radius<<endl;
                cout<<old_density<<endl;
                cout<<count_left<<" "<<count_right<<endl;
                cout<<error_left<<" "<<error_right<<endl;
                cout<<endl;*/
        }

//  (*density_count)=(*density_count)/radius;
  return 0;

}

int Octree::VoronoiDensityEstimator(const int search_index, const double search_x, const double search_y, const double search_z, const double radius, double* density_count){

  std::deque<int> nodelist;
  nodelist.push_back(0);

  double radius_sq = radius*radius;

  while(1) {

    if(nodelist.empty()) return 0;

    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }
    int n;
    double h=radius/2.0;
    double volume_weight_sum=0.0;
    double weight_sum=0.0;
    double weight;

  while(!nodelist.empty()) {

    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];

    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = - search_x + m_vCoordX[n];
      double dy = - search_y + m_vCoordY[n];
      double dz = - search_z + m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
//	cout<<search_index<<" "<<mass[n]<<" "<<mass[search_index]<<" "<<fabs(mass[n]/mass[search_index]-1)<<endl;
      if( (squared_distance < radius_sq) && (fabs(mass[n]/mass[search_index]-1)<0.1)) {
        double distance=sqrt(squared_distance);
        double u=distance/h;

	if (u<1)
		weight=(1-1.5*u*u+0.75*u*u*u);
	else
		weight=0.25*(2-u)*(2-u)*(2-u);

        weight_sum=weight_sum+weight;
        volume_weight_sum=volume_weight_sum+weight*VolumeVoronoi[n]/mass[n];
      }
    }
    nodelist.pop_front();
  }
  (*density_count)=weight_sum/volume_weight_sum;
//  cout<<search_index<<" "<<search_x<<" "<<search_y<<" "<<search_z<<" "<<VolumeVoronoi[search_index] <<" "<<weight_sum<<" "<<volume_weight_sum<<endl; 
  return 0;

}

int Octree::densityEstimator(const int search_index, const double search_x, const double search_y, const double search_z, const double radius, double* density_count, const double* volume, double* vmin, double* vmax){

  std::deque<int> nodelist;
  nodelist.push_back(0);
  (*density_count)=0;
//  int result_index = 0;

  double radius_sq = radius*radius;

  while(1) {

    if(nodelist.empty()) return 0;

    int inode = nodelist[0];

    int first_index = m_vFirstChildIndex[inode];
    if(first_index == -1) break;

    int last_index = first_index + m_vNumberOfChildren[inode];

    for(int i=first_index ; i < last_index ; i++) {
      if(is_node_intersect_search_region(m_vLowerLimitOfX[i], m_vUpperLimitOfX[i], m_vLowerLimitOfY[i], m_vUpperLimitOfY[i], m_vLowerLimitOfZ[i], m_vUpperLimitOfZ[i], search_x, search_y, search_z, radius)) nodelist.push_back(i);
    }
    nodelist.pop_front();
  }
    double h=radius/2;
    double c=10.0/7.0/M_PI/h/h;
//    double cl=3/M_PI/radius/radius/radius;
    int n;
    double count=0;

  while(!nodelist.empty()) {

    int index = nodelist[0];

    int first = m_vFirstParticleIndex[index];
    int end = first + m_vNumberOfContainedParticles[index];
    
    for (int i = first; i < end; i++) {
      n = m_vParticleKeyIndex[i].index;
      double dx = - search_x + m_vCoordX[n];
      double dy = - search_y + m_vCoordY[n];
      double dz = - search_z + m_vCoordZ[n];

      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius_sq) {
	double distance=sqrt(squared_distance);
	double u=distance/h;
	if(volume[n]>(*vmax))
		(*vmax)=volume[n];
	if(volume[n]<(*vmin))
		(*vmin)=volume[n];
	if (u<1){
//		(*density_count)=(*density_count)+c*(1-1.5*u*u+0.75*u*u*u);
		count=count+c*(1-1.5*u*u+0.75*u*u*u)*mass[n];
	}
	else{
//		(*density_count)=(*density_count)+c/4*(2-u)*(2-u)*(2-u);
		count=count+c/4*(2-u)*(2-u)*(2-u)*mass[n];
	}
//	(*density_count)=(*density_count)+cl*(radius-distance);
	
      }
	
    }
    nodelist.pop_front();
  }

	(*density_count)=count;

//  (*density_count)=(*density_count)/radius;
  return 0;

}
