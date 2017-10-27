CC     = g++
DEBUG  = -g
#OMP    = -fopenmp
OMP =
INCS = 
LIBS =  
#INCS   = -I ~/Software/lapack-3.4.2/lapacke/include
#LIBS   = -L ~/Software/lapack-3.4.2
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP)
LFLAGS = -Wall $(DEBUG) $(INCS) $(LIBS) $(OMP)
OBJS   = neighbour_searcher.o octree.o vtkfile_processor.o interpolation_main.o 

all: interpolation

interpolation: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o interpolation

neighbour_searcher.o: octree.h neighbour_searcher.h neighbour_searcher.cpp  
	$(CC) $(CFLAGS) neighbour_searcher.cpp

octree.o: octree.h octree.cpp
	$(CC) $(CFLAGS) octree.cpp

interpolation_main.o: vtkfile_processor.h interpolation_main.cpp
	$(CC) $(CFLAGS) interpolation_main.cpp

vtkfile_processor.o: octree.h neighbour_searcher.h vtkfile_processor.h vtkfile_processor.cpp
	$(CC) $(CFLAGS) vtkfile_processor.cpp




clean:
	rm *.o *~ interpolation output/* 

tar:
	tar cvzf reduce.gz *.h *.cpp makefile* 

