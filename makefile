CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
LAPACK_DIR=/home/xingyu/lapack-3.5.0
INCS   = -I $(LAPACK_DIR)/lapacke/include
LIBS   = -L $(LAPACK_DIR)
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP)
LFLAGS = -Wall $(DEBUG) $(INCS) $(LIBS) $(OMP)
OBJS   = boundary.o boundary_solid_shocktube.o boundary_solid_shocktube3d.o boundary_solid_tpshocktube.o boundary_solid_gresho.o boundary_rayleightaylor.o boundary_rayleightaylor_periodic.o boundary_rayleightaylor3d.o boundary_kelvinhelmholtz.o boundary_dambreak.o boundary_powder_target.o boundary_powder_target_3d.o eos.o geometry.o geometry_1d.o boundary_nozzle.o geometry_nozzle.o\
		 geometry_ballexp.o geometry_collision.o geometry_gresho.o geometry_powder_target.o geometry_powder_target_3d.o\
		 geometry_jet.o geometry_shocktube.o geometry_shocktube3d.o hexagonal_packing.o\
		 initializer.o lp_main.o lp_solver.o ls_solver.o\
         neighbour_searcher.o octree.o particle_data.o\
		 particle_viewer.o registrar.o state.o state_1d.o\
		 state_ballexp.o state_collision.o state_gresho.o state_powder_target.o state_powder_target_3d.o state_nozzle.o\
	     state_jet.o state_shocktube.o time_controller.o

all: lp

lp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o lp -lgomp -llapacke -llapack -lgfortran -lrefblas

boundary.o: boundary.h boundary.cpp
	$(CC) $(CFLAGS) boundary.cpp

boundary_solid_gresho.o: boundary.h boundary_solid_gresho.h boundary_solid_gresho.cpp
	$(CC) $(CFLAGS) boundary_solid_gresho.cpp

boundary_solid_shocktube.o: boundary.h boundary_solid_shocktube.h boundary_solid_shocktube.cpp
	$(CC) $(CFLAGS) boundary_solid_shocktube.cpp

boundary_solid_tpshocktube.o: boundary.h boundary_solid_tpshocktube.h boundary_solid_tpshocktube.cpp
	$(CC) $(CFLAGS) boundary_solid_tpshocktube.cpp

boundary_solid_shocktube3d.o: boundary.h boundary_solid_shocktube3d.h boundary_solid_shocktube3d.cpp
	$(CC) $(CFLAGS) boundary_solid_shocktube3d.cpp

boundary_rayleightaylor.o: boundary.h boundary_rayleightaylor.h boundary_rayleightaylor.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor.cpp

boundary_rayleightaylor_periodic.o: boundary.h boundary_rayleightaylor_periodic.h boundary_rayleightaylor_periodic.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor_periodic.cpp

boundary_rayleightaylor3d.o: boundary.h boundary_rayleightaylor3d.h boundary_rayleightaylor3d.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor3d.cpp

boundary_kelvinhelmholtz.o: boundary.h boundary_kelvinhelmholtz.h boundary_kelvinhelmholtz.cpp
	$(CC) $(CFLAGS) boundary_kelvinhelmholtz.cpp

boundary_dambreak.o: boundary.h boundary_dambreak.h boundary_dambreak.cpp
	$(CC) $(CFLAGS) boundary_dambreak.cpp

boundary_powder_target.o: boundary.h boundary_powder_target.h boundary_powder_target.cpp
	$(CC) $(CFLAGS) boundary_powder_target.cpp

boundary_powder_target_3d.o: boundary.h boundary_powder_target_3d.h boundary_powder_target_3d.cpp
	$(CC) $(CFLAGS) boundary_powder_target_3d.cpp

boundary_nozzle.o: boundary.h boundary_nozzle.h boundary_nozzle.cpp
	$(CC) $(CFLAGS) boundary_nozzle.cpp

eos.o: eos.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

geometry.o: geometry.h geometry.cpp
	$(CC) $(CFLAGS) geometry.cpp

geometry_1d.o: geometry.h geometry_1d.h geometry_1d.cpp
	$(CC) $(CFLAGS) geometry_1d.cpp

geometry_ballexp.o: geometry.h geometry_ballexp.h geometry_ballexp.cpp
	$(CC) $(CFLAGS) geometry_ballexp.cpp

geometry_collision.o: geometry.h geometry_collision.h geometry_collision.cpp
	$(CC) $(CFLAGS) geometry_collision.cpp

geometry_gresho.o: geometry.h geometry_gresho.h geometry_gresho.cpp
	$(CC) $(CFLAGS) geometry_gresho.cpp

geometry_jet.o: geometry.h geometry_jet.h geometry_jet.cpp
	$(CC) $(CFLAGS) geometry_jet.cpp

geometry_shocktube.o: geometry.h geometry_shocktube.h geometry_shocktube.cpp
	$(CC) $(CFLAGS) geometry_shocktube.cpp

geometry_shocktube3d.o: geometry.h geometry_shocktube3d.h geometry_shocktube3d.cpp
	$(CC) $(CFLAGS) geometry_shocktube3d.cpp

geometry_powder_target.o: geometry.h geometry_powder_target.h geometry_powder_target.cpp
	$(CC) $(CFLAGS) geometry_powder_target.cpp

geometry_powder_target_3d.o: geometry.h geometry_powder_target_3d.h geometry_powder_target_3d.cpp
	$(CC) $(CFLAGS) geometry_powder_target_3d.cpp

geometry_nozzle.o: geometry.h geometry_nozzle.h geometry_nozzle.cpp
	$(CC) $(CFLAGS) geometry_nozzle.cpp

hexagonal_packing.o: hexagonal_packing.h hexagonal_packing.cpp
	$(CC) $(CFLAGS) hexagonal_packing.cpp

initializer.o: initializer.h initializer.cpp geometry.h state.h eos.h hexagonal_packing.h\
				neighbour_searcher.h
	$(CC) $(CFLAGS) initializer.cpp

lp_main.o: lp_main.cpp initializer.h neighbour_searcher.h particle_data.h particle_viewer.h\
           lp_solver.h time_controller.h
	$(CC) $(CFLAGS) lp_main.cpp

lp_solver.o: lp_solver.h lp_solver.cpp neighbour_searcher.h eos.h particle_data.h\
           initializer.h ls_solver.h hexagonal_packing.h
	$(CC) $(CFLAGS) lp_solver.cpp

ls_solver.o: ls_solver.h ls_solver.cpp
	$(CC) $(CFLAGS) ls_solver.cpp

neighbour_searcher.o: neighbour_searcher.h neighbour_searcher.cpp octree.h initializer.h
	$(CC) $(CFLAGS) neighbour_searcher.cpp

octree.o: octree.h octree.cpp
	$(CC) $(CFLAGS) octree.cpp

particle_data.o: particle_data.h particle_data.cpp initializer.h
	$(CC) $(CFLAGS) particle_data.cpp

particle_viewer.o: particle_viewer.h particle_viewer.cpp particle_data.h
	$(CC) $(CFLAGS) particle_viewer.cpp

registrar.o: registrar.h registrar.cpp geometry.h state.h geometry_collision.h state_collision.h
	$(CC) $(CFLAGS) registrar.cpp

state.o: state.h state.cpp
	$(CC) $(CFLAGS) state.cpp

state_1d.o: state.h state_1d.h state_1d.cpp
	$(CC) $(CFLAGS) state_1d.cpp

state_ballexp.o: state.h state_ballexp.h state_ballexp.cpp
	$(CC) $(CFLAGS) state_ballexp.cpp

state_collision.o: state.h state_collision.h state_collision.cpp
	$(CC) $(CFLAGS) state_collision.cpp

state_gresho.o: state.h state_gresho.h state_gresho.cpp
	$(CC) $(CFLAGS) state_gresho.cpp

state_jet.o: state.h state_jet.h state_jet.cpp
	$(CC) $(CFLAGS) state_jet.cpp

state_shocktube.o: state.h state_shocktube.h state_shocktube.cpp
	$(CC) $(CFLAGS) state_shocktube.cpp

state_powder_target.o: state.h state_powder_target.h state_powder_target.cpp
	$(CC) $(CFLAGS) state_powder_target.cpp

state_powder_target_3d.o: state.h state_powder_target_3d.h state_powder_target_3d.cpp
	$(CC) $(CFLAGS) state_powder_target_3d.cpp

time_controller.o: time_controller.h time_controller.cpp lp_solver.h particle_viewer.h initializer.h
	$(CC) $(CFLAGS) time_controller.cpp

state_nozzle.o: state.h state_nozzle.h state_nozzle.cpp
	$(CC) $(CFLAGS) state_nozzle.cpp

clean:
	rm *.o *~ debug log save_init_param lp -r out

tar:
	tar cvzf lp_backup.tar.gz *.h *.cpp makefile* input* README* dox* run* data_processor/*.cpp data_processor/*.h data_processor/makefile

