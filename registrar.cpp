#include "registrar.h" 
#include "geometry_collision.h"
#include "state_collision.h"
#include "geometry_jet.h"
#include "state_jet.h"
#include "geometry_1d.h"
#include "state_1d.h"
#include "geometry_shocktube.h"
#include "geometry_shocktube3d.h"
#include "state_shocktube.h"
#include "geometry_gresho.h"
#include "state_gresho.h"
#include "boundary_solid_shocktube.h"
#include "boundary_solid_shocktube3d.h"
#include "boundary_solid_tpshocktube.h"
#include "boundary_rayleightaylor.h"
#include "boundary_rayleightaylor_periodic.h"
#include "boundary_rayleightaylor3d.h"
#include "boundary_kelvinhelmholtz.h"
#include "boundary_dambreak.h"
#include "boundary_solid_gresho.h"
#include "boundary_powder_target.h"
#include "boundary_powder_target_3d.h"
#include "boundary_nozzle.h"
#include "geometry_powder_target.h"
#include "geometry_powder_target_3d.h"
#include "geometry_nozzle.h"
#include "state_powder_target.h"
#include "state_powder_target_3d.h"
#include "state_nozzle.h"
#include "geometry_ballexp.h"
#include "state_ballexp.h"

namespace 
{
	
	GeometryRegistrar<Ball> r1("ball");
	GeometryRegistrar<Disk> r2("disk");	
	
	StateRegistrar<GaussianPressureState> s1("gauss_pressure");

	// 2d collision simulation
	GeometryRegistrar<DiskLeft> r3("disk_left");
	GeometryRegistrar<DiskRight> r4("disk_right");
	StateRegistrar<LeftUniformVelocityState> s2("left_uniform_velocity");
	StateRegistrar<RightUniformVelocityState> s3("right_uniform_velocity");

	StateRegistrar<UniformVelocityState> s4("uniform_velocity");


	// 1d
	GeometryRegistrar<Line> r5("line");
	StateRegistrar<GaussianPressure1DState> s5("1d_gauss_pressure");

	// 3d jet splash
	GeometryRegistrar<Jet3D> r6("jet3d");
	StateRegistrar<Jet3DState> s6("jet3dstate");
	
	// 1d jet splash
	GeometryRegistrar<Jet1D> r7("jet1d");
	StateRegistrar<Jet1DState> s7("jet1dstate");

	
	// 2d jet expansion
	GeometryRegistrar<Jet2DExpansion> r8("jet2dexp");
	StateRegistrar<Jet2DExpansionState> s8("jet2dexpstate");
	

	// 3d jet expansion
	GeometryRegistrar<Jet3DExpansion> r9("jet3dexp");
	StateRegistrar<Jet3DExpansionState> s9("jet3dexpstate");


	// 2d shock tube
	GeometryRegistrar<Shocktube2D> r10("shocktube2d");
	StateRegistrar<Shocktube2DState> s10("shocktube2dstate");
	BoundaryRegistrar<Shocktube2DSolidBoundary> b10("shocktube2dsolidboundary");

	// 2d Gresho test
	GeometryRegistrar<Gresho2D> r11("gresho2d");
	StateRegistrar<Gresho2DState> s11("gresho2dstate");
	BoundaryRegistrar<Gresho2DSolidBoundary> b11("gresho2dsolidboundary");

	// 3d gas ball expansion into vacuum test
	GeometryRegistrar<Ballexp3D> r12("ballexp3d");
	StateRegistrar<Ballexp3DState> s12("ballexp3dstate");

	// 2d jet merge: upper jet
	GeometryRegistrar<Jet2DMergeUpper> r13("jet2dmergeupper");
	StateRegistrar<Jet2DMergeUpperState> s13("jet2dmergeupperstate");

	// 2d jet merge: lower jet
	GeometryRegistrar<Jet2DMergeLower> r14("jet2dmergelower");
	StateRegistrar<Jet2DMergeLowerState> s14("jet2dmergelowerstate");

	// 2d jet collision
	GeometryRegistrar<Jet2DCollision> r15("jet2dcollision");
	StateRegistrar<Jet2DCollisionState> s15("jet2dcollisionstate");

	// 2d jet merge
	GeometryRegistrar<Jet2DMerge> r16("jet2dmerge");
	StateRegistrar<Jet2DMergeState> s16("jet2dmergestate");

	// 2d normal shock
	GeometryRegistrar<Shocktube2D> r17("shocktube2d");
	StateRegistrar<NormalShock2DState> s17("normalshock2dstate");
	BoundaryRegistrar<Shocktube2DSolidBoundary> b17("shocktube2dsolidboundary");

        // 2d simple wave
        GeometryRegistrar<Shocktube2D> r18("shocktube2d");
        StateRegistrar<SimpleWave2DState> s18("simplewave2dstate");
        BoundaryRegistrar<Shocktube2DSolidBoundary> b18("shocktube2dsolidboundary");

        // 2d ConvergentShock test
        GeometryRegistrar<Gresho2D> r19("gresho2d");
        StateRegistrar<ConvergentShock2DState> s19("convergentshock2dstate");
        BoundaryRegistrar<Gresho2DSolidBoundary> b19("gresho2dsolidboundary");

        // Sod's shocktube test: left
        GeometryRegistrar<Shocktube2DLeft> r20("sodleft");
        StateRegistrar<SodShocktube2DState> s20("sodstate");

        // Sod's shocktube test: right
        GeometryRegistrar<Shocktube2DRight> r21("sodright");

        GeometryRegistrar<BigShocktube2D> r22("bigshocktube2d");

        // 2d RayleighTaylor instability test
        GeometryRegistrar<RayleighTaylor2D> r23("rayleightaylor2d");
        StateRegistrar<RayleighTaylor2DState> s23("rayleightaylor2dstate");
        BoundaryRegistrar<RayleighTaylor2DBoundary> b23("rayleightaylor2dboundary");

	// 3d Shocktube
        GeometryRegistrar<Shocktube3D> r24("shocktube3d");
        BoundaryRegistrar<Shocktube3DSolidBoundary> b24("shocktube3dsolidboundary");

        // 3d RayleighTaylor instability test
        GeometryRegistrar<RayleighTaylor3D> r25("rayleightaylor3d");
        StateRegistrar<RayleighTaylor3DState> s25("rayleightaylor3dstate");
        BoundaryRegistrar<RayleighTaylor3DBoundary> b25("rayleightaylor3dboundary");

	// 2d Sod shocktube starting at t=0.5
        StateRegistrar<SodShocktube2DLaterState> s26("sodlaterstate");

        // 2d Dam Break test
        GeometryRegistrar<DamBreak2D> r27("dambreak2d");
        StateRegistrar<DamBreak2DState> s27("dambreak2dstate");
        BoundaryRegistrar<DamBreak2DBoundary> b27("dambreak2dboundary");

	// 2d mirror boundary test
        GeometryRegistrar<BoundaryTest2D> r28("boundarytest2d");
        StateRegistrar<BoundaryTest2DState> s28("boundarytest2dstate");

	// 2d powder traget test
	GeometryRegistrar<PowderTarget2D> r29("powdertarget2d");
	StateRegistrar<PowderTarget2DState> s29("powdertarget2dstate");
	BoundaryRegistrar<PowderTarget2DSolidBoundary> b29("powdertarget2dsolidboundary");

	//3d powder target test
        GeometryRegistrar<PowderTarget3D> r30("powdertarget3d");
        StateRegistrar<PowderTarget3DState> s30("powdertarget3dstate");
        BoundaryRegistrar<PowderTarget3DSolidBoundary> b30("powdertarget3dsolidboundary");

	//1d jet merge left
        GeometryRegistrar<Jet1DLeft> r31("jet1dleft");
        StateRegistrar<Jet1DLeftState> s31("jet1dleftstate");

	//1d jet merge right
        GeometryRegistrar<Jet1DRight> r32("jet1dright");
        StateRegistrar<Jet1DRightState> s32("jet1drightstate");

	//1d symmetric normal shock
        GeometryRegistrar<Jet1DCenter> r33("jet1dcenter");
        StateRegistrar<Jet1DCenterState> s33("jet1dcenterstate");

	//1d symmetric normal shock starting at a later time
        StateRegistrar<Jet1DLaterState> s34("jet1dlaterstate");

	//3d ball rotate
        StateRegistrar<Ballrotate3DState> s35("ballrotate3dstate");

	//3d ball pressure wave
        StateRegistrar<Ballpressurewave3DState> s36("ballpressurewave3dstate");

	//2d Yee vortex
        GeometryRegistrar<Yee2D> r37("yee2d");
        StateRegistrar<Yee2DState> s37("yee2dstate");
        BoundaryRegistrar<Yee2DSolidBoundary> b37("yee2dsolidboundary");

	//2d KelvinHelmholtz instability test
	GeometryRegistrar<KelvinHelmholtz2D> r38("kelvinhelmholtz2d");
        StateRegistrar<KelvinHelmholtz2DState> s38("kelvinhelmholtz2dstate");
        BoundaryRegistrar<KelvinHelmholtz2DBoundary> b38("kelvinhelmholtz2dboundary");

	//2d Sedov blastwave
        GeometryRegistrar<Sedov2D> r39("sedov2d");
        StateRegistrar<Sedov2DState> s39("sedov2dstate");
        BoundaryRegistrar<Sedov2DSolidBoundary> b39("sedov2dsolidboundary");

	//2d Noh problem
        StateRegistrar<Noh2DState> s40("noh2dstate");

	//2d three point shock test
        GeometryRegistrar<TPShocktube2D> r41("tpshocktube2d");
        StateRegistrar<TPShocktube2DState> s41("tpshocktube2dstate");
        BoundaryRegistrar<TPShocktube2DSolidBoundary> b41("tpshocktube2dsolidboundary");

        //2d RayleighTaylor instability test with periodic domain
        BoundaryRegistrar<RayleighTaylorPeriodic2DBoundary> b42("rayleightaylor2dperiodicboundary");

	//nozzle with inflow boundary condition
        GeometryRegistrar<Nozzle2DSimple> r43("nozzle2dsimple");
	GeometryRegistrar<Nozzle2D> r44("nozzle2d");
        GeometryRegistrar<Nozzle2DComplete> r45("nozzle2dcomplete");
        GeometryRegistrar<Nozzle3D> r46("nozzle3d");
        StateRegistrar<NozzleState> s43("nozzlestate");
	BoundaryRegistrar<NozzleInflowBoundary> b43("nozzleinflowboundary");
        BoundaryRegistrar<Nozzle2DSimpleSolidBoundary> b44("nozzle2dsimplesolidboundary");
	BoundaryRegistrar<Nozzle2DSolidBoundary> b45("nozzle2dsolidboundary");
        BoundaryRegistrar<Nozzle3DInflowBoundary> b46("nozzle3dinflowboundary");
        BoundaryRegistrar<NozzleInflowFixPressureBoundary> b47("nozzleinflowfixpressureboundary");
        BoundaryRegistrar<Nozzle3DInflowFixPressureBoundary> b48("nozzle3dinflowfixpressureboundary");
        BoundaryRegistrar<Nozzle3DSimpleSolidBoundary> b49("nozzle3dsimplesolidboundary");
        GeometryRegistrar<Nozzle2DRothe> r47("nozzle2drothe");
        GeometryRegistrar<Nozzle3DRothe> r48("nozzle3drothe");
        StateRegistrar<NozzleRotheState> s44("nozzlerothestate");
        BoundaryRegistrar<Nozzle3DSolidBoundary> b50("nozzle3dsolidboundary");
        BoundaryRegistrar<NozzleOutflowBoundary> b51("nozzleoutflowboundary");
        BoundaryRegistrar<Nozzle3DOutflowBoundary> b52("nozzle3doutflowboundary");
        BoundaryRegistrar<Nozzle2DBNLSolidBoundary> b53("nozzle2dbnlsolidboundary");
        BoundaryRegistrar<Nozzle3DBNLSolidBoundary> b54("nozzle3dbnlsolidboundary");
        BoundaryRegistrar<Nozzle3DSolidRightBoundary> b55("nozzle3dsolidrightboundary");


}
