#include "boundary_nozzle.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Start of Nozzle2DBoundary
////////////////////////////////////////////////////////////////////////////////////////

NozzleInflowBoundary::NozzleInflowBoundary():left(-1e-3), right(0), radius(0.75e-3), Uinflow(500), Pinflow(1e+5), Vinflow(12.5) {
}

int NozzleInflowBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *pressure = m_pParticleData->m_vPressure;
        double *velocityU = m_pParticleData->m_vVelocityU;
        double *velocityV = m_pParticleData->m_vVelocityV;
        double *volume = m_pParticleData->m_vVolume;
        double *volumeold = m_pParticleData->m_vVolumeOld;
	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
	double *mass = m_pParticleData->m_vMass;
	double *sound = m_pParticleData->m_vSoundSpeed;
	bool *LeftInflow = m_pParticleData->m_bLeftInflow;

	if(inflowEndIndex-fluidEndIndex==0) // no inflow particle, only happens at initialization
	{
		double tx=right-0.25*dx,ty=-radius+0.25*sqrt(3.0)*dx;
		while(ty<radius)
		{
			if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        x[inflowEndIndex]=tx-Uinflow*dt;
                        y[inflowEndIndex]=ty;
                        pressure[inflowEndIndex]=Pinflow;
                        velocityU[inflowEndIndex]=Uinflow;
                        velocityV[inflowEndIndex]=0;
                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
			localParSpacing[inflowEndIndex]=dx;
                        mass[inflowEndIndex]=sqrt(3)*0.5*dx*dx/Vinflow;
			sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                        inflowEndIndex++;
			ty+=sqrt(3.0)*dx;
		}
		tx-=0.5*dx;
		ty-=0.5*sqrt(3.0)*dx;
		if(ty>radius) ty-=sqrt(3.0)*dx;
                while(ty>-radius)
                {
                        if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        x[inflowEndIndex]=tx-Uinflow*dt;
                        y[inflowEndIndex]=ty;
                        pressure[inflowEndIndex]=Pinflow;
                        velocityU[inflowEndIndex]=Uinflow;
                        velocityV[inflowEndIndex]=0;
                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
                        localParSpacing[inflowEndIndex]=dx;
                        mass[inflowEndIndex]=sqrt(3)*0.5*dx*dx/Vinflow;
                        sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                        inflowEndIndex++;
                        ty-=sqrt(3.0)*dx;
                }
		
	}
		for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
		{
			x[index]+=Uinflow*dt;
                        if(LeftInflow[index]==false&&x[index]>left+dx)
                        {
                                if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                                x[inflowEndIndex]=x[index]-dx-Uinflow*dt;
                                y[inflowEndIndex]=y[index];
                                pressure[inflowEndIndex]=pressure[index];
                                velocityU[inflowEndIndex]=velocityU[index];
                                velocityV[inflowEndIndex]=velocityV[index];
                                volumeold[inflowEndIndex]=volume[inflowEndIndex]=volume[index];
                                localParSpacing[inflowEndIndex]=dx;
                                mass[inflowEndIndex]=mass[index];
                                sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                                LeftInflow[inflowEndIndex]=false;
                                inflowEndIndex++;
                                LeftInflow[index]=true;
                        }
			if(x[index]>right) 
			{
				if(index>fluidEndIndex)
					m_pParticleData->swap(index,fluidEndIndex);
				x[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
				y[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
				fluidEndIndex++;
//				for(size_t j=fluidEndIndex;j<=index;j++)
//				{
//					x[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//					y[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//				}
//				fluidEndIndex=index+1;//mark this particle (and those before it) as fluid particle
				
			}
		}
	m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
	m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
	m_pParticleData->m_iBoundaryNum=m_pParticleData->m_iInflowNum=inflowEndIndex-fluidEndIndex;
	m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
	m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;
}

NozzleInflowFixPressureBoundary::NozzleInflowFixPressureBoundary():left(-1e-3), right(0), average_rightlimit(0.1e-3),radius(0.75e-3), Uinflow(0), OldUinflow(0), Pinflow(1e+5), Vinflow(0), Pinitial(1e+5), Vinitial(12.5) {
//NozzleInflowFixPressureBoundary::NozzleInflowFixPressureBoundary():left(-3e-3), right(0), average_rightlimit(1e-3),radius(8.3e-3), Uinflow(0), OldUinflow(0), Pinflow(0), Vinflow(0), Pinitial(473.86), Vinitial(1.0/0.0053198) {
}
int NozzleInflowFixPressureBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *pressure = m_pParticleData->m_vPressure;
        double *velocityU = m_pParticleData->m_vVelocityU;
        double *velocityV = m_pParticleData->m_vVelocityV;
        double *volume = m_pParticleData->m_vVolume;
        double *volumeold = m_pParticleData->m_vVolumeOld;
	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
	double *mass = m_pParticleData->m_vMass;
	double *sound = m_pParticleData->m_vSoundSpeed;
	bool *LeftInflow = m_pParticleData->m_bLeftInflow;

	double Ufluid=0,Vfluid=0,Pfluid=0;
	double gamma=1.4;
	int count=0;
	for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
	{
		if(x[index]<average_rightlimit)
		{
//			Ufluid+=velocityU[index];
//			Vfluid+=volume[index];
			Pfluid+=pressure[index];
			count++;
		}
	}
//	Ufluid/=count;
//	Vfluid/=count;
	Pfluid/=count;
	OldUinflow=Uinflow;
//	Uinflow=Ufluid;
//	Vinflow=Vfluid*pow(Pfluid/Pinflow,1.0/gamma);
	Pinflow=Pfluid;
	Vinflow=Vinitial*pow(Pinitial/Pinflow,1.0/gamma);
	double v2=2.0*gamma/(gamma-1)*(Pinitial*Vinitial-Pinflow*Vinflow);
	if(v2>0)
		Uinflow=sqrt(v2);
	else
		Uinflow=0;

	if(inflowEndIndex-fluidEndIndex==0) // no inflow particle, only happens at initialization
	{
		double tx=right-0.25*dx,ty=-radius+0.25*sqrt(3.0)*dx;
		while(ty<radius)
		{
			if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        x[inflowEndIndex]=tx-0.5*(OldUinflow+Uinflow)*dt;
                        y[inflowEndIndex]=ty;
                        pressure[inflowEndIndex]=Pinflow;
                        velocityU[inflowEndIndex]=Uinflow;
                        velocityV[inflowEndIndex]=0;
                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
			localParSpacing[inflowEndIndex]=dx;
                        mass[inflowEndIndex]=sqrt(3)*0.5*dx*dx/Vinflow;
			sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                        inflowEndIndex++;
			ty+=sqrt(3.0)*dx;
		}
		tx-=0.5*dx;
		ty-=0.5*sqrt(3.0)*dx;
		if(ty>radius) ty-=sqrt(3.0)*dx;
                while(ty>-radius)
                {
                        if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        x[inflowEndIndex]=tx-0.5*(OldUinflow+Uinflow)*dt;
                        y[inflowEndIndex]=ty;
                        pressure[inflowEndIndex]=Pinflow;
                        velocityU[inflowEndIndex]=Uinflow;
                        velocityV[inflowEndIndex]=0;
                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
                        localParSpacing[inflowEndIndex]=dx;
                        mass[inflowEndIndex]=sqrt(3)*0.5*dx*dx/Vinflow;
                        sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                        inflowEndIndex++;
                        ty-=sqrt(3.0)*dx;
                }
		
	}
		for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
		{
			x[index]+=0.5*(OldUinflow+Uinflow)*dt;
                        velocityU[index]=Uinflow;
			pressure[index]=Pinflow;
                        volumeold[index]=volume[index]=Vinflow;
                        mass[index]=sqrt(3)*0.5*dx*dx/Vinflow;
                        sound[index]= m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
                        if(LeftInflow[index]==false&&x[index]>left+dx)
                        {
                                if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                                x[inflowEndIndex]=x[index]-dx-0.5*(OldUinflow+Uinflow)*dt;
                                y[inflowEndIndex]=y[index];
                                pressure[inflowEndIndex]=pressure[index];
                                velocityU[inflowEndIndex]=velocityU[index];
                                velocityV[inflowEndIndex]=velocityV[index];
                                volumeold[inflowEndIndex]=volume[inflowEndIndex]=volume[index];
                                localParSpacing[inflowEndIndex]=dx;
                                mass[inflowEndIndex]=mass[index];
                                sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                                LeftInflow[inflowEndIndex]=false;
                                inflowEndIndex++;
                                LeftInflow[index]=true;
                        }
			if(x[index]>right) 
			{
                                if(index>fluidEndIndex)
                                        m_pParticleData->swap(index,fluidEndIndex);
                                x[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                y[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                fluidEndIndex++;
//				for(size_t j=fluidEndIndex;j<=index;j++)
//				{
//					x[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//					y[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//				}
//				fluidEndIndex=index+1;//mark this particle (and those before it) as fluid particle
				
			}
		}
	m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
	m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
	m_pParticleData->m_iBoundaryNum=m_pParticleData->m_iInflowNum=inflowEndIndex-fluidEndIndex;
	m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
	m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;
}

Nozzle3DInflowBoundary::Nozzle3DInflowBoundary():left(-1e-3), right(0), radius(0.75e-3), Uinflow(1768), Pinflow(5.369e+4), Vinflow(12.5) {
}

int Nozzle3DInflowBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *z = m_pParticleData->m_vPositionZ;
        double *pressure = m_pParticleData->m_vPressure;
        double *velocityU = m_pParticleData->m_vVelocityU;
        double *velocityV = m_pParticleData->m_vVelocityV;
        double *velocityW = m_pParticleData->m_vVelocityW;
        double *volume = m_pParticleData->m_vVolume;
        double *volumeold = m_pParticleData->m_vVolumeOld;
	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
	double *mass = m_pParticleData->m_vMass;
	double *sound = m_pParticleData->m_vSoundSpeed;
	bool *LeftInflow = m_pParticleData->m_bLeftInflow;

	if(inflowEndIndex-fluidEndIndex==0) // no inflow particle, only happens at initialization
	{
		double tx=right-0.25*dx,ty=-radius+0.25*sqrt(3.0)/3.0*dx,tz=-radius+0.5*dx,cx,cy,cz;
		int ix=0,iy=0,iz=0;
		for(iz=0;1;iz++)
		{
			cz=tz+iz*1.0*sqrt(6.0)/3.0*dx;
			if(cz>radius) break;
			for(iy=0;1;iy++)
			{
				cy=ty+iy*0.5*sqrt(3.0)*dx+((iz+1)%2)*sqrt(3.0)/6.0*dx;
                                if(cy>radius) break;
				if(cy*cy+cz*cz>radius*radius) continue;
				for(ix=0;ix<1;ix++)
				{
					cx=tx-ix*dx-0.5*(iy%2)*dx;
					if(cx<left) break;
                        		if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        		x[inflowEndIndex]=cx;
                        		y[inflowEndIndex]=cy;
					z[inflowEndIndex]=cz;
        	        	        pressure[inflowEndIndex]=Pinflow;
                	        	velocityU[inflowEndIndex]=Uinflow;
                        		velocityV[inflowEndIndex]=0;
					velocityW[inflowEndIndex]=0;
		                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
		                        localParSpacing[inflowEndIndex]=dx;
		                        mass[inflowEndIndex]=sqrt(2)*0.5*dx*dx*dx/Vinflow;
		                        sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
		                        inflowEndIndex++;					
				}
			}

		}
	}
		for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
		{
			x[index]+=Uinflow*dt;
                        if(LeftInflow[index]==false&&x[index]>left+dx)
                        {
                                if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                                x[inflowEndIndex]=x[index]-dx-Uinflow*dt;
                                y[inflowEndIndex]=y[index];
                                z[inflowEndIndex]=z[index];
                                pressure[inflowEndIndex]=pressure[index];
                                velocityU[inflowEndIndex]=velocityU[index];
                                velocityV[inflowEndIndex]=velocityV[index];
                                velocityW[inflowEndIndex]=velocityW[index];
                                volumeold[inflowEndIndex]=volume[inflowEndIndex]=volume[index];
                                localParSpacing[inflowEndIndex]=dx;
                                mass[inflowEndIndex]=mass[index];
                                sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                                LeftInflow[inflowEndIndex]=false;
                                inflowEndIndex++;
                                LeftInflow[index]=true;
                        }
			if(x[index]>right) 
			{
                                if(index>fluidEndIndex)
                                        m_pParticleData->swap(index,fluidEndIndex);
                                x[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                y[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
				z[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                fluidEndIndex++;
//				for(size_t j=fluidEndIndex;j<=index;j++)
//				{
//					x[j]+=0*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//					y[j]+=0*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//                                        z[j]+=0*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//				}
//				fluidEndIndex=index+1;//mark this particle (and those before it) as fluid particle
				
			}
		}
	m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
	m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
	m_pParticleData->m_iBoundaryNum=m_pParticleData->m_iInflowNum=inflowEndIndex-fluidEndIndex;
	m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
	m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;
}

Nozzle3DInflowFixPressureBoundary::Nozzle3DInflowFixPressureBoundary():left(-0.5e-3), right(0), average_rightlimit(0.2e-3),radius(0.75e-3), Uinflow(0), OldUinflow(0), Pinflow(1e+5), Vinflow(0), Pinitial(1e+5), Vinitial(12.5) {
//Nozzle3DInflowFixPressureBoundary::Nozzle3DInflowFixPressureBoundary():left(-1e-3), right(0), average_rightlimit(1e-3),radius(8.3e-3), Uinflow(0), OldUinflow(0), Pinflow(0), Vinflow(0), Pinitial(473.86), Vinitial(1.0/0.0053198) {
}
int Nozzle3DInflowFixPressureBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
	size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        double *z = m_pParticleData->m_vPositionZ;
        double *pressure = m_pParticleData->m_vPressure;
        double *velocityU = m_pParticleData->m_vVelocityU;
        double *velocityV = m_pParticleData->m_vVelocityV;
        double *velocityW = m_pParticleData->m_vVelocityW;
        double *volume = m_pParticleData->m_vVolume;
        double *volumeold = m_pParticleData->m_vVolumeOld;
	double *localParSpacing = m_pParticleData->m_vLocalParSpacing;
	double *mass = m_pParticleData->m_vMass;
	double *sound = m_pParticleData->m_vSoundSpeed;
	bool *LeftInflow = m_pParticleData->m_bLeftInflow;

        double Ufluid=0,Vfluid=0,Pfluid=0;
        double gamma=1.4;
        int count=0;

        for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
        {
                if(x[index]<average_rightlimit)
                {
//                        Ufluid+=velocityU[index];
//                        Vfluid+=volume[index];
                        Pfluid+=pressure[index];
                        count++;
                }
        }
	if(count==0)
		std::cout<<"Warning: count==0"<<std::endl;
//        Ufluid/=count;
//        Vfluid/=count;
        Pfluid/=count;
        OldUinflow=Uinflow;
//        Uinflow=Ufluid;
//        Vinflow=Vfluid*pow(Pfluid/Pinflow,1.0/gamma);
        Pinflow=Pfluid;
        Vinflow=Vinitial*pow(Pinitial/Pinflow,1.0/gamma);
        double v2=2.0*gamma/(gamma-1)*(Pinitial*Vinitial-Pinflow*Vinflow);
        if(v2>0)
                Uinflow=sqrt(v2);
        else
                Uinflow=0;

	if(inflowEndIndex-fluidEndIndex==0) // no inflow particle, only happens at initialization
	{
		double tx=right-0.25*dx,ty=-radius+0.25*sqrt(3.0)/3.0*dx,tz=-radius+0.5*dx,cx,cy,cz;
		int ix=0,iy=0,iz=0;
		for(iz=0;1;iz++)
		{
			cz=tz+iz*1.0*sqrt(6.0)/3.0*dx;
			if(cz>radius) break;
			for(iy=0;1;iy++)
			{
				cy=ty+iy*0.5*sqrt(3.0)*dx+((iz+1)%2)*sqrt(3.0)/6.0*dx;
                                if(cy>radius) break;
				if(cy*cy+cz*cz>radius*radius) continue;
				for(ix=0;ix<1;ix++)
				{
					cx=tx-ix*dx-0.5*(iy%2)*dx;
					if(cx<left) break;
                        		if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                        		x[inflowEndIndex]=cx;
                        		y[inflowEndIndex]=cy;
					z[inflowEndIndex]=cz;
        	        	        pressure[inflowEndIndex]=Pinflow;
                	        	velocityU[inflowEndIndex]=Uinflow;
                        		velocityV[inflowEndIndex]=0;
					velocityW[inflowEndIndex]=0;
		                        volumeold[inflowEndIndex]=volume[inflowEndIndex]=Vinflow;
		                        localParSpacing[inflowEndIndex]=dx;
		                        mass[inflowEndIndex]=sqrt(2)*0.5*dx*dx*dx/Vinflow;
		                        sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
		                        inflowEndIndex++;					
				}
			}

		}
	}
		for(size_t index=fluidEndIndex;index<inflowEndIndex;index++)
		{
                        x[index]+=0.5*(OldUinflow+Uinflow)*dt;
                        velocityU[index]=Uinflow;
                        volumeold[index]=volume[index]=Vinflow;
                        mass[index]=sqrt(2)*0.5*dx*dx*dx/Vinflow;
			pressure[index]=Pinflow;
                        sound[index]= m_pEOS->getSoundSpeed(pressure[index],1./volume[index]);
                        if(LeftInflow[index]==false&&x[index]>left+dx)
                        {
                                if(inflowEndIndex>=m_pParticleData->m_iCapacity) return 1;//too many
                                x[inflowEndIndex]=x[index]-dx-0.5*(OldUinflow+Uinflow)*dt;
                                y[inflowEndIndex]=y[index];
                                z[inflowEndIndex]=z[index];
                                pressure[inflowEndIndex]=pressure[index];
                                velocityU[inflowEndIndex]=velocityU[index];
                                velocityV[inflowEndIndex]=velocityV[index];
                                velocityW[inflowEndIndex]=velocityW[index];
                                volumeold[inflowEndIndex]=volume[inflowEndIndex]=volume[index];
                                localParSpacing[inflowEndIndex]=dx;
                                mass[inflowEndIndex]=mass[index];
                                sound[inflowEndIndex]= m_pEOS->getSoundSpeed(pressure[inflowEndIndex],1./volume[inflowEndIndex]);
                                LeftInflow[inflowEndIndex]=false;
                                inflowEndIndex++;
                                LeftInflow[index]=true;
                        }
			if(x[index]>right) 
			{
                                if(index>fluidEndIndex)
                                        m_pParticleData->swap(index,fluidEndIndex);
                                x[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                y[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                z[fluidEndIndex]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
                                fluidEndIndex++;
//				for(size_t j=fluidEndIndex;j<=index;j++)
//				{
//					x[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//					y[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//                                        z[j]+=1*0.25*dx*(2*((double)rand()/(double)RAND_MAX)-1);
//				}
//				fluidEndIndex=index+1;//mark this particle (and those before it) as fluid particle
				
			}
		}
	m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
	m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
	m_pParticleData->m_iBoundaryNum=m_pParticleData->m_iInflowNum=inflowEndIndex-fluidEndIndex;
	m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
	m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;
}

Nozzle2DSimpleSolidBoundary::Nozzle2DSimpleSolidBoundary():radius(8.3e-3), thickness(2e-3) {
        bo = radius-thickness;
}

int Nozzle2DSimpleSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	if((y<bo && y>-bo) || y>radius || y<-radius) return 0;

	xb.push_back(x);

	if(y>0) yb.push_back(2*radius-y);
	else	yb.push_back(-2*radius-y);

        pressureb.push_back(pressure);

	if(vy*y>0)
	{
		vxb.push_back(vx);
		vyb.push_back(-vy);
	}
	else
	{
		vxb.push_back(vx);
		vyb.push_back(vy);
	}
        return 1;
}

Nozzle2DBNLSolidBoundary::Nozzle2DBNLSolidBoundary():x1(0e-3),r1(0.75e-3),x2(4e-3),r2(0.25e-3),x3(7e-3),r3(0.5e-3),x4(8e-3),r4(1.65e-3), thickness(0.2e-3){}
int Nozzle2DBNLSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	int count=0;
	double l,lr,r,rr;
	int sign=1;
        if((x<x2) && ((y<r1 && y>r1-thickness) || (y>-r1 && y<-r1+thickness)))
	{
		count++;
	        xb.push_back(x);
	
        	if(y>0) yb.push_back(2*r1-y);
        	else    yb.push_back(-2*r1-y);

        	pressureb.push_back(pressure);

        	if(vy*y>0)
        	{
                	vxb.push_back(vx);
                	vyb.push_back(-vy);
        	}
        	else
        	{
                	vxb.push_back(vx);
                	vyb.push_back(vy);
        	}
	}
	if(x<x2&&x>x2-thickness && fabs(y)>r2 && fabs(y)<r1)
	{
                count++;
                xb.push_back(2*x2-x);
		yb.push_back(y);

                pressureb.push_back(pressure);

                if(vx>0)
                {
                        vxb.push_back(-vx);
                        vyb.push_back(vy);
                }
                else
                {
                        vxb.push_back(vx);
                        vyb.push_back(vy);
                }
	}
        if(count==2)
        {
                count++;
                xb.push_back(2*x2-x);

                if(y>0) yb.push_back(2*r1-y);
                else    yb.push_back(-2*r1-y);

                pressureb.push_back(pressure);

                if(vx+vy*y/fabs(y)>0)
                {
                        vxb.push_back(-vx);
                        vyb.push_back(-vy);
                }
                else
                {
                        vxb.push_back(vx);
                        vyb.push_back(vy);
                }
        }
	if(x<x2)
		return count;
        if(x<x3)
        {
                l=x2,r=x3,lr=r2,rr=r3;
        }
	else if(x<x4)
	{
		l=x3,r=x4,lr=r3,rr=r4;
	}
	else return count;
	if(y<0)
	{
		y=-y;
		vy=-vy;
		sign=-1;
	}
	double dis=((rr-lr)*x-(r-l)*y+r*lr-rr*l)/sqrt((rr-lr)*(rr-lr)+(r-l)*(r-l));
	if(dis<0 || dis>thickness) return count;
	double k=-(r-l)/(rr-lr);
	double dx=sqrt(4*dis*dis/(1+k*k));
	if(k<0) dx=-dx;
	if(x+dx<x2) return count;
	xb.push_back(x+dx);
	yb.push_back(sign*(y+dx*k));
	pressureb.push_back(pressure);
	double normal_vy=fabs(k)/sqrt(1+k*k);
	double normal_vx=normal_vy/k;
	double ip = ip2d(vx,vy,normal_vx,normal_vy);
	if(ip<=0) {//leaving boundary
		vxb.push_back(vx);
		vyb.push_back(vy);
        }
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                vyb.push_back(sign*(vy-2.*ip*normal_vy));
        }
	return count+1;
}

Nozzle2DSolidBoundary::Nozzle2DSolidBoundary():x0(-2e-3), r0(0.751e-3), x1(0+0e-3),r1(0.75e-3),x2(4e-3+0e-3),r2(0.25e-3),x3(7e-3+0e-3),r3(0.5e-3),x4(8e-3+0e-3),r4(1.65e-3), thickness(0.1e-3){}
//Nozzle2DSolidBoundary::Nozzle2DSolidBoundary():x0(0), r0(8.301e-3), x1(0+10e-3),r1(8.3e-3),x2(9.9593e-3+10e-3),r2(2.55e-3),x3(60.65e-3+10e-3),r3(21e-3),x4(60.65e-3+10e-3),r4(21e-3), thickness(2e-3){}
int Nozzle2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	double l,lr,r,rr;
	int sign=1;
	if(x<x1)
	{
		l=x0,r=x1,lr=r0,rr=r1;
	}
	else if(x<x2)
	{
                l=x1,r=x2,lr=r1,rr=r2;
	}
        else if(x<x3)
        {
                l=x2,r=x3,lr=r2,rr=r3;
        }
	else if(x<x4)
	{
		l=x3,r=x4,lr=r3,rr=r4;
	}
	else return 0;
	if(y<0)
	{
		y=-y;
		vy=-vy;
		sign=-1;
	}
	double dis=((rr-lr)*x-(r-l)*y+r*lr-rr*l)/sqrt((rr-lr)*(rr-lr)+(r-l)*(r-l));
	if(dis<0 || dis>thickness) return 0;
	double k=-(r-l)/(rr-lr);
	double dx=sqrt(4*dis*dis/(1+k*k));
	if(k<0) dx=-dx;
	xb.push_back(x+dx);
	yb.push_back(sign*(y+dx*k));
	pressureb.push_back(pressure);
	double normal_vy=fabs(k)/sqrt(1+k*k);
	double normal_vx=normal_vy/k;
	double ip = ip2d(vx,vy,normal_vx,normal_vy);
	if(ip<=0) {//leaving boundary
		vxb.push_back(vx);
		vyb.push_back(vy);
//                vxb.push_back(vx-ip*normal_vx);
//                vyb.push_back(sign*(vy-ip*normal_vy));       
        }
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                vyb.push_back(sign*(vy-2.*ip*normal_vy));
        }
	return 1;
}
//Nozzle3DSimpleSolidBoundary::Nozzle3DSimpleSolidBoundary():radius(0.75e-3), thickness(0.2e-3) {
Nozzle3DSimpleSolidBoundary::Nozzle3DSimpleSolidBoundary():radius(8.3e-3), thickness(2e-3) {
        bo = radius-thickness;
}

int Nozzle3DSimpleSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	double r=sqrt(y*y+z*z);
        if(r<bo || r>radius) return 0;

	double normaly=y/r,normalz=z/r;
	double ip=ip2d(vy,vz,normaly,normalz);
        xb.push_back(x);

	yb.push_back(normaly*(2*radius-r));
	zb.push_back(normalz*(2*radius-r));

        pressureb.push_back(pressure);

        if(ip<=0)
        {
                vxb.push_back(vx);
                vyb.push_back(vy);
		vzb.push_back(vz);
        }
        else
        {
                vxb.push_back(vx);
                vyb.push_back(vy-2*ip*normaly);
		vzb.push_back(vz-2*ip*normalz);
        }
        return 1;
}

Nozzle3DSolidRightBoundary::Nozzle3DSolidRightBoundary():x0(20e-3),r0(0.501e-3),thickness(0.4e-3){}

int Nozzle3DSolidRightBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {
	if(x<x0-thickness || x>x0)
		return 0;
	xb.push_back(2*x0-x);
	yb.push_back(y);
	zb.push_back(z);
	pressureb.push_back(pressure);
	vxb.push_back(vx);
	vyb.push_back(vy);
	vzb.push_back(vz);
	double yz=sqrt(y*y+z*z);
	if(yz<r0-thickness || yz>r0)
		return 1;
        double normaly=y/yz,normalz=z/yz;
        double ip=ip2d(vy,vz,normaly,normalz);
        xb.push_back(2*x0-x);

        yb.push_back(normaly*(2*r0-yz));
        zb.push_back(normalz*(2*r0-yz));

        pressureb.push_back(pressure);

        if(ip<=0)
        {
                vxb.push_back(vx);
                vyb.push_back(vy);
                vzb.push_back(vz);
        }
        else
        {
                vxb.push_back(vx);
                vyb.push_back(vy-2*ip*normaly);
                vzb.push_back(vz-2*ip*normalz);
        }
        return 2;	
}

//Nozzle3DSolidBoundary::Nozzle3DSolidBoundary():x0(-2e-3), r0(0.751e-3), x1(0+0e-3),r1(0.75e-3),x2(4e-3+0e-3),r2(0.25e-3),x3(7e-3+0e-3),r3(0.5e-3),x4(30e-3+0e-3),r4(0.501e-3), thickness(0.4e-3){}
Nozzle3DSolidBoundary::Nozzle3DSolidBoundary():x0(-2e-3), r0(0.751e-3), x1(0+0e-3),r1(0.75e-3),x2(4e-3+0e-3),r2(0.25e-3),x3(7e-3+0e-3),r3(0.5e-3),x4(8e-3+0e-3),r4(1.65e-3), thickness(0.4e-3){}
//Nozzle3DSolidBoundary::Nozzle3DSolidBoundary():x0(0.1), r0(8.301e-3), x1(0+10e-3),r1(8.3e-3),x2(9.9593e-3+10e-3),r2(2.55e-3),x3(60.65e-3+10e-3),r3(21e-3),x4(60.65e-3+10e-3),r4(21e-3), thickness(2e-3){}
int Nozzle3DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	double l,lr,r,rr;
	int sign=1;
        double yz=sqrt(y*y+z*z);
	if(x<x1)
	{
		l=x0,r=x1,lr=r0,rr=r1;
	}
	else if(x<x2)
	{
                l=x1,r=x2,lr=r1,rr=r2;
	}
        else if(x<x3 || yz<r3-(x-x3)*1.973)
//        else if(x<x3)
        {
                l=x2,r=x3,lr=r2,rr=r3;
        }
	else if(x<x4)
	{
		l=x3,r=x4,lr=r3,rr=r4;
	}
	else return 0;
	{

        double normaly=y/yz,normalz=z/yz;
	double ty=yz;

        double tvy=ip2d(vy,vz,normaly,normalz);
	double tvz=ip2d(vy,vz,-normalz,normaly);

	double dis=((rr-lr)*x-(r-l)*ty+r*lr-rr*l)/sqrt((rr-lr)*(rr-lr)+(r-l)*(r-l));
	if(dis<0 || dis>thickness) return 0;
	double k=-(r-l)/(rr-lr);
	double dx=sqrt(4*dis*dis/(1+k*k));
	if(k<0) dx=-dx;
	xb.push_back(x+dx);
	double tyb=sign*(ty+dx*k);
	yb.push_back(tyb*normaly);
	zb.push_back(tyb*normalz);
	pressureb.push_back(pressure);
	double normal_vy=fabs(k)/sqrt(1+k*k);
	double normal_vx=normal_vy/k;
	double ip = ip2d(vx,tvy,normal_vx,normal_vy);
	double tvyb;
	if(ip<=0) {//leaving boundary
		vxb.push_back(vx);
		tvyb=tvy;
//                vxb.push_back(vx-ip*normal_vx);
//                vyb.push_back(sign*(vy-ip*normal_vy));       
        }
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                tvyb=sign*(tvy-2.*ip*normal_vy);
        }
        vyb.push_back(normaly*tvyb-normalz*tvz);
        vzb.push_back(normalz*tvyb+normaly*tvz);
	}
	if(x>x1)
		return 1;
	else
	{
		l=x1,r=x2,lr=r1,rr=r2;
		sign=1;
        double yz=sqrt(y*y+z*z);

        double normaly=y/yz,normalz=z/yz;
        double ty=yz;

        double tvy=ip2d(vy,vz,normaly,normalz);
        double tvz=ip2d(vy,vz,-normalz,normaly);

        double dis=((rr-lr)*x-(r-l)*ty+r*lr-rr*l)/sqrt((rr-lr)*(rr-lr)+(r-l)*(r-l));
        if(dis<0 || dis>thickness) return 1;
        double k=-(r-l)/(rr-lr);
        double dx=sqrt(4*dis*dis/(1+k*k));
        if(k<0) dx=-dx;
	if(x+dx<x1)
		return 1;
        xb.push_back(x+dx);
        double tyb=sign*(ty+dx*k);
        yb.push_back(tyb*normaly);
        zb.push_back(tyb*normalz);
        pressureb.push_back(pressure);
        double normal_vy=fabs(k)/sqrt(1+k*k);
        double normal_vx=normal_vy/k;
        double ip = ip2d(vx,tvy,normal_vx,normal_vy);
        double tvyb;
        if(ip<=0) {//leaving boundary
                vxb.push_back(vx);
                tvyb=tvy;
	}
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                tvyb=sign*(tvy-2.*ip*normal_vy);
        }
        vyb.push_back(normaly*tvyb-normalz*tvz);
        vzb.push_back(normalz*tvyb+normaly*tvz);
	return 2;
	}
}

Nozzle3DBNLSolidBoundary::Nozzle3DBNLSolidBoundary():x1(0e-3),r1(0.75e-3),x2(4e-3),r2(0.25e-3),x3(7e-3),r3(0.5e-3),x4(8e-3),r4(1.65e-3), thickness(0.4e-3){}
int Nozzle3DBNLSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

	int count=0;
	double l,lr,r,rr;
	int sign=1;
        double yz=sqrt(y*y+z*z);

        double normaly=y/yz,normalz=z/yz;
        double ty=yz;

        double tvy=ip2d(vy,vz,normaly,normalz);
        double tvz=ip2d(vy,vz,-normalz,normaly);
        if(x<x2 && yz<r1 && yz>r1-thickness)
	{
		count++;
	        xb.push_back(x);
	        yb.push_back(normaly*(2*r1-yz));
        	zb.push_back(normalz*(2*r1-yz));	
        	pressureb.push_back(pressure);
	        if(tvy<=0)
        	{
                	vxb.push_back(vx);
                	vyb.push_back(vy);
                	vzb.push_back(vz);
        	}
        	else
        	{
                	vxb.push_back(vx);
	                vyb.push_back(vy-2*tvy*normaly);
        	        vzb.push_back(vz-2*tvy*normalz);
        	}
	}
	if(x<x2&&x>x2-thickness && yz>r2 && yz<r1 && (x2-x)/(yz-r2)<0.847)
//        if(x<x2&&x>x2-thickness && yz>r2 && yz<r1)
	{
                count++;
                xb.push_back(2*x2-x);
		yb.push_back(y);
		zb.push_back(z);
                pressureb.push_back(pressure);
                if(vx>0)
                {
                        vxb.push_back(-vx);
                        vyb.push_back(vy);
			vzb.push_back(vz);
                }
                else
                {
                        vxb.push_back(vx);
                        vyb.push_back(vy);
			vzb.push_back(vz);
                }
	}
        if(count==2)
        {
                count++;
                xb.push_back(2*x2-x);
		yb.push_back(normaly*(2*r1-yz));
                zb.push_back(normalz*(2*r1-yz));

                pressureb.push_back(pressure);

                if(vx+tvy>0)
                {
                        vxb.push_back(-vx);
                        vyb.push_back(vy-2*tvy*normaly);
                        vzb.push_back(vz-2*tvy*normalz);
                }
                else
                {
                        vxb.push_back(vx);
                        vyb.push_back(vy);
			vzb.push_back(vz);
                }
        }
	if(x<x2)
		return count;
        else if(x<x3 || yz<r3-(x-x3)*1.973)
//        else if(x<x3)
        {
                l=x2,r=x3,lr=r2,rr=r3;
        }
	else if(x<x4)
	{
		l=x3,r=x4,lr=r3,rr=r4;
	}
	else return count;

	double dis=((rr-lr)*x-(r-l)*ty+r*lr-rr*l)/sqrt((rr-lr)*(rr-lr)+(r-l)*(r-l));
	if(dis<0 || dis>thickness) return count;
	double k=-(r-l)/(rr-lr);
	double dx=sqrt(4*dis*dis/(1+k*k));
	if(k<0) dx=-dx;
	if(x+dx<x2) return count;
	if(x<x3 && (x+dx-x2)/(fabs(ty+dx*k)-r2)<0.847) return count;
	xb.push_back(x+dx);
	double tyb=sign*(ty+dx*k);
	yb.push_back(tyb*normaly);
	zb.push_back(tyb*normalz);
	pressureb.push_back(pressure);
	double normal_vy=fabs(k)/sqrt(1+k*k);
	double normal_vx=normal_vy/k;
	double ip = ip2d(vx,tvy,normal_vx,normal_vy);
	double tvyb;
	if(ip<=0) {//leaving boundary
		vxb.push_back(vx);
		tvyb=tvy;
//                vxb.push_back(vx-ip*normal_vx);
//                vyb.push_back(sign*(vy-ip*normal_vy));       
        }
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                tvyb=sign*(tvy-2.*ip*normal_vy);
        }
        vyb.push_back(normaly*tvyb-normalz*tvz);
        vzb.push_back(normalz*tvyb+normaly*tvz);
	return count+1;
}

NozzleOutflowBoundary::NozzleOutflowBoundary():xmin(-0.01),xmax(0.020),ymin(-0.010),ymax(0.010) {
}

int NozzleOutflowBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
        size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
        for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
	{
		if(x[index]<xmin || x[index]>xmax || y[index]<ymin || y[index]>ymax)
		{
//			std::cout<<"Delete:"<<index<<" "<<x[index]<<" "<<y[index]<<std::endl;
			if(index+1<fluidEndIndex)
				m_pParticleData->swap(index,fluidEndIndex-1);
			if(fluidEndIndex<inflowEndIndex)
				m_pParticleData->swap(fluidEndIndex-1,inflowEndIndex-1);
			fluidEndIndex--;
			inflowEndIndex--;
		}
	}
        m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
        m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
        m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
        m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;	
}


Nozzle3DOutflowBoundary::Nozzle3DOutflowBoundary():xmin(-0.01),xmax(0.020),rmax(0.010) {
}

int Nozzle3DOutflowBoundary::UpdateInflowBoundary(ParticleData *m_pParticleData, EOS* m_pEOS, double dt, double dx) {
        size_t fluidStartIndex = m_pParticleData->m_iFluidStartIndex;
        size_t fluidEndIndex = m_pParticleData->m_iFluidStartIndex + m_pParticleData->m_iFluidNum;
        size_t inflowEndIndex = fluidEndIndex + m_pParticleData->m_iInflowNum;
        double *x = m_pParticleData->m_vPositionX;
        double *y = m_pParticleData->m_vPositionY;
	double *z = m_pParticleData->m_vPositionZ;
        for(size_t index=fluidStartIndex;index<fluidEndIndex;index++)
	{
		double r2=y[index]*y[index]+z[index]*z[index];
		if(x[index]<xmin || x[index]>xmax || r2>rmax*rmax)
		{
//			std::cout<<"Delete:"<<index<<" "<<x[index]<<" "<<y[index]<<" "<<z[index]<<std::endl;
			if(index+1<fluidEndIndex)
			{
				m_pParticleData->swap(index,fluidEndIndex-1);
//				m_pParticleData->makezero(fluidEndIndex-1);
			}
			if(fluidEndIndex<inflowEndIndex)
			{
				m_pParticleData->swap(fluidEndIndex-1,inflowEndIndex-1);
//				m_pParticleData->makezero(inflowEndIndex-1);
			}
			fluidEndIndex--;
			inflowEndIndex--;
		}
	}
        m_pParticleData->m_iFluidNum=fluidEndIndex-fluidStartIndex;
        m_pParticleData->m_iBoundaryStartIndex=fluidEndIndex;
        m_pParticleData->m_iGhostStartIndex=inflowEndIndex;
        m_pParticleData->m_iTotalNum=inflowEndIndex-fluidStartIndex;
	return 0;	
}
