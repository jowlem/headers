/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas
 */

/** \file
 * A fluid constrained between a hot bottom wall (no-slip for the velocity) and a cold
 * top wall (no-slip for the velocity). The lateral walls are periodic. Under the
 * influence of gravity, convection rolls are formed. Thermal effects are modelled
 * by means of a Boussinesq approximation: the fluid is incompressible, and the influence
 * of the temperature is visible only through a body-force term, representing buoyancy
 * effects. The temperature field obeys an advection-diffusion equation.
 *
 * The simulation is first created in a fully symmetric manner. The symmetry is therefore
 * not spontaneously broken; while the temperature drops linearly between the hot and
 * and cold wall, the convection rolls fail to appear at this point. In a second stage, a
 * random noise is added to trigger the instability.
 *
 * This application is technically a bit more advanced than the other ones, because it
 * illustrates the concept of data processors. In the present case, they are used to
 * create the initial condition, and to trigger the instability.
 **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include "headersJo.hh"
#include "headersJo.h"

#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor

#define ADYNAMICS AdvectionDiffusionRLBdynamics
#define NSDYNAMICS GuoExternalForceCompleteRegularizedBGKdynamics

/// Initialization of the volume fraction field.
template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
struct IniVolFracProcessor3D : public BoxProcessingFunctional3D_L<T,adDescriptor>
{
    IniVolFracProcessor3D(RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,adDescriptor>& adLattice)
    {
        Dot3D absoluteOffset = adLattice.getLocation();


      T nz=parameters.getNz();
      T nx=parameters.getNx();
      T ny=parameters.getNy();
      T up=0.135;
      T low=0.25;

      T kx;
      T ky;
            for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                      plint absoluteX = absoluteOffset.x + iX;
                      plint absoluteY = absoluteOffset.y + iY;
                        plint absoluteZ = absoluteOffset.z + iZ;

                        	T VolFrac;
                          T rand_val=(double)rand()/RAND_MAX;
                          T amplitude=0.5;

                          kx=(2*3.1416)/(nx/1);
                          ky=(2*3.1416)/(ny/2);


    			if (absoluteZ==(nz-1)-(int)((up/(up+low))*nz)) {
            VolFrac=1-(rand_val*amplitude*(1+cos(absoluteX*kx)*cos(absoluteY*ky)));
          }
          else
    			if(absoluteZ<=(nz-1)-(int)((up/(up+low))*nz)-1)
    				VolFrac=0.0;
    			else VolFrac=1.;



                    Array<T,adDescriptor<T>::d> jEq(0., 0., 0.);
                    adLattice.get(iX,iY,iZ).defineDensity(VolFrac);
                    iniCellAtEquilibrium(adLattice.get(iX,iY,iZ), VolFrac, jEq);
                }
            }
        }
    }

    virtual IniVolFracProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniVolFracProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
private :
    RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

/// Initialization of the density field.
template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
struct IniDensityProcessor3D : public BoxProcessingFunctional3D_L<T,adDescriptor>
{
    IniDensityProcessor3D(RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,adDescriptor>& deLattice)
    {
        Dot3D absoluteOffset = deLattice.getLocation();


        T nz=parameters.getNz();

        T up=0.135;
        T low=0.25;

              for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                  for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                      for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                          plint absoluteZ = absoluteOffset.z + iZ;

                          	T dens;


      			if(absoluteZ<=(nz-1)-(int)((up/(up+low))*nz))
      				dens= 1.0084;
      			else dens=1.;



                    Array<T,adDescriptor<T>::d> jEq(0., 0., 0.);
                    deLattice.get(iX,iY,iZ).defineDensity(dens);
                    iniCellAtEquilibrium(deLattice.get(iX,iY,iZ), dens, jEq);
                }
            }
        }
    }
    virtual IniDensityProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniDensityProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
private :
    RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters;
};


void ExpSetup (
        MultiBlockLattice3D<T, NSDESCRIPTOR>& nsLattice,
        MultiBlockLattice3D<T, ADESCRIPTOR>& adLattice,
        MultiBlockLattice3D<T, ADESCRIPTOR>& deLattice,
        OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& nsBoundaryCondition,
        /*OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>& adBoundaryCondition,
        OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>& deBoundaryCondition,*/
        RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> &parameters )
{
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    Box3D upper(1, nx, 1, ny, 2*(nz-1)/3+1, nz-1);
    Box3D lower(1, nx, 1, ny, 0, 2*(nz-1)/3);

    Box3D bottom(0,nx-1,0,ny-1,0,0);
    Box3D top(0,nx-1,0,ny-1,nz-1,nz-1);

    Box3D front(nx-1,nx-1,0,ny-1,1,nz-2);
    Box3D back(0,0,0,ny-1,1,nz-2);

    Box3D left(0,nx-1,0,0,1,nz-2);
    Box3D right(0,nx-1,ny-1,ny-1,1,nz-2);

    T rho0f = 0.;


    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, top, boundary::dirichlet );
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, bottom, boundary::dirichlet );
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, front, boundary::dirichlet );
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, back, boundary::dirichlet );
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, left, boundary::dirichlet );
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries (nsLattice, right, boundary::dirichlet );


    //defineDynamics(nsLattice, nsLattice.getBoundingBox(), new BounceBack<T, NSDESCRIPTOR> (rho0f) );
  //  defineDynamics(adLattice, adLattice.getBoundingBox(), new BounceBack<T, ADESCRIPTOR> (rho0f) );
  //  defineDynamics(deLattice, deLattice.getBoundingBox(), new BounceBack<T, ADESCRIPTOR> (rho0f) );


    //defineDynamics(nsLattice, top, new BounceBack<T, NSDESCRIPTOR>(rho0f) );
    //defineDynamics(nsLattice, bottom, new BounceBack<T, NSDESCRIPTOR>(rho0f) );
    //defineDynamics(nsLattice, right, new BounceBack<T, NSDESCRIPTOR>(rho0f) );
    //defineDynamics(nsLattice, left, new BounceBack<T, NSDESCRIPTOR>(rho0f) );
    //defineDynamics(nsLattice, front, new BounceBack<T, NSDESCRIPTOR>(rho0f) );
    //defineDynamics(nsLattice, back, new BounceBack<T, NSDESCRIPTOR>(rho0f) );


    defineDynamics(adLattice, top, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(adLattice, bottom, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(adLattice, right, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(adLattice, left, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(adLattice, front, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(adLattice, back, new BounceBack<T, ADESCRIPTOR>(rho0f) );

    defineDynamics(deLattice, top, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(deLattice, bottom, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(deLattice, right, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(deLattice, left, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(deLattice, front, new BounceBack<T, ADESCRIPTOR>(rho0f) );
    defineDynamics(deLattice, back, new BounceBack<T, ADESCRIPTOR>(rho0f) );





      //adBoundaryCondition.setTemperatureConditionOnBlockBoundaries(adLattice, adLattice.getBoundingBox(), boundary::neumann );
      //adBoundaryCondition.setTemperatureConditionOnBlockBoundaries(adLattice, bottom,  boundary::neumann );

      //setBoundaryDensity(adLattice, adLattice.getBoundingBox(), 1e-4 );



    // Boundary conditions for DE

      //deBoundaryCondition.addTemperatureBoundary0N (front, deLattice, boundary::freeslip);
      //deBoundaryCondition.addTemperatureBoundary0P (back, deLattice, boundary::freeslip);
      //deBoundaryCondition.addTemperatureBoundary1N (left, deLattice, boundary::freeslip);
      //deBoundaryCondition.addTemperatureBoundary1P (right, deLattice, boundary::freeslip);
      //deBoundaryCondition.addTemperatureBoundary2N (bottom, deLattice, boundary::freeslip);
      //deBoundaryCondition.addTemperatureBoundary2P (top, deLattice, boundary::freeslip);


      //deBoundaryCondition.setTemperatureConditionOnBlockBoundaries(deLattice, deLattice.getBoundingBox(), boundary::neumann );
      //deBoundaryCondition.setTemperatureConditionOnBlockBoundaries(deLattice, bottom, boundary::neumann );
      //setBoundaryDensity(deLattice, deLattice.getBoundingBox(), 1. );


    initializeAtEquilibrium(nsLattice, nsLattice.getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    //initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    //initializeAtEquilibrium(deLattice, deLattice.getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );


    applyProcessingFunctional(
            new IniVolFracProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters),
            adLattice.getBoundingBox(), adLattice );

    applyProcessingFunctional(
            new IniDensityProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters),
            deLattice.getBoundingBox(), deLattice );

    deLattice.initialize();
    adLattice.initialize();
    nsLattice.initialize();

}



void writeVTK(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& adLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& deLattice,
              RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(nsLattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(nsLattice), "velocity", dx/dt);
    // Temperature is the order-0 moment of the advection-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    vtkOut.writeData<float>(*computeDensity(adLattice), "VolFrac", (T)1);
    //vtkOut.writeData<3, float>(*computeVelocity(adLattice), "ADVel", dx/dt);
}

void writeGif(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& adLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& deLattice, int iT)
{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    const plint nz = nsLattice.getNz();
    Box3D slice((nx-1)/2, (nx-1)/2, 0, ny-1, 0, nz-1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(nsLattice, slice),
                               imSize, imSize);
    // Temperature is the order-0 moment of the adveadction-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(createFileName("VolFrac", iT, 6),
                               *computeDensity(adLattice, slice),
                               imSize, imSize);

    imageWriter.writeScaledGif(createFileName("Density", iT, 6),
                                *computeDensity(deLattice, slice),
                                imSize, imSize);
}

plint particleTimeFactor = 1;
T particleProbabilityPerCell = 1;   // Probability of injecting a particle into an injection cell at each time step.
T cutOffSpeedSqr = 1.e-3; // Criterion to eliminate particles with very small velocity.


class BoxInjection {
public:
    BoxInjection ()
    { }
    bool operator()(Array<T,3> const& pos) const {
        return (true);
    }

};

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::timer("simTime").start();




    const T lx  = 0.075;
    const T ly  = 0.303;
    const T lz  = 0.385;
    //const T uLb  = 7.8987e-6;
    const T uCar = 16e-3;
    const T Gr = 1e6;
    //const T Kappa = 1e-6;
    const T Ri = 11.1;

    const plint resolution = 75;

    global::directories().setOutputDir("./tmp/");



    RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters (
            Ri,
            Gr,
            uCar,
            resolution,
            lx,
            ly,
            lz );


    T rho0=1.0;
    T rhoP=2.55;
    T tunedg=9.81*parameters.getLatticeGravity();



    writeLogFile(parameters,"palabos.log");


    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    Box3D upper(1, nx, 1, ny, 2*(nz-1)/3+1, nz);
    Box3D absorb(1,nx,1,ny,1,1);

    T nsOmega = parameters.getSolventOmega();
    T adOmega = parameters.getTemperatureOmega();


    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice  (
            nx,ny,nz,new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega) );


    MultiBlockLattice3D<T, ADESCRIPTOR> adLattice (
            nx,ny,nz,new ADYNAMICS<T, ADESCRIPTOR>(adOmega) );

    OnLatticeBoundaryCondition3D<T, NSDESCRIPTOR>* nsBoundaryCondition = createLocalBoundaryCondition3D<T, NSDESCRIPTOR>();


    MultiBlockLattice3D<T, ADESCRIPTOR> deLattice (
            nx,ny,nz,new ADYNAMICS<T, ADESCRIPTOR>(adOmega) );



    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);
    deLattice.toggleInternalStatistics(false);

    ExpSetup(nsLattice, adLattice, deLattice, *nsBoundaryCondition, /*adBoundaryCondition, *deBoundaryCondition,*/ parameters);

      Array<T,NSDESCRIPTOR<T>::d> forceOrientation(T(),T(),(T)1);



      integrateProcessingFunctional (
              new CouplingProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR> (),
                  nsLattice.getBoundingBox(),
                  nsLattice, deLattice, 0);

      integrateProcessingFunctional (
            new DVFCouplingProcessor3D<T,ADESCRIPTOR,ADESCRIPTOR> (
                rhoP, (parameters.getDeltaT()/parameters.getDeltaX())),
                deLattice.getBoundingBox(),
                deLattice, adLattice, 1 );


    integrateProcessingFunctional (
            new BuoyanTermProcessor3Dv2<T,NSDESCRIPTOR,ADESCRIPTOR> (
                tunedg, rho0,
                rhoP, forceOrientation),
            nsLattice.getBoundingBox(),
            nsLattice, adLattice, 1 );

            //particles


                      MultiParticleField3D<DenseParticleField3D<T,ADESCRIPTOR> >* particles=0;
                      particles = new MultiParticleField3D<DenseParticleField3D<T,ADESCRIPTOR> > (
                          adLattice.getMultiBlockManagement(),
                          defaultMultiBlockPolicy3D().getCombinedStatistics() );

                      std::vector<MultiBlock3D*> particleArg;
                      particleArg.push_back(particles);

                      std::vector<MultiBlock3D*> particleFluidArg;
                      particleFluidArg.push_back(particles);
                      particleFluidArg.push_back(&adLattice);

                      // Functional that advances the particles to their new position at each
                      //   predefined time step.
                      integrateProcessingFunctional (
                              new AdvanceParticlesEveryWhereFunctional3D<T,ADESCRIPTOR>(cutOffSpeedSqr),
                              adLattice.getBoundingBox(), particleArg, 0);
                      // Functional that assigns the particle velocity according to the particle's
                      //   position in the fluid.
                      integrateProcessingFunctional (
                              new FluidToParticleCoupling3D<T,ADESCRIPTOR>((T)particleTimeFactor),
                              adLattice.getBoundingBox(), particleFluidArg, 1 );

                    // Definition of simple mass-less particles.
                    Particle3D<T,ADESCRIPTOR>* particleTemplate=0;
                    particleTemplate = new PointParticle3D<T,ADESCRIPTOR>(0, Array<T,3>(0.,0.,0.), Array<T,3>(0.,0.,0.));


                    integrateProcessingFunctional (
                            new AnalyticalInjectRandomParticlesFunctional3D<T,ADESCRIPTOR,BoxInjection> (
                                particleTemplate, particleProbabilityPerCell, BoxInjection() ),
                            upper, particleArg, 0 );

                    integrateProcessingFunctional (
                                    new AbsorbParticlesFunctional3D<T,ADESCRIPTOR>,
                                    absorb, particleArg, 0 );

                    particles->executeInternalProcessors();




    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for ExpSetup:" << tIni << endl;
    global::timer("s100imTime").start();

    plint evalTime =5000;
    plint iT = 0;
    plint maxT = 5000000;
    //plint statIter = 1;
    plint saveIter = 5000;
    util::ValueTracer<T> converge((T)1,(T)100,1.0e-3);
    //bool convergedOnce = false;

    // Main loop over time iterations.
    for (iT = 0; iT <= maxT; ++iT)
    {
        if (iT == (evalTime))
        {
            T tEval = global::timer("simTime").stop();
            T remainTime = (tEval - tIni) / (T)evalTime * (T)maxT/(T)3600;
            global::timer("simTime").start();
            pcout << "Remaining " << (plint)remainTime << " hours, and ";
            pcout << (plint)((T)60*(remainTime - (T)((plint)remainTime))+0.5) << " minutes." << endl;
        }
//        if (iT == statIter)
//        {
//
//		Kappa=0;
//            int zDirection = 2;
//            T nusselt = computeNusseltNumber (
//                            nsLattice, adLattice,
//                            nsLattice.getBoundingBox(),
//                            zDirection, parameters.getDeltaX(),
//                            param125eters.getLatticeKappa(), parameters.getDeltaTemperature());
//            converge.takeValue(nusselt,true);
//        }
//        if (converge.hasConverged())
//        {
//            if (!convergedOnce)
//            {
//                convergedOnce = true;
//                converge.resetValues();
//                converge.setEpsilon(1.0e-14);
//                applyProcessingFunctional(
//                    new PerturbVolFracityProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters),
//                    adLattice.getBoundingBox(),adLattice);
//                pcout << "Intermetiate convergence.\n";
//            }
//            else
//            {
//                pcout << "Simulation is over.\n";
//                break;
//            }
//        }


        if (iT % saveIter == 0)
        {

          pcout << "Number of particles in the tank: "
                << countParticles(*particles, particles->getBoundingBox()) << std::endl;
            pcout << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
            writeVTK(nsLattice, adLattice,deLattice, parameters, iT);
            //pcout << "The NS average energy is " << computeAverageEnergy(nsLattice) << endl;
            //pcout << "The AD average energy is " << computeAverageEnergy(adLattice) << endl;

            pcout << iT << " : Writing gif." << endl;
            writeGif(nsLattice,adLattice,deLattice, iT);
            pcout << "Write visualization files." << std::endl;
            VtkImageOutput3D<T> vtkOut("volume", 1.);
            pcout << "Write particle output file." << endl;
            writeAsciiParticlePos(*particles, "particle_positions.dat");
            writeParticleVtk(*particles, "particles.vtk");
        }

        // Lattice Boltzmann iteration step.
        adLattice.collideAndStream();
        deLattice.collideAndStream();
        nsLattice.collideAndStream();
    }

    writeGif(nsLattice,adLattice,deLattice, iT);

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd-tIni;
    T nx100 = nsLattice.getNx()/(T)100;
    T ny100 = nsLattice.getNy()/(T)100;
    T nz100 = nsLattice.getNz()/(T)100;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << nx100*ny100*nz100*(T)iT/totalTime << endl;

    delete particles;
    return 0;
}
