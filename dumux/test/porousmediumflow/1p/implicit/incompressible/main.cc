// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OnePTests
 * \brief Test for the one-phase CC model
 */

#include <config.h>

#include <ctime>
#include <iostream>

// Support for quad precision has to be included before any other Dune module:
#include <dumux/common/quad.hh>

#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/pointcloudvtkwriter.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/mpfa/scvgradients.hh>

#include <dumux/assembly/fvassembler.hh>

#include "problem.hh"

//! Function to write out the scv-wise velocities (overload for mpfa)
template<class FVGridGeometry, class GridVariables, class Sol,
         std::enable_if_t<FVGridGeometry::discMethod == Dumux::DiscretizationMethod::ccmpfa, int> = 0>
void writeMpfaVelocities(const FVGridGeometry& fvGridGeometry,
                         const GridVariables& gridVariables,
                         const Sol& x)
{
    using Scalar = typename GridVariables::Scalar;
    using GlobalPos = typename FVGridGeometry::SubControlVolume::GlobalPosition;

    const auto velocities = Dumux::CCMpfaScvGradients::computeVelocities(fvGridGeometry, gridVariables, x, /*phaseIdx*/0);
    Dumux::PointCloudVtkWriter<Scalar, GlobalPos> writer(velocities.first);
    writer.addPointData(velocities.second, "velocity (m/s)");
    writer.write("mpfa_scv_velocities");
}

//! Function to write out the scv-wise velocities (overload for NOT mpfa)
template<class FVGridGeometry, class GridVariables, class Sol,
         std::enable_if_t<FVGridGeometry::discMethod != Dumux::DiscretizationMethod::ccmpfa, int> = 0>
void writeMpfaVelocities(const FVGridGeometry& fvGridGeometry,
                         const GridVariables& gridVariables,
                         const Sol& x)
{}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // make assemble and attach linear system
    using Assembler = FVAssembler<TypeTag, NUMDIFFMETHOD>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables);
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
    assembler->setLinearSystem(A, r);

    Dune::Timer timer;
    // assemble the local jacobian and the residual
    Dune::Timer assemblyTimer;
    if (mpiHelper.rank() == 0) std::cout << "Assembling linear system ..." << std::flush;
    assembler->assembleJacobianAndResidual(x);
    assemblyTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;

    // we solve Ax = -r to save update and copy
    (*r) *= -1.0;

    // solve the linear system
    Dune::Timer solverTimer;
    using LinearSolver = SSORCGBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    if (mpiHelper.rank() == 0) std::cout << "Solving linear system using " + linearSolver->name() + "..." << std::flush;
    linearSolver->solve(*A, x, *r);
    solverTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // the grid variables need to be up to date for subsequent output
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    gridVariables->update(x);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;

    // output result to vtk
    vtkWriter.write(1.0);

    timer.stop();

    const bool checkIsConstantVelocity = getParam<bool>("Problem.CheckIsConstantVelocity", false);
    if(checkIsConstantVelocity)
    {
        // instantiate the velocity output
        VelocityOutput velocityOutput(*gridVariables);
        using VelocityVector = typename VelocityOutput::VelocityVector;
        VelocityVector velocity;

        constexpr bool isBox = FVGridGeometry::discMethod == Dumux::DiscretizationMethod::box;
        constexpr int dimWorld = FVGridGeometry::GridView::dimensionworld;
        const auto numCells = leafGridView.size(0);
        const auto numDofs = fvGridGeometry->numDofs();
        auto numVelocities = (isBox && dimWorld == 1) ? numCells : numDofs;

        velocity.resize(numVelocities);

        const auto exactVel = problem->velocity();

        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);

            auto fvGeometry = localView(*fvGridGeometry);
            auto elemVolVars = localView(gridVariables->curGridVolVars());

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, x);

            velocityOutput.calculateVelocity(velocity, elemVolVars, fvGeometry, element, 0);

            using Scalar = Grid::ctype;
            // the y-component of the velocity should be exactly reproduced
            // the x-component should be zero
            // use a relative comparison for the y-component and an absolute one for the x-component
            if(Dune::FloatCmp::ne(velocity[eIdx][dimWorld-1], exactVel[dimWorld-1], /*eps*/1e-8) ||
               Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(velocity[eIdx][0], exactVel[0], /*eps*/1e-10))
                    DUNE_THROW(Dune::InvalidStateException, "Velocity is not exactly reproduced");
        }
    }

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    // For the mpfa test, write out the gradients in the scv centers
    if (getParam<bool>("IO.WriteMpfaVelocities", false))
        writeMpfaVelocities(*fvGridGeometry, *gridVariables, x);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

}
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
