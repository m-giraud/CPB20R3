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
 * \ingroup FacetTests
 * \brief The problem for the lower-dimensional domain in the single-phase facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIMPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIMPROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePLowDimProblem;

namespace Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimBox { using InheritsFrom = std::tuple<OnePLowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = Dune::FoamGrid<1, 2>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = OnePSpatialParams< GetPropType<TypeTag, Properties::FVGridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};
} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief The lower-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePLowDimProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

public:
    OnePLowDimProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                      std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManagerPtr,
                      const std::string& paramGroup = "LowDim")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    , aperture_(getParam<Scalar>("Problem.FractureAperture"))
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Evaluates the exact solution at a given position.
    Scalar exact(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::cosh;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        return cos(x)*cosh(y);
    }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exact(globalPos)); }

    //! Sets the aperture as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return aperture_; }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar aperture_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
