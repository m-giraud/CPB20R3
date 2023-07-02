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
  * \ingroup FreeflowNCModel
  * \copydoc Dumux::FreeflowNCFluxVariablesImpl
  */
#ifndef DUMUX_FREEFLOW_NC_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_FREEFLOW_NC_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FreeflowNCFluxVariablesImpl;

/*!
 * \ingroup FreeflowNCModel
 * \brief The flux variables class for the multi-component free-flow model using the staggered grid discretization.
 */
template<class TypeTag>
class FreeflowNCFluxVariablesImpl<TypeTag, DiscretizationMethod::staggered>
: public NavierStokesFluxVariables<TypeTag>
{
    using ParentType = NavierStokesFluxVariables<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    static constexpr auto numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    template<class ElementVolumeVariables, class ElementFaceVariables, class FluxVariablesCache>
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace& scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto diffusiveFluxes = MolecularDiffusionType::flux(problem, element, fvGeometry, elemVolVars, scvf);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            auto upwindTerm = [compIdx](const auto& volVars)
            {
                const auto density = useMoles ? volVars.molarDensity() : volVars.density();
                const auto fraction =  useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
                return density * fraction;
            };

            flux[compIdx] = ParentType::advectiveFluxForCellCenter(problem, elemVolVars, elemFaceVars, scvf, upwindTerm);

            flux[compIdx] += useMoles ? diffusiveFluxes[compIdx] : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
        }

        // in case one balance is substituted by the total mass balance
        if (ModelTraits::replaceCompEqIdx() < numComponents)
        {
            flux[ModelTraits::replaceCompEqIdx()] = std::accumulate(flux.begin(), flux.end(), 0.0);
        }

        return flux;
    }
};

} // end namespace

#endif
