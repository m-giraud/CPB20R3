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
 * \ingroup KEpsilonModel
 * \copydoc Dumux::KEpsilonFluxVariablesImpl
 */
#ifndef DUMUX_KEPSILON_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_KEPSILON_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/fluxvariables.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseFluxVariables, DiscretizationMethod discMethod>
class KEpsilonFluxVariablesImpl;

/*!
 * \ingroup KEpsilonModel
 * \brief The flux variables class for the k-epsilon model using the staggered grid discretization.
 */
template<class TypeTag, class BaseFluxVariables>
class KEpsilonFluxVariablesImpl<TypeTag, BaseFluxVariables, DiscretizationMethod::staggered>
: public BaseFluxVariables
{
    using ParentType = BaseFluxVariables;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;

    static constexpr int turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx - ModelTraits::dim();
    static constexpr int dissipationEqIdx = Indices::dissipationEqIdx - ModelTraits::dim();

public:

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element &element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace &scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux = ParentType::computeMassFlux(problem, element, fvGeometry,
                                                                      elemVolVars, elemFaceVars, scvf, fluxVarsCache);

        // calculate advective flux
        auto upwindTermK = [](const auto& volVars)
        {
            return volVars.turbulentKineticEnergy();
        };
        auto upwindTermEpsilon = [](const auto& volVars)
        {
            return volVars.dissipation();
        };

        flux[turbulentKineticEnergyEqIdx]
            = ParentType::advectiveFluxForCellCenter(problem, elemVolVars, elemFaceVars, scvf, upwindTermK);
        flux[dissipationEqIdx ]
            = ParentType::advectiveFluxForCellCenter(problem, elemVolVars, elemFaceVars, scvf, upwindTermEpsilon);

        // calculate diffusive flux
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // effective diffusion coefficients
        Scalar insideCoeff_k = insideVolVars.kinematicEddyViscosity() / insideVolVars.sigmaK();
        Scalar outsideCoeff_k = outsideVolVars.kinematicEddyViscosity() / outsideVolVars.sigmaK();
        Scalar insideCoeff_e = insideVolVars.kinematicEddyViscosity() / insideVolVars.sigmaEpsilon();
        Scalar outsideCoeff_e = outsideVolVars.kinematicEddyViscosity() / outsideVolVars.sigmaEpsilon();
        static const auto kEpsilonEnableKinematicViscosity_
            = getParamFromGroup<bool>(problem.paramGroup(), "KEpsilon.EnableKinematicViscosity", true);
        if (kEpsilonEnableKinematicViscosity_)
        {
            insideCoeff_k += insideVolVars.kinematicViscosity();
            outsideCoeff_k += outsideVolVars.kinematicViscosity();
            insideCoeff_e += insideVolVars.kinematicViscosity();
            outsideCoeff_e += outsideVolVars.kinematicViscosity();
        }

        // scale by extrusion factor
        insideCoeff_k *= insideVolVars.extrusionFactor();
        outsideCoeff_k *= outsideVolVars.extrusionFactor();
        insideCoeff_e *= insideVolVars.extrusionFactor();
        outsideCoeff_e *= outsideVolVars.extrusionFactor();

        // average and distance
        // is more stable with simple/unweighted arithmetic mean
        Scalar coeff_k = arithmeticMean(insideCoeff_k, outsideCoeff_k);
        Scalar coeff_e = arithmeticMean(insideCoeff_e, outsideCoeff_e);
        Scalar distance = 0.0;

        // adapt boundary handling
        if (scvf.boundary())
        {
            distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            coeff_k = insideCoeff_k;
            coeff_e = insideCoeff_e;
        }
        else
        {
            distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();
        }

        const auto bcTypes = problem.boundaryTypes(element, scvf);
        if (!(scvf.boundary() && (bcTypes.isOutflow(Indices::turbulentKineticEnergyEqIdx)
                                  || bcTypes.isSymmetry()
                                  || problem.isOnWall(scvf))))
        {
            if (!(insideVolVars.isMatchingPoint() && outsideVolVars.isMatchingPoint())
                || !(insideVolVars.isMatchingPoint() && outsideVolVars.inNearWallRegion())
                || !(insideVolVars.inNearWallRegion() && outsideVolVars.isMatchingPoint()))
            {
                flux[turbulentKineticEnergyEqIdx]
                    += coeff_k / distance
                       * (insideVolVars.turbulentKineticEnergy() - outsideVolVars.turbulentKineticEnergy())
                       * scvf.area();
            }
        }
        if (!(scvf.boundary() && (bcTypes.isOutflow(Indices::dissipationEqIdx)
                                  || bcTypes.isSymmetry())))
        {
            flux[dissipationEqIdx]
                += coeff_e / distance
                   * (insideVolVars.dissipation() - outsideVolVars.dissipation())
                   * scvf.area();
        }
        return flux;
    }

    /*!
     * \brief Returns the momentum flux over all staggered faces.
     */
    FacePrimaryVariables computeMomentumFlux(const Problem& problem,
                                             const Element& element,
                                             const SubControlVolumeFace& scvf,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        return ParentType::computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars)
               + ParentType::computeLateralMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars)
               + 2.0 / ModelTraits::dim() * insideVolVars.density() * insideVolVars.turbulentKineticEnergy()
                 * scvf.area() * scvf.directionSign() * insideVolVars.extrusionFactor();
    }
};

} // end namespace

#endif
