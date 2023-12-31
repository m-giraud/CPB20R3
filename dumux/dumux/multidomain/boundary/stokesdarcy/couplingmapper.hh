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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <memory>
#include <unordered_map>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension.
 */
template<class MDTraits>
class StokesDarcyCouplingMapper
{
    using Scalar = typename MDTraits::Scalar;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto cellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto faceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<2>::Index();


private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;

    struct ElementMapInfo
    {
        std::size_t eIdx;
        std::size_t scvfIdx;
        std::size_t flipScvfIdx;
    };

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    using CouplingManager = GetPropType<StokesTypeTag, Properties::CouplingManager>;

    static_assert(GetPropType<SubDomainTypeTag<stokesIdx>, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::staggered,
                  "The free flow domain must use the staggered discretization");

    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::cctpfa,
                  "The Darcy domain must use the CCTpfa discretization");
public:

    /*!
     * \brief Constructor
     */
    StokesDarcyCouplingMapper(const CouplingManager& couplingManager) : couplingManager_(couplingManager) {}

    /*!
     * \brief Main update routine
     */
    template<class Stencils>
    void computeCouplingMapsAndStencils(Stencils& darcyToStokesCellCenterStencils,
                                        Stencils& darcyToStokesFaceStencils,
                                        Stencils& stokesCellCenterToDarcyStencils,
                                        Stencils& stokesFaceToDarcyStencils)
    {
        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.fvGridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.fvGridGeometry();

        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.numScvf(), false);

        auto darcyFvGeometry = localView(darcyFvGridGeometry);
        auto stokesFvGeometry = localView(stokesFvGridGeometry);
        const auto& stokesGridView = stokesFvGridGeometry.gridView();

        for(const auto& stokesElement : elements(stokesGridView))
        {
            stokesFvGeometry.bindElement(stokesElement);

            for(const auto& scvf : scvfs(stokesFvGeometry))
            {
                // skip the DOF if it is not on the boundary
                if(!scvf.boundary())
                    continue;

                // get element intersecting with the scvf center
                // for robustness, add epsilon in unit outer normal direction
                const auto eps = (scvf.center() - stokesElement.geometry().center()).two_norm()*1e-8;
                auto globalPos = scvf.center(); globalPos.axpy(eps, scvf.unitOuterNormal());
                const auto darcyElementIdx = intersectingEntities(globalPos, darcyFvGridGeometry.boundingBoxTree());

                // skip if no intersection was found
                if(darcyElementIdx.empty())
                    continue;

                // sanity check
                if(darcyElementIdx.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");

                const auto stokesElementIdx = stokesFvGridGeometry.elementMapper().index(stokesElement);

                const auto darcyDofIdx = darcyElementIdx[0];


                stokesFaceToDarcyStencils[scvf.dofIndex()].push_back(darcyDofIdx);
                stokesCellCenterToDarcyStencils[stokesElementIdx].push_back(darcyDofIdx);

                darcyToStokesFaceStencils[darcyElementIdx[0]].push_back(scvf.dofIndex());
                darcyToStokesCellCenterStencils[darcyElementIdx[0]].push_back(stokesElementIdx);

                const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx[0]);
                darcyFvGeometry.bindElement(darcyElement);

                // find the corresponding Darcy sub control volume face
                for(const auto& darcyScvf : scvfs(darcyFvGeometry))
                {
                    const Scalar distance = (darcyScvf.center() - scvf.center()).two_norm();

                    if(distance < eps)
                    {
                        isCoupledDarcyScvf_[darcyScvf.index()] = true;
                        darcyElementToStokesElementMap_[darcyElementIdx[0]].push_back({stokesElementIdx, scvf.index(), darcyScvf.index()});
                        stokesElementToDarcyElementMap_[stokesElementIdx].push_back({darcyElementIdx[0], darcyScvf.index(), scvf.index()});
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a Darcy scvf is coupled to the other domain
     */
    bool isCoupledDarcyScvf(std::size_t darcyScvfIdx) const
    {
        return isCoupledDarcyScvf_[darcyScvfIdx];
    }

    /*!
     * \brief A map that returns all Stokes elements coupled to a Darcy element
     */
    const auto& darcyElementToStokesElementMap() const
    {
        return darcyElementToStokesElementMap_;
    }

    /*!
     * \brief A map that returns all Darcy elements coupled to a Stokes element
     */
    const auto& stokesElementToDarcyElementMap() const
    {
        return stokesElementToDarcyElementMap_;
    }

private:
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> darcyElementToStokesElementMap_;
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> stokesElementToDarcyElementMap_;

    std::vector<bool> isCoupledDarcyScvf_;

    const CouplingManager& couplingManager_;
};

} // end namespace Dumux

#endif
