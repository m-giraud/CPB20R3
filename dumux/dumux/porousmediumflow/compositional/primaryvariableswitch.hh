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
 * \ingroup PorousmediumflowModels
 * \brief The primary variable switch base class for compositional models.
 */

#ifndef DUMUX_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Empty class for models without pri var switch.
 */
class NoPrimaryVariableSwitch
{
public:
    template<typename... Args>
    NoPrimaryVariableSwitch(Args&&...) {}

    template<typename... Args> void reset(Args&&...) {}
    template<typename... Args> bool wasSwitched(Args&&...) const { return false; }
    template<typename... Args> bool update(Args&&...) { return false; }
    template<typename... Args> void updateSwitchedVolVars(Args&&...) {}
    template<typename... Args> void updateSwitchedFluxVarsCache(Args&&...) {}
};

/*!
 * \ingroup PorousmediumflowModels
 * \brief The primary variable switch controlling the phase presence state variable.
 */
template<class Implementation>
class PrimaryVariableSwitch
{
public:
    PrimaryVariableSwitch(int verbosity = 1)
    : verbosity_(verbosity)
    {}

    //! If the primary variables were recently switched
    bool wasSwitched(std::size_t dofIdxGlobal) const
    {
        return wasSwitched_[dofIdxGlobal];
    }

    //! Reset all flags
    void reset(const std::size_t numDofs)
    {
        wasSwitched_.resize(numDofs, false);
    }

    /*!
     * \brief Updates the variable switch / phase presence.
     *
     * \param curSol The current solution to be updated / modified
     * \param gridVariables The secondary variables on the grid
     * \param problem The problem
     * \param fvGridGeometry The finite-volume grid geometry
     */
    template<class SolutionVector, class GridVariables, class Problem>
    bool update(SolutionVector& curSol,
                GridVariables& gridVariables,
                const Problem& problem,
                const typename GridVariables::GridGeometry& fvGridGeometry)
    {
        bool switched = false;
        visited_.assign(wasSwitched_.size(), false);
        std::size_t countSwitched = 0;
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            // make sure FVElementGeometry is bound to the element
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            const auto curElemSol = elementSolution(element, curSol, fvGridGeometry);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                if (!visited_[dofIdxGlobal])
                {
                    // Note this implies that volume variables don't differ
                    // in any sub control volume associated with the dof!
                    visited_[dofIdxGlobal] = true;
                    // Compute volVars on which grounds we decide
                    // if we need to switch the primary variables
                    auto& volVars = getVolVarAccess(gridVariables.curGridVolVars(), elemVolVars, scv);
                    volVars.update(curElemSol, problem, element, scv);

                    if (asImp_().update_(curSol[dofIdxGlobal], volVars, dofIdxGlobal, scv.dofPosition()))
                    {
                        switched = true;
                        ++countSwitched;
                    }
                }
            }
        }

        if (verbosity_ > 0 && countSwitched > 0)
            std::cout << "Switched primary variables at " << countSwitched << " dof locations on processor "
                      << fvGridGeometry.gridView().comm().rank() << "." << std::endl;

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag for our partition.
        if (fvGridGeometry.gridView().comm().size() > 1)
            switched = fvGridGeometry.gridView().comm().max(switched);

        return switched;
    }

    /*!
     * \brief Updates the volume variables whose primary variables were
     *        switched.
     *
     * Required when volume variables are cached globally.
     */
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<GridVariables::GridVolumeVariables::cachingEnabled, int> = 0>
    void updateSwitchedVolVars(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector& sol)
    {
        // make sure FVElementGeometry is bound to the element
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bindElement(element);

        // update the secondary variables if global caching is enabled
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdxGlobal = scv.dofIndex();
            if (asImp_().wasSwitched(dofIdxGlobal))
            {
                const auto elemSol = elementSolution(element, sol, fvGridGeometry);
                auto& volVars = gridVariables.curGridVolVars().volVars(scv);
                volVars.update(elemSol, problem, element, scv);
            }
        }
    }

    /*!
     * \brief Updates the fluxVars cache for dof whose primary variables were
     *        switched.
     *
     * Required when flux variables are cached globally (not for box method).
     */
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<(GridVariables::GridFluxVariablesCache::cachingEnabled &&
                               GridVariables::GridGeometry::discMethod != DiscretizationMethod::box), int> = 0>
    void updateSwitchedFluxVarsCache(const Problem& problem,
                                     const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                                     const typename GridVariables::GridGeometry& fvGridGeometry,
                                     GridVariables& gridVariables,
                                     const SolutionVector& sol)
    {
        // update the flux variables if global caching is enabled
        const auto dofIdxGlobal = fvGridGeometry.dofMapper().index(element);

        if (asImp_().wasSwitched(dofIdxGlobal))
        {
            // make sure FVElementGeometry and the volume variables are bound
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bind(element);
            auto curElemVolVars = localView(gridVariables.curGridVolVars());
            curElemVolVars.bind(element, fvGeometry, sol);
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
        }
    }

    /*!
     * \brief Updates the the primary variables state at the boundary.
     *
     * Required when a Dirichlet BC differes from the initial conditon (only for box method).
     */
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<(GridVariables::GridGeometry::discMethod == DiscretizationMethod::box), int> = 0>
    void updateBoundary(const Problem& problem,
                        const typename GridVariables::GridGeometry& fvGridGeometry,
                        GridVariables& gridVariables,
                        SolutionVector& sol)
    {
        std::vector<bool> stateChanged(sol.size(), false);
        std::size_t countChanged = 0;

        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            // skip if the element is not at a boundary
            if (!fvGeometry.hasBoundaryScvf())
                continue;

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (const auto& scv : scvs(fvGeometry))
            {
                // this implies that state is set equal for all scvs associated with the dof
                const auto dofIdx = scv.dofIndex();
                if (!fvGridGeometry.dofOnBoundary(dofIdx) || stateChanged[dofIdx])
                    continue;

                const auto bcTypes = problem.boundaryTypes(element, scv);
                if (bcTypes.hasDirichlet())
                {
                    const auto dirichletValues = problem.dirichlet(element, scv);

                    if (sol[dofIdx].state() != dirichletValues.state())
                    {
                        if (verbosity() > 1)
                            std::cout << "Changing primary variable state at boundary (" << sol[dofIdx].state()
                                      << ") to the one given by the Dirichlet condition (" << dirichletValues.state() << ") at dof " << dofIdx
                                      << ", coordinates: " << scv.dofPosition()
                                      << std::endl;

                        // make sure the solution vector has the right state (given by the Dirichlet BC)
                        sol[dofIdx].setState(dirichletValues.state());
                        stateChanged[dofIdx] = true;
                        ++countChanged;

                        // overwrite initial with Dirichlet values
                        for (int eqIdx = 0; eqIdx < SolutionVector::block_type::dimension; ++eqIdx)
                        {
                            if (bcTypes.isDirichlet(eqIdx))
                            {
                                const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                                sol[dofIdx][pvIdx] = dirichletValues[pvIdx];
                            }
                        }
                    }
                }
            }

            // update the volVars if caching is enabled
            if (GridVariables::GridVolumeVariables::cachingEnabled && countChanged > 0)
            {
                const auto curElemSol = elementSolution(element, sol, fvGridGeometry);
                for (const auto& scv : scvs(fvGeometry))
                {
                    if (stateChanged[scv.dofIndex()])
                    {
                        auto& volVars = getVolVarAccess(gridVariables.curGridVolVars(), elemVolVars, scv);
                        volVars.update(curElemSol, problem, element, scv);
                    }
                }
            }
        }

        if (verbosity_ > 0 && countChanged > 0)
            std::cout << "Changed primary variable states and solution values at boundary to Dirichlet states and values at " << countChanged << " dof locations on processor "
                      << fvGridGeometry.gridView().comm().rank() << "." << std::endl;
    }

    //! Do nothing when volume variables are not cached globally.
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<!GridVariables::GridVolumeVariables::cachingEnabled, int> = 0>
    void updateSwitchedVolVars(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector &uCurrentIter) const {}

    //! Do nothing when flux variables are not cached globally or the box method is used.
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<(!GridVariables::GridFluxVariablesCache::cachingEnabled ||
                               GridVariables::GridGeometry::discMethod == DiscretizationMethod::box), int> = 0>
    void updateSwitchedFluxVarsCache(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector& sol) const {}

    //! Do nothing when the box method is not used.
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<(GridVariables::GridGeometry::discMethod != DiscretizationMethod::box), int> = 0>
    void updateBoundary(const Problem& problem,
                        const typename GridVariables::GridGeometry& fvGridGeometry,
                        GridVariables& gridVariables,
                        SolutionVector& sol) const {}

    //! The verbosity level
    int verbosity() const
    { return verbosity_; }

protected:

    //! Return actual implementation (static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    //! Return actual implementation (static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // Perform variable switch at a degree of freedom location
    template<class VolumeVariables, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 std::size_t dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate if the primary variable switch would switch
        // to be implemented by the deriving class
        DUNE_THROW(Dune::NotImplemented, "This model seems to use a primary variable switch but none is implemented!");
    }

    std::vector<bool> wasSwitched_;
    std::vector<bool> visited_;

private:
    template<class GridVolumeVariables, class ElementVolumeVariables, class SubControlVolume>
    static auto getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    -> std::enable_if_t<!GridVolumeVariables::cachingEnabled, decltype(elemVolVars[scv])>
    { return elemVolVars[scv]; }

    template<class GridVolumeVariables, class ElementVolumeVariables, class SubControlVolume>
    static auto getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    -> std::enable_if_t<GridVolumeVariables::cachingEnabled, decltype(gridVolVars.volVars(scv))>
    { return gridVolVars.volVars(scv); }

    int verbosity_; //!< The verbosity level of the primary variable switch
};

} // end namespace dumux

#endif
