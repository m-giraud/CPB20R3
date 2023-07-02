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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of linear elastic geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_FV_SPATIAL_PARAMS_HH
#define DUMUX_GEOMECHANICS_ELASTIC_FV_SPATIAL_PARAMS_HH

#include <memory>

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a lameParamsAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasLameParamsAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.lameParamsAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Detail
#endif

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of linear elastic geomechanical problems
 */
template<class Scalar, class FVGridGeometry, class Implementation>
class FVSpatialParamsElastic
{
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! The constructor
    FVSpatialParamsElastic(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : fvGridGeometry_(fvGridGeometry)
    {}

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        static_assert(SolidSystem::isInert(), "Elastic model can only be used with inert solid systems");

        // when there is only one component, the volume fraction must be one
        if (SolidSystem::numInertComponents == 1)
            return 1.0;

        // otherwise we require the user to define the solid composition
        return asImp_().template inertVolumeFractionAtPos<SolidSystem>(scv.center(), compIdx);
    }

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param globalPos The global position
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "The spatial parameters do not provide inertVolumeFractionAtPos() method."); }

    /*!
     * \brief Define the Lame parameters
     * \note  These are possibly solution dependent and are evaluated
     *        for an integration point inside the element. Therefore,
     *        a flux variables cache object is passed to this function
     *        containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param fvGeometry The local finite volume geometry
     * \param elemVolVars Primary/Secondary variables inside the element
     * \param fluxVarsCache Contains data on shape functions at the integration point
     * \return lame parameters
     */
    template<class ElemVolVars, class FluxVarsCache>
    decltype(auto) lameParams(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElemVolVars& elemVolVars,
                              const FluxVarsCache& fluxVarsCache) const
    {
        static_assert(decltype(isValid(Detail::hasLameParamsAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const LameParams& lameParams(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemVolVars& elemVolVars,\n"
        "                                      const FluxVarsCache& fluxVarsCache) const\n\n");

        return asImp_().lameParamsAtPos(fluxVarsCache.ipGlobal());
    }

    //! The finite volume grid geometry
    const FVGridGeometry& fvGridGeometry() const { return *fvGridGeometry_; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;
};
} // end namespace Dumuxs
#endif
