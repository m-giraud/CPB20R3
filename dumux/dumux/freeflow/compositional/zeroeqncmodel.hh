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
 *
 * \brief A single-phase, multi-component Reynolds-Averaged Navier-Stokes 0-Eq. model
 *
 * \copydoc dumux/freeflow/compositional/navierstokesncmodel.hh
 */

#ifndef DUMUX_ZEROEQ_NC_MODEL_HH
#define DUMUX_ZEROEQ_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/iofields.hh>
#include <dumux/freeflow/rans/zeroeq/model.hh>

#include "iofields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component ZeroEq model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component isothermal ZeroEq model
struct ZeroEqNC { using InheritsFrom = std::tuple<NavierStokesNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup ZeroEqModel
 * \brief Traits for the Reynolds-averaged Navier-Stokes 0-Eq. model
 */
template<int dimension, int nComp, bool useM, int replaceCompEqIdx>
struct ZeroEqNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, useM, replaceCompEqIdx>
{
    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }
};

//! The model traits of the isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ZeroEqNC>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::FVGridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
public:
    using type = ZeroEqNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ZeroEqNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using CompositionalVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = ZeroEqVolumeVariables<Traits, CompositionalVolVars>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::ZeroEqNC> { using type = FreeflowNCIOFields<RANSIOFields, true/*turbulenceModel*/>; };

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component ZeroEq model
//////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component non-isothermal ZeroEq models
struct ZeroEqNCNI { using InheritsFrom = std::tuple<NavierStokesNCNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ZeroEqNCNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::FVGridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
    using IsothermalModelTraits = ZeroEqNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ZeroEqNCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = ZeroEqVolumeVariables<Traits, NCVolVars>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::ZeroEqNCNI>
{
private:
    using IsothermalIOFields = FreeflowNCIOFields<RANSIOFields, true/*turbulenceModel*/>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalIOFields, true/*turbulenceModel*/>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
