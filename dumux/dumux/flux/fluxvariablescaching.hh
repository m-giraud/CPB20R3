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
 * \ingroup Flux
 * \brief Classes related to flux variables caching
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH
#define DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH

namespace Dumux {
namespace FluxVariablesCaching {

#ifndef DOXYGEN // hide the empty caches from doxygen

//! The empty filler class corresponding to EmptyCache
struct EmptyCacheFiller
{
    template<typename... Args>
    static void fill(Args&&... args) {}
};

// an empty cache filler
// \note Never use the _EmptyCache directly as it lead to ambiguous definitions
struct _EmptyCache
{ using Filler = EmptyCacheFiller; };

#endif // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief Empty caches to use in a constitutive flux law/process, e.g. Darcy's law
 */
class EmptyAdvectionCache : public _EmptyCache {};
class EmptyDiffusionCache : public _EmptyCache {};
class EmptyHeatConductionCache : public _EmptyCache {};

} // end namespace FluxVariablesCaching
} // end namespace Dumux

#endif
