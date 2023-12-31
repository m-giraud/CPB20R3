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
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 */
#ifndef VAN_GENUCHTEN_PARAMS_HH
#define VAN_GENUCHTEN_PARAMS_HH

#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 *
 *       In this implementation setting either the \f$\mathrm{n}\f$ or \f$\mathrm{m}\f$ shape parameter
 *       automatically calculates the other. I.e. they cannot be set independently.
 */
template<class ScalarT>
class VanGenuchtenParams
{
public:
    using Scalar = ScalarT;

    VanGenuchtenParams()
    {}

    VanGenuchtenParams(Scalar vgAlpha, Scalar vgn)
    {
        setVgAlpha(vgAlpha);
        setVgn(vgn);
    }

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        return Dune::FloatCmp::eq(vgAlpha_, otherParams.vgAlpha(), /*eps*/1e-6*vgAlpha_)
               && Dune::FloatCmp::eq(vgn_, otherParams.vgn(), /*eps*/1e-6*vgn_);
    }

    /*!
     * \brief Return the \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

    /*!
     * \brief Set the \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$ of van Genuchten's
     *        curve.
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$\mathrm{m}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgm() const
    { return vgm_; }

    /*!
     * \brief Set the \f$\mathrm{m}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$\mathrm{n = \frac{1}{1 - m}}\f$
     */
    void setVgm(Scalar m)
    { vgm_ = m; vgn_ = 1/(1 - vgm_); }

    /*!
     * \brief Return the \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgn() const
    { return vgn_; }

    /*!
     * \brief Set the \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$\mathrm{m = 1 - \frac{1}{n}}\f$
     */
    void setVgn(Scalar n)
    { vgn_ = n; vgm_ = 1 - 1/vgn_; }

private:
    Scalar vgAlpha_;
    Scalar vgm_;
    Scalar vgn_;
};
} // namespace Dumux

#endif
