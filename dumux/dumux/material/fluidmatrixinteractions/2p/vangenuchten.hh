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
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations according to van Genuchten.
 */
#ifndef VAN_GENUCHTEN_HH
#define VAN_GENUCHTEN_HH

#include "vangenuchtenparams.hh"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the van Genuchten capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 * For general info: EffToAbsLaw
 *
 * \see VanGenuchtenParams
 */
template <class ScalarT, class ParamsT = VanGenuchtenParams<ScalarT> >
class VanGenuchten
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f$\mathrm{
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     }\f$
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input. Note that for pc(swe = 0.0) = inf, have a look at RegularizedVanGenuchten if this is a problem.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar pc = pow(pow(swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
        return pc;
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
     \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     }\f$
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          The effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       i.e. sw(pc < 0.0) = 0.0, by clamping the input to the physical bounds.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar sw = pow(pow(params.vgAlpha()*pc, params.vgn()) + 1, -params.vgm());
        return sw;
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar endPointPc(const Params &params)
    { return 0.0; }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to van Genuchten.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
     \overline{S}_w^{-1/m} / \overline{S}_w / m
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if swe is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar powSwe = pow(swe, -1/params.vgm());
        return - 1.0/params.vgAlpha() * pow(powSwe - 1, 1.0/params.vgn() - 1)/params.vgn()
                                      * powSwe/swe/params.vgm();
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure according to van Genuchten.
     *
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        using std::pow;
        using std::max;

        pc = max(pc, 0.0); // the equation below is undefined for negative pcs

        const Scalar powAlphaPc = pow(params.vgAlpha()*pc, params.vgn());
        return -pow(powAlphaPc + 1, -params.vgm()-1)*params.vgm()*powAlphaPc/pc*params.vgn();
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::sqrt;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar r = 1.0 - pow(1.0 - pow(swe, 1.0/params.vgm()), params.vgm());
        return sqrt(swe)*r*r;
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrw_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::sqrt;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = 1.0 - pow(swe, 1.0/params.vgm());
        const Scalar xToM = pow(x, params.vgm());
        return (1.0 - xToM)/sqrt(swe) * ( (1.0 - xToM)/2 + 2*xToM*(1.0-x)/x );
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        return pow(1 - swe, 1.0/3) * pow(1 - pow(swe, 1.0/params.vgm()), 2*params.vgm());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     *
     * \note Instead of undefined behaviour if pc is not in the valid range, we return a valid number,
     *       by clamping the input.
     */
    static Scalar dkrn_dswe(const Params &params, Scalar swe)
    {
        using std::pow;
        using std::min;
        using std::max;

        swe = min(max(swe, 0.0), 1.0); // the equation below is only defined for 0.0 <= sw <= 1.0

        const Scalar x = pow(swe, 1.0/params.vgm());
        return -pow(1.0 - x, 2*params.vgm()) * pow(1.0 - swe, -2.0/3) * (1.0/3 + 2*x/swe);
    }

};

} // end namespace Dumux

#endif
