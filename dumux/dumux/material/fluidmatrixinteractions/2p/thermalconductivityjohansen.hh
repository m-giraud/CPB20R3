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
 * \brief   Relation for the saturation-dependent effective thermal conductivity
 */
#ifndef THERMALCONDUCTIVITY_JOHANSEN_HH
#define THERMALCONDUCTIVITY_JOHANSEN_HH

#include <cmath>
#include <algorithm>

namespace Dumux {

struct JohansenIndices
{
    static const int wPhaseIdx = 0;
    static const int nPhaseIdx = 1;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the saturation-dependent effective thermal conductivity
 *
 * The Johansen method (Johansen 1975 \cite johansen1977 ) computes the thermal conductivity of dry and the
 * wet soil material and uses a root function of the wetting saturation to compute the
 * effective thermal conductivity for a two-phase fluidsystem. The individual thermal
 * conductivities are calculated as geometric mean of the thermal conductivity of the porous
 * material and of the respective fluid phase.
 * The material law is:
 * \f$\mathrm{[
 \lambda_\text{eff} = \lambda_{\text{dry}} + \sqrt{(S_w)} \left(\lambda_\text{wet} - \lambda_\text{dry}\right)
 }\f$
 *
 * with
 * \f$\mathrm{
 \lambda_\text{wet} = \lambda_{solid}^{\left(1-\phi\right)}*\lambda_w^\phi
 }\f$
 * and the semi-empirical relation
 *
 * \f$\mathrm{
 \lambda_\text{dry} = \frac{0.135*\rho_s*\phi + 64.7}{\rho_s - 0.947 \rho_s*\phi}.
 }\f$
 *
 * Source: Phdthesis (Johansen1975) Johansen, O. Thermal conductivity of soils Norw. Univ. of Sci. Technol., Trondheim, Norway, 1975 \cite johansen1977
 */
template<class Scalar, class Indices = JohansenIndices>
class ThermalConductivityJohansen
{
public:
    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Johansen (1975) \cite johansen1977 .
     *
     * \param volVars volume variables
     * \param spatialParams spatial parameters
     * \param element element (to be passed to spatialParams)
     * \param fvGeometry fvGeometry (to be passed to spatialParams)
     * \param scv the sub control volume (to be passed to spatialParams)
     *
     * \return Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Johansen (1975) \cite johansen1977 <BR>
     *
     * This formulation is semi-empirical and fitted to quartz sand.
     * This gives an interpolation of the effective thermal conductivities of a porous medium
     * filled with the non-wetting phase and a porous medium filled with the wetting phase.
     * These two effective conductivities are computed as geometric mean of the solid and the
     * fluid conductivities and interpolated with the Kersten number.<br>
     * Johansen, O. 1975. Thermal conductivity of soils. Ph.D. diss. Norwegian Univ.
     *                    of Sci. and Technol., Trondheim. (Draft Transl. 637. 1977. U.S. Army
     *                    Corps of Eng., Cold Regions Res. and Eng. Lab., Hanover, NH.) \cite johansen1977
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry, class SubControlVolume>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars,
                                               const SpatialParams& spatialParams,
                                               const Element& element,
                                               const FVGeometry& fvGeometry,
                                               const SubControlVolume& scv)
    {
        const Scalar sw = volVars.saturation(Indices::wPhaseIdx);
        const Scalar lambdaW = volVars.fluidThermalConductivity(Indices::wPhaseIdx);
        const Scalar lambdaN = volVars.fluidThermalConductivity(Indices::nPhaseIdx);
        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();
        const Scalar rhoSolid = volVars.solidDensity();

        return effectiveThermalConductivity(sw, lambdaW, lambdaN, lambdaSolid, porosity, rhoSolid);
    }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Johansen (1975) \cite johansen1977 .
     *
     * \param Sw The saturation of the wetting phase
     * \param lambdaW The thermal conductivity of the wetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaN The thermal conductivity of the non-wetting phase in \f$\mathrm{[W/(m K)]}\f$
     * \param lambdaSolid The thermal conductivity of the solid phase in \f$\mathrm{[W/(m K)]}\f$
     * \param porosity The porosity
     * \param rhoSolid The density of solid phase in \f$\mathrm{[kg/m^3]}\f$
     *
     * \return Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$ after Johansen (1975) \cite johansen1977
     */
    static Scalar effectiveThermalConductivity(const Scalar Sw,
                                               const Scalar lambdaW,
                                               const Scalar lambdaN,
                                               const Scalar lambdaSolid,
                                               const Scalar porosity,
                                               const Scalar rhoSolid)
    {
        using std::max;
        const Scalar satW = max<Scalar>(0.0, Sw);

        const Scalar kappa = 15.6; // fitted to medium quartz sand
        const Scalar rhoBulk = rhoSolid*porosity;
        // lambdaSolid^(1-porosity) * lambdaW^porosity =

        using std::pow;
        const Scalar lSat = lambdaSolid * pow(lambdaW / lambdaSolid, porosity);
        const Scalar lDry = (0.135*rhoBulk + 64.7)/(rhoSolid - 0.947*rhoBulk);
        const Scalar Ke = (kappa*satW)/(1+(kappa-1)*satW);// Kersten number, equation 13

        return lDry + Ke * (lSat - lDry); // equation 14
    }
};
} // end namespace Dumux
#endif
