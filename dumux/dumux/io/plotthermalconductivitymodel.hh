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
 * \ingroup InputOutput
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
#define DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH

#include <string>
#include <vector>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux {

// forward declaration
template<class Scalar> class GnuplotInterface;

/*!
 * \ingroup InputOutput
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
template<class Scalar, class ThermalConductivityModel, class FluidSystem>
class PlotThermalConductivityModel
{
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;

    // phase indices
    enum
    {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx
    };

public:
    /*!
     * \brief Constructor
     *
     * Initializes the fluid system.
     *
     * \param temperature temperature in \f$\mathrm{[K]}\f$
     * \param pressure reference pressure in \f$\mathrm{[Pa]}\f$
     */
    PlotThermalConductivityModel(Scalar temperature = 283.15,
                                 Scalar pressure = 1e5)
    : numIntervals_(1000)
    {
        FluidState fluidstate;
        fluidstate.setTemperature(temperature);
        fluidstate.setPressure(phase0Idx, pressure);
        fluidstate.setPressure(phase1Idx, pressure);
        lambdaW_ = FluidSystem::thermalConductivity(fluidstate, phase0Idx);
        lambdaN_ = FluidSystem::thermalConductivity(fluidstate, phase1Idx);
    }

    /*!
     * \brief Add a effective thermal conductivity-saturation curve to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param porosity The porosity
     * \param rhoSolid The density of the solid phase
     * \param lambdaSolid The conductivity of the solid phase
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addlambdaeffcurve(GnuplotInterface<Scalar> &gnuplot,
                           Scalar porosity,
                           Scalar rhoSolid,
                           Scalar lambdaSolid,
                           Scalar lowerSat = 0.0,
                           Scalar upperSat = 1.0,
                           std::string curveName = "lambdaeff",
                           std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> lambda(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            lambda[i] = ThermalConductivityModel::effectiveThermalConductivity(sw[i], lambdaW_,
                                                                               lambdaN_, lambdaSolid,
                                                                               porosity, rhoSolid);
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("effective thermal conductivity [W/(m K)]");
        gnuplot.addDataSetToPlot(sw, lambda, curveName, curveOptions);
    }

private:
    int numIntervals_;
    Scalar lambdaN_;
    Scalar lambdaW_;
};

} // end namespace Dumux

#endif // DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH