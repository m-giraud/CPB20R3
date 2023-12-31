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
 *
 * \brief Test for the adaptive 2p2c model in 2d
 */
#include <config.h>

#if HAVE_UG

#include "test_adaptive2p2c2dproblem.hh"
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                       "\t-TimeManager.TEnd      End of the simulation [s] \n"
                                       "\t-TimeManager.DtInitial Initial timestep size [s] \n"
                                       "\t-Grid.Cells   Resolution in x- and y-direction [-]\n"
                                       "\t-Grid.UpperRight      Length and height of the domain [m]\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}
// The main function using the standard start procedure
int main(int argc, char** argv)
{
    using TypeTag = TTAG(Adaptive2p2c2d);
    return Dumux::start<TypeTag>(argc, argv, usage);
}

#else

#include <iostream>

int main()
{
#warning You need to have UGGrid to run this test
    std::cerr << "You need to have UGGrid to run this test\n";
    return 77;
}
#endif
