#ifndef DUNE_GRIDGLUE_COMMON_AREAWRITER_HH
#define DUNE_GRIDGLUE_COMMON_AREAWRITER_HH

#include <ostream>
#include <string>

namespace Dune {
namespace GridGlue {

template<int side, typename Glue>
void write_glue_area_vtk(const Glue& glue, std::ostream& out);

template<int side, typename Glue>
void write_glue_area_vtk(const Glue& glue, const std::string& filename);

/**
 * Write area covered by grid glue.
 *
 * Write a VTK file for each side that indicate where the grid glue
 * intersections are defined.  A cell data field is provided that gives
 * the relative area that is covered by glue intersections: if the value
 * is zero no intersections are defined on the element, if the value is one
 * the entire element should be covered (provided the glue is injective).
 * Note that also values greater than one can be reached if the mapping is
 * not injective.
 *
 * This method is intended to be used for debugging purposes only.
 *
 * \param glue GridGlue for which the coverage by intersections should
 *             be computed
 * \param base prefix of filenames to use.  The generated files are named
 *             <code><em>base</em>-inside.vtk</code> and
 *             <code><em>base</em>-outside.vtk</code>.
 */
template<typename Glue>
void write_glue_areas_vtk(const Glue& glue, const std::string& base);

} /* namespace GridGlue */
} /* namespace Dune */

#include "areawriter_impl.hh"

#endif
