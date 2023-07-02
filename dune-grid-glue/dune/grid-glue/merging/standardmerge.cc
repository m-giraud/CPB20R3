// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include "standardmerge.hh"

namespace Dune {
namespace GridGlue {

#define DECL
#define STANDARD_MERGE_INSTANTIATE(T,A,B,C) \
  DECL template \
  void StandardMerge<T,A,B,C>::build(const std::vector<Dune::FieldVector<T,C> >& grid1_coords, \
                                     const std::vector<unsigned int>& grid1_elements, \
                                     const std::vector<Dune::GeometryType>& grid1_element_types, \
                                     const std::vector<Dune::FieldVector<T,C> >& grid2_coords, \
                                     const std::vector<unsigned int>& grid2_elements, \
                                     const std::vector<Dune::GeometryType>& grid2_element_types \
                                     )

STANDARD_MERGE_INSTANTIATE(double,1,1,1);
STANDARD_MERGE_INSTANTIATE(double,2,2,2);
STANDARD_MERGE_INSTANTIATE(double,3,3,3);
#undef STANDARD_MERGE_INSTANTIATE
#undef DECL

} /* namespace GridGlue */
} /* namespace Dune */
