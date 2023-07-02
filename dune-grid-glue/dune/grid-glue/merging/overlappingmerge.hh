// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_MERGING_OVERLAPPINGMERGE_HH
#define DUNE_GRIDGLUE_MERGING_OVERLAPPINGMERGE_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/standardmerge.hh>
#include <dune/grid-glue/merging/computeintersection.hh>

namespace Dune {
namespace GridGlue {

/** \brief Computing overlapping grid intersections for grids of different dimensions

   \tparam dim1 Grid dimension of grid 1
   \tparam dim2 Grid dimension of grid 2
   \tparam dimworld World dimension
   \tparam T Type used for coordinates
 */
template<int dim1, int dim2, int dimworld, typename T = double>
class OverlappingMerge
      : public StandardMerge<T,dim1,dim2,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  //typedef Dune::FieldVector<T, dim>  LocalCoords;

  OverlappingMerge()
  {}

protected:
  typedef typename StandardMerge<T,dim1,dim2,dimworld>::SimplicialIntersection SimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

   The result is a set of simplices.

   \param grid1ElementType Type of the first element to be intersected
   \param grid1ElementCorners World coordinates of the corners of the first element

   \param grid2ElementType Type of the second element to be intersected
   \param grid2ElementCorners World coordinates of the corners of the second element

   */
  void computeIntersections(const Dune::GeometryType& grid1ElementType,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                           std::bitset<(1<<dim1)>& neighborIntersects1,
                           unsigned int grid1Index,
                           const Dune::GeometryType& grid2ElementType,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                           std::bitset<(1<<dim2)>& neighborIntersects2,
                           unsigned int grid2Index,
                           std::vector<SimplicialIntersection>& intersections);

private:
  bool inPlane(std::vector<FieldVector<T,dimworld> >& points);

};

} /* namespace Dune::GridGlue */
} /* namespace Dune */

#include "overlappingmerge.cc"


#endif // DUNE_GRIDGLUE_MERGING_OVERLAPPINGMERGE_HH
