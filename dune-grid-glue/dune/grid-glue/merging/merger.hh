// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_MERGING_MERGER_HH
#define DUNE_GRIDGLUE_MERGING_MERGER_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dune/grid-glue/merging/intersectionlist.hh>

namespace Dune {
namespace GridGlue {

/** \brief Abstract base for all classes that take extracted grids and build sets of intersections

   \tparam ctype The type used for coordinates (assumed to be the same for both grids)
   \tparam grid1Dim Dimension of the grid1 grid
   \tparam grid2Dim Dimension of the grid2 grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
class Merger
{
public:

  /// @brief the local coordinate type for the grid1 coordinates
  typedef Dune::FieldVector<ctype, grid1Dim>  Grid1Coords;

  /// @brief the local coordinate type for the grid2 coordinates
  typedef Dune::FieldVector<ctype, grid2Dim>  Grid2Coords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<ctype, dimworld>  WorldCoords;

  using IntersectionList = Dune::GridGlue::IntersectionList<Grid1Coords, Grid2Coords>;

  /**
   * @brief builds the merged grid
   *
   * Note that the indices are used consequently throughout the whole class interface just like they are
   * introduced here.
   *
   * @param grid1_coords the grid1 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid1_elements array with all grid1 elements represented as corner indices into @c grid1_coords
   * @param grid1_element_types array with the GeometryType of the elements listed grid1_elements
   * @param grid2_coords the grid2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid2_elements just like with the grid1_elements and grid1_coords
   * @param grid2_element_types array with the GeometryType of the elements listed grid2_elements
   */
  virtual void build(const std::vector<Dune::FieldVector<ctype,dimworld> >& grid1_coords,
                     const std::vector<unsigned int>& grid1_elements,
                     const std::vector<Dune::GeometryType>& grid1_element_types,
                     const std::vector<Dune::FieldVector<ctype,dimworld> >& grid2_coords,
                     const std::vector<unsigned int>& grid2_elements,
                     const std::vector<Dune::GeometryType>& grid2_element_types) = 0;

  /** @brief get the number of simplices in the merged grid
      The indices are then in 0..nSimplices()-1
   */
  unsigned int nSimplices() const
    { return intersectionList()->size(); }

  virtual void clear() = 0;

  /**
   * list of intersections
   *
   * \note only valid after `build()` was called
   */
  virtual std::shared_ptr<IntersectionList> intersectionList() const = 0;

  /**
    * doc me
    */
  template<int n>
  unsigned int parents(unsigned int idx) const {
    return intersectionList()->template parents<n>(idx);
  }

  /**
   * @brief get index of grid-n's parent simplex for given merged grid simplex
   * @tparam n specify which grid
   * @param idx index of the merged grid simplex
   * @return index of the parent simplex
   */
  template<int n>
  unsigned int parent(unsigned int idx, unsigned int parId = 0) const
  {
    return intersectionList()->template parent<n>(idx, parId);
  }

  /**
   * @brief get the grid-n parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "parent<n>")
   * @tparam n specify which grid
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in grid-n grid1
   */
  template<int n>
  auto parentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const
  {
    return intersectionList()->template corner<n>(idx, corner, parId);
  }

  /** \brief Counts the number of times the computeIntersection method has been called
   *
   * Used temporarily to speed up the implementation
   */
  unsigned int counter;
};

} /* namespace GridGlue */
} /* namespace Dune */

#endif
