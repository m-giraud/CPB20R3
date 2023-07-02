// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    conformingmerge.hh
 *  Version:     1.0
 *  Created on:  Sep 14, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: implementation of the Merger concept for conforming interfaces
 *
 */
/**
 * @file
 * @brief
 * Implementation of the Merger concept for conforming interfaces
 */

#ifndef DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH
#define DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH

#include <iomanip>
#include <vector>
#include <algorithm>
#include <bitset>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid-glue/merging/standardmerge.hh>

namespace Dune {

  namespace GridGlue {

/** \brief Implementation of the Merger concept for conforming interfaces

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class ConformingMerge
  : public StandardMerge<T,dim,dim,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  T tolerance_;

  typedef typename StandardMerge<T,dim,dim,dimworld>::SimplicialIntersection SimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersections(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects1,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects2,
                                   unsigned int grid2Index,
                                   std::vector<SimplicialIntersection>& intersections);

public:

  static constexpr T default_tolerance = 1e-4;

  ConformingMerge(T tolerance = default_tolerance) :
    tolerance_(tolerance)
  {}
};

template<int dim, int dimworld, typename T>
constexpr T ConformingMerge<dim, dimworld, T>::default_tolerance;

template<int dim, int dimworld, typename T>
void ConformingMerge<dim, dimworld, T>::computeIntersections(const Dune::GeometryType& grid1ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                                            std::bitset<(1<<dim)>& neighborIntersects1,
                                                            unsigned int grid1Index,
                                                            const Dune::GeometryType& grid2ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                                            std::bitset<(1<<dim)>& neighborIntersects2,
                                                            unsigned int grid2Index,
                                                            std::vector<SimplicialIntersection>& intersections)
{
  this->counter++;

  // A few consistency checks
  assert((unsigned int)(Dune::ReferenceElements<T,dim>::general(grid1ElementType).size(dim)) == grid1ElementCorners.size());
  assert((unsigned int)(Dune::ReferenceElements<T,dim>::general(grid2ElementType).size(dim)) == grid2ElementCorners.size());
  // any intersection we may find will be the entire elements.
  neighborIntersects1.reset();
  neighborIntersects2.reset();

  // the intersection is either conforming or empty, hence the GeometryTypes have to match
  if (grid1ElementType != grid2ElementType)
    return;

  // ////////////////////////////////////////////////////////////
  //   Find correspondences between the different corners
  // ////////////////////////////////////////////////////////////
  std::vector<int> other(grid1ElementCorners.size(), -1);

  for (unsigned int i=0; i<grid1ElementCorners.size(); i++) {

    for (unsigned int j=0; j<grid2ElementCorners.size(); j++) {

      if ( (grid1ElementCorners[i]-grid2ElementCorners[j]).two_norm() < tolerance_ ) {

        other[i] = j;
        break;

      }

    }

    // No corresponding grid2 vertex found for this grid1 vertex
    if (other[i] == -1)
      return;

  }

  // ////////////////////////////////////////////////////////////
  //   Set up the new remote intersection
  // ////////////////////////////////////////////////////////////

  const auto& refElement = Dune::ReferenceElements<T,dim>::general(grid1ElementType);

  /** \todo Currently the Intersections have to be simplices */
  if (grid1ElementType.isSimplex()) {

    intersections.emplace_back(grid1Index, grid2Index);

    for (int i=0; i<refElement.size(dim); i++) {
      intersections.back().corners0[0][i] = refElement.position(i,dim);
      intersections.back().corners1[0][i] = refElement.position(other[i],dim);
    }

  } else if (dim == 2 && grid1ElementType.isQuadrilateral()) {

    // split the quadrilateral into two triangles
    const unsigned int subVertices[2][3] = {{0,1,3}, {0,3,2}};

    for (int i=0; i<2; i++) {

      SimplicialIntersection newSimplicialIntersection(grid1Index, grid2Index);

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.corners0[0][j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.corners1[0][j] = refElement.position(subVertices[i][other[j]],dim);
      }

      intersections.push_back(newSimplicialIntersection);

    }

  } else if (grid1ElementType.isHexahedron()) {

    // split the hexahedron into five tetrahedra
    // This can be removed if ever we allow Intersections that are not simplices
    const unsigned int subVertices[5][4] = {{0,1,3,5}, {0,3,2,6}, {4,5,0,6}, {6,7,6,3}, {6,0,5,3}};

    for (int i=0; i<5; i++) {

      SimplicialIntersection newSimplicialIntersection(grid1Index, grid2Index);

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.corners0[0][j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.corners1[0][j] = refElement.position(subVertices[i][other[j]],dim);
      }

      intersections.push_back(newSimplicialIntersection);

    }

  } else
    DUNE_THROW(Dune::GridError, "Unsupported element type");

}

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH
