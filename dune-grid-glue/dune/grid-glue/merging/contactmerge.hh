// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
  @file
  @brief Merge two grid boundary surfaces that may be a positive distance apart
 */

#ifndef DUNE_GRIDGLUE_MERGING_CONTACTMERGE_HH
#define DUNE_GRIDGLUE_MERGING_CONTACTMERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>
#include <memory>
#include <functional>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/deprecated.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/standardmerge.hh>
#include <dune/grid-glue/gridglue.hh>

namespace Dune {
namespace GridGlue {

/** \brief Merge two codimension-1 surfaces that may be a positive distance apart

  \tparam dimworld  Dimension of the world coordinates.
  \tparam T Type used for coordinates
 */
template<int dimworld, typename T = double>
class ContactMerge
: public StandardMerge<T,dimworld-1,dimworld-1,dimworld>
{
    static constexpr int dim = dimworld-1;

    static_assert( dim==1 || dim==2,
            "ContactMerge yet only handles the cases dim==1 and dim==2!");

    typedef StandardMerge<T,dim,dim,dimworld> Base;
public:

    /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

    /// @brief the numeric type used in this interface
    typedef T ctype;

    /// @brief the coordinate type used in this interface
    typedef Dune::FieldVector<T, dimworld>  WorldCoords;

    /// @brief the coordinate type used in this interface
    typedef Dune::FieldVector<T, dim>  LocalCoords;

    /// @brief Type of the projection, closest point or outer normal projection
    enum ProjectionType {OUTER_NORMAL, CLOSEST_POINT};
    /**
     * @brief Construct merger given overlap and possible projection directions.
     *
     * @param allowedOverlap Allowed overlap of the surfaces
     * @param domainDirections Projection direction field for the first surface that differ from the default normal field
     * @param targetDirections Projection direction field for the second surface that differ from the default normal field
     */
    ContactMerge(const T allowedOverlap=T(0),
                 std::function<WorldCoords(WorldCoords)> domainDirections=nullptr,
                 std::function<WorldCoords(WorldCoords)> targetDirections=nullptr,
                 ProjectionType type = OUTER_NORMAL)
        : domainDirections_(domainDirections), targetDirections_(targetDirections),
          overlap_(allowedOverlap), type_(type)
    {}

    /**
     * @brief Construct merger given overlap and type of the projection
     * @param allowedOverlap Allowed overlap of the surfacs
     * @param type Type of the projection
     */
    ContactMerge(const T allowedOverlap, ProjectionType type)
        : overlap_(allowedOverlap),
          type_(type)
    {}

    /**
     * @brief Set surface direction functions
     *
     * The matching of the geometries offers the possibility to specify a function for
     * the exact evaluation of domain surface normals. If no such function is specified
     * (default) normals are interpolated.
     * @param value the new function (or nullptr to unset the function)
     */
    inline
    void setSurfaceDirections(std::function<WorldCoords(WorldCoords)> domainDirections,
                              std::function<WorldCoords(WorldCoords)> targetDirections)
    {
        domainDirections_ = domainDirections;
        targetDirections_ = targetDirections;
        this->valid = false;
    }

    //! Set the allowed overlap of the surfaces.
    void setOverlap(T overlap)
    {
        overlap_ = overlap;
    }

    //!Get the allowed overlap of the surfaces.
    T getOverlap() const
    {
        return overlap_;
    }

    /**
     * \brief set minimum angle in radians between normals at <code>x</code> and <code>Φ(x)</code>
     */
    void minNormalAngle(T angle)
    {
        using std::cos;
        maxNormalProduct_ = cos(angle);
    }

    /**
     * \brief get minimum angle in radians between normals at <code>x</code> and <code>Φ(x)</code>
     */
    T minNormalAngle() const
    {
        using std::acos;
        return acos(maxNormalProduct_);
    }

protected:
  typedef typename StandardMerge<T,dimworld-1,dimworld-1,dimworld>::SimplicialIntersection SimplicialIntersection;

private:
    /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
     */
    std::function<WorldCoords(WorldCoords)> domainDirections_;
    std::vector<WorldCoords> nodalDomainDirections_;

    /** \brief Vector field on the target surface which prescribes a 'forward'
      direction.

      We use the normals of the target side to increase projection
      robustness.  If these cannot be computed from the surface directly
      (e.g. because it is not properly oriented), they can be given
      explicitly through the targetDirections field.
     */
    std::function<WorldCoords(WorldCoords)> targetDirections_;
    std::vector<WorldCoords> nodalTargetDirections_;

    //! Allow some overlap, i.e. also look in the negative projection directions
    T overlap_;

    //! The type of the projection, i.e. closest point projection or outer normal
    ProjectionType type_;

    /**
     * See Projection::m_max_normal_product
     */
    T maxNormalProduct_ = T(-0.1);

    /** \brief Compute the intersection between two overlapping elements
     *
     *   The result is a set of simplices.
     */
    void computeIntersections(const Dune::GeometryType& grid1ElementType,
            const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
            std::bitset<(1<<dim)>& neighborIntersects1,
            unsigned int grid1Index,
            const Dune::GeometryType& grid2ElementType,
            const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
            std::bitset<(1<<dim)>& neighborIntersects2,
            unsigned int grid2Index,
            std::vector<SimplicialIntersection>& intersections) override;

    /**
      * @copydoc StandardMerge<T,grid1Dim,grid2Dim,dimworld>::build
     */
protected:
    void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
        const std::vector<unsigned int>& grid1Elements,
        const std::vector<Dune::GeometryType>& grid1ElementTypes,
        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
        const std::vector<unsigned int>& grid2Elements,
        const std::vector<Dune::GeometryType>& grid2ElementTypes) override
    {
        std::cout<<"ContactMerge building grid!\n";
        // setup the nodal direction vectors
        setupNodalDirections(grid1Coords, grid1Elements, grid1ElementTypes,
                grid2Coords, grid2Elements, grid2ElementTypes);

        Base::build(grid1Coords, grid1Elements, grid1ElementTypes,
                   grid2Coords, grid2Elements, grid2ElementTypes);

    }

private:

    //! Compute local coordinates of a corner
    static LocalCoords localCornerCoords(int i, const Dune::GeometryType& gt)
    {
        const auto& ref = Dune::ReferenceElements<T,dim>::general(gt);
        return ref.position(i,dim);
    }

protected:

    //! Order the corners of the intersection polytope in cyclic order
    void computeCyclicOrder(const std::vector<std::array<LocalCoords,2> >& polytopeCorners,
            const LocalCoords& center, std::vector<int>& ordering) const;

    //! Setup the direction vectors containing the directions for each vertex
    void setupNodalDirections(const std::vector<WorldCoords>& coords1,
            const std::vector<unsigned int>& elements1,
            const std::vector<Dune::GeometryType>& elementTypes1,
            const std::vector<WorldCoords>& coords2,
            const std::vector<unsigned int>& elements2,
            const std::vector<Dune::GeometryType>& elementTypes2);

    //! If no direction field was specified compute the outer normal field
    void computeOuterNormalField(const std::vector<WorldCoords>& coords,
            const std::vector<unsigned int>& elements,
            const std::vector<Dune::GeometryType>& elementTypes,
            std::vector<WorldCoords>& normals);

    //! Remove all multiples
    void removeDoubles(std::vector<std::array<LocalCoords,2> >& polytopeCorners);
};

} /* namespace GridGlue */
} /* namespace Dune */

#include "contactmerge.cc"

#endif // DUNE_GRIDGLUE_MERGING_CONTACTMERGE_HH
