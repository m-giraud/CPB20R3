// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * @brief Central component of the module implementing the coupling of two grids.
 * @author Gerrit Buse, Christian Engwer
 */


#ifndef DUNE_GRIDGLUE_GRIDGLUE_HH
#define DUNE_GRIDGLUE_GRIDGLUE_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/stdstreams.hh>

#include "adapter/gridgluecommunicate.hh"
#include <dune/grid-glue/merging/merger.hh>

#include <dune/common/parallel/mpitraits.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/plocalindex.hh>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/common/parallel/interface.hh>

namespace Dune {
namespace GridGlue {

// forward declarations
template<typename P0, typename P1>
class GridGlue;

template<typename P0, typename P1>
class IntersectionData;

template<typename P0, typename P1, int inside, int outside>
class Intersection;

template<typename P0, typename P1, int inside, int outside>
class IntersectionIterator;

template<typename P0, typename P1>
class IntersectionIndexSet;

/**
 * @class GridGlue
 * @brief sequential adapter to couple two grids at specified close together boundaries
 *
 * @tparam P0 patch (extractor) to use for grid 0
 * @tparam P1 patch (extractor) to use for grid 1
 *
 * \todo adapt member names according to style guide
 */
template<typename P0, typename P1>
class GridGlue
{
private:

  /*   F R I E N D S   */

  friend class IntersectionData<P0,P1>;
  friend class Intersection<P0,P1,0,1>;
  friend class Intersection<P0,P1,1,0>;
  friend class IntersectionIterator<P0,P1,0,1>;
  friend class IntersectionIterator<P0,P1,1,0>;
  friend class IntersectionIndexSet<P0,P1>;

  /*   P R I V A T E   T Y P E S   */

  /** \brief GlobalId type of an intersection (used for communication) */
  typedef ::Dune::GridGlue::GlobalId GlobalId;

  /** \brief LocalIndex type of an intersection (used for communication) */
  typedef Dune::ParallelLocalIndex <Dune::PartitionType> LocalIndex;

  /** \brief ParallelIndexSet type (used for communication communication) */
  typedef Dune::ParallelIndexSet <GlobalId, LocalIndex> PIndexSet;

public:

  /*   P U B L I C   T Y P E S   A N D   C O N S T A N T S   */

  /** coupling patch of grid <code>side</code> */
  template<int side>
  using GridPatch = std::conditional_t<side == 0, P0, std::conditional_t<side == 1, P1, void>>;

  /** GridView of grid <code>side</code> */
  template<int side>
  using GridView = typename GridPatch<side>::GridView;

  /** Grid type of grid <code>side</code> */
  template<int side>
  using Grid = typename GridView<side>::Grid;

  /** \brief GridView of grid 0 (aka domain grid) */
  using Grid0View DUNE_DEPRECATED_MSG("please use GridView<0> instead") = GridView<0>;

  /** \brief Grid 0 type */
  using Grid0 DUNE_DEPRECATED_MSG("please use Grid<0> instead") = Grid<0>;

  /** \brief Coupling patch of grid 0 */
  using Grid0Patch DUNE_DEPRECATED_MSG("please use GridPatch<0> instead") = GridPatch<0>;

  template<int side>
  static constexpr auto griddim()
    { return GridPatch<side>::dim; }

  template<int side>
  static constexpr auto griddimworld()
    { return GridPatch<side>::dimworld; }

  /** \brief dimension of the grid 0 extractor */
  DUNE_DEPRECATED_MSG("please use griddim<0>() instead")
  static constexpr auto grid0dim = griddim<0>();

  DUNE_DEPRECATED_MSG("please use griddim<0>() instead")
  static constexpr auto domdim = griddim<0>();

  /** \brief world dimension of the grid 0 extractor */
  DUNE_DEPRECATED_MSG("please use griddimworld<0>() instead")
  static constexpr auto grid0dimworld = griddimworld<0>();

  DUNE_DEPRECATED_MSG("please use griddimworld<0>() instead")
  static constexpr auto domdimworld = griddimworld<0>();

  /** \brief GridView of grid 1 (aka target grid) */
  using Grid1View DUNE_DEPRECATED_MSG("please use GridView<0> instead") = GridView<1>;

  /** \brief Grid 1 type */
  using Grid1 DUNE_DEPRECATED_MSG("please use Grid<1> instead") = Grid<1>;

  /** \brief Coupling patch of grid 1 */
  using Grid1Patch DUNE_DEPRECATED_MSG("please use GridPatch<1> instead") = GridPatch<1>;

  /** \todo */
  typedef unsigned int IndexType;

  /** \brief dimension of the grid 1 extractor */
  DUNE_DEPRECATED_MSG("please use griddim<1>() instead")
  static constexpr auto grid1dim = griddim<1>();

  DUNE_DEPRECATED_MSG("please use griddim<1>() instead")
  static constexpr auto tardim = griddim<1>();

  /** \brief world dimension of the grid 1 extractor */
  DUNE_DEPRECATED_MSG("please use griddimworld<1>() instead")
  static constexpr auto grid1dimworld = griddimworld<1>();

  DUNE_DEPRECATED_MSG("please use griddimworld<1>() instead")
  static constexpr auto tardimworld = griddimworld<1>();

  /** \brief export the world dimension
   * This is the maximum of the extractors' world dimensions.
   */
  static constexpr int dimworld = (int)griddimworld<0>() > (int)griddimworld<1>() ? (int)griddimworld<0>() : (int)griddimworld<1>();

  /** \brief The type used for coordinates
   */
  typedef typename PromotionTraits<typename GridView<0>::ctype,
                                   typename GridView<1>::ctype>::PromotedType ctype;

  /** \brief The type used for coordinate vectors */
  typedef Dune::FieldVector<ctype, dimworld>                   Coords;

  /** \brief type of grid elements on side <code>side</code> */
  template<int side>
  using GridElement = typename GridView<side>::Traits::template Codim<0>::Entity;

  /** \brief type of grid vertices on side <code>side</code> */
  template<int side>
  using GridVertex = typename GridView<side>::Traits::template Codim<Grid<side>::dimension>::Entity;

  /** \brief The type of the Grid0 elements */
  using Grid0Element DUNE_DEPRECATED_MSG("please use GridElement<0> instead") = GridElement<0>;

  /** \brief The type of the Grid0 vertices */
  using Grid0Vertex DUNE_DEPRECATED_MSG("please use GridVertex<0> instead") = GridVertex<0>;

  /** \brief The type of the Grid1 elements */
  using Grid1Element DUNE_DEPRECATED_MSG("please use GridElement<1> instead") = GridElement<1>;

  /** \brief The type of the Grid1 vertices */
  using Grid1Vertex DUNE_DEPRECATED_MSG("please use GridVertex<1> instead") = GridVertex<1>;

  /** \brief Instance of a Merger */
  typedef Dune::GridGlue::Merger<ctype,
      Grid<0>::dimension - GridPatch<0>::codim,
      Grid<1>::dimension - GridPatch<1>::codim,
      dimworld> Merger;

  /** \brief Type of remote intersection objects */
  typedef Dune::GridGlue::Intersection<P0,P1,0,1> Intersection;

  /** \brief Type of remote intersection indexSet */
  typedef Dune::GridGlue::IntersectionIndexSet<P0,P1> IndexSet;

  /** \brief Type of the iterator that iterates over remove intersections */
  template<int side>
  using IntersectionIterator = Dune::GridGlue::IntersectionIterator<P0, P1, side, (side+1) % 2>;

  /** \todo Please doc me! */
  using Grid0IntersectionIterator DUNE_DEPRECATED_MSG("please use IntersectionIterator<0> instead") = IntersectionIterator<0>;

  /** \todo Please doc me! */
  using Grid1IntersectionIterator DUNE_DEPRECATED_MSG("please use IntersectionIterator<1> instead") = IntersectionIterator<1>;

private:

  /*   M E M B E R   V A R I A B L E S   */

  using GridPatches = std::tuple<
      const std::shared_ptr< const GridPatch<0> >,
      const std::shared_ptr< const GridPatch<1> >
    >;

  const GridPatches patches_;

  /// @brief the surface merging utility
  const std::shared_ptr<Merger> merger_;

  /// @brief number of intersections
  IndexType index__sz = 0;

#if HAVE_MPI
  /// @brief MPI_Comm which this GridGlue is working on
  MPI_Comm mpicomm_;

  /// @brief parallel indexSet for the intersections with a local grid0 entity
  PIndexSet patch0_is_;

  /// @brief parallel indexSet for the intersections with a local grid1 entity
  PIndexSet patch1_is_;

  /// @brief keeps information about which process has which intersection
  Dune::RemoteIndices<PIndexSet> remoteIndices_;
#endif // HAVE_MPI

  /// \todo
  typedef Dune::GridGlue::IntersectionData<P0,P1> IntersectionData;

  /// @brief a vector with intersection elements
  mutable std::vector<IntersectionData>   intersections_;

protected:

  /**
   * @brief after building the merged grid the intersection can be updated
   * through this method (for internal use)
   *
   * @param patch0coords the patch0 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param patch0entities array with all patch0 entities represented as corner indices into @c patch0coords. Free of (potentially heterogeneous) block structure, the last component of one entity is immediately followed by the first component of the next entity here.
   * @param patch0types array with all patch0 entities types
   * @param patch0rank  rank of the process where patch0 was extracted
   *
   * @param patch1coords the patch2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param patch1entities just like with the patch0entities and patch0corners
   * @param patch1types array with all patch1 entities types
   * @param patch1rank  rank of the process where patch1 was extracted
   *
   */
  void mergePatches(const std::vector<Dune::FieldVector<ctype,dimworld> >& patch0coords,
                    const std::vector<unsigned int>& patch0entities,
                    const std::vector<Dune::GeometryType>& patch0types,
                    const int patch0rank,
                    const std::vector<Dune::FieldVector<ctype,dimworld> >& patch1coords,
                    const std::vector<unsigned int>& patch1entities,
                    const std::vector<Dune::GeometryType>& patch1types,
                    const int patch1rank);


  template<typename Extractor>
  void extractGrid (const Extractor & extractor,
                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                    std::vector<unsigned int> & faces,
                    std::vector<Dune::GeometryType>& geometryTypes) const;

public:

  /*   C O N S T R U C T O R S   A N D   D E S T R U C T O R S   */

  /**
   * @brief constructor
   *
   * Initializes components but does not "glue" the surfaces. The surfaces
   * are extracted from the grids here though.
   * @param gp0 the grid0 patch
   * @param gp1 the grid1 patch
   * @param merger The merger object that is used to compute the merged grid. This class has
   * to be a model of the SurfaceMergeConcept.
   */
  GridGlue(const std::shared_ptr< const GridPatch<0> >& gp0, const std::shared_ptr< const GridPatch<1> >& gp1, const std::shared_ptr<Merger>& merger);

  /*   G E T T E R S   */

  /** \todo Please doc me! */
  template<int P>
  const GridPatch<P>& patch() const
  {
    return *std::get<P>(patches_);
  }

  /**
   * @brief getter for the GridView of patch P
   * @return the object
   */
  template<int P>
  const GridView<P>& gridView() const
  {
    return std::get<P>(patches_)->gridView();
  }


  /*   F U N C T I O N A L I T Y   */

  void build();

  /*   I N T E R S E C T I O N S   A N D   I N T E R S E C T I O N   I T E R A T O R S   */

  /**
   * @brief gets an iterator over all remote intersections in the merged grid between grid0 and grid1
   * @tparam I select inside grid I=0 or I=1
   *
   * @return the iterator
   */
  template<int I = 0>
  IntersectionIterator<I> ibegin() const
  {
    return {this, 0};
  }


  /**
   * @brief gets the (general) end-iterator for grid glue iterations
   * @tparam I select inside grid I=0 or I=1
   *
   * @return the iterator
   */
  template<int I = 0>
  IntersectionIterator<I> iend() const
  {
    return {this, index__sz};
  }


  /*! \brief Communicate information on the MergedGrid of a GridGlue

     Template parameter is a model of Dune::GridGlue::CommDataHandle

     \param data GridGlueDataHandle
     \param iftype Interface for which the Communication should take place
     \param dir Communication direction (Forward means grid0 to grid1, Backward is the reverse)

     \todo fix mixed communication: seq->par use commSeq, par->seq use commPar
     \todo add directed version communicate<FROM,TO, DH,DT>(data,iftype,dir)
   */
  template<class DataHandleImp, class DataTypeImp>
  void communicate (Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> & data,
                    Dune::InterfaceType iftype, Dune::CommunicationDirection dir) const;

  /*
   * @brief return an IndexSet mapping from Intersection to IndexType
   */
  IndexSet indexSet() const
  {
    return IndexSet(this);
  }

#if QUICKHACK_INDEX
  // indexset size
  size_t indexSet_size() const
  {
    return index__sz;
  }

#endif

  Intersection getIntersection(int i) const
  {
    return Intersection(this, & intersections_[i]);
  }

  size_t size() const
  {
    return index__sz;
  }

};

} // end namespace GridGlue
} // end namespace Dune

#include "adapter/rangegenerators.hh"

#include "adapter/gridglue.cc"

#include "adapter/intersection.hh"
#include "adapter/intersectioniterator.hh"
#include "adapter/intersectionindexset.hh"

#endif // DUNE_GRIDGLUE_GRIDGLUE_HH
