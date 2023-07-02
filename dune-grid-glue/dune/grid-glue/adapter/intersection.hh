// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file
   @author Christian Engwer
   @brief Model of the Intersection concept provided by GridGlue.
 */

#ifndef DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH
#define DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH

#include <algorithm>
#include <memory>
#include <tuple>

#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid-glue/gridglue.hh>

#define ONLY_SIMPLEX_INTERSECTIONS

namespace Dune {
  namespace GridGlue {

    // forward declaration
    template<typename P0, typename P1>
    class IntersectionIndexSet;

    /**
       @brief storage class for Dune::GridGlue::Intersection related data
     */
    template<typename P0, typename P1>
    class IntersectionData
    {
    public:
      typedef ::Dune::GridGlue::GridGlue<P0, P1> GridGlue;

      typedef typename GridGlue::IndexType IndexType;

      /** \brief Dimension of the world space of the intersection */
      static constexpr int coorddim = GridGlue::dimworld;

    private:
      // intermediate quantities
      template<int side>
      static constexpr int dim() { return GridGlue::template GridView<side>::Grid::dimension - GridGlue::template GridPatch<side>::codim; }

    public:
      /** \brief Dimension of the intersection */
      static constexpr int mydim = dim<0>() < dim<1>() ? dim<0>() : dim<1>();

      template<int side>
      using GridLocalGeometry = AffineGeometry<
        typename GridGlue::template GridView<side>::ctype, mydim, GridGlue::template GridView<side>::dimension>;

      using Grid0LocalGeometry DUNE_DEPRECATED_MSG("please use GridLocalGeometry<0> instead") = GridLocalGeometry<0>;
      using Grid1LocalGeometry DUNE_DEPRECATED_MSG("please use GridLocalGeometry<1> instead") = GridLocalGeometry<1>;

      template<int side>
      using GridGeometry = AffineGeometry<
        typename GridGlue::template GridView<side>::ctype, mydim, GridGlue::template GridView<side>::dimensionworld>;

      using Grid0Geometry DUNE_DEPRECATED_MSG("please use GridGeometry<0> instead") = GridGeometry<0>;
      using Grid1Geometry DUNE_DEPRECATED_MSG("please use GridGeometry<1> instead") = GridGeometry<1>;

      template<int side>
      using GridIndexType = typename GridGlue::template GridView<side>::IndexSet::IndexType;

      using Grid0IndexType DUNE_DEPRECATED_MSG("please use GridIndexType<0> instead") = GridIndexType<0>;
      using Grid1IndexType DUNE_DEPRECATED_MSG("please use GridIndexType<1> instead") = GridIndexType<1>;

      /** \brief Constructor the n'th IntersectionData of a given GridGlue */
      IntersectionData(const GridGlue& glue, unsigned int mergeindex, unsigned int offset, bool grid0local, bool grid1local);

      /** \brief Default Constructor */
      IntersectionData()
        {
          std::get<0>(sideData_).gridlocal = false;
          std::get<1>(sideData_).gridlocal = false;
        }

      /* Accessor Functions */

      template<int side>
      const GridLocalGeometry<side>& localGeometry(unsigned int parentId = 0) const
        { return *std::get<side>(sideData_).gridlocalgeom[parentId]; }

      template<int side>
      const GridGeometry<side>& geometry() const
        { return *std::get<side>(sideData_).gridgeom; }

      template<int side>
      bool local() const
        { return std::get<side>(sideData_).gridlocal; }

      template<int side>
      IndexType index(unsigned int parentId = 0) const
        { return std::get<side>(sideData_).gridindices[parentId]; }

      template<int side>
      IndexType parents() const
        { return std::get<side>(sideData_).gridindices.size(); }

    private:
      template<int side>
      void initializeGeometry(const GridGlue& glue, unsigned mergeindex);

      /*   M E M B E R   V A R I A B L E S   */

    public:
      /// @brief index of this intersection after GridGlue interface
      IndexType index_;

    private:
      template<int side>
      struct SideData {
        /** \brief true if the associated grid entity is local */
        bool gridlocal;

        /** \brief indices of the associated local grid entity */
        std::vector< GridIndexType<side> > gridindices;

        /** \brief embedding of intersection into local grid entity coordinates */
        /* TODO [C++17 or DUNE-2.6]: use std::optional */
        std::vector< std::unique_ptr< GridLocalGeometry<side> > > gridlocalgeom;

        /**
         * global intersection geometry on grid `side` side.
         *
         * This is the same as g∘i for the first embedding i into a grid
         * entity as stored in gridlocalgeom and that entities global
         * geometry g.
         */
        /* TODO [C++17 or DUNE-2.6]: use std::optional */
        std::unique_ptr< GridGeometry<side> > gridgeom;
      };

      std::tuple< SideData<0>, SideData<1> > sideData_;
    };

    template<typename P0, typename P1>
    template<int side>
    void IntersectionData<P0, P1>::initializeGeometry(const GridGlue& glue, unsigned mergeindex)
    {
      auto& data = std::get<side>(sideData_);

      const unsigned n_parents = glue.merger_->template parents<side>(mergeindex);

      // init containers
      data.gridindices.resize(n_parents);
      data.gridlocalgeom.resize(n_parents);

      // default values
      data.gridindices[0] = 0;

      static constexpr int nSimplexCorners = mydim + 1;
      using ctype = typename GridGlue::ctype;

      // initialize the local and the global geometries of grid `side`

      // compute the coordinates of the subface's corners in codim 0 entity local coordinates
      static constexpr int elementdim = GridGlue::template GridView<side>::template Codim<0>::Geometry::mydimension;

      // coordinates within the subentity that contains the remote intersection
      std::array<Dune::FieldVector< ctype, dim<side>() >, nSimplexCorners> corners_subEntity_local;

      for (unsigned int par = 0; par < n_parents; ++par) {
        for (int i = 0; i < nSimplexCorners; ++i)
          corners_subEntity_local[i] = glue.merger_->template parentLocal<side>(mergeindex, i, par);

        // Coordinates of the remote intersection corners wrt the element coordinate system
        std::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

        if (data.gridlocal) {
          data.gridindices[par] = glue.merger_->template parent<side>(mergeindex,par);

          typename GridGlue::template GridPatch<side>::LocalGeometry
            gridLocalGeometry = glue.template patch<side>().geometryLocal(data.gridindices[par]);
          for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
            corners_element_local[i] = gridLocalGeometry.global(corners_subEntity_local[i]);
          }

          // set the corners of the local geometry
#ifdef ONLY_SIMPLEX_INTERSECTIONS
#  if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
          const Dune::GeometryType type = Dune::GeometryTypes::simplex(mydim);
#  else
          const Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
#  endif
#else
#error Not Implemented
#endif
          data.gridlocalgeom[par] = std::make_unique< GridLocalGeometry<side> >(type, corners_element_local);

          // Add world geometry only for 0th parent
          if (par == 0) {
            typename GridGlue::template GridPatch<side>::Geometry
              gridWorldGeometry = glue.template patch<side>().geometry(data.gridindices[par]);

            // world coordinates of the remote intersection corners
            std::array<Dune::FieldVector<ctype, GridGlue::template GridView<side>::dimensionworld>, nSimplexCorners> corners_global;

            for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
              corners_global[i]        = gridWorldGeometry.global(corners_subEntity_local[i]);
            }

            data.gridgeom = std::make_unique< GridGeometry<side> >(type, corners_global);
          }
        }
      }
    }

    //! \todo move this functionality to GridGlue
    template<typename P0, typename P1>
    IntersectionData<P0, P1>::IntersectionData(const GridGlue& glue, unsigned int mergeindex, unsigned int offset,
                                               bool grid0local, bool grid1local)
      : index_(mergeindex+offset)
    {
      // if an invalid index is given do not proceed!
      // (happens when the parent GridGlue initializes the "end"-Intersection)
      assert (0 <= mergeindex || mergeindex < glue.index__sz);

      std::get<0>(sideData_).gridlocal = grid0local;
      std::get<1>(sideData_).gridlocal = grid1local;

      initializeGeometry<0>(glue, mergeindex);
      initializeGeometry<1>(glue, mergeindex);
    }

    /**
       @brief
       @todo doc me
     */
    template<typename P0, typename P1, int inside, int outside>
    struct IntersectionTraits
    {
      using GridGlue = ::Dune::GridGlue::GridGlue<P0, P1>;
      using IntersectionData = Dune::GridGlue::IntersectionData<P0, P1>;

      using InsideGridView = typename GridGlue::template GridView<inside>;
      using OutsideGridView = typename GridGlue::template GridView<outside>;

      using InsideLocalGeometry = typename IntersectionData::template GridLocalGeometry<inside>;
      using OutsideLocalGeometry = typename IntersectionData::template GridLocalGeometry<outside>;

      using Geometry = typename IntersectionData::template GridGeometry<inside>;
      using OutsideGeometry = typename IntersectionData::template GridGeometry<outside>;

      static constexpr auto coorddim = IntersectionData::coorddim;
      static constexpr auto mydim = IntersectionData::mydim;
      static constexpr int insidePatch = inside;
      static constexpr int outsidePatch = outside;

      using ctype = typename GridGlue::ctype;
      using LocalCoordinate = Dune::FieldVector<ctype, mydim>;
      using GlobalCoordinate = Dune::FieldVector<ctype, coorddim>;
    };

    /** @brief The intersection of two entities of the two patches of a GridGlue
     */
    template<typename P0, typename P1, int I, int O>
    class Intersection
    {

    public:

      typedef IntersectionTraits<P0,P1,I,O> Traits;

      typedef typename Traits::GridGlue GridGlue;
      typedef typename Traits::IntersectionData IntersectionData;


      typedef typename Traits::InsideGridView InsideGridView;
      typedef typename Traits::InsideLocalGeometry InsideLocalGeometry;

      typedef typename Traits::OutsideGridView OutsideGridView;
      typedef typename Traits::OutsideLocalGeometry OutsideLocalGeometry;
      typedef typename Traits::OutsideGeometry OutsideGeometry;

      typedef typename Traits::Geometry Geometry;
      typedef typename Traits::ctype ctype;

      typedef typename InsideGridView::Traits::template Codim<0>::Entity InsideEntity;
      typedef typename OutsideGridView::Traits::template Codim<0>::Entity OutsideEntity;

      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      /** \brief dimension of the world space of the intersection */
      static constexpr auto coorddim = Traits::coorddim;

      /** \brief dimension of the intersection */
      static constexpr auto mydim = Traits::mydim;

      /** \brief inside patch */
      static constexpr int insidePatch = Traits::insidePatch;

      /** \brief outside patch */
      static constexpr int outsidePatch = Traits::outsidePatch;

      // typedef unsigned int IndexType;

    private:
      /**
       * \brief codimension of the intersection with respect to the world dimension
       */
      static constexpr int codimensionWorld = coorddim - mydim;

    public:
      /*   C O N S T R U C T O R S   */

      /** \brief Constructor for a given Dataset */
      Intersection(const GridGlue* glue, const IntersectionData*  i) :
        glue_(glue), i_(i) {}

      /*   F U N C T I O N A L I T Y   */

      /** \brief Return element on the inside of this intersection.
       */
      InsideEntity
      inside(unsigned int parentId = 0) const
      {
        assert(self());
        return glue_->template patch<I>().element(i_->template index<I>(parentId));
      }

      /** \brief Return element on the outside of this intersection.
       */
      OutsideEntity
      outside(unsigned int parentId = 0) const
      {
        assert(neighbor());
        return glue_->template patch<O>().element(i_->template index<O>(parentId));
      }

      /** \brief Return true if intersection is conforming */
      bool conforming() const
      {
        throw Dune::NotImplemented();
      }

      /** \brief Geometric information about this intersection in local coordinates of the inside() element.
       */
      const InsideLocalGeometry& geometryInInside(unsigned int parentId = 0) const
      {
        return i_->template localGeometry<I>(parentId);
      }

      /** \brief Geometric information about this intersection in local coordinates of the outside() element.
       */
      const OutsideLocalGeometry& geometryInOutside(unsigned int parentId = 0) const
      {
        return i_->template localGeometry<O>(parentId);
      }

      /** \brief Geometric information about this intersection as part of the inside grid.
       *
       * This is the same geometry as the application of the first
       * embedding into the "inside" entity and then this entities
       * global geometry.
       */
      const Geometry& geometry() const
      {
        return i_->template geometry<I>();
      }

      /** \brief Geometric information about this intersection as part of the outside grid.
       *
       * This is the same geometry as the application of the first
       * embedding into the "outside" entity and then this entities
       * global geometry.
       */
      const OutsideGeometry& geometryOutside() const // DUNE_DEPRECATED
      {
        return i_->template geometry<O>();
      }

      /** \brief Type of reference element for this intersection */
      Dune::GeometryType type() const
      {
        #ifdef ONLY_SIMPLEX_INTERSECTIONS
        #  if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
        return Dune::GeometryTypes::simplex(mydim);
        #  else
        static const Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
        return type;
        #  endif
        #else
        #error Not Implemented
        #endif
      }


      /** \brief For parallel computations: Return true if inside() entity exists locally */
      bool self() const
      {
        return i_->template local<I>();
      }

      /** \brief Return number of embeddings into local grid0 (grid1) entities. */
      size_t neighbor(unsigned int g = 0) const
      {
          if (g == 0 && i_->template local<O>()) {
            return i_->template parents<O>();
          } else if (g == 1 && i_->template local<I>()) {
            return i_->template parents<I>();
          }
          return 0;
      }

      /** \brief Local number of codim 1 entity in the inside() Entity where intersection is contained in. */
      int indexInInside(unsigned int parentId = 0) const
      {
        assert(self());
        return glue_->template patch<I>().indexInInside(i_->template index<I>(parentId));
      }

      /** \brief Local number of codim 1 entity in outside() Entity where intersection is contained in. */
      int indexInOutside(unsigned int parentId = 0) const
      {
        assert(neighbor());
        return glue_->template patch<O>().indexInInside(i_->template index<O>(parentId));
      }

      /** \brief Return an outer normal (length not necessarily 1)
       *
       * The outer normal is given with respect to the \ref geometry().
       */
      GlobalCoordinate outerNormal(const LocalCoordinate &local) const
      {
        GlobalCoordinate normal;

        if (codimensionWorld == 0)
          DUNE_THROW(Dune::Exception, "There is no normal vector to a full-dimensional intersection");
        else if (codimensionWorld == 1) {
          /* TODO: Implement the general n-ary cross product here */
          const auto jacobianTransposed = geometry().jacobianTransposed(local);
          if (mydim==1) {
            normal[0] = - jacobianTransposed[0][1];
            normal[1] =   jacobianTransposed[0][0];
          } else if (mydim==2) {
            normal[0] =   (jacobianTransposed[0][1] * jacobianTransposed[1][2] - jacobianTransposed[0][2] * jacobianTransposed[1][1]);
            normal[1] = - (jacobianTransposed[0][0] * jacobianTransposed[1][2] - jacobianTransposed[0][2] * jacobianTransposed[1][0]);
            normal[2] =   (jacobianTransposed[0][0] * jacobianTransposed[1][1] - jacobianTransposed[0][1] * jacobianTransposed[1][0]);
          } else
            DUNE_THROW(Dune::NotImplemented, "Remote intersections don't implement the 'outerNormal' method for " << mydim << "-dimensional intersections yet");
        } else
          DUNE_THROW(Dune::NotImplemented, "Remote intersections don't implement the 'outerNormal' method for intersections with codim >= 2 yet");

        return normal;
      }

      /** \brief Return a unit outer normal
       *
       * The outer normal is given with respect to the \ref geometry().
       */
      GlobalCoordinate unitOuterNormal(const LocalCoordinate &local) const
      {
        GlobalCoordinate normal = outerNormal(local);
        normal /= normal.two_norm();
        return normal;
      }

      /** \brief Return an outer normal with the length of the integration element
       *
       * The outer normal is given with respect to the \ref geometry().
       */
      GlobalCoordinate integrationOuterNormal(const LocalCoordinate &local) const
      {
        return (unitOuterNormal(local) *= geometry().integrationElement(local));
      }

      /** \brief Unit outer normal at the center of the intersection
       *
       * The outer normal is given with respect to the \ref geometry().
       */
      GlobalCoordinate centerUnitOuterNormal () const
      {
        return unitOuterNormal(ReferenceElements<ctype,mydim>::general(type()).position(0,0));
      }

      /**
       * \brief Return a copy of the intersection with inside and outside switched.
       */
      Intersection<P0,P1,O,I> flip() const
      {
        return Intersection<P0,P1,O,I>(glue_,i_);
      }

#ifdef QUICKHACK_INDEX
      typedef typename GridGlue::IndexType IndexType;

      IndexType index() const
      {
        return i_->index_;
      }

#endif

    private:

      friend class IntersectionIndexSet<P0,P1>;

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief the grid glue entity this is built on
      const GridGlue*       glue_;

      /// @brief the underlying remote intersection
      const IntersectionData* i_;
    };


  } // end namespace GridGlue
} // end namespace Dune

#endif // DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH
