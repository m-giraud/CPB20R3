#ifndef DUNE_GRIDGLUE_ADAPTER_INTERSECTIONINDEXSET_HH
#define DUNE_GRIDGLUE_ADAPTER_INTERSECTIONINDEXSET_HH

#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/adapter/intersection.hh>

#ifndef ONLY_SIMPLEX_INTERSECTIONS
// we currently support only one intersection type. If we want to support more,
// we have to think about the semantics of our IndexSet
#error Not Implemented
#endif

namespace Dune {
    namespace GridGlue {

        template<typename P0, typename P1>
        class IntersectionIndexSet
        {
            friend class ::Dune::GridGlue::GridGlue<P0,P1>;
            typedef ::Dune::GridGlue::GridGlue<P0,P1> GridGlue;

        public:

            /** \brief The type used for the indices */
            typedef typename GridGlue::IndexType IndexType;
            /** \brief The type used for the size */
            typedef size_t SizeType;

            /** @brief Map Dune::GridGlue::Intersection to index.
            */
            template<int I, int O>
            IndexType index (const Intersection<P0,P1,I,O> & i) const
            {
                return i.i_->index_;
            }

            /** @brief Return total number of intersections.
             */
            SizeType size () const
            {
                return glue_->size();
            }

  private:

            /** construct from a given GridGlue object */
            IntersectionIndexSet(const GridGlue * g) :
                glue_(g) {}

            const GridGlue * glue_;
  };

    } // end namespace GridGlue
} // end namespace Dune

#endif // DUNE_GRIDGLUE_ADAPTER_INTERSECTIONINDEXSET_HH
