// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file

   @brief Implement iterators over GridGlue intersections
   @author Christian Engwer
 */

#ifndef DUNE_GRIDGLUE_ADAPTER_INTERSECTIONITERATOR_HH
#define DUNE_GRIDGLUE_ADAPTER_INTERSECTIONITERATOR_HH

#include <dune/grid-glue/gridglue.hh>

namespace Dune {
  namespace GridGlue {

    /** @todo documentation */
    template<typename P0, typename P1, int inside, int outside>
    class IntersectionIterator :
      public Dune::ForwardIteratorFacade< IntersectionIterator<P0,P1,inside,outside>,
          const Intersection<P0,P1,inside,outside> >
    {
    public:

      typedef ::Dune::GridGlue::GridGlue<P0, P1> GridGlue;
      typedef ::Dune::GridGlue::Intersection<P0,P1,inside,outside> Intersection;

      IntersectionIterator(const GridGlue * glue, unsigned int i)
        : glue_(glue),
          index_(i),
          intersection_(glue_, & glue_->intersections_[index_])
      {}

      const Intersection& dereference() const
      {
        assert(("never dereference the end iterator" &&
                index_ != glue_->index__sz));
        return intersection_;
      }

      void increment()
      {
        intersection_ = Intersection(glue_, & glue_->intersections_[++index_]);
      }

      bool equals(const IntersectionIterator& iter) const
      {
        return iter.index_ == index_;
      }

    private:

      const GridGlue*   glue_;
      unsigned int index_;

      Intersection intersection_;
    };

  } // end namespace GridGlue
} // end namespace Dune

#endif // DUNE_GRIDGLUE_ADAPTER_INTERSECTIONITERATOR_HH
