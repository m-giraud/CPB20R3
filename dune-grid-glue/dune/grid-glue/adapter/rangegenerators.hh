#ifndef DUNE_GRIDGLUE_ADAPTER_RANGEGENERATORS_HH
#define DUNE_GRIDGLUE_ADAPTER_RANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>

namespace Dune {
namespace GridGlue {

/**
 * Static tag representing reversal of in- and outside of intersecions.
 */
template<bool reverse>
struct Reverse
  : std::integral_constant<bool, reverse>
{
  typedef Reverse type;

  constexpr
  Reverse<!reverse> operator!() const
    { return {}; }
};

#ifdef DOXYGEN

/**
 * Static tag representing reversal of in- and outside of intersections.
 * \relates Reverse
 */
const Reverse<true> reversed;

/**
 * \brief Iterate over all intersections of a GridGlue.
 *
 * This function returns an object representing the range of intersections
 * with respect to the GridGlue glue. Its main purpose is to enable iteration
 * over these intersections by means of a range-based for loop:
 *
 * \code
 * // Iterate over all intersections of a GridGlue in various ways
 * using Dune::GridGlue::GridGlue;
 * using Dune::GridGlue::Reverse;
 * using Dune::GridGlue::reversed;
 *
 * GridGlue<...> glue;
 * for (const auto& in : intersections(glue)) { ... }
 * for (const auto& in : intersections(glue, reversed)) { ... }
 * for (const auto& in : intersections(glue, !reversed)) { ... }
 * for (const auto& in : intersections(glue, Reversed<true>())) { ... }
 * \endcode
 *
 * The in- and outside of the intersection can be reversed by passing
 * `reversed` as the second argument. The fourth form can be used in
 * case a template parameter for reversal is required.
 *
 * \since dune-common 2.4
 * \relatesalso Dune::GridGlue::GridGlue
 * \param glue GridGlue to obtain the intersections from
 * \param reverse Tag to indicate reversal of in- and outside of intersections
 * \returns an unspecified object that is guaranteed to fulfill the interface
 *          of IteratorRange and that can be iterated over using a range-based
 *          for loop.
 * \see Dune::GridGlue::Intersection
 */
template<...>
IteratorRange<...>
intersections(const GridGlue<...>& glue, const Reverse<...>& reverse = !reversed);

#else

namespace {
const Reverse<true> reversed = {};
} /* namespace */

template<typename P0, typename P1, bool reverse = false>
IteratorRange< typename GridGlue<P0, P1>::template IntersectionIterator<reverse ? 1 : 0> >
intersections(const GridGlue<P0, P1>& glue, const Reverse<reverse>& = {})
{
  const static int side = reverse ? 1 : 0;
  return {glue.template ibegin<side>(), glue.template iend<side>()};
}

#endif // DOXYGEN

} /* namespace GridGlue */
} /* namespace Dune */

#endif
