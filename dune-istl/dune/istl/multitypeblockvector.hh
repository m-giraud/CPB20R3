// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_MULTITYPEBLOCKVECTOR_HH
#define DUNE_ISTL_MULTITYPEBLOCKVECTOR_HH

#include <cmath>
#include <iostream>
#include <tuple>

#include <dune/common/dotproduct.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>

#include "istlexception.hh"

// forward declaration
namespace Dune {
  template < typename... Args >
  class MultiTypeBlockVector;
}

#include "gsetc.hh"

namespace Dune {

  /**
      @addtogroup ISTL_SPMV
      @{
   */




  /** @addtogroup DenseMatVec
      @{
   */
  template <typename Arg0, typename... Args>
  struct FieldTraits< MultiTypeBlockVector<Arg0, Args...> >
  {
    typedef typename FieldTraits<Arg0>::field_type field_type;
    typedef typename FieldTraits<Arg0>::real_type real_type;
  };
  /**
      @}
   */

  /**
      @brief A Vector class to support different block types

      This vector class combines elements of different types known at compile-time.
   */
  template < typename... Args >
  class MultiTypeBlockVector
  : public std::tuple<Args...>
  {
    /** \brief Helper type */
    typedef std::tuple<Args...> TupleType;
  public:

    /**
     * \brief Get the constructors from tuple
     */
    using TupleType::TupleType;

    /**
     * own class' type
     */
    typedef MultiTypeBlockVector<Args...> type;

    /** \brief The type used for scalars
     *
     * The current code hardwires it to 'double', which is far from nice.
     * On the other hand, it is not clear what the correct type is.  If the MultiTypeBlockVector class
     * is instantiated with several vectors of different field_types, what should the resulting
     * field_type be?
     */
    typedef double field_type;

    /** \brief Return the number of vector entries */
    static constexpr std::size_t size()
    {
      return sizeof...(Args);
    }

    /**
     * number of elements
     */
    int count() const
    {
      return sizeof...(Args);
    }

    /** \brief Random-access operator
     *
     * This method mimicks the behavior of normal vector access with square brackets like, e.g., v[5] = 1.
     * The problem is that the return type is different for each value of the argument in the brackets.
     * Therefore we implement a trick using std::integral_constant.  To access the first entry of
     * a MultiTypeBlockVector named v write
     * \code
     *  MultiTypeBlockVector<A,B,C,D> v;
     *  std::integral_constant<std::size_t,0> _0;
     *  v[_0] = ...
     * \endcode
     * The name '_0' used here as a static replacement of the integer number zero is arbitrary.
     * Any other variable name can be used.  If you don't like the separate variable, you can writee
     * \code
     *  MultiTypeBlockVector<A,B,C,D> v;
     *  v[std::integral_constant<std::size_t,0>()] = ...
     * \endcode
     */
    template< std::size_t index >
    typename std::tuple_element<index,TupleType>::type&
    operator[] ( const std::integral_constant< std::size_t, index > indexVariable )
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }

    /** \brief Const random-access operator
     *
     * This is the const version of the random-access operator.  See the non-const version for a full
     * explanation of how to use it.
     */
    template< std::size_t index >
    const typename std::tuple_element<index,TupleType>::type&
    operator[] ( const std::integral_constant< std::size_t, index > indexVariable ) const
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }

    /** \brief Assignment operator
     */
    template<typename T>
    void operator= (const T& newval) {
      Dune::Hybrid::forEach(*this, [&](auto&& entry) {
        entry = newval;
      });
    }

    /**
     * operator for MultiTypeBlockVector += MultiTypeBlockVector operations
     */
    void operator+= (const type& newv) {
      using namespace Dune::Hybrid;
      forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) {
        (*this)[i] += newv[i];
      });
    }

    /**
     * operator for MultiTypeBlockVector -= MultiTypeBlockVector operations
     */
    void operator-= (const type& newv) {
      using namespace Dune::Hybrid;
      forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) {
        (*this)[i] -= newv[i];
      });
    }

    // Once we have the IsNumber traits class the following
    // three implementations could be replaced by:
    //
    //    template<class T,
    //      std::enable_if_t< IsNumber<T>::value, int> = 0>
    //    void operator*= (const T& w) {
    //      Dune::Hybrid::forEach(*this, [&](auto&& entry) {
    //        entry *= w;
    //      });
    //    }

    void operator*= (const int& w) {
      Dune::Hybrid::forEach(*this, [&](auto&& entry) {
        entry *= w;
      });
    }

    void operator*= (const float& w) {
      Dune::Hybrid::forEach(*this, [&](auto&& entry) {
        entry *= w;
      });
    }

    void operator*= (const double& w) {
      Dune::Hybrid::forEach(*this, [&](auto&& entry) {
        entry *= w;
      });
    }

    field_type operator* (const type& newv) const {
      using namespace Dune::Hybrid;
      return accumulate(integralRange(Hybrid::size(*this)), field_type(0), [&](auto&& a, auto&& i) {
        return a + (*this)[i]*newv[i];
      });
    }

    field_type dot (const type& newv) const {
      using namespace Dune::Hybrid;
      return accumulate(integralRange(Hybrid::size(*this)), field_type(0), [&](auto&& a, auto&& i) {
        return a + (*this)[i].dot(newv[i]);
      });
    }

    /** \brief Compute the squared Euclidean norm
     */
    typename FieldTraits<field_type>::real_type two_norm2() const {
      using namespace Dune::Hybrid;
      return accumulate(*this, typename FieldTraits<field_type>::real_type(0), [&](auto&& a, auto&& entry) {
        return a + entry.two_norm2();
      });
    }

    /** \brief Compute the Euclidean norm
     */
    typename FieldTraits<field_type>::real_type two_norm() const {return sqrt(this->two_norm2());}

    /** \brief Compute the maximum norm
     */
    typename FieldTraits<field_type>::real_type infinity_norm() const
    {
      using namespace Dune::Hybrid;
      using std::max;
      using real_type = typename FieldTraits<field_type>::real_type;

      real_type result = 0.0;
      // Compute max norm tracking appearing nan values
      // if the field type supports nan.
      ifElse(has_nan<field_type>(), [&](auto&& id) {
        // This variable will preserve any nan value
        real_type nanTracker = 1.0;
        forEach(*this, [&](auto&& entry) {
          real_type entryNorm = entry.infinity_norm();
          result = max(entryNorm, result);
          nanTracker += entryNorm;
        });
        // Incorporate possible nan value into result
        result *= (nanTracker / nanTracker);
      }, [&](auto&& id) {
        forEach(*this, [&](auto&& entry) {
          result = max(entry.infinity_norm(), result);
        });
      });
      return result;
    }

    /** \brief Axpy operation on this vector (*this += a * y)
     *
     * \tparam Ta Type of the scalar 'a'
     */
    template<typename Ta>
    void axpy (const Ta& a, const type& y) {
      using namespace Dune::Hybrid;
      forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) {
        (*this)[i].axpy(a, y[i]);
      });
    }

  };



  /** \brief Send MultiTypeBlockVector to an outstream
   */
  template <typename... Args>
  std::ostream& operator<< (std::ostream& s, const MultiTypeBlockVector<Args...>& v) {
    using namespace Dune::Hybrid;
    forEach(integralRange(Dune::Hybrid::size(v)), [&](auto&& i) {
      s << "\t(" << i << "):\n" << v[i] << "\n";
    });
    return s;
  }



} // end namespace

#endif
