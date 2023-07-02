#ifndef DUNE_GRIDGLUE_MERGING_INTERSECTIONLIST_HH
#define DUNE_GRIDGLUE_MERGING_INTERSECTIONLIST_HH 1

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

namespace Dune {
namespace GridGlue {

/**
 * \tparam L0 type of local coordinates in the first grid
 * \tparam L1 type of local coordinates in the second grid
 */
template<typename L0, typename L1>
class IntersectionListProvider
{
public:

  /**
   * type of local coordinates in the first grid
   */
  using Local0 = L0;

  /**
   * type of local coordinates in the second grid
   */
  using Local1 = L1;

  /**
   * type used to index intersections
   */
  using Index = unsigned int;

  /**
   * number of intersections in the merged grid
   */
  virtual std::size_t size() const = 0;

  /**
   * number of embeddings of the `intersection`-th intersection into the first grid
   *
   * \param intersection number of the intersection
   */
  virtual std::size_t parents0(Index intersection) const = 0;

  /**
   * number of embeddings of the `intersection`-th intersection into the second grid
   *
   * \param intersection number of the intersection
   */
  virtual std::size_t parents1(Index intersection) const = 0;

  /**
   * parent entity of an embedding of an intersection in the first grid
   *
   * \param intersection number of the intersection
   * \param index number of the embedding of the intersection
   */
  virtual Index parent0(Index intersection, unsigned index) const = 0;

  /**
   * parent entity of an embedding of an intersection in the second grid
   *
   * \param intersection number of the intersection
   * \param index number of the embedding of the intersection
   */
  virtual Index parent1(Index intersection, unsigned index) const = 0;

  /**
   * corner local coordinates of an embedding of an intersection in the first grid
   *
   * \param intersection number of the intersection
   * \param corner number of the corner
   * \param index number of the embedding of the intersection
   */
  virtual Local0 corner0(Index intersection, unsigned corner, unsigned index) const = 0;

  /**
   * corner local coordinates of an embedding of an intersection in the second grid
   *
   * \param intersection number of the intersection
   * \param corner number of the corner
   * \param index number of the embedding of the intersection
   */
  virtual Local1 corner1(Index intersection, unsigned corner, unsigned index) const = 0;
};

namespace Impl {

template<typename P, int I>
struct IntersectionListLocal
{};

template<typename P>
struct IntersectionListLocal<P, 0>
{
  static std::size_t parents(const P& p, typename P::Index intersection)
    { return p.parents0(intersection); }

  static typename P::Index parent(const P& p, typename P::Index intersection, unsigned index)
    { return p.parent0(intersection, index); }

  static typename P::Local0 corner(const P& p, typename P::Index intersection, unsigned corner, unsigned index)
    { return p.corner0(intersection, corner, index); }
};

template<typename P>
struct IntersectionListLocal<P, 1>
{
  static std::size_t parents(const P& p, typename P::Index intersection)
    { return p.parents1(intersection); }

  static typename P::Index parent(const P& p, typename P::Index intersection, unsigned index)
    { return p.parent1(intersection, index); }

  static typename P::Local1 corner(const P& p, typename P::Index intersection, unsigned corner, unsigned index)
    { return p.corner1(intersection, corner, index); }
};

} /* namespace Impl */

/**
 * \tparam L0 type of local coordinates in the first grid
 * \tparam L1 type of local coordinates in the second grid
 */
template<typename Local0, typename Local1>
class IntersectionList
{
public:
  using Provider = IntersectionListProvider<Local0, Local1>;
  using Index = typename Provider::Index;

  IntersectionList(const std::shared_ptr<Provider>& provider)
    : impl_(provider)
    {}

  /**
   * number of intersections in the merged grid
   */
  std::size_t size() const
    { return impl_->size(); }

  /**
   * number of embeddings of the `intersection`-th intersection into the `I`-th grid
   *
   * \tparam I number of the grid (0 or 1)
   * \param intersection number of the intersection
   */
  template<int I>
  std::size_t parents(Index intersection) const
    {
      static_assert(I == 0 or I == 1, "I must be 0 or 1");
      // TODO [C++17]: use `if constexpr` instead of indirection
      return Impl::IntersectionListLocal<Provider, I>::parents(*impl_, intersection);
    }

  /**
   * parent entity of an embedding of an intersection in the `I`-th grid
   *
   * \tparam I number of the grid (0 or 1)
   * \param intersection number of the intersection
   * \param index number of the embedding of the intersection
   */
  template<int I>
  Index parent(Index intersection, unsigned index = 0) const
    {
      static_assert(I == 0 or I == 1, "I must be 0 or 1");
      // TODO [C++17]: use `if constexpr` instead of indirection
      return Impl::IntersectionListLocal<Provider, I>::parent(*impl_, intersection, index);
    }

  /**
   * corner local coordinates of an embedding of an intersection in the `I`-th grid
   *
   * \tparam I number of the grid (0 or 1)
   * \param intersection number of the intersection
   * \param corner number of the corner
   * \param index number of the embedding of the intersection
   */
  template<int I>
  auto corner(Index intersection, unsigned corner, unsigned index = 0) const
    {
      static_assert(I == 0 or I == 1, "I must be 0 or 1");
      // TODO [C++17]: use `if constexpr` instead of indirection
      return Impl::IntersectionListLocal<Provider, I>::corner(*impl_, intersection, corner, index);
    }

private:
  std::shared_ptr<Provider> impl_;
};

/**
 * intersection list provider storing simplicial intersections
 *
 * \tparam dim0 dimension of the first grid
 * \tparam dim1 dimension of the second grid
 */
template<int dim0, int dim1>
class SimplicialIntersectionListProvider final
  : public IntersectionListProvider< FieldVector<double, dim0>, FieldVector<double, dim1> >
{
  using Base = IntersectionListProvider< FieldVector<double, dim0>, FieldVector<double, dim1> >;

public:
  using Index = typename Base::Index;
  using Local0 = FieldVector<double, dim0>;
  using Local1 = FieldVector<double, dim1>;

  template<int I>
  using Local = std::conditional_t< I == 0, Local0, Local1 >;

  /**
   * simplical intersection
   */
  struct SimplicialIntersection
  {
  private:
    static constexpr int intersectionDim = dim0 < dim1 ? dim0 : dim1;
    static constexpr int nVertices = intersectionDim + 1;

  public:
    SimplicialIntersection() = default;
    SimplicialIntersection(Index parent0, Index parent1)
      : parents0{parent0}
      , parents1{parent1}
      {}

    /**
     * type of a set of corner local coordinates for an embedding
     */
    template<int I>
    using Corners = std::array<Local<I>, nVertices>;

    /**
     * list of sets of corner local coordinates for embeddings into the first grid
     */
    std::vector< Corners<0> > corners0 = std::vector< Corners<0> >(1);

    /**
     * list of parent entities for embeddings into the first grid
     */
    std::vector< Index > parents0 = std::vector< Index >(1);

    /**
     * list of sets of corner local coordinates for embeddings into the second grid
     */
    std::vector< Corners<1> > corners1 = std::vector< Corners<1> >(1);

    /**
     * list of parent entities for embeddings into the second grid
     */
    std::vector< Index > parents1 = std::vector< Index >(1);
  };

  SimplicialIntersectionListProvider() = default;
  SimplicialIntersectionListProvider(std::vector<SimplicialIntersection>&& intersections)
    : intersections_(std::move(intersections))
    {}

  auto& intersections()
    { return intersections_; }

  std::size_t size() const override
    { return intersections_.size(); }

  std::size_t parents0(Index intersection) const override
    { return intersections_[intersection].parents0.size(); }

  std::size_t parents1(Index intersection) const override
    { return intersections_[intersection].parents1.size(); }

  Index parent0(Index intersection, unsigned index) const override
    { return intersections_[intersection].parents0[index]; }

  Index parent1(Index intersection, unsigned index) const override
    { return intersections_[intersection].parents1[index]; }

  Local0 corner0(Index intersection, unsigned corner, unsigned index) const override
    { return intersections_[intersection].corners0[index][corner]; }

  Local1 corner1(Index intersection, unsigned corner, unsigned index) const override
    { return intersections_[intersection].corners1[index][corner]; }

  void clear()
    {
      intersections_.clear();
    }

private:
  std::vector<SimplicialIntersection> intersections_;
};

} /* namespace GridGlue */
} /* namespace Dune */

#endif
