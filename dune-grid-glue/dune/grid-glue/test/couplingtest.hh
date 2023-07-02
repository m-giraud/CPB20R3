// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_TEST_COUPLINGTEST_HH
#define DUNE_GRIDGLUE_TEST_COUPLINGTEST_HH

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/grid-glue/gridglue.hh>

template <class IntersectionIt>
bool testIntersection(const IntersectionIt & rIIt, double eps)
{
  bool success = true;

  typedef typename IntersectionIt::value_type Intersection;
  // Dimension of the intersection
  const int dim = Intersection::mydim;

  // Dimension of world coordinates
  const int coorddim = Intersection::coorddim;

  // Create a set of test points
  const Dune::QuadratureRule<double, dim>& quad = Dune::QuadratureRules<double, dim>::rule(rIIt->type(), 3);

  for (unsigned int l=0; l<quad.size(); l++) {
    const auto inside = rIIt->inside();
    const auto outside = rIIt->outside();

    Dune::FieldVector<double, dim> quadPos = quad[l].position();

    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> localGrid0Pos =
      inside.geometry().global(rIIt->geometryInInside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> globalGrid0Pos =
      rIIt->geometry().global(quadPos);

    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> localGrid1Pos =
      outside.geometry().global(rIIt->geometryInOutside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> globalGrid1Pos =
      rIIt->geometryOutside().global(quadPos);

    // Test whether local grid0 position is consistent with global grid0 position
    if ( (localGrid0Pos-globalGrid0Pos).two_norm() >= eps )
    {
      std::cout << __FILE__ << ":" << __LINE__ << ": error: assert( (localGrid0Pos-globalGrid0Pos).two_norm() < eps ) failed\n";
      std::cerr << "localGrid0Pos  = " << localGrid0Pos << "\n";
      std::cerr << "globalGrid0Pos = " << globalGrid0Pos << "\n";
      success = false;
    }

    // Test whether local grid1 position is consistent with global grid1 position
    if ( (localGrid1Pos-globalGrid1Pos).two_norm() >= eps )
    {
      std::cout << __FILE__ << ":" << __LINE__ << ": error: assert( (localGrid1Pos-globalGrid1Pos).two_norm() < eps ) failed\n";
      std::cerr << "localGrid1Pos  = " << localGrid1Pos << "\n";
      std::cerr << "globalGrid1Pos = " << globalGrid1Pos << "\n";
      success = false;
    }

    // Here we assume that the two interfaces match geometrically:
    if ( (globalGrid0Pos-globalGrid1Pos).two_norm() >= eps )
    {
      std::cout << __FILE__ << ":" << __LINE__ << ": error: assert( (globalGrid0Pos-globalGrid1Pos).two_norm() < eps ) failed\n";
      std::cerr << "localGrid0Pos  = " << localGrid0Pos << "\n";
      std::cerr << "globalGrid0Pos = " << globalGrid0Pos << "\n";
      std::cerr << "localGrid1Pos  = " << localGrid1Pos << "\n";
      std::cerr << "globalGrid1Pos = " << globalGrid1Pos << "\n";
      success = false;
    }

    // Test the normal vector methods.  At least test whether they don't crash
    if (coorddim - dim == 1) // only test for codim 1
    {
      rIIt->outerNormal(quadPos);
      rIIt->unitOuterNormal(quadPos);
      rIIt->integrationOuterNormal(quadPos);
      rIIt->centerUnitOuterNormal();
    }
  }

  return success;
}


template <class GlueType>
void testCoupling(const GlueType& glue, double eps = 1e-12)
{
  bool success = true;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 6)
  using View0Mapper = Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::template GridView<0> >;
  using View1Mapper = Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::template GridView<1> >;
  View0Mapper view0mapper(glue.template gridView<0>(), Dune::mcmgElementLayout());
  View1Mapper view1mapper(glue.template gridView<1>(), Dune::mcmgElementLayout());
#else
  using View0Mapper = Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::template GridView<0>, Dune::MCMGElementLayout >;
  using View1Mapper = Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::template GridView<1>, Dune::MCMGElementLayout >;
  View0Mapper view0mapper(glue.template gridView<0>());
  View1Mapper view1mapper(glue.template gridView<1>());
#endif

  std::vector<unsigned int> countInside0(view0mapper.size());
  std::vector<unsigned int> countOutside1(view1mapper.size());
  std::vector<unsigned int> countInside1(view1mapper.size(), 0);
  std::vector<unsigned int> countOutside0(view0mapper.size(), 0);

  // ///////////////////////////////////////
  //   IndexSet
  // ///////////////////////////////////////

  {
    size_t count = 0;
    for (auto rIIt = glue.template ibegin<0>(); rIIt != glue.template iend<0>(); ++rIIt) count ++;
    typename GlueType::IndexSet is = glue.indexSet();
    if(is.size() != glue.size())
      DUNE_THROW(Dune::Exception,
        "Inconsistent size information: indexSet.size() " << is.size() << " != GridGlue.size() " << glue.size());
    if(is.size() != count)
      DUNE_THROW(Dune::Exception,
        "Inconsistent size information: indexSet.size() " << is.size() << " != iterator count " << count);
    std::vector<bool> visited(count, false);
    for (auto rIIt = glue.template ibegin<0>(); rIIt != glue.template iend<0>(); ++rIIt) {
      size_t idx = is.index(*rIIt);
      if(idx >= count)
        DUNE_THROW(Dune::Exception,
          "Inconsistent IndexSet: index " << idx << " out of range, size is " << count);
      if(visited[idx] != false)
        DUNE_THROW(Dune::Exception,
          "Inconsistent IndexSet: visited index " << idx << " twice");
      visited[idx] = true;
    }
    // make sure that we have a consecutive zero starting index set
    for (size_t i = 0; i<count; i++)
    {
      if (visited[i] != true)
        DUNE_THROW(Dune::Exception,
          "Non-consective IndexSet: " << i << " missing.");
    }
  }


  // ///////////////////////////////////////
  //   MergedGrid centric Grid0->Grid1
  // ///////////////////////////////////////

  {
    for (auto rIIt = glue.template ibegin<0>(); rIIt != glue.template iend<0>(); ++rIIt)
    {
      if (rIIt->self() && rIIt->neighbor())
      {
        const auto index0 = view0mapper.index(rIIt->inside());
        const auto index1 = view1mapper.index(rIIt->outside());

        countInside0[index0]++;
        countOutside1[index1]++;
        success = success && testIntersection(rIIt, eps);
      }
    }
  }

  // ///////////////////////////////////////
  //   MergedGrid centric Grid1->Grid0
  // ///////////////////////////////////////

  {
    for (auto rIIt = glue.template ibegin<1>(); rIIt != glue.template iend<1>(); ++rIIt)
    {
      if (rIIt->self() && rIIt->neighbor())
      {
        const auto index1 = view1mapper.index(rIIt->inside());
        const auto index0 = view0mapper.index(rIIt->outside());

        countInside1[index1]++;
        countOutside0[index0]++;
        success = success && testIntersection(rIIt, eps);
      }
    }
  }

  if (! success)
    DUNE_THROW(Dune::Exception, "Test failed, see above for details.");
}

#endif // DUNE_GRIDGLUE_TEST_COUPLINGTEST_HH
