#include "config.h"

#include <array>
#include <iostream>
#include <vector>

#include <dune/grid-glue/merging/contactmerge.hh>

const static int dim = 3;
typedef double Real;

struct MyMerger
  : public Dune::GridGlue::ContactMerge<3, Real>
{
  typedef Dune::GridGlue::ContactMerge<3, Real> Base;
  using Base::computeCyclicOrder;
};

typedef Dune::FieldVector<Real, dim-1> Local;
typedef std::array<Local, 2> Corner;
typedef std::vector<Corner> Corners;

bool
testComputeCyclicOrder()
{
  bool pass = true;

  Corners corners;
  for (unsigned i = 0; i < 4; ++i) {
    Local local;
    switch (i) {
    case 0:
      local[0] = -0.25; local[1] = 1;
      break;
    case 1:
      local[0] = 1; local[1] = 0;
      break;
    case 2:
      local[0] = 0; local[1] = 0;
      break;
    case 3:
      local[0] = 0.75; local[1] = 1;
      break;
    }
    corners.push_back(Corner{{local, local}});
  }

  Local center(0);
  for (const auto& c : corners)
    center += c[0];
  center /= corners.size();

  std::vector<int> ordering;

  MyMerger merger;
  merger.computeCyclicOrder(corners, center, ordering);

  /* expected cycle is 0 -> 3 -> 1 -> 2, but we don't know where it starts */
  std::vector<int> next{3, 2, 0, 1};
  for (std::size_t i = 0; i < ordering.size(); ++i) {
    const int current = ordering[i];
    const int expected = next[current];
    const int got = ordering[(i+1) % ordering.size()];
    if (got != expected) {
      std::cerr << "FAIL: computeCyclicOrder: " << got << " follows "
                << i << ", but expected " << expected << std::endl;
      pass = false;
    }
  }

  return pass;
}

int
main()
{
  bool pass = true;
  pass &= testComputeCyclicOrder();
  return pass ? 0 : 1;
}
