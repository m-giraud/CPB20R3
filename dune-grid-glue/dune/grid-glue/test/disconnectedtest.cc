#include "config.h"

/** \file
 * \brief Tests a disconnected coupling boundary
 *
 * StandardMerge contains some extra logic to make sure it finds the complete coupling boundary
 * if that boundary has more than one connected component.  Since StandardMerge is an advancing
 * front-type algorithm this does not work automatically.  This test therefore constructs
 * a setting where the coupling boundary consists of two connected components, and tests whether
 * ContactMerge (which is based on StandardMerge) finds both components.
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid-glue/adapter/gridgluevtkwriter.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/merger.hh>
#include <dune/grid-glue/merging/contactmerge.hh>

const int dim = 3;
const int codim = 1;
typedef Dune::UGGrid<dim> Grid;
typedef Grid::LeafGridView GridView;
typedef std::shared_ptr<Grid> GridPtr;
typedef Dune::GridFactory<Grid> GridFactory;
typedef Dune::FieldVector<double, dim> Vector;
typedef Dune::GridGlue::Merger<Grid::ctype, dim-codim, dim-codim, dim> MyMerger;
using Element = Grid::Codim<0>::Entity;

std::function<bool(const Element&, unsigned int)>
make_z_plane_predicate(double z)
{
  return [z](const Element& element, unsigned int subentity) -> bool {
    using std::abs;
    const auto geometry = element.template subEntity<1>(subentity).geometry();
    const auto global = geometry.center();
    const auto epsilon = std::numeric_limits<Element::Geometry::ctype>::epsilon();
    return abs(global[GridView::dimensionworld-1] - z) < epsilon;
  };
}

void insertCube(GridFactory& factory, unsigned int& vertex, const Vector& offset)
{
  std::vector<unsigned int> vertices;

  for (unsigned i(0); i < 1<<dim; ++i) {
    Vector v(offset);
    for (unsigned j(0); j < dim; ++j) {
      if (i & (1 << j))
        v[j] += 1;
    }
    factory.insertVertex(v);
    vertices.push_back(vertex++);
  }

  factory.insertElement(Dune::GeometryTypes::cube(dim), vertices);
}

bool testDisconnected(const std::string& name, std::shared_ptr<MyMerger> merger)
{
  std::cout << "TEST: " << name << std::endl
            << "===================" << std::endl;

  bool pass(true);

  std::array<GridPtr, 2> grids;
  for (unsigned i(0); i < grids.size(); ++i) {
    unsigned int vertex(0);
    GridFactory factory;
    Vector offset(0); offset[dim-1] = i;
    insertCube(factory, vertex, offset);
    offset[0] = 2;
    insertCube(factory, vertex, offset);
    grids[i] = GridPtr(factory.createGrid());
  }

  typedef Dune::GridGlue::Codim1Extractor<GridView> Extractor;
  const Extractor::Predicate predicate = make_z_plane_predicate(1.0);
  std::array<std::shared_ptr<Extractor>, 2> extractors{{
      std::make_shared<Extractor>(grids[0]->leafGridView(), predicate),
      std::make_shared<Extractor>(grids[1]->leafGridView(), predicate),
  }};
  for (unsigned i(0); i < extractors.size(); ++i) {
    auto n = extractors[i]->nCoords();
    if (n != 8) {
      std::cerr << "FAIL: Extracted patch on grid " << i << " has " << n << " coordinates (expected 8)" << std::endl;
      pass = false;
    }
  }

  typedef Dune::GridGlue::GridGlue<Extractor, Extractor> Glue;
  Glue glue(extractors[0], extractors[1], merger);
  glue.build();
  if (glue.size() != 4) {
    std::cerr << "FAIL: " << glue.size() << " remote intersections found (expected 4)." << std::endl;
    pass = false;
  }

  if (!pass) {
    std::string filename("disconnected-");
    filename += name;
    Dune::GridGlue::GridGlueVtkWriter::write<Glue>(glue, filename);
  }

  return pass;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  bool pass(true);

  {
    auto merger = std::make_shared< Dune::GridGlue::ContactMerge<dim> >(0.0);
    merger->enableFallback(true);
    pass &= testDisconnected("ContactMerge", merger);
  }

  return pass ? 0 : 1;
}
