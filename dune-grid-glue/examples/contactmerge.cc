#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid-glue/adapter/gridgluevtkwriter.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/contactmerge.hh>

const unsigned dim = 3;
using Coordinates = Dune::EquidistantOffsetCoordinates<double, dim>;
using Grid = Dune::YaspGrid<dim, Coordinates>;
using Element = Grid::Codim<0>::Entity;
using Extractor = Dune::GridGlue::Codim1Extractor<Grid::LeafGridView>;
using GridGlue = Dune::GridGlue::GridGlue<Extractor, Extractor>;
using ContactMerge = Dune::GridGlue::ContactMerge<dim, Grid::ctype>;

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Grid grid0{{0., 0., 0.}, {1., 1., 1.}, {10, 10, 10}};
  Grid grid1{{.12, 0.23, 1.05}, {1.12, 1.23, 2.05}, {10, 10, 10}};

  auto truePredicate = [](const Element&, unsigned int) { return true; };

  auto extractor0 = std::make_shared<Extractor>(grid0.leafGridView(), truePredicate);
  auto extractor1 = std::make_shared<Extractor>(grid1.leafGridView(), truePredicate);

  auto merger = std::make_shared<ContactMerge>();

  GridGlue glue(extractor0, extractor1, merger);
  glue.build();

  Dune::GridGlue::GridGlueVtkWriter::write(glue, "contactmerge");

  return 0;
}
