// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

/** \brief Returns always true */
template<typename GridView>
typename Codim0Extractor<GridView>::Predicate
makeTruePredicate()
{
  using Element = typename GridView::Traits::template Codim<0>::Entity;
  auto predicate = [](const Element&, unsigned int) -> bool {
    return true;
  };
  return predicate;
}

/** \brief trafo from dim to dim+1 */
template<int dim, int dimw, class ctype>
class MixedDimTrafo
  : public AnalyticalCoordFunction< ctype, dim, dimw, MixedDimTrafo<dim,dimw,ctype> >
{
  static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  static_assert(dim==1, "MixedDimTrafo currently assumes dim==1");
public:

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dimw> &y ) const
  {
    y[0] = x[0]+0.2;
    y[1] = x[0]+0.1;
  }
};


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  const static int dim0 = 2;
  const static int dim1 = 1;
  const static int dimworld = dim0;

  // /////////////////////////////////////////////////////////
  //   Make a 2d unit cube grid and a 1d grid embedded in 2d
  // /////////////////////////////////////////////////////////

  using Grid0 = Dune::YaspGrid<dim0>;

  std::array<int, dim0> elements0;
  elements0.fill(2);
  FieldVector<double, dim0> upper0(1);

  Grid0 grid0(upper0, elements0);

  using Grid1 = Dune::YaspGrid<dim1>;

  std::array<int, dim1> elements1;
  elements1.fill(2);
  FieldVector<double, dim1> upper1(1);

  typedef MixedDimTrafo<dim1, dimworld, double> Transformation;
  typedef GeometryGrid<Grid1, Transformation> LiftedGrid;

  Grid1 cubeGrid1_in(upper1, elements1);

  Transformation trafo;   // transform dim-1 to dim
  LiftedGrid grid1(cubeGrid1_in, trafo);

  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename Grid0::LeafGridView DomGridView;
  typedef typename LiftedGrid::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(grid0.leafGridView(), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(grid1.leafGridView(), tardesc);
  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  // The following code is out-commented, because the test functionality
  // doesn't actually work yet.
  auto merger = std::make_shared< OverlappingMerge<dim0, dim1, dimworld> >();
  GlueType glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}
