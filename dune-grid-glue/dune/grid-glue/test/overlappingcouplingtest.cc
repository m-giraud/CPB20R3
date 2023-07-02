// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <array>
#include <memory>

#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_HYBRIDTESTGRIDS
#  include <doc/grids/gridfactory/hybridtestgrids.hh>
#endif

#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/gridglue.hh>

#include <dune/grid-glue/merging/conformingmerge.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

/** \brief Returns always true */
template<typename GridView>
typename Dune::GridGlue::Codim0Extractor<GridView>::Predicate
makeTruePredicate()
{
  using Element = typename GridView::Traits::template Codim<0>::Entity;
  auto predicate = [](const Element&, unsigned int) -> bool {
    return true;
  };
  return predicate;
}

template <int dim>
void testCubeGrids(std::shared_ptr< Merger<double,dim,dim,dim> > merger, const FieldVector<double,dim>& gridOffset)
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  using GridType = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;

  std::array<int, dim> elements;
  elements.fill(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(lower, upper, elements);

  lower += gridOffset;
  upper += gridOffset;

  GridType grid1(lower, upper, elements);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(grid0.leafGridView(), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(grid1.leafGridView(), tardesc);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  GlueType glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim>
void testSimplexGrids(std::shared_ptr< Merger<double,dim,dim,dim> > merger, const FieldVector<double,dim>& gridOffset)
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  using GridType = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;

  std::array<int, dim> elements;
  elements.fill(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(lower, upper, elements);

  lower += gridOffset;
  upper += gridOffset;

  GridType grid1(lower, upper, elements);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(grid0.leafGridView(), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(grid1.leafGridView(), tardesc);

  Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


#if HAVE_UG
template <int dim>
void testSimplexGridsUG(std::shared_ptr< Merger<double,dim,dim,dim> > merger, const FieldVector<double,dim>& gridOffset)
{
  // /////////////////////////////////////////////////////////////////
  //   Make two triangle grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lowerLeft(0);
  FieldVector<double,dim> upperRight(1);
  std::array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 10);

  StructuredGridFactory<GridType> factory;
  std::shared_ptr<GridType> grid0 = factory.createSimplexGrid(lowerLeft, upperRight, elements);

  lowerLeft  += gridOffset;
  upperRight += gridOffset;
  std::shared_ptr<GridType> grid1 = factory.createSimplexGrid(lowerLeft, upperRight, elements);

  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(grid0->leafGridView(), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(grid1->leafGridView(), tardesc);

  Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}
#endif

#if HAVE_UG && HAVE_HYBRIDTESTGRIDS
template <int dim>
void testHybridGridsUG(std::shared_ptr< Merger<double,dim,dim,dim> > merger, const FieldVector<double,dim>& gridOffset)
{
  // /////////////////////////////////////////////////////////////////////////
  //   Create the hybrid test grid from dune-grid twice and shift it once
  //   wrt the other grid.
  // /////////////////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  std::unique_ptr<Dune::UGGrid<dim> > grid0(make2DHybridTestGrid<Dune::UGGrid<dim> >());
  std::unique_ptr<Dune::UGGrid<dim> > grid1(make2DHybridTestGrid<Dune::UGGrid<dim> >());


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(grid0->leafGridView(), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(grid1->leafGridView(), tardesc);

  Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}
#endif


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  // //////////////////////////////////////////////////////////
  //   Test with the OverlappingMerge implementation
  // //////////////////////////////////////////////////////////
  auto overlappingMerge1d = std::make_shared< OverlappingMerge<1,1,1,double> >();
  auto overlappingMerge2d = std::make_shared< OverlappingMerge<2,2,2,double> >();

  testCubeGrids<1>(overlappingMerge1d, FieldVector<double,1>(0.05));
  testCubeGrids<2>(overlappingMerge2d, FieldVector<double,2>(0.05));

  testSimplexGrids<1>(overlappingMerge1d, FieldVector<double,1>(0.05));
#if HAVE_UG
  testSimplexGridsUG<2>(overlappingMerge2d, FieldVector<double,2>(0.05));
#endif
#if HAVE_UG && HAVE_HYBRIDTESTGRIDS
  testHybridGridsUG<2>(overlappingMerge2d, FieldVector<double,2>(0.05));
#endif

  // //////////////////////////////////////////////////////////
  //   Test with the ConformingMerge implementation
  // //////////////////////////////////////////////////////////

  auto conformingMerge1d = std::make_shared< ConformingMerge<1,1,double> >();
  auto conformingMerge2d = std::make_shared< ConformingMerge<2,2,double> >();
  auto conformingMerge3d = std::make_shared< ConformingMerge<3,3,double> >();

  testCubeGrids<1>(conformingMerge1d, FieldVector<double,1>(0));
  testCubeGrids<2>(conformingMerge2d, FieldVector<double,2>(0));

  testSimplexGrids<1>(conformingMerge1d, FieldVector<double,1>(0));
#if HAVE_UG
  testSimplexGridsUG<2>(conformingMerge2d, FieldVector<double,2>(0));
#endif
#if HAVE_UG && HAVE_HYBRIDTESTGRIDS
  testHybridGridsUG<2>(conformingMerge2d, FieldVector<double,2>(0));
#endif
}
