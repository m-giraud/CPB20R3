// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

template<typename GridView>
typename Dune::GridGlue::Codim1Extractor<GridView>::Predicate
makeHorizontalFacePredicate(double sliceCoord)
{
  using Element = typename GridView::Traits::template Codim<0>::Entity;
  auto predicate = [sliceCoord](const Element& element, unsigned int face) -> bool {
    const int dim = GridView::dimension;
    const auto& refElement = Dune::ReferenceElements<double, dim>::general(element.type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(element.geometry().corner(refElement.subEntity(face,1,i,dim))[1] - sliceCoord) > 1e-6 )
        return false;

    return true;
  };
  return predicate;
};

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

/** \brief trafo from dim to dim+1 */
template<int dim, int dimw, class ctype>
class MixedDimTrafo
  : public AnalyticalCoordFunction< ctype, dim, dimw, MixedDimTrafo<dim,dimw,ctype> >
{
  static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  double yOffset_;
public:
  MixedDimTrafo(double yOffset) : yOffset_(yOffset) {}

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dimw> &y ) const
  {
    y = yOffset_;
    y[0] = x[0];
    for (int i=2; i<dimw; i++)
      y[i] = x[i-1];
  }
};

template<int dim, int dimw, typename ctype = double>
struct Embedding
  : AnalyticalCoordFunction<ctype, dim, dimw, Embedding<dim, dimw, ctype> >
{
  static_assert(dimw >= dim, "Embeddings are only possible into a higher-dimensional space");

  void evaluate(const Dune::FieldVector<ctype, dim>& x, Dune::FieldVector<ctype, dimw>& y) const
  {
    y = ctype(0);
    for (unsigned i = 0; i < dim; ++i)
      y[i] = x[i];
  }
};

template <int dim>
void test1d2dCouplingMatchingDimworld()
{
  double slice = 0.0;

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  using GridType2d = Dune::YaspGrid<dim>;

  std::array<int, dim> elements;
  elements.fill(1);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid0(upper, elements);

  using GridType1d_ = Dune::YaspGrid<dim-1>;

  std::array<int, dim-1> elements1d;
  elements1d.fill(1);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d_ cubeGrid1_(upper1d, elements1d);

  using Embedding1 = Embedding<dim-1, dim>;
  using GridType1d = Dune::GeometryGrid<GridType1d_, Embedding1>;
  Embedding1 embedding1;
  GridType1d cubeGrid1(cubeGrid1_, embedding1);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  typedef typename GridType1d::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeHorizontalFacePredicate<DomGridView>(0);
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelGridView(0), tardesc);
  tarEx->positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  auto merger = std::make_shared< Dune::GridGlue::ContactMerge<dim, double> >();
  merger->minNormalAngle(0.0);

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
void test2d1dCouplingMatchingDimworld()
{
  double slice = 0.0;

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  using GridType1d_ = Dune::YaspGrid<dim-1>;

  std::array<int, dim-1> elements1d;
  elements1d.fill(1);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d_ cubeGrid0_(upper1d, elements1d);

  using Embedding0 = Embedding<dim-1, dim>;
  using GridType1d = Dune::GeometryGrid<GridType1d_, Embedding0>;
  Embedding0 embedding0;
  GridType1d cubeGrid0(cubeGrid0_, embedding0);

  using GridType2d = Dune::YaspGrid<dim>;

  std::array<int, dim> elements;
  elements.fill(1);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid1(upper, elements);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType1d::LevelGridView DomGridView;
  typedef typename GridType2d::LevelGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeHorizontalFacePredicate<TarGridView>(0);

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelGridView(0), tardesc);
  domEx->positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  auto merger = std::make_shared< Dune::GridGlue::ContactMerge<dim, double> >();
  merger->minNormalAngle(0.0);

  GlueType glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}


template <int dim, bool par=false>
void test1d2dCoupling(double slice=0.0)
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  using GridType2d = Dune::YaspGrid<dim>;

  std::array<int, dim> elements;
  elements.fill(1);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid0(upper, elements);

  using GridType1d = Dune::YaspGrid<dim-1>;

  std::array<int, dim-1> elements1d;
  elements1d.fill(1);
  FieldVector<double,dim-1> upper1d(1);

  typedef GeometryGrid<GridType1d, MixedDimTrafo<dim-1,dim,double> > LiftedGridType;

  GridType1d cubeGrid1_in(upper1d, elements1d);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim

  LiftedGridType cubeGrid1(cubeGrid1_in, trafo);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  // typedef typename GridType1d::LevelGridView TarGridView;
  typedef typename LiftedGridType::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeHorizontalFacePredicate<DomGridView>(slice);
  const typename TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelView(0), tardesc);
  tarEx->positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  auto merger = std::make_shared< Dune::GridGlue::ContactMerge<dim, double> >();
  merger->minNormalAngle(0.0);

  GlueType glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}


template <int dim, bool par=false>
void test2d1dCoupling(double slice=0.0)
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  using GridType1d = Dune::YaspGrid<dim-1>;

  std::array<int, dim-1> elements1d;
  elements1d.fill(1);
  FieldVector<double,dim-1> upper1d(1);

  typedef GeometryGrid<GridType1d, MixedDimTrafo<dim-1,dim,double> > LiftedGridType;

  GridType1d cubeGrid0_in(upper1d, elements1d);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim

  LiftedGridType cubeGrid0(cubeGrid0_in, trafo);

  using GridType2d = Dune::YaspGrid<dim>;

  std::array<int, dim> elements;
  elements.fill(1);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid1(upper, elements);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename LiftedGridType::LevelGridView DomGridView;
  typedef typename GridType2d::LevelGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const typename TarExtractor::Predicate tardesc = makeHorizontalFacePredicate<TarGridView>(slice);

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelGridView(0), tardesc);
  domEx->positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  auto merger = std::make_shared< Dune::GridGlue::ContactMerge<dim, double> >();
  merger->minNormalAngle(0.0);

  GlueType glue(domEx, tarEx, merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  // /////////////////////////////////////////////////////////////
  //   First set of tests: the grid have different dimensions,
  //   but the world dimension is the same for both of them.
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d == matching =============================\n";
  test1d2dCouplingMatchingDimworld<2>();
  std::cout << "====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == matching =============================\n";
  test2d1dCouplingMatchingDimworld<2>();
  std::cout << "====================================================\n";

#define _3DTEST 0
#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == matching =============================\n";
  test2d1dCouplingMatchingDimworld<3>();
  std::cout << "====================================================\n";
#endif

  // /////////////////////////////////////////////////////////////
  //   Second set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==0.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d === nonmatching ==========================\n";
  // test1d2dCoupling<2>();
  // test1d2dCoupling<2, true>();
  std::cout << "=====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching ==========================\n";
  test2d1dCoupling<2>();
  test2d1dCoupling<2, true>();
  std::cout << "====================================================\n";

#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching ==========================\n";
  test2d1dCoupling<3>();
  test2d1dCoupling<3, true>();
  std::cout << "====================================================\n";
#endif

  // /////////////////////////////////////////////////////////////
  //   Third set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==1.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d == nonmatching top ======================\n";
  // test1d2dCoupling<2>(1.0);
  // test1d2dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching top ======================\n";
  test2d1dCoupling<2>(1.0);
  test2d1dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching top ======================\n";
  test2d1dCoupling<3>(1.0);
  test2d1dCoupling<3, true>(1.0);
  std::cout << "====================================================\n";
#endif

}
