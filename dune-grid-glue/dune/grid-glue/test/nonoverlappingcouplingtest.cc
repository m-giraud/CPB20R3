// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <array>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid-glue/extractors/codim1extractor.hh>

#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>
#include <dune/grid-glue/test/communicationtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

template<typename GridView>
typename Dune::GridGlue::Codim1Extractor<GridView>::Predicate
makeVerticalFacePredicate(double sliceCoord)
{
  using Element = typename GridView::Traits::template Codim<0>::Entity;
  auto predicate = [sliceCoord](const Element& element, unsigned int face) -> bool {
    const int dim = GridView::dimension;
    const auto& refElement = Dune::ReferenceElements<double, dim>::general(element.type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(element.geometry().corner(refElement.subEntity(face,1,i,dim))[0] - sliceCoord) > 1e-6 )
        return false;

    return true;
  };
  return predicate;
}

/** \brief trafo used for yaspgrids */
template<int dim, typename ctype>
class ShiftTrafo
  : public AnalyticalCoordFunction< ctype, dim, dim, ShiftTrafo<dim,ctype> >
{
  double shift;
public:
  ShiftTrafo(double x) : shift(x) {};

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dim> &y ) const
  {
    y = x;
    y[0] += shift;
  }
};

template <int dim>
void testMatchingCubeGrids()
{

  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  using GridType = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;

  std::array<int, dim> elements;
  elements.fill(8);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType cubeGrid0(lower, upper, elements);

  lower[0] += 1;
  upper[0] += 1;

  GridType cubeGrid1(lower, upper, elements);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType::LevelGridView DomGridView;
  typedef typename GridType::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeVerticalFacePredicate<DomGridView>(1);
  const typename TarExtractor::Predicate tardesc = makeVerticalFacePredicate<TarGridView>(1);

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelGridView(0), tardesc);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  auto contactMerger = std::make_shared<ContactMergeImpl>(0.01);
  GlueType contactGlue(domEx, tarEx, contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
}


template <int dim>
void testNonMatchingCubeGrids()
{

  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  using GridType = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;

  std::array<int, dim> elements;
  elements.fill(8);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType cubeGrid0(lower, upper, elements);

  elements.fill(10);
  lower[0] += 1;
  upper[0] += 1;

  GridType cubeGrid1(lower, upper, elements);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType::LevelGridView DomGridView;
  typedef typename GridType::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeVerticalFacePredicate<DomGridView>(1);
  const typename TarExtractor::Predicate tardesc = makeVerticalFacePredicate<TarGridView>(1);

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0.levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1.levelGridView(0), tardesc);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  auto contactMerger = std::make_shared<ContactMergeImpl>(0.01);
  GlueType contactGlue(domEx, tarEx, contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
}


template<int dim, bool par>
class MeshGenerator
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  using GridType = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> >;

  std::shared_ptr<GridType> generate()
  {
    std::array<int, dim> elements;
    elements.fill(8);
    FieldVector<double,dim> lower(0);
    FieldVector<double,dim> upper(1);

    if (tar)
    {
      elements.fill(10);
      lower[0] += 1;
      upper[0] += 1;
    }

    return std::make_shared<GridType>(lower, upper, elements);
  }
};


template<int dim>
class MeshGenerator<dim, true>
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  typedef YaspGrid<dim> HostGridType;
  typedef GeometryGrid<HostGridType, ShiftTrafo<dim,double> > GridType;

  std::shared_ptr<GridType> generate()
  {
    std::array<int,dim> elements;
    std::fill(elements.begin(), elements.end(), 8);
    std::bitset<dim> periodic(0);
    FieldVector<double,dim> size(1);
    int overlap = 1;
    double shift = 0.0;

    if (tar)
    {
      std::fill(elements.begin(), elements.end(), 10);
      shift = 1.0;
    }

    HostGridType * hostgridp = new HostGridType(
      size, elements, periodic, overlap
#if HAVE_MPI
      , MPI_COMM_WORLD
#endif
      );
    ShiftTrafo<dim,double> * trafop = new ShiftTrafo<dim,double>(shift);
    return std::make_shared<GridType>(*hostgridp, *trafop);
  }
};


template <int dim, class DomGen, class TarGen>
void testParallelCubeGrids()
{
  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  typedef typename DomGen::GridType GridType0;
  typedef typename TarGen::GridType GridType1;

  DomGen domGen(0);
  TarGen tarGen(1);

  double slice = 1.0;

  std::shared_ptr<GridType0> cubeGrid0 = domGen.generate();
  std::shared_ptr<GridType1> cubeGrid1 = tarGen.generate();

  // ////////////////////////////////////////
  //   Set up Traits
  // ////////////////////////////////////////

  typedef typename GridType0::LevelGridView DomGridView;
  typedef typename GridType1::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  const typename DomExtractor::Predicate domdesc = makeVerticalFacePredicate<DomGridView>(slice);
  const typename TarExtractor::Predicate tardesc = makeVerticalFacePredicate<TarGridView>(slice);

  auto domEx = std::make_shared<DomExtractor>(cubeGrid0->levelGridView(0), domdesc);
  auto tarEx = std::make_shared<TarExtractor>(cubeGrid1->levelGridView(0), tardesc);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  auto contactMerger = std::make_shared<ContactMergeImpl>(0.01);
  GlueType contactGlue(domEx, tarEx, contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
}

#if HAVE_MPI
void eh( MPI_Comm *comm, int *err, ... )
{
  int len = 1024;
  char error_txt[len];

  MPI_Error_string(*err, error_txt, &len);
  assert(len <= 1024);
  DUNE_THROW(Dune::Exception, "MPI ERROR -- " << error_txt);
}
#endif // HAVE_MPI

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::dinfo.attach(std::cout);

#if HAVE_MPI
  MPI_Errhandler errhandler;
  MPI_Comm_create_errhandler(eh, &errhandler);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
#endif

  // 2d Tests
  typedef MeshGenerator<2,false>  Seq;
  typedef MeshGenerator<2,true>   Par;

  // Test two unit squares
  std::cout << "==== 2D hybrid =============================================\n";
  testMatchingCubeGrids<2>();
  std::cout << "============================================================\n";
  testNonMatchingCubeGrids<2>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Seq,Seq>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Par,Seq>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Seq,Par>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Par,Par>();
  std::cout << "============================================================\n";

  // 3d Tests
#if ! HAVE_MPI
  typedef MeshGenerator<3,false>  Seq3d;
  typedef MeshGenerator<3,true>   Par3d;

  // Test two unit cubes
  std::cout << "==== 3D hybrid =============================================\n";
  testMatchingCubeGrids<3>();
  std::cout << "============================================================\n";
  testNonMatchingCubeGrids<3>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Seq3d,Seq3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Par3d,Seq3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Seq3d,Par3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Par3d,Par3d>();
  std::cout << "============================================================\n";
#endif // HAVE_MPI

  return 0;
}
