#include <fstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {
namespace GridGlue {

namespace AreaWriterImplementation {

template<int dimgrid>
struct FacetLayout
{
  bool contains(Dune::GeometryType gt) const
    {
      return gt.dim() == dimgrid - 1;
    }
};

template<typename GridView>
void write_facet_geometry(const GridView& gv, std::ostream& out)
{
  using Coordinate = Dune::FieldVector<double, 3>;

  std::vector<Coordinate> corners;
  for (const auto& facet : facets(gv)) {
    const auto geometry = facet.geometry();
    for (int i = 0; i < geometry.corners(); ++i) {
      /* VTK always needs 3-dim coordinates... */
      const auto c0 = geometry.corner(i);
      Coordinate c1;
      for (int d = 0; d < GridView::dimensionworld; ++d)
        c1[d] = c0[d];
      for (int d = GridView::dimensionworld; d < Coordinate::dimension; ++d)
        c1[d] = double(0);
      corners.push_back(c1);
    }
  }

  {
    out << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << corners.size() << " double\n";
    for (const auto& c : corners)
      out << c << "\n";
  }
  {
    out << "CELLS " << gv.size(1) << " " << (gv.size(1) + corners.size()) << "\n";
    std::size_t c = 0;
    for (const auto& facet : facets(gv)) {
      const auto geometry = facet.geometry();
      out << geometry.corners();
      for (int i = 0; i < geometry.corners(); ++i, ++c)
        out << " " << c;
      out << "\n";
    }
  }
  {
    out << "CELL_TYPES " << gv.size(1) << "\n";
    for (const auto& facet : facets(gv)) {
      const auto type = facet.type();
      if (type.isVertex())
        out << "1\n";
      else if (type.isLine())
        out << "2\n";
      else if (type.isTriangle())
        out << "5\n";
      else if (type.isQuadrilateral())
        out << "9\n";
      else if (type.isTetrahedron())
        out << "10\n";
      else
        DUNE_THROW(Dune::Exception, "Unhandled geometry type");
    }
  }
}

} /* namespace AreaWriterImplementation */

template<int side, typename Glue>
void write_glue_area_vtk(const Glue& glue, std::ostream& out)
{
  using GridView = typename std::decay< decltype(glue.template gridView<side>()) >::type;
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView, AreaWriterImplementation::FacetLayout>;
  using ctype = typename GridView::ctype;

  const GridView gv = glue.template gridView<side>();
  Mapper mapper(gv);
  std::vector<ctype> coveredArea(mapper.size(), ctype(0));
  std::vector<ctype> totalArea(mapper.size(), ctype(1));

  for (const auto& in : intersections(glue, Reverse<side == 1>())) {
    const auto element = in.inside();
    const auto index = mapper.subIndex(element, in.indexInInside(), 1);
    coveredArea[index] += in.geometryInInside().volume();

    const auto& refElement = Dune::ReferenceElements<ctype, GridView::dimension>::general(element.type());
    const auto& subGeometry = refElement.template geometry<1>(in.indexInInside());
    totalArea[index] = subGeometry.volume();
  }

  for (std::size_t i = 0; i < coveredArea.size(); ++i)
    coveredArea[i] /= totalArea[i];

  out << "# vtk DataFile Version 2.0\n"
      << "Filename: Glue Area\n"
      << "ASCII\n";

  AreaWriterImplementation::write_facet_geometry(gv, out);

  out << "CELL_DATA " << coveredArea.size() << "\n"
      << "SCALARS CoveredArea double 1\n"
      << "LOOKUP_TABLE default\n";
  for (const auto& value : coveredArea)
    out << value << "\n";
}

template<int side, typename Glue>
void write_glue_area_vtk(const Glue& glue, const std::string& filename)
{
  std::ofstream out(filename.c_str());
  write_glue_area_vtk<side>(glue, out);
}

template<typename Glue>
void write_glue_areas_vtk(const Glue& glue, const std::string& base)
{
  {
    std::string filename = base;
    filename += "-inside.vtk";
    write_glue_area_vtk<0>(glue, filename);
  }
  {
    std::string filename = base;
    filename += "-outside.vtk";
    write_glue_area_vtk<1>(glue, filename);
  }
}

} /* namespace GridGlue */
} /* namespace Dune */
