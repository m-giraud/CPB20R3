#include <fstream>

namespace Dune {
namespace GridGlue {

namespace ProjectionWriterImplementation {

template<unsigned side, typename Coordinate, typename Corners>
void write_points(const Projection<Coordinate>& projection, const Corners& corners, std::ostream& out)
{
  using namespace ProjectionImplementation;
  using std::get;
  const unsigned other_side = 1 - side;

  for (const auto& c : get<side>(corners))
    out << c << "\n";

  for (const auto& i : get<side>(projection.images())) {
    const auto global = interpolate(i, get<other_side>(corners));
    out << global << "\n";
  }
}

template<unsigned side, typename Coordinate, typename Normals>
void write_normals(const Projection<Coordinate>& projection, const Normals& normals, std::ostream& out)
{
  using namespace ProjectionImplementation;
  using std::get;
  const unsigned other_side = 1 - side;

  for (const auto& n : get<side>(normals))
    out << n << "\n";

  for (const auto& x : get<side>(projection.images())) {
    const auto n = interpolate_unit_normals(x, get<other_side>(normals));
    out << n << "\n";
  }
}

template<typename Coordinate, typename Corners>
void write_edge_intersection_points(const Projection<Coordinate>& projection, const Corners& corners, std::ostream& out)
{
  using namespace ProjectionImplementation;
  using std::get;

  for (std::size_t i = 0; i < projection.numberOfEdgeIntersections(); ++i) {
    const auto& local = projection.edgeIntersections()[i].local;
    out << interpolate(local[0], get<0>(corners)) << "\n"
        << interpolate(local[1], get<1>(corners)) << "\n";
  }
}

template<typename Coordinate, typename Normals>
void write_edge_intersection_normals(const Projection<Coordinate>& projection, const Normals& normals, std::ostream& out)
{
  using namespace ProjectionImplementation;
  using std::get;

  for (std::size_t i = 0; i < projection.numberOfEdgeIntersections(); ++i) {
    const auto& local = projection.edgeIntersections()[i].local;
    const auto n0 = interpolate_unit_normals(local[0], get<0>(normals));
    const auto n1 = interpolate_unit_normals(local[1], get<1>(normals));

    out << n0 << "\n"
        << n1 << "\n";
  }
}

template<unsigned side, typename Coordinate>
void write_success(const Projection<Coordinate>& projection, std::ostream& out)
{
  using std::get;

  out << side << "\n";

  const auto& success = get<side>(projection.success());
  for (std::size_t i = 0; i < success.size(); ++i)
    out << (success[i] ? "1\n" : "0\n");
}

} /* namespace ProjectionWriterImplementation */

template<typename Coordinate, typename Corners, typename Normals>
void write(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals,
           std::ostream& out)
{
  using namespace ProjectionWriterImplementation;

  const auto numberOfEdgeIntersections = projection.numberOfEdgeIntersections();
  const auto nPoints = 12 + 2 * numberOfEdgeIntersections;

  out << "# vtk DataFile Version2.0\n"
      << "Filename: projection\n"
      << "ASCII\n"
      << "DATASET UNSTRUCTURED_GRID\n"
      << "POINTS " << nPoints << " double\n";
  write_points<0>(projection, corners, out);
  write_points<1>(projection, corners, out);
  write_edge_intersection_points(projection, corners, out);
  out << "CELLS " << (8 + numberOfEdgeIntersections) << " " << (26 + 3 * numberOfEdgeIntersections) << "\n";
  out << "3 0 1 2\n" "2 0 3\n" "2 1 4\n" "2 2 5\n"
      << "3 6 7 8\n" "2 6 9\n" "2 7 10\n" "2 8 11\n";
  for (std::size_t i = 0; i < numberOfEdgeIntersections; ++i)
    out << "2 " << (12 + 2*i) << " " << (12 + 2*i + 1) << "\n";
  out << "CELL_TYPES " << (8 + numberOfEdgeIntersections) << "\n" "5\n3\n3\n3\n" "5\n3\n3\n3\n";
  for (std::size_t i = 0; i < numberOfEdgeIntersections; ++i)
    out << "3\n";
  out << "CELL_DATA " << (8 + numberOfEdgeIntersections) << "\n";
  out << "SCALARS success int 1\n"
      << "LOOKUP_TABLE success\n";
  write_success<0>(projection, out);
  write_success<1>(projection, out);
  for (std::size_t i = 0; i < numberOfEdgeIntersections; ++i)
    out << "2\n";
  out << "LOOKUP_TABLE success 2\n"
      << "1.0 0.0 0.0 1.0\n"
      << "0.0 1.0 0.0 1.0\n";
  out << "POINT_DATA " << nPoints << "\n"
      << "NORMALS normals double\n";
  write_normals<0>(projection, normals, out);
  write_normals<1>(projection, normals, out);
  write_edge_intersection_normals(projection, normals, out);
}

template<typename Coordinate, typename Corners, typename Normals>
void write(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals,
           const std::string& filename)
{
  std::ofstream out(filename.c_str());
  write(projection, corners, normals, out);
}

template<typename Coordinate, typename Corners, typename Normals>
void print(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals)
{
  using namespace ProjectionWriterImplementation;

  std::cout << "Side 0 corners and images:\n";
  write_points<0>(projection, corners, std::cout);
  std::cout << "Side 0 success:\n";
  write_success<0>(projection, std::cout);
  std::cout << "Side 1 corners and images:\n";
  write_points<1>(projection, corners, std::cout);
  std::cout << "Side 1 success:\n";
  write_success<1>(projection, std::cout);
  std::cout << "Side 0 normals and projected normals:\n";
  write_normals<0>(projection, normals, std::cout);
  std::cout << "Side 1 normals and projected normals:\n";
  write_normals<1>(projection, normals, std::cout);
  std::cout << projection.numberOfEdgeIntersections() << " edge intersections:\n";
  write_edge_intersection_points(projection, corners, std::cout);
}

} /* namespace GridGlue */
} /* namespace Dune */
