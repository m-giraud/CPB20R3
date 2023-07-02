/**
 * @file
 * @brief Write all remote intersections to a AmiraMesh file
 */

#ifndef DUNE_GRIDGLUE_ADAPTER_GRIDGLUEAMIRAWRITER_HH
#define DUNE_GRIDGLUE_ADAPTER_GRIDGLUEAMIRAWRITER_HH

#include <fstream>
#include <sstream>
#include <type_traits>

namespace Dune {
namespace GridGlue {

/** \brief Write remote intersections to a AmiraMesh file for debugging purposes
 */
class GridGlueAmiraWriter
{

  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
  */
  template <class Glue, int side>
  static void writeIntersections(const Glue& glue, const std::string& filename)
  {
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fgrid;

    fgrid.open(filename.c_str());

    using GridView = typename Glue::template GridView<side>;
    const int dim = GridView::dimension;
    const int domdimw = GridView::dimensionworld;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
        coordinatePadding += " 0";

    int overlaps = glue.size();

    if (dim==3) {

      fgrid << "# HyperSurface 0.1 ASCII \n" << std::endl;
      fgrid<<"\n";
      fgrid<<"Parameters {\n";
      fgrid<<"    Materials {\n";
      fgrid<<"        outside {\n";
      fgrid<<"            Id 0\n";
      fgrid<<"        }\n";
      fgrid<<"        inside {\n";
      fgrid<<"            Id 1\n";
      fgrid<<"        }\n";
      fgrid<<"    }\n";
      fgrid<<"\n";
      fgrid<<"}\n";

      // ////////////////////////////////////////////
      //   Write vertices
      // ////////////////////////////////////////////

      //use dim and not dim+1
      fgrid<<"\nVertices "<< overlaps*(dim)<<"\n";
      auto isEnd = glue.template iend<side>();
      for (auto isIt = glue.template ibegin<side>(); isIt != isEnd; ++isIt)
      {
        const auto& geometry = isIt->geometry();
        for (int i = 0; i < geometry.corners(); ++i)
          fgrid << geometry.corner(i) << coordinatePadding << std::endl;
      }

      // ////////////////////////////////////////////
      //   Write triangles
      // ////////////////////////////////////////////

      fgrid<<"NBranchingPoints 0\n";
      fgrid<<"NVerticesOnCurves 0\n";
      fgrid<<"BoundaryCurves 0\n";
      fgrid<<"Patches 1\n";
      fgrid<<"{\n";
      fgrid<<"InnerRegion inside\n";
      fgrid<<"OuterRegion outside\n";
      fgrid<<"BoundaryID 0\n";
      fgrid<<"BranchingPoints 0";
      fgrid<<"\n";

      fgrid<<"Triangles "<<overlaps<<std::endl;

      for (int i=0;i<overlaps; i++)
        fgrid<<i*dim+1<<" "<<i*dim+2<<" "<<i*dim+3<<std::endl;
      fgrid<<"}\n";

    } else if (dim==2) {

      fgrid << "# AmiraMesh 3D ASCII 2.0 \n";
      fgrid<<"\n";
      fgrid<<"define Lines "<<3*overlaps<<"\n";
      fgrid<<"nVertices "<<2*overlaps<<"\n";
      fgrid<<"\n";
      fgrid<<"Parameters {\n";
      fgrid<<"    ContentType \"HxLineSet\" \n";
      fgrid<<"}\n";
      fgrid<<"\n";
      fgrid<<"Lines { int LineIdx } @1\n";
      fgrid<<"Vertices { float[3] Coordinates } @2\n";
      fgrid<<"\n";
      fgrid<<"# Data section follows\n";
      fgrid<<"@1 \n";
      for (int i=0; i<overlaps;i++)
        fgrid<<2*i<<"\n"<<2*i+1<<"\n"<<-1<<"\n";
      fgrid<<"\n";
      fgrid<<"@2 \n";

      auto isEnd = glue.template iend<side>();
      for (auto isIt = glue.template ibegin<side>(); isIt != isEnd; ++isIt) {
        const auto& geometry = isIt->geometry();
        for (int i = 0; i <2; ++i)
          fgrid << geometry.corner(i) <<" "<<0<<"\n";
      }
    }

    fgrid.close();
  }

public:
  template<typename Glue>
  static void write(const Glue& glue, const std::string& path, int appendix=1)
  {
    std::ostringstream name0;
    name0 << path;
    name0 << "/domain.surf" << std::setw(3) << std::setfill('0') << appendix;

    // Write extracted grid and remote intersection on the grid1-side
    writeIntersections<Glue,0>(glue,name0.str());

    std::ostringstream name1;
    name1 << path;
    name1 << "/target.surf" << std::setw(3) << std::setfill('0') << appendix;

    writeIntersections<Glue,1>(glue, name1.str());
  }

};

} // namespace GridGlue
} // namespace Dune

#endif // DUNE_GRIDGLUE_ADAPTER_GRIDGLUEAMIRAWRITER_HH
