// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueVtkWriter.hh
 *  Version:     1.0
 *  Created on:  Mar 5, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Class thought to make graphical debugging of couplings easier.
 *
 */
/**
 * @file
 * @brief Write all remote intersections to a vtk file for debugging
 */

#ifndef DUNE_GRIDGLUE_ADAPTER_GRIDGLUEVTKWRITER_HH
#define DUNE_GRIDGLUE_ADAPTER_GRIDGLUEVTKWRITER_HH


#include <fstream>
#include <iomanip>
#include <type_traits>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {
namespace GridGlue {

/** \brief Write remote intersections to a vtk file for debugging purposes
 */
class GridGlueVtkWriter
{

  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeExtractedPart(const Glue& glue, const std::string& filename)
  {
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fgrid;

    fgrid.open(filename.c_str());

    using GridView = typename Glue::template GridView<side>;
    using Extractor = typename Glue::template GridPatch<side>;

    typedef typename GridView::ctype ctype;

    const int domdimw = GridView::dimensionworld;
    const int patchDim = Extractor::dim - Extractor::codim;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
      coordinatePadding += " 0";

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << filename << "\nASCII" << std::endl;

    // WRITE POINTS
    // ----------------
    std::vector<typename Extractor::Coords> coords;
    glue.template patch<side>().getCoords(coords);

    fgrid << ((patchDim==3) ? "DATASET UNSTRUCTURED_GRID" : "DATASET POLYDATA") << std::endl;
    fgrid << "POINTS " << coords.size() << " " << Dune::className<ctype>() << std::endl;

    for (size_t i=0; i<coords.size(); i++)
      fgrid << coords[i] << coordinatePadding << std::endl;

    fgrid << std::endl;

    // WRITE POLYGONS
    // ----------------

    std::vector<typename Extractor::VertexVector> faces;
    std::vector<Dune::GeometryType> geometryTypes;
    glue.template patch<side>().getFaces(faces);
    glue.template patch<side>().getGeometryTypes(geometryTypes);

    unsigned int faceCornerCount = 0;
    for (size_t i=0; i<faces.size(); i++)
      faceCornerCount += faces[i].size();

    fgrid << ((patchDim==3) ? "CELLS " : "POLYGONS ")
          << geometryTypes.size() << " " << geometryTypes.size() + faceCornerCount << std::endl;

    for (size_t i=0; i<faces.size(); i++) {

      fgrid << faces[i].size();

      // vtk expects the vertices to by cyclically ordered
      // therefore unfortunately we have to deal with several element types on a case-by-case basis
      if (geometryTypes[i].isSimplex()) {
        for (int j=0; j<patchDim+1; j++)
          fgrid << " " << faces[i][j];

      } else if (geometryTypes[i].isQuadrilateral()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2];

      } else if (geometryTypes[i].isPyramid()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2] << " " << faces[i][4];

      } else if (geometryTypes[i].isPrism()) {
        fgrid << " " << faces[i][0] << " " << faces[i][2] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][5] << " " << faces[i][4];

      } else if (geometryTypes[i].isHexahedron()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2]
              << " " << faces[i][4] << " " << faces[i][5]
              << " " << faces[i][7] << " " << faces[i][6];

      } else {
        DUNE_THROW(Dune::NotImplemented, "Geometry type " << geometryTypes[i] << " not supported yet");
      }

      fgrid << std::endl;
    }

    fgrid << std::endl;

    // 3d VTK files need an extra section specifying the CELL_TYPES aka GeometryTypes
    if (patchDim==3) {

      fgrid << "CELL_TYPES " << geometryTypes.size() << std::endl;

      for (size_t i=0; i<geometryTypes.size(); i++) {
        if (geometryTypes[i].isSimplex())
          fgrid << "10" << std::endl;
        else if (geometryTypes[i].isHexahedron())
          fgrid << "12" << std::endl;
        else if (geometryTypes[i].isPrism())
          fgrid << "13" << std::endl;
        else if (geometryTypes[i].isPyramid())
          fgrid << "14" << std::endl;
        else
          DUNE_THROW(Dune::NotImplemented, "Geometry type " << geometryTypes[i] << " not supported yet");

      }

    }

#if 0
    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (gridSubEntityData.size()-1);

    fgrid << "CELL_DATA " << gridSubEntityData.size() << std::endl;
    fgrid << "SCALARS property_coding " << Dune::className<ctype>() << " 1" << std::endl;
    fgrid << "LOOKUP_TABLE default" << std::endl;

    for (typename GridSubEntityData::const_iterator sEIt = gridSubEntityData.begin();
         sEIt != gridSubEntityData.end();
         ++sEIt, accum += delta)
    {
      // "encode" the parent with one color...
      fgrid << accum << std::endl;
    }
#endif
    fgrid.close();
  }


  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeIntersections(const Glue& glue, const std::string& filename)
  {
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fmerged;

    fmerged.open(filename.c_str());

    using GridView = typename Glue::template GridView<side>;
    using RemoteIntersectionIterator = typename Glue::template IntersectionIterator<side>;

    typedef typename GridView::ctype ctype;

    const int domdimw = GridView::dimensionworld;
    const int intersectionDim = Glue::Intersection::mydim;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
      coordinatePadding += " 0";

    int overlaps = glue.size();

    // WRITE POINTS
    // ----------------
    using Extractor = typename Glue::template GridPatch<0>;
    std::vector<typename Extractor::Coords> coords;
    glue.template patch<side>().getCoords(coords);

    // the merged grid (i.e. the set of remote intersections
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << filename << "\nASCII" << std::endl;
    fmerged << ((intersectionDim==3) ? "DATASET UNSTRUCTURED_GRID" : "DATASET POLYDATA") << std::endl;
    fmerged << "POINTS " << overlaps*(intersectionDim+1) << " " << Dune::className<ctype>() << std::endl;

    for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
         isIt != glue.template iend<side>();
         ++isIt)
    {
      for (int i = 0; i < isIt->geometry().corners(); ++i)
        fmerged << isIt->geometry().corner(i) << coordinatePadding << std::endl;
    }

    // WRITE POLYGONS
    // ----------------

    std::vector<typename Extractor::VertexVector> faces;
    std::vector<Dune::GeometryType> geometryTypes;
    glue.template patch<side>().getFaces(faces);
    glue.template patch<side>().getGeometryTypes(geometryTypes);

    unsigned int faceCornerCount = 0;
    for (size_t i=0; i<faces.size(); i++)
      faceCornerCount += faces[i].size();

    int grid0SimplexCorners = intersectionDim+1;
    fmerged << ((intersectionDim==3) ? "CELLS " : "POLYGONS ")
            << overlaps << " " << (grid0SimplexCorners+1)*overlaps << std::endl;

    for (int i = 0; i < overlaps; ++i) {
      fmerged << grid0SimplexCorners;
      for (int j=0; j<grid0SimplexCorners; j++)
        fmerged << " " << grid0SimplexCorners*i+j;
      fmerged << std::endl;
    }

    // 3d VTK files need an extra section specifying the CELL_TYPES aka GeometryTypes
    if (intersectionDim==3) {

      fmerged << "CELL_TYPES " << overlaps << std::endl;

      for (int i = 0; i < overlaps; i++)
        fmerged << "10" << std::endl;

    }

#if 0
    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (gridSubEntityData.size()-1);

    fmerged << "CELL_DATA " << overlaps << std::endl;
    fmerged << "SCALARS property_coding " << Dune::className<ctype>() << " 1" << std::endl;
    fmerged << "LOOKUP_TABLE default" << std::endl;

    for (typename GridSubEntityData::const_iterator sEIt = gridSubEntityData.begin();
         sEIt != gridSubEntityData.end();
         ++sEIt, accum += delta)
    {
      // ...and mark all of its merged grid parts with the same color
      for (int j = 0; j < sEIt->first.second; ++j)
        fmerged << accum << std::endl;
    }
#endif
    fmerged.close();
  }

public:
  template<typename Glue>
  static void write(const Glue& glue, const std::string& filenameTrunk)
  {

    // Write extracted grid and remote intersection on the grid0-side
    writeExtractedPart<Glue,0>(glue,
                               filenameTrunk + "-grid0.vtk");

    writeIntersections<Glue,0>(glue,
                               filenameTrunk + "-intersections-grid0.vtk");

    // Write extracted grid and remote intersection on the grid1-side
    writeExtractedPart<Glue,1>(glue,
                               filenameTrunk + "-grid1.vtk");

    writeIntersections<Glue,1>(glue,
                               filenameTrunk + "-intersections-grid1.vtk");

  }

};

} /* namespace GridGlue */
} /* namespace Dune */

#endif // DUNE_GRIDGLUE_ADAPTER_GRIDGLUEVTKWRITER_HH
