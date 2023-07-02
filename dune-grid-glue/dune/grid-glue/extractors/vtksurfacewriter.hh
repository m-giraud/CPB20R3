// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    VtkSurfaceWriter.hh
 *  Version:     1.0
 *  Created on:  Jan 16, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: helper class for graphical output of grids in generic representation
 *
 */
/**
 * @file
 * @brief helper class for graphical output of grids in generic representation
 */

#ifndef DUNE_GRIDGLUE_EXTRACTORS_VTKSURFACEWRITER_HH
#define DUNE_GRIDGLUE_EXTRACTORS_VTKSURFACEWRITER_HH

#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>

#include "../adapter/gridgluevtkwriter.hh"

namespace Dune {

  namespace GridGlue {

class VtkSurfaceWriter
{
public:


  VtkSurfaceWriter(const char* filename) : filename_(filename)
  {}

  ~VtkSurfaceWriter()
  {}

  void setFilename(const char* name)
  {
    if (std::strlen(name) > 0)
      this->filename_ = name;
  }


  template<typename K>
  void writeSurface(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, int dim)
  {
    std::ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->filename_);
    fos.open(buffer);
    fos << std::setprecision(8) << std::setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << std::endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    fos.close();
  }


  template<typename K, typename T>
  void writeSurfaceElementData(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, const std::vector<T>& data, const char* dataname, int dim)
  {
    std::ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->filename_);
    fos.open(buffer);
    fos << std::setprecision(8) << std::setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << std::endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    this->writeCellData(data, dataname, dim, fos);
    fos.close();
  }


  template<typename K, typename T>
  void writeSurfaceVertexData(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, const std::vector<T>& data, const char* dataname, int dim)
  {
    std::ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->filename_);
    fos.open(buffer);
    fos << std::setprecision(8) << std::setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << std::endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    this->writePointData(data, dataname, dim, fos);
    fos.close();
  }

protected:

  template<typename K>
  void writePoints(const std::vector<K>& coords, int dim, std::ofstream& fos)
  {
    fos << "DATASET POLYDATA\nPOINTS " << coords.size() << " " << TypeNames[Nametraits<K>::nameidx] << std::endl;
    for (unsigned int i = 0; i < coords.size(); ++i)
    {
      fos << coords[i][0];
      if (dim == 2)
        fos << " " << coords[i][1] << " 0 \n" << coords[i][0] << " " << coords[i][1] << " 0.01" << std::endl;
      else       // dim == 3
        fos << " " << coords[i][1] << " "  << coords[i][2] << std::endl;
    }
  }

  void writePolygons(const std::vector<unsigned int>& indices, const int* corners, int ncorners, int dim, std::ofstream& fos)
  {
    if (dim == 2)
    {
      fos << "POLYGONS " << indices.size()/2 << " " << 5*(indices.size() / 2) << std::endl;
      for (unsigned int i = 0; i < indices.size(); i += 2)
        fos << "4 " << 2*indices[i] << " " << 2*indices[i+1] << " " << 2*indices[i+1]+1 << " "<< 2*indices[i]+1 << std::endl;

      // arbitrary shapes - ignored here!
      //                      int sum = ncorners;
      //                      for (int i = 0; i < ncorners; ++i)
      //                              sum += (corners[i] > 2 ? corners[i] : 3);
      //
      //                      fos << "POLYGONS " << ncorners << " " << sum << std::endl;
      //                      int index = 0;
      //                      for (int i = 0; i < ncorners; ++i)
      //                      {
      //                              // write the first index twice if it is an egde
      //                              // => triangle instead of edge - paraview can display it then
      //                              if (corners[i] > 2)
      //                                      fos << corners[i];
      //                              else
      //                                      fos << "3 " << indices[index];
      //
      //                              for (int j = 0; j < corners[i]; ++j)
      //                                      fos << " " << indices[index++];
      //                              fos << std::endl;
      //                       }
    }
    else
    {
      int sum = ncorners;
      for (int i = 0; i < ncorners; ++i)
        sum += corners[i];
      fos << "POLYGONS " << ncorners << " " << sum << std::endl;
      int index = 0;
      for (int i = 0; i < ncorners; ++i)
      {
        fos << corners[i];
        for (int j = 0; j < corners[i]; ++j)
          fos << " " << indices[index++];
        fos << std::endl;
      }
    }
  }

  template<typename T>
  void writeCellData(const std::vector<T>& data, const char* dataname, int dim, std::ofstream& fos)
  {
    fos << "CELL_DATA " << data.size()*(dim == 2 ? 2 : 1) << std::endl;
    fos << "SCALARS " << dataname << " " << TypeNames[Nametraits<T>::nameidx] << " 1" << std::endl;
    fos << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      fos << data[i] << std::endl;
      if (dim == 2)
        fos << data[i] << std::endl;
    }
  }

  template<typename T>
  void writePointData(const std::vector<T>& data, const char* dataname, int dim, std::ofstream& fos)
  {
    fos << "POINT_DATA " << data.size()*(dim == 2 ? 2 : 1) << std::endl;
    fos << "SCALARS " << dataname << " " << TypeNames[Nametraits<T>::nameidx] << " 1" << std::endl;
    fos << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      fos << data[i] << std::endl;
      if (dim == 2)
        fos << data[i] << std::endl;
    }
  }


private:
  const char*  filename_;
};

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_EXTRACTORS_VTKSURFACEWRITER_HH
