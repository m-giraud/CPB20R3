// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    codim1extractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Oliver Sander, Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: class for grid extractors extracting surface grids
 *
 */
/**
 * @file
 * @brief Grid extractor class for codim 1 subgrids
 */

#ifndef DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH
#define DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH

#include "extractor.hh"

#include <array>
#include <deque>
#include <functional>

#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>
#include <dune/grid-glue/common/crossproduct.hh>

namespace Dune {

  namespace GridGlue {

/**
 * extractor for codim-1 entities (facets)
 */
template<typename GV>
class Codim1Extractor : public Extractor<GV,1>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  using Extractor<GV,1>::dimworld;
  using Extractor<GV,1>::dim;
  using Extractor<GV,1>::codim;
  using Extractor<GV,1>::cube_corners;
  typedef typename Extractor<GV,1>::IndexType IndexType;

  /// @brief compile time number of corners of surface simplices
  static constexpr int simplex_corners = dim;

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;

  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef std::function<bool(const Element&, unsigned int subentity)> Predicate;

  // import typedefs from base class
  typedef typename Extractor<GV,1>::SubEntityInfo SubEntityInfo;
  typedef typename Extractor<GV,1>::ElementInfo ElementInfo;
  typedef typename Extractor<GV,1>::VertexInfo VertexInfo;
  typedef typename Extractor<GV,1>::CoordinateInfo CoordinateInfo;
  typedef typename Extractor<GV,1>::VertexInfoMap VertexInfoMap;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   * @param predicate a predicate to mark entities for extraction (unary functor returning bool)
   */
  Codim1Extractor(const GV& gv, const Predicate& predicate)
    :  Extractor<GV,1>(gv)
  {
    std::cout << "This is Codim1Extractor on a <" << dim
              << "," << dimworld << "> grid!"
              << std::endl;
    update(predicate);
  }

private:

  /**
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a predicate class.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * @param predicate a predicate that "selects" the faces to add to the surface
   */
  void update(const Predicate& predicate);

};


template<typename GV>
void Codim1Extractor<GV>::update(const Predicate& predicate)
{
  // free everything there is in this object
  this->clear();

  // In this first pass iterate over all entities of codim 0.
  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.
  {
    // several counter for consecutive indexing are needed
    int simplex_index = 0;
    int vertex_index = 0;
    IndexType eindex = 0;     // supress warning

    // needed later for insertion into a std::set which only
    // works with const references

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<SubEntityInfo> temp_faces;

    // iterate over interior codim 0 elements on the grid
    for (const auto& elmt : elements(this->gv_, Partitions::interior))
    {
      Dune::GeometryType gt = elmt.type();

      // if some face is part of the surface add it!
      if (elmt.hasBoundaryIntersections())
      {
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
        eindex = this->cellMapper_.index(elmt);
        this->elmtInfo_.emplace(eindex, ElementInfo(simplex_index, elmt, 0));

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (const auto& in : intersections(this->gv_, elmt))
        {
          // Stop only at selected boundary faces
          if (!in.boundary() or !predicate(elmt, in.indexInInside()))
            continue;

          const auto& refElement = Dune::ReferenceElements<ctype, dim>::general(gt);
          // get the corner count of this face
          const int face_corners = refElement.size(in.indexInInside(), 1, dim);

          // now we only have to care about the 3D case, i.e. a triangle face can be
          // inserted directly whereas a quadrilateral face has to be divided into two triangles
          switch (face_corners)
          {
          case 2 :
          case 3:
          {
            // we have a simplex here

            // register the additional face(s)
            this->elmtInfo_.at(eindex).faces++;

            // add a new face to the temporary collection
            temp_faces.emplace_back(eindex, in.indexInInside(),
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
                                    Dune::GeometryTypes::simplex(dim-codim)
#else
                                    Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)
#endif
                                    );

            std::vector<FieldVector<ctype,dimworld> > cornerCoords(face_corners);

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < face_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = refElement.subEntity(in.indexInInside(), 1, i, dim);

              // get the vertex pointer and the index from the index set
              const Vertex vertex = elmt.template subEntity<dim>(vertex_number);
              cornerCoords[i] = vertex.geometry().corner(0);

              IndexType vindex = this->gv_.indexSet().template index<dim>(vertex);

              // remember the vertex' number in parent element's vertices
              temp_faces.back().corners[i].num = vertex_number;

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->vtxInfo_.find(vindex);
              if (vimit == this->vtxInfo_.end())
              {
                // insert into the map
                this->vtxInfo_.emplace(vindex, VertexInfo(vertex_index, vertex));
                // remember the vertex as a corner of the current face in temp_faces
                temp_faces.back().corners[i].idx = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only insert the index into the simplices array
                temp_faces.back().corners[i].idx = vimit->second.idx;
              }
            }

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want them to be oriented such that all boundary edges
            // point in the counterclockwise direction.  Therefore, we check the orientation of the
            // new face and possibly switch the two vertices.
            FieldVector<ctype,dimworld> realNormal = in.centerUnitOuterNormal();

            // Compute segment normal
            FieldVector<ctype,dimworld> reconstructedNormal;
            if (dim==2)  // boundary face is a line segment
            {
              reconstructedNormal[0] =  cornerCoords[1][1] - cornerCoords[0][1];
              reconstructedNormal[1] =  cornerCoords[0][0] - cornerCoords[1][0];
            } else {    // boundary face is a triangle
              FieldVector<ctype,dimworld> segment1 = cornerCoords[1] - cornerCoords[0];
              FieldVector<ctype,dimworld> segment2 = cornerCoords[2] - cornerCoords[0];
              reconstructedNormal = crossProduct(segment1, segment2);
            }
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);

            // now increase the current face index
            simplex_index++;
            break;
          }
          case 4 :
          {
            assert(dim == 3 && cube_corners == 4);
            // we have a quadrilateral here
            std::array<unsigned int, 4> vertex_indices;
            std::array<unsigned int, 4> vertex_numbers;

            // register the additional face(s) (2 simplices)
            this->elmtInfo_.at(eindex).faces += 2;

            std::array<FieldVector<ctype,dimworld>, 4> cornerCoords;

            // get the vertex pointers for the quadrilateral's corner vertices
            // and try for each of them whether it is already inserted or not
            for (int i = 0; i < cube_corners; ++i)
            {
              // get the number of the vertex in the parent element
              vertex_numbers[i] = refElement.subEntity(in.indexInInside(), 1, i, dim);

              // get the vertex pointer and the index from the index set
              const Vertex vertex = elmt.template subEntity<dim>(vertex_numbers[i]);
              cornerCoords[i] = vertex.geometry().corner(0);

              IndexType vindex = this->gv_.indexSet().template index<dim>(vertex);

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->vtxInfo_.find(vindex);
              if (vimit == this->vtxInfo_.end())
              {
                // insert into the map
                this->vtxInfo_.emplace(vindex, VertexInfo(vertex_index, vertex));
                // remember this vertex' index
                vertex_indices[i] = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only remember the vertex' index
                vertex_indices[i] = vimit->second.idx;
              }
            }

            // now introduce the two triangles subdividing the quadrilateral
            // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
            // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1

            // add a new face to the temporary collection for the first tri
            temp_faces.emplace_back(eindex, in.indexInInside(),
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
                                    Dune::GeometryTypes::simplex(dim-codim)
#else
                                    Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)
#endif
                                    );
            temp_faces.back().corners[0].idx = vertex_indices[0];
            temp_faces.back().corners[1].idx = vertex_indices[1];
            temp_faces.back().corners[2].idx = vertex_indices[2];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[0];
            temp_faces.back().corners[1].num = vertex_numbers[1];
            temp_faces.back().corners[2].num = vertex_numbers[2];

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want the triangle vertices on counterclockwise order,
            // when viewed from the outside of the grid. Therefore, we check the orientation of the
            // new face and possibly switch two vertices.
            FieldVector<ctype,dimworld> realNormal = in.centerUnitOuterNormal();

            // Compute segment normal
            FieldVector<ctype,dimworld> reconstructedNormal = crossProduct(cornerCoords[1] - cornerCoords[0],
                                                                           cornerCoords[2] - cornerCoords[0]);
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);


            // add a new face to the temporary collection for the second tri
            temp_faces.emplace_back(eindex, in.indexInInside(),
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 6)
                                    Dune::GeometryTypes::simplex(dim-codim)
#else
                                    Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)
#endif
                                    );
            temp_faces.back().corners[0].idx = vertex_indices[3];
            temp_faces.back().corners[1].idx = vertex_indices[2];
            temp_faces.back().corners[2].idx = vertex_indices[1];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[3];
            temp_faces.back().corners[1].num = vertex_numbers[2];
            temp_faces.back().corners[2].num = vertex_numbers[1];

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want the triangle vertices on counterclockwise order,
            // when viewed from the outside of the grid. Therefore, we check the orientation of the
            // new face and possibly switch two vertices.
            // Compute segment normal
            reconstructedNormal = crossProduct(cornerCoords[2] - cornerCoords[3],
                                               cornerCoords[1] - cornerCoords[3]);
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);

            simplex_index+=2;
            break;
          }
          default :
            DUNE_THROW(Dune::NotImplemented, "the extractor does only work for triangle and quadrilateral faces (" << face_corners << " corners)");
            break;
          }
        }         // end loop over found surface parts
      }
    }     // end loop over elements

    std::cout << "added " << simplex_index << " subfaces\n";

    // allocate the array for the face specific information...
    this->subEntities_.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->subEntities_.begin());
  }


  // now first write the array with the coordinates...
  this->coords_.resize(this->vtxInfo_.size());
  for (const auto& vinfo : this->vtxInfo_)
  {
    // get a pointer to the associated info object
    CoordinateInfo* current = &this->coords_[vinfo.second.idx];
    // store this coordinates index // NEEDED?
    current->index = vinfo.second.idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = vinfo.first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    const auto vtx = this->grid().entity(vinfo.second.p);
    current->coord = vtx.geometry().corner(0);
  }

}

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH
