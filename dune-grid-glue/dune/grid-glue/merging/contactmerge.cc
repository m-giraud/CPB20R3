// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/grid-glue/common/crossproduct.hh>
#include <dune/grid-glue/common/projection.hh>

namespace Dune {
namespace GridGlue {

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeIntersections(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects1,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects2,
                                   unsigned int grid2Index,
                                   std::vector<SimplicialIntersection>& intersections)
{
    using std::get;

    std::vector<std::array<LocalCoords,2> > polytopeCorners;

    // Initialize
    neighborIntersects1.reset();
    neighborIntersects2.reset();

    const int nCorners1 = grid1ElementCorners.size();
    const int nCorners2 = grid2ElementCorners.size();

    if (nCorners1 != dimworld)
      DUNE_THROW(Dune::Exception, "element1 must have " << dimworld << " corners, but has " << nCorners1);
    if (nCorners2 != dimworld)
      DUNE_THROW(Dune::Exception, "element2 must have " << dimworld << " corners, but has " << nCorners2);

    // The grid1 projection directions
    std::vector<WorldCoords> directions1(nCorners1);
    for (size_t i=0; i<directions1.size(); i++)
        directions1[i] = nodalDomainDirections_[this->grid1ElementCorners_[grid1Index][i]];

    // The grid2 projection directions
    std::vector<WorldCoords> directions2(nCorners2);
    for (size_t i=0; i<directions2.size(); i++)
        directions2[i] = nodalTargetDirections_[this->grid2ElementCorners_[grid2Index][i]];

    // The difference between the closest point projection and the normal projection is just the ordering
    // of the involved surfaces. The closest point projection is done along the outer normal field of grid2
    // (due to being a best approximation) and the outer normal projection is using the outer normal field
    // of grid1 instead.
    std::array<decltype(std::cref(grid1ElementCorners)),2> cornersRef ={std::cref(grid1ElementCorners), std::cref(grid2ElementCorners)};
    std::array<decltype(std::cref(directions1)),2> directionsRef ={std::cref(directions1), std::cref(directions2)};
    std::array<Dune::GeometryType,2> elementTypes = {grid1ElementType, grid2ElementType};

    // Determine which is the grid we use for outer normal projection
    const size_t domGrid = (type_==ProjectionType::OUTER_NORMAL) ? 0 : 1;
    const size_t tarGrid = (type_==ProjectionType::OUTER_NORMAL) ? 1 : 0;

    /////////////////////////////////////////////////////
    //  Compute all corners of the intersection polytope
    /////////////////////////////////////////////////////

    const auto corners = std::tie(cornersRef[domGrid].get(),cornersRef[tarGrid].get());
    const auto normals = std::tie(directionsRef[domGrid].get(), directionsRef[tarGrid].get());
    Projection<WorldCoords> p(overlap_, maxNormalProduct_);
    p.project(corners, normals);

    /* projection */
    {
      const auto& success = get<0>(p.success());
      const auto& images = get<0>(p.images());
      for (unsigned i = 0; i < dimworld; ++i) {
        if (success[i]) {
          std::array<LocalCoords, 2> corner;
          corner[domGrid] = localCornerCoords(i, elementTypes[domGrid]);
          for (unsigned j = 0; j < dim; ++j)
            corner[tarGrid][j] = images[i][j];
          polytopeCorners.push_back(corner);
        }
      }
    }

    /* inverse projection */
    {
      const auto& success = get<1>(p.success());
      const auto& preimages = get<1>(p.images());
      for (unsigned i = 0; i < dimworld; ++i) {
        if (success[i]) {
          std::array<LocalCoords, 2> corner;
          for (unsigned j = 0; j < dim; ++j)
            corner[domGrid][j] = preimages[i][j];
          corner[tarGrid] = localCornerCoords(i, elementTypes[tarGrid]);
          polytopeCorners.push_back(corner);
        }
      }
    }

    /* edge intersections */
    {
      for (unsigned i = 0; i < p.numberOfEdgeIntersections(); ++i) {
        std::array<LocalCoords, 2> corner;
        const auto& local = p.edgeIntersections()[i].local;
        for (unsigned j = 0; j < dim; ++j) {
          corner[domGrid][j] = local[0][j];
          corner[tarGrid][j] = local[1][j];
        }
        polytopeCorners.push_back(corner);
      }
    }

    // check which neighbors might also intersect
    std::array<decltype(std::ref(neighborIntersects1)),2> neighborIntersectsRef = {std::ref(neighborIntersects1), std::ref(neighborIntersects2)};
    const auto& refTar = Dune::ReferenceElements<T,dim>::general(elementTypes[tarGrid]);
    for (int i=0; i<refTar.size(1); i++) {

        // if all face corners hit the other element then
        // the neighbor might also intersect

        bool intersects(true);
        for (int k=0; k<refTar.size(i,1,dim); k++)
            intersects &= get<1>(p.success())[refTar.subEntity(i,1,k,dim)];

        if (intersects)
            neighborIntersectsRef[tarGrid].get()[i] = true;
    }

    const auto& refDom = Dune::ReferenceElements<T,dim>::general(elementTypes[domGrid]);
    for (int i=0; i<refDom.size(1); i++) {

        // if all face corners hit the other element then
        // the neighbor might also intersect

        bool intersects(true);
        for (int k=0; k<refDom.size(i,1,dim); k++)
            intersects &= get<0>(p.success())[refDom.subEntity(i,1,k,dim)];

        if (intersects)
            neighborIntersectsRef[domGrid].get()[i] = true;
    }

    // Compute the edge intersections
    for (unsigned i = 0; i < p.numberOfEdgeIntersections(); ++i) {
      const auto& edge = p.edgeIntersections()[i].edge;
      neighborIntersects1[edge[domGrid]] = true;
      neighborIntersects2[edge[tarGrid]] = true;
    }

    // remove possible doubles
    removeDoubles(polytopeCorners);

    // Compute an interior point of the polytope
    int nPolyCorners = polytopeCorners.size();

    // If the polytope is degenerated then there is no intersection
    if (nPolyCorners<dimworld)
        return;

    // If the polytope is a simplex return it
    if (nPolyCorners==dim+1) {

     //   std::cout<<"Add intersection: 1\n";
        typename Base::SimplicialIntersection intersect(grid1Index, grid2Index);

        for (int j=0;j<dim+1; j++) {
            intersect.corners0[0][j]=polytopeCorners[j][0];
            intersect.corners1[0][j]=polytopeCorners[j][1];
        }
        intersections.push_back(intersect);

        return;
    }

    // At this point we must have dimworld>=3

    ///////////////////////////////////////////////////////////////////////////////
    //  Compute a point in the middle of the polytope and order all corners cyclic
    //////////////////////////////////////////////////////////////////////////////

    std::array<LocalCoords,2> center;
    center[0] = 0; center[1] = 0;
    for (int i=0; i<nPolyCorners; i++) {
        center[0].axpy(1.0/nPolyCorners,polytopeCorners[i][0]);
        center[1].axpy(1.0/nPolyCorners,polytopeCorners[i][1]);
    }

    // Order cyclic
    std::vector<int> ordering;
    computeCyclicOrder(polytopeCorners,center[0],ordering);

    //////////////////////////////////////
    // Add intersections
    ////////////////////////////////

    for (size_t i=1; i<polytopeCorners.size()-1; i++) {

        typename Base::SimplicialIntersection intersect(grid1Index, grid2Index);

        for (int j=0;j<dim; j++) {
            intersect.corners0[0][j]=polytopeCorners[ordering[i+j]][0];
            intersect.corners1[0][j]=polytopeCorners[ordering[i+j]][1];
        }

        // last corner is the first for all intersections
        intersect.corners0[0][dim]=polytopeCorners[ordering[0]][0];
        intersect.corners1[0][dim]=polytopeCorners[ordering[0]][1];

        intersections.push_back(intersect);
    }
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeCyclicOrder(const std::vector<std::array<LocalCoords,2> >& polytopeCorners,
                        const LocalCoords& center, std::vector<int>& ordering) const
{
    ordering.resize(polytopeCorners.size());

    for (size_t k=0; k<ordering.size(); k++)
        ordering[k] = k;

    //TODO Do I have to order triangles to get some correct orientation?
    if (polytopeCorners.size()<=3)
        return;

    // compute angles inside the polygon plane w.r.t to this axis
    LocalCoords  edge0 = polytopeCorners[1][0] - polytopeCorners[0][0];

    // Compute a vector that is perpendicular to the edge but lies in the polytope plane
    // So we have a unique ordering
    LocalCoords  edge1 = polytopeCorners[2][0] - polytopeCorners[0][0];
    LocalCoords normal0 = edge1;
    normal0.axpy(-(edge0*edge1),edge0);

    std::vector<T> angles(polytopeCorners.size());

    for (size_t i=0; i<polytopeCorners.size(); i++) {

        LocalCoords  edge = polytopeCorners[i][0] - center;

        T x(edge*edge0);
        T y(edge*normal0);

        angles[i] = std::atan2(y, x);
        if (angles[i]<0)
            angles[i] += 2*M_PI;
    }

    // bubblesort

    for (int i=polytopeCorners.size(); i>1; i--){
        bool swapped = false;

        for (int j=0; j<i-1; j++){

            if (angles[j] > angles[j+1]){
                swapped = true;
                std::swap(angles[j], angles[j+1]);
                std::swap(ordering[j], ordering[j+1]);
            }
        }

        if (!swapped)
            break;
    }
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::setupNodalDirections(const std::vector<WorldCoords>& coords1,
                              const std::vector<unsigned int>& elements1,
                              const std::vector<Dune::GeometryType>& elementTypes1,
                              const std::vector<WorldCoords>& coords2,
                              const std::vector<unsigned int>& elements2,
                              const std::vector<Dune::GeometryType>& elementTypes2)
{
    if (domainDirections_) {

        // Sample the provided analytical contact direction field
        nodalDomainDirections_.resize(coords1.size());
        for (size_t i=0; i<coords1.size(); i++)
            nodalDomainDirections_[i] = domainDirections_(coords1[i]);
    } else
        computeOuterNormalField(coords1,elements1,elementTypes1, nodalDomainDirections_);

    if (targetDirections_) {

        // Sample the provided analytical target direction field
        nodalTargetDirections_.resize(coords2.size());
        for (size_t i=0; i<coords2.size(); i++)
            nodalTargetDirections_[i] = targetDirections_(coords2[i]);
    } else
        computeOuterNormalField(coords2,elements2,elementTypes2, nodalTargetDirections_);
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeOuterNormalField(const std::vector<WorldCoords>& coords,
                                const std::vector<unsigned int>& elements,
                                const std::vector<Dune::GeometryType>& elementTypes,
                                std::vector<WorldCoords>& normals)
{
    normals.assign(coords.size(),WorldCoords(0));


    int offset = 0;

    for (size_t i=0; i<elementTypes.size(); i++) {

        int nCorners = Dune::ReferenceElements<T,dim>::general(elementTypes[i]).size(dim);

        // For segments 1, for triangles or quadrilaterals take the first 2
        std::array<WorldCoords, dim> edges;
        for (int j=1; j<=dim; j++)
            edges[j-1] = coords[elements[offset + j]] - coords[elements[offset]];

        WorldCoords elementNormal;

        if (dim==1) {
            elementNormal[0] = edges[0][1]; elementNormal[1] = -edges[0][0];
        } else
            elementNormal = crossProduct(edges[0], edges[1]);

        elementNormal /= elementNormal.two_norm();

        for (int j=0; j<nCorners;j++)
            normals[elements[offset + j]] += elementNormal;

        offset += nCorners;
    }

    for (size_t i=0; i<coords.size(); i++)
        normals[i] /= normals[i].two_norm();
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::removeDoubles(std::vector<std::array<LocalCoords,2> >& polytopeCorners)
{

    size_t counter(1);
    for (size_t i=1; i<polytopeCorners.size(); i++) {
        bool contained = false;
        for (size_t j=0; j<counter; j++)
            if ( (polytopeCorners[j][0]-polytopeCorners[i][0]).two_norm()<1e-10) {
                assert((polytopeCorners[j][1]-polytopeCorners[i][1]).two_norm()<1e-10);
                contained = true;
                break;
            }

        if (!contained) {
            if (counter < i)
                polytopeCorners[counter] = polytopeCorners[i];
            counter++;
        }
    }
    polytopeCorners.resize(counter);
}

} /* namespace GridGlue */
} /* namespace Dune */
