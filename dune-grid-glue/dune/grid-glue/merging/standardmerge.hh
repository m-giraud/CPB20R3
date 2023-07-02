// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect
 */

#ifndef DUNE_GRIDGLUE_MERGING_STANDARDMERGE_HH
#define DUNE_GRIDGLUE_MERGING_STANDARDMERGE_HH


#include <iostream>
#include <iomanip>
#include <vector>
#include <stack>
#include <set>
#include <utility>
#include <map>
#include <memory>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/timer.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/intersectionlist.hh>
#include <dune/grid-glue/merging/merger.hh>
#include <dune/grid-glue/merging/computeintersection.hh>

namespace Dune {
namespace GridGlue {

/** \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect

   Many merger algorithms consist of two parts: on the one hand there is a mechanism that produces pairs of
   elements that may intersect.  On the other hand there is an algorithm that computes the intersection of two
   given elements.  For the pairs-producing algorithm there appears to be a canonical choice, namely the algorithm
   by Gander and Japhet described in 'An Algorithm for Non-Matching Grid Projections with Linear Complexity,
   M.J. Gander and C. Japhet, Domain Decomposition Methods in Science and Engineering XVIII,
   pp. 185--192, Springer-Verlag, 2009.'  This class implements this algorithm, calling a pure virtual function
   computeIntersection() to compute the intersection between two elements.  Actual merger implementations
   can derive from this class and only implement computeIntersection().

   \tparam T The type used for coordinates (assumed to be the same for both grids)
   \tparam grid1Dim Dimension of the grid1 grid
   \tparam grid2Dim Dimension of the grid2 grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template<class T, int grid1Dim, int grid2Dim, int dimworld>
class StandardMerge
  : public Merger<T,grid1Dim,grid2Dim,dimworld>
{
  using Base = Merger<T, grid1Dim, grid2Dim, dimworld>;

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief Type used for local coordinates on the grid1 side
  using Grid1Coords = typename Base::Grid1Coords;

  /// @brief Type used for local coordinates on the grid2 side
  using Grid2Coords = typename Base::Grid2Coords;

  /// @brief the coordinate type used in this interface
  using WorldCoords = typename Base::WorldCoords;

  using IntersectionList = typename Base::IntersectionList;

protected:

  /** \brief The computed intersections */
  using IntersectionListProvider = SimplicialIntersectionListProvider<grid1Dim, grid2Dim>;
  using SimplicialIntersection = typename IntersectionListProvider::SimplicialIntersection;
  using RemoteSimplicialIntersection = SimplicialIntersection;

  bool valid = false;

  StandardMerge()
    : intersectionListProvider_(std::make_shared<IntersectionListProvider>())
    , intersectionList_(std::make_shared<IntersectionList>(intersectionListProvider_))
    {}

  virtual ~StandardMerge() = default;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices stored in the vector intersections.
   */
  virtual void computeIntersections(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                                   unsigned int grid2Index,
                                   std::vector<SimplicialIntersection>& intersections) = 0;

  /** \brief Compute the intersection between two overlapping elements
   * \return true if at least one intersection point was found
   */
  bool computeIntersection(unsigned int candidate0, unsigned int candidate1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                           const std::vector<Dune::GeometryType>& grid1_element_types,
                           std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                           const std::vector<Dune::GeometryType>& grid2_element_types,
                           std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                           bool insert = true);

  /*   M E M B E R   V A R I A B L E S   */

  std::shared_ptr<IntersectionListProvider> intersectionListProvider_;
  std::shared_ptr<IntersectionList> intersectionList_;

  /** \brief Temporary internal data */
  std::vector<std::vector<unsigned int> > grid1ElementCorners_;
  std::vector<std::vector<unsigned int> > grid2ElementCorners_;

  std::vector<std::vector<int> > elementNeighbors1_;
  std::vector<std::vector<int> > elementNeighbors2_;

public:

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @copydoc Merger<T,grid1Dim,grid2Dim,dimworld>::build
   */
  void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1_Coords,
             const std::vector<unsigned int>& grid1_elements,
             const std::vector<Dune::GeometryType>& grid1_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& grid2_coords,
             const std::vector<unsigned int>& grid2_elements,
             const std::vector<Dune::GeometryType>& grid2_element_types) override;


  /*   P R O B I N G   T H E   M E R G E D   G R I D   */

  void clear() override
  {
    // Delete old internal data, from a possible previous run
    intersectionListProvider_->clear();
    purge(grid1ElementCorners_);
    purge(grid2ElementCorners_);

    valid = false;
  }

  std::shared_ptr<IntersectionList> intersectionList() const final
  {
    assert(valid);
    return intersectionList_;
  }

  void enableFallback(bool fallback)
  {
      m_enableFallback = fallback;
  }

  void enableBruteForce(bool bruteForce)
  {
      m_enableBruteForce = bruteForce;
  }

private:
  /**
   * Enable fallback in case the advancing-front algorithm does not find an intersection.
   */
  bool m_enableFallback = false;

  /**
   * Enable brute force implementation instead of advancing-front algorithm.
   */
  bool m_enableBruteForce = false;

  auto& intersections()
    { return intersectionListProvider_->intersections(); }

  /** clear arbitrary containers */
  template<typename V>
  static void purge(V & v)
  {
    v.clear();
    V v2(v);
    v.swap(v2);
  }

  /**
   * Do a brute-force search to find one pair of intersecting elements
   * to start or continue the advancing-front type algorithm.
   */
  void generateSeed(std::vector<int>& seeds,
                    Dune::BitSetVector<1>& isHandled2,
                    std::stack<unsigned>& candidates2,
                    const std::vector<Dune::FieldVector<T, dimworld> >& grid1Coords,
                    const std::vector<Dune::GeometryType>& grid1_element_types,
                    const std::vector<Dune::FieldVector<T, dimworld> >& grid2Coords,
                    const std::vector<Dune::GeometryType>& grid2_element_types);

  /**
    * Insert intersections into this->intersection_ and return index
    */
  int insertIntersections(unsigned int candidate1, unsigned int candidate2,std::vector<SimplicialIntersection>& intersections);

  /**
   * Find a grid2 element intersecting the candidate1 grid1 element by brute force search
   */
  int bruteForceSearch(int candidate1,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                       const std::vector<Dune::GeometryType>& grid1_element_types,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                       const std::vector<Dune::GeometryType>& grid2_element_types);

  /**
   * Get the index of the intersection in intersections_ (= size if it is a new intersection)
   */
  std::pair<bool, unsigned int>
  intersectionIndex(unsigned int grid1Index, unsigned int grid2Index,
                                 SimplicialIntersection& intersection);

  /**
   * get the neighbor relations between the given elements
   */
  template <int gridDim>
  void computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                                  const std::vector<std::vector<unsigned int> >& gridElementCorners,
                                  std::vector<std::vector<int> >& elementNeighbors);

  void buildAdvancingFront(
    const std::vector<Dune::FieldVector<T,dimworld> >& grid1_Coords,
    const std::vector<unsigned int>& grid1_elements,
    const std::vector<Dune::GeometryType>& grid1_element_types,
    const std::vector<Dune::FieldVector<T,dimworld> >& grid2_coords,
    const std::vector<unsigned int>& grid2_elements,
    const std::vector<Dune::GeometryType>& grid2_element_types
    );

  void buildBruteForce(
    const std::vector<Dune::FieldVector<T,dimworld> >& grid1_Coords,
    const std::vector<unsigned int>& grid1_elements,
    const std::vector<Dune::GeometryType>& grid1_element_types,
    const std::vector<Dune::FieldVector<T,dimworld> >& grid2_coords,
    const std::vector<unsigned int>& grid2_elements,
    const std::vector<Dune::GeometryType>& grid2_element_types
    );
};


/* IMPLEMENTATION */

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::computeIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                      const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                      std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                      const std::vector<Dune::GeometryType>& grid2_element_types,
                                                                      std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                                                                      bool insert)
{
  // Select vertices of the grid1 element
  int grid1NumVertices = grid1ElementCorners_[candidate0].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid1ElementCorners(grid1NumVertices);
  for (int i=0; i<grid1NumVertices; i++)
    grid1ElementCorners[i] = grid1Coords[grid1ElementCorners_[candidate0][i]];

  // Select vertices of the grid2 element
  int grid2NumVertices = grid2ElementCorners_[candidate1].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid2ElementCorners(grid2NumVertices);
  for (int i=0; i<grid2NumVertices; i++)
    grid2ElementCorners[i] = grid2Coords[grid2ElementCorners_[candidate1][i]];

  // ///////////////////////////////////////////////////////
  //   Compute the intersection between the two elements
  // ///////////////////////////////////////////////////////

  std::vector<SimplicialIntersection> intersections(0);

  // compute the intersections
  computeIntersections(grid1_element_types[candidate0], grid1ElementCorners,
                       neighborIntersects1, candidate0,
                       grid2_element_types[candidate1], grid2ElementCorners,
                       neighborIntersects2, candidate1,
                       intersections);

  // insert intersections if needed
  if(insert && !intersections.empty())
      insertIntersections(candidate0,candidate1,intersections);

  // Have we found an intersection?
  return !intersections.empty() || neighborIntersects1.any() || neighborIntersects2.any();

}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::bruteForceSearch(int candidate1,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                  const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                  const std::vector<Dune::GeometryType>& grid2_element_types)
{
  std::bitset<(1<<grid1Dim)> neighborIntersects1;
  std::bitset<(1<<grid2Dim)> neighborIntersects2;
  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    bool intersectionFound = computeIntersection(i, candidate1,
                                                 grid1Coords, grid1_element_types, neighborIntersects1,
                                                 grid2Coords, grid2_element_types, neighborIntersects2,
                                                 false);

    // if there is an intersection, i is our new seed candidate on the grid1 side
    if (intersectionFound)
      return i;

  }

  return -1;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
template<int gridDim>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::
computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                           const std::vector<std::vector<unsigned int> >& gridElementCorners,
                           std::vector<std::vector<int> >& elementNeighbors)
{
  typedef std::vector<unsigned int> FaceType;
  typedef std::map<FaceType, std::pair<unsigned int, unsigned int> > FaceSetType;

  ///////////////////////////////////////////////////////////////////////////////////////
  //  First: grid 1
  ///////////////////////////////////////////////////////////////////////////////////////
  FaceSetType faces;
  elementNeighbors.resize(gridElementTypes.size());

  for (size_t i=0; i<gridElementTypes.size(); i++)
    elementNeighbors[i].resize(Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]).size(1), -1);

  for (size_t i=0; i<gridElementTypes.size(); i++) { //iterate over all elements
    const auto& refElement = Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]);

    for (size_t j=0; j<(size_t)refElement.size(1); j++) { // iterate over all faces of the element

      FaceType face;
      // extract element face
      for (size_t k=0; k<(size_t)refElement.size(j,1,gridDim); k++)
        face.push_back(gridElementCorners[i][refElement.subEntity(j,1,k,gridDim)]);

      // sort the face vertices to get rid of twists and other permutations
      std::sort(face.begin(), face.end());

      typename FaceSetType::iterator faceHandle = faces.find(face);

      if (faceHandle == faces.end()) {

        // face has not been visited before
        faces.insert(std::make_pair(face, std::make_pair(i,j)));

      } else {

        // face has been visited before: store the mutual neighbor information
        elementNeighbors[i][j] = faceHandle->second.first;
        elementNeighbors[faceHandle->second.first][faceHandle->second.second] = i;

        faces.erase(faceHandle);

      }

    }

  }
}

// /////////////////////////////////////////////////////////////////////
//   Compute the intersection of all pairs of elements
//   Linear algorithm by Gander and Japhet, Proc. of DD18
// /////////////////////////////////////////////////////////////////////

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                        const std::vector<unsigned int>& grid1_elements,
                                                        const std::vector<Dune::GeometryType>& grid1_element_types,
                                                        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                        const std::vector<unsigned int>& grid2_elements,
                                                        const std::vector<Dune::GeometryType>& grid2_element_types
                                                        )
{

  std::cout << "StandardMerge building merged grid..." << std::endl;
  Dune::Timer watch;

  clear();
  // clear global intersection list
  intersectionListProvider_->clear();
  this->counter = 0;

  // /////////////////////////////////////////////////////////////////////
  //   Copy element corners into a data structure with block-structure.
  //   This is not as efficient but a lot easier to use.
  //   We may think about efficiency later.
  // /////////////////////////////////////////////////////////////////////

  // first the grid1 side
  grid1ElementCorners_.resize(grid1_element_types.size());

  unsigned int grid1CornerCounter = 0;

  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    // Select vertices of the grid1 element
    int numVertices = Dune::ReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
    grid1ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid1ElementCorners_[i][j] = grid1_elements[grid1CornerCounter++];

  }

  // then the grid2 side
  grid2ElementCorners_.resize(grid2_element_types.size());

  unsigned int grid2CornerCounter = 0;

  for (std::size_t i=0; i<grid2_element_types.size(); i++) {

    // Select vertices of the grid2 element
    int numVertices = Dune::ReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
    grid2ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid2ElementCorners_[i][j] = grid2_elements[grid2CornerCounter++];

  }

  ////////////////////////////////////////////////////////////////////////
  //  Compute the face neighbors for each element
  ////////////////////////////////////////////////////////////////////////

  computeNeighborsPerElement<grid1Dim>(grid1_element_types, grid1ElementCorners_, elementNeighbors1_);
  computeNeighborsPerElement<grid2Dim>(grid2_element_types, grid2ElementCorners_, elementNeighbors2_);

  std::cout << "setup took " << watch.elapsed() << " seconds." << std::endl;

  if (m_enableBruteForce)
    buildBruteForce(grid1Coords, grid1_elements, grid1_element_types, grid2Coords, grid2_elements, grid2_element_types);
  else
    buildAdvancingFront(grid1Coords, grid1_elements, grid1_element_types, grid2Coords, grid2_elements, grid2_element_types);

  valid = true;
  std::cout << "intersection construction took " << watch.elapsed() << " seconds." << std::endl;
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::buildAdvancingFront(
  const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
  const std::vector<unsigned int>& grid1_elements,
  const std::vector<Dune::GeometryType>& grid1_element_types,
  const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
  const std::vector<unsigned int>& grid2_elements,
  const std::vector<Dune::GeometryType>& grid2_element_types
  )
{
  ////////////////////////////////////////////////////////////////////////
  //   Data structures for the advancing-front algorithm
  ////////////////////////////////////////////////////////////////////////

  std::stack<unsigned int> candidates1;
  std::stack<unsigned int> candidates2;

  std::vector<int> seeds(grid2_element_types.size(), -1);

  // /////////////////////////////////////////////////////////////////////
  //   Do a brute-force search to find one pair of intersecting elements
  //   to start the advancing-front type algorithm with.
  // /////////////////////////////////////////////////////////////////////

  // Set flag if element has been handled
  Dune::BitSetVector<1> isHandled2(grid2_element_types.size());

  // Set flag if the element has been entered in the queue
  Dune::BitSetVector<1> isCandidate2(grid2_element_types.size());

  generateSeed(seeds, isHandled2, candidates2, grid1Coords, grid1_element_types, grid2Coords, grid2_element_types);

  // /////////////////////////////////////////////////////////////////////
  //   Main loop
  // /////////////////////////////////////////////////////////////////////

  std::set<unsigned int> isHandled1;

  std::set<unsigned int> isCandidate1;

  while (!candidates2.empty()) {

    // Get the next element on the grid2 side
    unsigned int currentCandidate2 = candidates2.top();
    int seed = seeds[currentCandidate2];
    assert(seed >= 0);

    candidates2.pop();
    isHandled2[currentCandidate2] = true;

    // Start advancing front algorithm on the grid1 side from the 'seed' element that
    // we stored along with the current grid2 element
    candidates1.push(seed);

    isHandled1.clear();
    isCandidate1.clear();

    while (!candidates1.empty()) {

      unsigned int currentCandidate1 = candidates1.top();
      candidates1.pop();
      isHandled1.insert(currentCandidate1);

      // Test whether there is an intersection between currentCandidate0 and currentCandidate1
      std::bitset<(1<<grid1Dim)> neighborIntersects1;
      std::bitset<(1<<grid2Dim)> neighborIntersects2;
      bool intersectionFound = computeIntersection(currentCandidate1, currentCandidate2,
                                                   grid1Coords,grid1_element_types, neighborIntersects1,
                                                   grid2Coords,grid2_element_types, neighborIntersects2);

      for (size_t i=0; i<neighborIntersects2.size(); i++)
        if (neighborIntersects2[i] && elementNeighbors2_[currentCandidate2][i] != -1)
          seeds[elementNeighbors2_[currentCandidate2][i]] = currentCandidate1;

      // add neighbors of candidate0 to the list of elements to be checked
      if (intersectionFound) {

        for (size_t i=0; i<elementNeighbors1_[currentCandidate1].size(); i++) {

          int neighbor = elementNeighbors1_[currentCandidate1][i];

          if (neighbor == -1)            // do nothing at the grid boundary
            continue;

          if (isHandled1.find(neighbor) == isHandled1.end()
              && isCandidate1.find(neighbor) == isCandidate1.end()) {
            candidates1.push(neighbor);
            isCandidate1.insert(neighbor);
          }

        }

      }

    }

    // We have now found all intersections of elements in the grid1 side with currentCandidate2
    // Now we add all neighbors of currentCandidate2 that have not been treated yet as new
    // candidates.

    // Do we have an unhandled neighbor with a seed?
    bool seedFound = !candidates2.empty();
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      // Add all unhandled intersecting neighbors to the queue
      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0] && seeds[neighbor]>-1) {

        isCandidate2[neighbor][0] = true;
        candidates2.push(neighbor);
        seedFound = true;
      }
    }

    if (seedFound || !m_enableFallback)
      continue;

    // There is no neighbor with a seed, so we need to be a bit more aggressive...
    // get all neighbors of currentCandidate2, but not currentCandidate2 itself
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0]) {

        // Get a seed element for the new grid2 element
        // Look for an element on the grid1 side that intersects the new grid2 element.
        int seed = -1;

        // Look among the ones that have been tested during the last iteration.
        for (typename std::set<unsigned int>::iterator seedIt = isHandled1.begin();
             seedIt != isHandled1.end(); ++seedIt) {

          std::bitset<(1<<grid1Dim)> neighborIntersects1;
          std::bitset<(1<<grid2Dim)> neighborIntersects2;
          bool intersectionFound = computeIntersection(*seedIt, neighbor,
                                                       grid1Coords, grid1_element_types, neighborIntersects1,
                                                       grid2Coords, grid2_element_types, neighborIntersects2,
                                                       false);

          // if the intersection is nonempty, *seedIt is our new seed candidate on the grid1 side
          if (intersectionFound) {
            seed = *seedIt;
            Dune::dwarn << "Algorithm entered first fallback method and found a new seed in the build algorithm." <<
                     "Probably, the neighborIntersects bitsets computed in computeIntersection specialization is wrong." << std::endl;
            break;
          }

        }

        if (seed < 0) {
          // The fast method didn't find a grid1 element that intersects with
          // the new grid2 candidate.  We have to do a brute-force search.
          seed = bruteForceSearch(neighbor,
                                  grid1Coords,grid1_element_types,
                                  grid2Coords,grid2_element_types);
          Dune::dwarn << "Algorithm entered second fallback method. This probably should not happen." << std::endl;

        }

        // We have tried all we could: the candidate is 'handled' now
        isCandidate2[neighbor] = true;

        // still no seed?  Then the new grid2 candidate isn't overlapped by anything
        if (seed < 0)
          continue;

        // we have a seed now
        candidates2.push(neighbor);
        seeds[neighbor] = seed;
        seedFound = true;

      }

    }

    /* Do a brute-force search if there is still no seed:
     * There might still be a disconnected region out there.
     */
    if (!seedFound && candidates2.empty()) {
      generateSeed(seeds, isHandled2, candidates2, grid1Coords, grid1_element_types, grid2Coords, grid2_element_types);
    }
  }
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::buildBruteForce(
  const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
  const std::vector<unsigned int>& grid1_elements,
  const std::vector<Dune::GeometryType>& grid1_element_types,
  const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
  const std::vector<unsigned int>& grid2_elements,
  const std::vector<Dune::GeometryType>& grid2_element_types
  )
{
  std::bitset<(1<<grid1Dim)> neighborIntersects1;
  std::bitset<(1<<grid2Dim)> neighborIntersects2;

  for (unsigned i = 0; i < grid1_element_types.size(); ++i) {
    for (unsigned j = 0; j < grid2_element_types.size(); ++j) {
      (void) computeIntersection(i, j,
                                 grid1Coords, grid1_element_types, neighborIntersects1,
                                 grid2Coords, grid2_element_types, neighborIntersects2);
    }
  }
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::generateSeed(std::vector<int>& seeds, Dune::BitSetVector<1>& isHandled2, std::stack<unsigned>& candidates2, const std::vector<Dune::FieldVector<T, dimworld> >& grid1Coords, const std::vector<Dune::GeometryType>& grid1_element_types, const std::vector<Dune::FieldVector<T, dimworld> >& grid2Coords, const std::vector<Dune::GeometryType>& grid2_element_types)
{
  for (std::size_t j=0; j<grid2_element_types.size(); j++) {

    if (seeds[j] > 0 || isHandled2[j][0])
      continue;

    int seed = bruteForceSearch(j,grid1Coords,grid1_element_types,grid2Coords,grid2_element_types);

    if (seed >= 0) {
      candidates2.push(j);        // the candidate and a seed for the candidate
      seeds[j] = seed;
      break;
    } else // If the brute force search did not find any intersection we can skip this element
      isHandled2[j] = true;
  }
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::insertIntersections(unsigned int candidate1, unsigned int candidate2,
                                                                        std::vector<SimplicialIntersection>& intersections)
{
    typedef typename std::vector<SimplicialIntersection>::size_type size_t;
    int count = 0;

    for (size_t i = 0; i < intersections.size(); ++i) {
        // get the intersection index of the current intersection from intersections in this->intersections
      bool found;
      unsigned int index;
      std::tie(found, index) = intersectionIndex(candidate1,candidate2,intersections[i]);

        if (found && index >= this->intersections().size()) { //the intersection is not yet contained in this->intersections
            this->intersections().push_back(intersections[i]);   // insert

            ++count;
        } else if (found) {
            auto& intersection = this->intersections()[index];

            // insert each grid1 element and local representation of intersections[i] with parent candidate1
            for (size_t j = 0; j < intersections[i].parents0.size(); ++j) {
                intersection.parents0.push_back(candidate1);
                intersection.corners0.push_back(intersections[i].corners0[j]);
            }

            // insert each grid2 element and local representation of intersections[i] with parent candidate2
            for (size_t j = 0; j < intersections[i].parents1.size(); ++j) {
                intersection.parents1.push_back(candidate2);
                intersection.corners1.push_back(intersections[i].corners1[j]);
            }

            ++count;
        } else {
            Dune::dwarn << "Computed the same intersection twice!" << std::endl;
        }
    }
    return count;
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
std::pair<bool, unsigned int>
StandardMerge<T,grid1Dim,grid2Dim,dimworld>::intersectionIndex(unsigned int grid1Index, unsigned int grid2Index,
                                                                            SimplicialIntersection& intersection) {


    // return index in intersections_ if at least one local representation of a Simplicial Intersection (SI)
    // of intersections_ is equal to the local representation of one element in intersections

    std::size_t n_intersections = this->intersections().size();
    if (grid1Dim == grid2Dim)
      return {true, n_intersections};

    T eps = 1e-10;

    for (std::size_t i = 0; i < n_intersections; ++i) {

        // compare the local representation of the subelements of the SI
        for (std::size_t ei = 0; ei < this->intersections()[i].parents0.size(); ++ei) // merger subelement
        {
            if (this->intersections()[i].parents0[ei] == grid1Index)
            {
                for (std::size_t er = 0; er < intersection.parents0.size(); ++er) // list subelement
                {
                    bool found_all = true;
                    // compare the local coordinate representations
                    for (std::size_t ci = 0; ci < this->intersections()[i].corners0[ei].size(); ++ci)
                    {
                        Dune::FieldVector<T,grid1Dim> ni = this->intersections()[i].corners0[ei][ci];
                        bool found_ni = false;
                        for (std::size_t cr = 0; cr < intersection.corners0[er].size(); ++cr)
                        {
                            Dune::FieldVector<T,grid1Dim> nr = intersection.corners0[er][cr];

                            found_ni = found_ni || ((ni-nr).infinity_norm() < eps);
                            if (found_ni)
                                break;
                        }
                        found_all = found_all && found_ni;

                        if (!found_ni)
                            break;
                    }

                    if (found_all && (this->intersections()[i].parents1[ei] != grid2Index))
                      return {true, i};
                    else if (found_all)
                      return {false, 0};
                }
            }
        }

        // compare the local representation of the subelements of the SI
        for (std::size_t ei = 0; ei < this->intersections()[i].parents1.size(); ++ei) // merger subelement
        {
            if (this->intersections()[i].parents1[ei] == grid2Index)
            {
                for (std::size_t er = 0; er < intersection.parents1.size(); ++er) // list subelement
                {
                    bool found_all = true;
                    // compare the local coordinate representations
                    for (std::size_t ci = 0; ci < this->intersections()[i].corners1[ei].size(); ++ci)
                    {
                        Dune::FieldVector<T,grid2Dim> ni = this->intersections()[i].corners1[ei][ci];
                        bool found_ni = false;
                        for (std::size_t cr = 0; cr < intersection.corners1[er].size(); ++cr)
                        {
                            Dune::FieldVector<T,grid2Dim> nr = intersection.corners1[er][cr];
                            found_ni = found_ni || ((ni-nr).infinity_norm() < eps);

                            if (found_ni)
                                break;
                        }
                        found_all = found_all && found_ni;

                        if (!found_ni)
                            break;
                    }

                    if (found_all && (this->intersections()[i].parents0[ei] != grid1Index))
                      return {true, i};
                    else if (found_all)
                      return {false, 0};
                }
            }
        }
    }

    return {true, n_intersections};
}

#define DECL extern
#define STANDARD_MERGE_INSTANTIATE(T,A,B,C) \
  DECL template \
  void StandardMerge<T,A,B,C>::build(const std::vector<Dune::FieldVector<T,C> >& grid1Coords, \
                                     const std::vector<unsigned int>& grid1_elements, \
                                     const std::vector<Dune::GeometryType>& grid1_element_types, \
                                     const std::vector<Dune::FieldVector<T,C> >& grid2Coords, \
                                     const std::vector<unsigned int>& grid2_elements, \
                                     const std::vector<Dune::GeometryType>& grid2_element_types \
                                     )

STANDARD_MERGE_INSTANTIATE(double,1,1,1);
STANDARD_MERGE_INSTANTIATE(double,2,2,2);
STANDARD_MERGE_INSTANTIATE(double,3,3,3);
#undef STANDARD_MERGE_INSTANTIATE
#undef DECL

} /* namespace GridGlue */
} /* namespace Dune */

#endif // DUNE_GRIDGLUE_MERGING_STANDARDMERGE_HH
