#ifndef DUNE_GRIDGLUE_MERGING_COMPUTEINTERSECTION_HH
#define DUNE_GRIDGLUE_MERGING_COMPUTEINTERSECTION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune {
namespace GridGlue {

template<int dimWorld, int dim1, int dim2, typename T = double>
class ComputationMethod {
public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = dim1;
    static const int grid2Dimension = dim2;
    static const int intersectionDimension = (dim1 < dim2)?(dim1):(dim2);

    static bool computeIntersectionPoints(const std::vector<Vector> X,
                                                  const std::vector<Vector> Y,
                                                  std::vector<std::vector<int> >& SX,
                                                  std::vector<std::vector<int> >& SY,
                                                  std::vector<Vector>& P);
    static void grid1_subdivisions(const std::vector<Vector> elementCorners,
                                     std::vector<std::vector<unsigned int> >& subElements,
                                     std::vector<std::vector<int> >& faceIds);
    static void grid2_subdivisions(const std::vector<Vector> elementCorners,
                                     std::vector<std::vector<unsigned int> >& subElements,
                                     std::vector<std::vector<int> >& faceIds);
};

/**
 *  @brief Intersection computation method for two elements of arbitrary dimension
 *  @tparam dimWorld The world dimension
 *  @tparam T the field type
 */
template<class CM>
class IntersectionComputation {
private:
    typedef typename CM::Vector V;
    const int dimWorld = V::dimension;
    const int dim1 = CM::grid1Dimension;
    const int dim2 = CM::grid2Dimension;
public:
    /**
     * @brief Compute the intersection of two elements X and Y
     * Compute the intersection of two elements X and Y, where X is of dimension dim1
     * and Y is of dimension dim2 and return a vector P containing the corner points of
     * the intersection polyhedron.
     *
     * @param X ordered vector of corners defining the element with dimension dim1
     * @param Y ordered vector of corners defining the element with dimension dim2
     * @param SX indices of points in P on the index i face of X, where i denotes the vector index in SX
     * @param SY indices of points in P on the index i face of Y, where i denotes the vector index in SY
     * @param P the intersection points
     * @return true if enough intersection points were found to define at least one intersection element
     */
    static bool computeIntersection(const std::vector<V>& X,
                                    const std::vector<V>& Y,
                                    std::vector<std::vector<int> >& SX,
                                    std::vector<std::vector<int> >& SY,
                                    std::vector<V>& P);

    /**
     * @brief Order Points in the point list P face-wise such that a subsimplex subdivision can be constructed
     * @tparam isDim the intersection geometry dimension
     * @param centroid the center of all points
     * @param SX indices of points in P on the index i face of X, where i denotes the vector index in SX
     * @param SY indices of points in P on the index i face of Y, where i denotes the vector index in SY
     * @param P the point list (remains unchanged)
     * @param H ordered list of P indices
     */
    template<int isDim, int dW>
    static void orderPoints(const V& centroid,
                            const std::vector<std::vector<int> >& SX,
                            const std::vector<std::vector<int> >& SY,
                            const std::vector<V>& P,
                            std::vector<std::vector<int> >& H)
    {
        if (isDim > 1)
            orderPoints_(std::integral_constant<int,isDim>(),std::integral_constant<int,dW>(),
                    centroid, SX, SY, P,H);
    }

private:
    static void orderPoints_(std::integral_constant<int,1>,
                             std::integral_constant<int,1>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> >& H) {}
    static void orderPoints_(std::integral_constant<int,1>,
                             std::integral_constant<int,2>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> >& H) {}
    static void orderPoints_(std::integral_constant<int,1>,
                             std::integral_constant<int,3>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> >& H) {}
    static void orderPoints_(std::integral_constant<int,2>,
                             std::integral_constant<int,2>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> >& H);
    static void orderPoints_(std::integral_constant<int,2>,
                             std::integral_constant<int,3>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> >& H);
    static void orderPoints_(std::integral_constant<int,3>,
                             std::integral_constant<int,3>,
                             const V& centroid,
                             const std::vector<std::vector<int> >& SX,
                             const std::vector<std::vector<int> >& SY,
                             const std::vector<V>& P,
                             std::vector<std::vector<int> > & H);

    /**
     * @brief Order points counterclockwise
     * @tparam isDim the intersection geometry dimension
     * @param centroid the center of all points
     * @param id list of counterclockwise point indices in P
     * @param P the point list
     */
    static void orderPointsCC(std::integral_constant<int,2>,
                              const V& centroid,
                              std::vector<int> &id,
                              const std::vector<V>& P);
    static void orderPointsCC(std::integral_constant<int,3>,
                              const V& centroid,
                              std::vector<int> &id,
                              const std::vector<V>& P);

    /**
     * @brief Removes duplicate entries from the vector p.
     * @param the list duplicat entries are removed from
     */
    static void removeDuplicates( std::vector<int> & p);

    /**
     * @brief Checks if index set is contained already in H
     * b=NewFace(H,id) checks if a permutation of the vector id is contained in a row of the matrix H
     * @param id
     * @param H
     * @return true if the index set is contained in H
     */
    static bool newFace3D(const std::vector<int>& id,
                          const std::vector<std::vector<int> >& H);
};

template<class V>
inline int insertPoint(const V p, std::vector<V>& P)
{
    double eps= 1e-8;     // tolerance for identical nodes
    std::size_t k=0 ;

    if (P.size()>0) {

        while ((k<P.size())&&
               ((p - P[k]).infinity_norm()>eps*(P[k].infinity_norm()) &&
                (p - P[k]).infinity_norm()>eps*(p.infinity_norm())) &&
                !(p.infinity_norm() < eps && P[k].infinity_norm() <eps &&
                   (p - P[k]).infinity_norm() <  eps))
            k++ ;

        if (k>=P.size())
            P.push_back(p) ;        //  new node is not contained in P

    }
    else
        P.push_back(p);

    return k ;
}


} /* namespace Dune::GridGlue */
} /* namespace Dune */

#include "simplexintersection.cc"
#include "computeintersection.cc"

#endif
