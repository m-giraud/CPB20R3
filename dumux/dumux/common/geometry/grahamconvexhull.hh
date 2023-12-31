/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief A function to compute the convex hull of a point cloud
 *        and a function to triangulate the polygon area spanned by the convex hull
 */
#ifndef DUMUX_GRAHAM_CONVEX_HULL_HH
#define DUMUX_GRAHAM_CONVEX_HULL_HH

#include <vector>
#include <array>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Returns the orientation of a sequence a-->b-->c in one plane (defined by normal vector)
 * \return -1   if a-->b-->c forms a counter-clockwise turn (given the normal vector)
 *         +1   for a clockwise turn,
 *          0   if they are on one line (colinear)
 */
template<class ctype>
int getOrientation(const Dune::FieldVector<ctype, 3>& a,
                   const Dune::FieldVector<ctype, 3>& b,
                   const Dune::FieldVector<ctype, 3>& c,
                   const Dune::FieldVector<ctype, 3>& normal)
{
    const auto d = b-a;
    const auto e = c-b;
    const auto f = Dumux::crossProduct(d, e);
    const auto area = f*normal;
    return Dumux::sign(-area);
}

/*!
 * \ingroup Geometry
 * \brief Compute the points making up the convex hull around the given set of unordered points
 * \note We assume that all points are coplanar and there are no indentical points in the list
 */
template<class ctype>
std::vector<Dune::FieldVector<ctype, 3>>
grahamConvexHull2d3d(const std::vector<Dune::FieldVector<ctype, 3>>& points)
{
    auto copyPoints = points;
    return grahamConvexHull2d3d(copyPoints);
}

/*!
 * \ingroup Geometry
 * \brief Compute the points making up the convex hull around the given set of unordered points
 * \note We assume that all points are coplanar and there are no indentical points in the list
 * \note This algorithm changes the order of the given points a bit
 *       as they are unordered anyway this shouldn't matter too much
 */
template<class ctype>
std::vector<Dune::FieldVector<ctype, 3>>
grahamConvexHull2d3d(std::vector<Dune::FieldVector<ctype, 3>>& points)
{
    using Point = Dune::FieldVector<ctype, 3>;
    std::vector<Point> convexHull;

    // return empty convex hull
    if (points.size() < 3)
        return convexHull;

    // return the points (already just one triangle)
    if (points.size() == 3)
        return points;

    // try to compute the normal vector of the plane
    const auto a = points[1] - points[0];
    auto b = points[2] - points[0];
    auto normal = Dumux::crossProduct(a, b);

    // make sure the normal vector has non-zero length
    std::size_t k = 2;
    auto norm = normal.two_norm();
    while (norm == 0.0 && k < points.size()-1)
    {
        b = points[++k];
        normal = Dumux::crossProduct(a, b);
        norm = normal.two_norm();
    }

    // if all given points are colinear -> return empty convex hull
    if (norm == 0.0)
        return convexHull;

    using std::sqrt;
    const auto eps = 1e-7*sqrt(norm);
    normal /= norm;

    // find the element with the smallest x coordinate (if x is the same, smallest y coordinate, and so on...)
    auto minIt = std::min_element(points.begin(), points.end(), [&eps](const auto& a, const auto& b)
    {
        using std::abs;
        return (abs(a[0]-b[0]) > eps ? a[0] < b[0] : (abs(a[1]-b[1]) > eps ? a[1] < b[1] : (a[2] < b[2])));
    });

    // swap the smallest element to the front
    std::iter_swap(minIt, points.begin());

    // choose the first (min element) as the pivot point
    // sort in counter-clockwise order around pivot point
    const auto pivot = points[0];
    std::sort(points.begin()+1, points.end(), [&](const auto& a, const auto& b)
    {
        const int order = getOrientation(pivot, a, b, normal);
        if (order == 0)
             return (a-pivot).two_norm() < (b-pivot).two_norm();
        else
             return (order == -1);
    });

    // push the first three points
    convexHull.reserve(50);
    convexHull.push_back(points[0]);
    convexHull.push_back(points[1]);
    convexHull.push_back(points[2]);

    // This is the heart of the algorithm
    // pop_back until the last point in the queue forms a counter-clockwise oriented line
    // with the next two vertices. Then add these points to the queue.
    for (std::size_t i = 3; i < points.size(); ++i)
    {
        Point p = convexHull.back();
        convexHull.pop_back();
        // keep popping until the orientation a->b->currentp is counter-clockwise
        while (getOrientation(convexHull.back(), p, points[i], normal) != -1)
        {
            // make sure the queue doesn't get empty
            if (convexHull.size() == 1)
            {
                // before we reach i=size-1 there has to be a good candidate
                // as not all points are colinear (a non-zero plane normal exists)
                assert(i < points.size()-1);
                p = points[i++];
            }
            else
            {
                p = convexHull.back();
                convexHull.pop_back();
            }
        }

        // add back the last popped point and this point
        convexHull.emplace_back(std::move(p));
        convexHull.push_back(points[i]);
    }

    return convexHull;
}

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points of the convex hull
 * \note Assumes all points of the convex hull are coplanar
 * \note This inserts a mid point and connects all corners with that point to triangles
 */
template<class ctype>
std::vector<std::array<Dune::FieldVector<ctype, 3>, 3> >
triangulateConvexHull(const std::vector<Dune::FieldVector<ctype, 3>>& convexHull)
{
    using Point = Dune::FieldVector<ctype, 3>;
    using Triangle = std::array<Point, 3>;

    if (convexHull.size() < 3)
        DUNE_THROW(Dune::InvalidStateException, "Try to triangulate point cloud with less than 3 points!");

    if (convexHull.size() == 3)
        return std::vector<Triangle>(1, {convexHull[0], convexHull[1], convexHull[2]});

    Point midPoint(0.0);
    for (const auto p : convexHull)
        midPoint += p;
    midPoint /= convexHull.size();

    std::vector<Triangle> triangulation;
    triangulation.reserve(convexHull.size());

    for (std::size_t i = 0; i < convexHull.size()-1; ++i)
        triangulation.emplace_back(Triangle{midPoint, convexHull[i], convexHull[i+1]});

    triangulation.emplace_back(Triangle{midPoint, convexHull[convexHull.size()-1], convexHull[0]});

    return triangulation;
}

} // end namespace Dumux

# endif
