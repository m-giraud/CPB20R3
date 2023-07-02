#ifndef DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH
#define DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER2_HH

#include <array>
#include <bitset>
#include <tuple>

namespace Dune {
namespace GridGlue {

/**
 * \brief Projection of a line (triangle) on another line (triangle).
 *
 * This class implements methods to project a line (2d) or triangle (3d) on
 * another line (triangle) along normal field given by values at the corners.
 */
template<typename Coordinate>
class Projection
{
public:
  /**
   * \brief Intersection between two edges of a triangle.
   *
   * See also \ref Projection<Coordinate>::edgeIntersections()
   */
  struct EdgeIntersection
  {
    /**
     * \brief Edge numbers in image and preimage triangle.
     */
    std::array<unsigned, 2> edge;

    /**
     * \brief Local coordinates of intersection and distance along normals.
     *
     * Local coordinate of intersection point in barycentric coordinates with
     * respect to image and preimage triangle.
     */
    std::array<Coordinate, 2> local;
  };

  /**
   * \brief dimension of coordinates
   */
  constexpr static unsigned dim = Coordinate::dimension;

  /**
   * \brief maximum number of edge-edge intersections
   *
   * See also \seealso edgeIntersections()
   */
  constexpr static unsigned maxEdgeIntersections = dim == 3 ? 9 : 0;

  static_assert(dim == 2 || dim == 3, "Projection only implemented for dim=2 or dim=3");

  /**
   * \brief Scalar type.
   */
  typedef typename Coordinate::field_type Field;

  /**
   * \brief List of corner images.
   *
   * This type is used to return the list of images Φ(xᵢ) of the corners xᵢ
   * in barycentric coordinates with respect to the image simplex.
   * The last entry is used to return the (signed) distance along the normal.
   */
  typedef std::array<Coordinate, dim> Images;

  /**
   * List of corner preimages.
   *
   * This is used as \ref Images, but for the preimages Φ⁻¹(yᵢ) of the corners
   * yᵢ of the image simplex.
   */
  typedef Images Preimages;

private:
  /**
   * \brief Overlap allowed for the projection to be considered valid.
   */
  const Field m_overlap;

  /**
   * \brief Maximum value for scalar product ν(x)·ν(Φ(x)) of normals
   *
   * The normals at <code>x</code> and <code>Φ(x)</code> are expected
   * to be opposing to some degree.  This value is used to indicate
   * how much they are allowed to deviate from this by ensuring that
   * <code>ν(x)·ν(Φ(x)) ≤ m_max_normal_product</code>.
   */
  const Field m_max_normal_product;

  /**
   * \brief epsilon used for floating-point comparisons.
   *
   * See also \seealso epsilon(Field)
   */
  Field m_epsilon = Field(1e-12);

  /** \copydoc images() */
  std::tuple<Images, Preimages> m_images;

  /** \copydoc success() */
  std::tuple<std::bitset<dim>, std::bitset<dim> > m_success;

  /** \copydoc numberOfEdgeIntersections() */
  unsigned m_number_of_edge_intersections;

  /** \copydoc edgeIntersections() */
  std::array<EdgeIntersection, maxEdgeIntersections> m_edge_intersections;

  /**
   * \brief Forward projection successful for all corners <code>xᵢ</code>
   *
   * If <code>true</code>, the forward projection was successful, that is
   * Φ(xᵢ) could be computed for all xᵢ.
   *
   * \warning Note that this only means Φ(xᵢ) lie in the plane spanned by the
   *          image simplex which is required to compute the inverse
   *          projection Φ⁻¹(yᵢ). The bitset \ref m_success should be used to
   *          check whether the projection is feasible.
   */
  bool m_projection_valid;

  /**
   * \brief Compute forward projection Φ(xᵢ) for all xᵢ.
   *
   * \copydetails project
   */
  template<typename Corners, typename Normals>
  void doProjection(const std::tuple<Corners&, Corners&>& corners, const std::tuple<Normals&, Normals&>& normals);

  /**
   * \brief Compute inverse projection Φ⁻¹(yᵢ) for all yᵢ.
   *
   * \note This requires the forward projection was already computed by
   *       \ref doProjection.
   *
   * \copydetails project
   */
  template<typename Corners, typename Normals>
  void doInverseProjection(const std::tuple<Corners&, Corners&>& corners, const std::tuple<Normals&, Normals&>& normals);

  /**
   * \brief Compute intersections between projected edges and edges of the image simplex.
   *
   * \note This requires the forward and inverse projections were already
   *       computed by \ref doProjection and \ref doInverseProjection.
   *
   * \copydetails project
   */
  template<typename Corners, typename Normals>
  void doEdgeIntersection(const std::tuple<Corners&, Corners&>& corners, const std::tuple<Normals&, Normals&>& normals);

  /**
   * \brief Check if projection is feasible.
   *
   * Given a point <code>x</code>, its image <code>px</code> in barycentric
   * coordinates together with the signed distance along the normal at
   * <code>x</code> in the last entry of <code>px</code> and the corners and
   * normals of the image simplex given in <code>corners</code> and
   * <code>normals</code>, this method checks that the projection is feasible.
   * This means:
   *
   * <ul>
   *   <li><code>px</code> is inside the image simplex</li>
   *   <li>The signed distance given is not smaller than <code>-\ref m_overlap</code></li>
   *   <li>The signed distance along the normal at <code>px</code> is not smaller than <code>-\ref m_overlap</li>
   *   <li>The angle between the normals at <code>x</code> and <code>px</code> is at least \ref m_minimum_angle_between_normals
   * </ul>
   *
   * \param x euclidean coordinate of point to project
   * \param nx outer normal ν(x) at <code>x</code>
   * \param px barycentric coordinates of projected point;
   *           last entry is distance along normal
   * \param corners corners of image simplex
   * \param normals normals of image simplex
   * \return <code>true</code> if the projection is feasible, <code>false</code> otherwise.
   */
  template<typename Corners, typename Normals>
  inline bool projectionFeasible(const Coordinate& x, const Coordinate& nx, const Coordinate& px, const Corners& corners, const Normals& normals) const;

public:
  /**
   * \param overlap allowed overlap
   * \param max_normal_product maximum value for scalar product ν(x)·ν(Φ(x))
   */
  Projection(const Field overlap = Field(0), const Field max_normal_product = Field(-0.1));

  /**
   * \brief Set epsilon used for floating-point comparisons.
   *
   * \param epsilon new epsilon used for floating-point comaprisons
   */
  void epsilon(const Field epsilon);

  /**
   * \brief Do the actual projection.
   *
   * \param corners euclidean coordinates of corners of preimage and image
   * \param normals normals at corners of preimage and image
   * \tparam Corners list of corner coordinates, should be
   *                 <code>std::vector<Coordinate></code> or
   *                 <code>std::array<Coordinate, n></code>
   * \tparam Normals list of corner normals, should be
   *                 <code>std::vector<Coordinate></code> or
   *                 <code>std::array<Coordinate, n></code>
   */
  template<typename Corners, typename Normals>
  void project(const std::tuple<Corners&, Corners&>& corners, const std::tuple<Normals&, Normals&>& normals);

  /**
   * \brief Images and preimages of corners.
   *
   * Returns a pair of arrays. The first array contains the images
   * <code>Φ(xᵢ)</code> of the corners <code>xᵢ</code>. The second
   * array contains the preimages <code>Φ⁻¹(yⱼ)</code> of the
   * corners <code>yⱼ</code>.
   *
   * The first d-1 values are the barycentric coordinates with respect
   * to the corners of the (pre)image, the last value is the signed
   * distance between the projected point and its (pre)image along the
   * normal at the projected preimage corner or the inverse projected
   * image corner.
   *
   * \note \ref project() must be called before this method can be used.
   *
   * \returns pair of arrays giving <code>((Φ(xᵢ))ᵢ, (Φ⁻¹(yⱼ))ⱼ)</code> in barycentric coordinates
   *
   * \ref success()
   */
  const std::tuple<Images, Preimages>& images() const
    { return m_images; }

  /**
   * \brief Indicate whether projection (inverse projection) is valid for each corner or not.
   *
   * Returns a pair of bitsets. The first bitset indicates if the projection
   * <code>Φ(xᵢ)</code> is valid for each corner <code>xᵢ</code>, that is
   * that <code>Φ(xᵢ)</code> could be computed and lies in the image simplex.
   * The second bitset indicates the same for the inverse projection
   * <code>Φ⁻¹(yⱼ)</code> for the corners <code>yⱼ</code>.
   *
   * \note \ref project() must be called before this method can be used.
   *
   * \returns pair of bitsets indicating success of (inverse) projection at
   *          corners <code>xᵢ</code> (<code>yⱼ</code>)
   */
  const std::tuple<std::bitset<dim>, std::bitset<dim> >& success() const
    { return m_success; }

  /**
   * \brief Number of edge intersections.
   *
   * \note \ref project() must be called before this method can be used.
   *
   * \ref edgeIntersections()
   */
  unsigned numberOfEdgeIntersections() const
    { return m_number_of_edge_intersections; }

  /**
   * \brief Edge-edge intersections.
   *
   * \note \ref project() must be called before this method can be used.
   *
   * \warning Only the first \ref numberOfEdgeIntersections() entries are valid
   *          edge intersections.
   */
  const std::array<EdgeIntersection, maxEdgeIntersections>& edgeIntersections() const
    { return m_edge_intersections; }
};

} /* namespace GridGlue */
} /* namespace Dune */

#include "projection_impl.hh"

#endif
