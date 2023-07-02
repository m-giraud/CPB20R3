#include "config.h"

#include <array>
#include <cmath>
#include <iostream>
#include <tuple>

#include <dune/grid-glue/common/projection.hh>
#include <dune/grid-glue/common/projectionwriter.hh>

using Dune::GridGlue::Projection;

bool
test_project_simple()
{
  using std::get;

  bool pass = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
  const Corners c1{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
  const auto& corners = std::tie(c0, c1);

  const Normals n0{{{0, 0, 1}, {0, 0, 1}, {0, 0, 1}}};
  const Normals n1{{{0, 0, -1}, {0, 0, -1}, {0, 0, -1}}};
  const auto& normals = std::tie(n0, n1);

  P p;
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success = get<0>(p.success());
    if (!success.all()) {
      std::cout << "ERROR: test_project_simple: not all Phi(x_i) are inside the image triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
    const auto& images = get<0>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_simple: corner " << i
                  << " was projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (!success.all()) {
      std::cout << "ERROR: test_project_simple: not all Phi^{-1}(y_i) are inside the preimage triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
    const auto& images = get<1>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_simple: corner " << i
                  << " was inverse-projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* No edge intersections. */
  if (p.numberOfEdgeIntersections() != 0) {
    std::cout << "ERROR: test_project_simple: there were unexpected edge intersections" << std::endl;
    pass = false;
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project_simple.vtk");

  return pass;
}

bool
test_project_simple2()
{
  using std::get;
  using std::sqrt;

  bool pass = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
  const Corners c1{{{0, 0, 1}, {2, 0, 1}, {0, 2, 1}}};
  const auto& corners = std::tie(c0, c1);

  Normals n0{{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
  for (auto& n : n0)
    n /= n.two_norm();
  Normals n1{{{0, 0, -1}, {-1, 0, -1}, {0, -1, -1}}};
  for (auto& n : n1)
    n /= n.two_norm();
  const auto& normals = std::tie(n0, n1);

  P p;
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success = get<0>(p.success());
    if (!success.all()) {
      std::cout << "ERROR: test_project_simple2: not all Phi(x_i) are inside the image triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, sqrt(2.)}, {0, 1, sqrt(2.)}}};
    const auto& images = get<0>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_simple2: corner " << i
                  << " was projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (!success.all()) {
      std::cout << "ERROR: test_project_simple: not all Phi^{-1}(y_i) are inside the preimage triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, sqrt(2.)}, {0, 1, sqrt(2.)}}};
    const auto& images = get<1>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_simple2: corner " << i
                  << " was inverse-projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* No edge intersections. */
  if (p.numberOfEdgeIntersections() != 0) {
    std::cout << "ERROR: test_project_simple2: there were unexpected edge intersections" << std::endl;
    pass = false;
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project_simple2.vtk");

  return pass;
}

bool
test_project3()
{
  /* Test with intersting numbers.
   * The result was not checked, but is given to prevent accidential changes... Hopefully it is right.
   */
  using std::get;

  bool pass = true;

  using V = Dune::FieldVector<double, 3>;
  using P = Projection<V>;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{2.85, 4., 3.}, {2.85, 5., 3.}, {1.9, 5., 3.}}};
  const Corners c1{{{2., 4.38091, 3.61115}, {-0.6, 4.38091, 3.61115}, {2., 5.08885, 4.31909}}};
  const auto& corners = std::tie(c0, c1);

  const Normals n0{{{0., 0., 1.}, {0., 0., 1.}, {0., 0., 1.}}};
  const Normals n1{{{0., 0.587785, -0.809017}, {0., 0.62962, -0.776903}, {0., 0.809017, -0.587785}}};
  const auto& normals = std::tie(n0, n1);

  P p;
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success = get<0>(p.success());
    if (success[0] || success[1] || !success[2]) {
      std::cout << "ERROR: test_project3: unexpected Phi(x_i) are in the image triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{
        {-0.326923, -0.538054, 0.23024},
        {-0.326923, 0.874495, 1.23024},
        {0.0384615, 0.874495, 1.23024}
    }};
    const auto& images = get<0>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: corner " << i
                  << " was projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (success.any()) {
      std::cout << "ERROR: test_project3: not all Phi^{-1}(y_i) are outside the preimage triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{
      {-0.513827, 0.894737, 0.494431},
      {-3.25067, 3.63158, 0.474804},
      {0.194113, 0.894737, 0.775341}
    }};
    const auto& images = get<1>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: corner " << i
                  << " was inverse-projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* Two edge intersections. */
  if (p.numberOfEdgeIntersections() != 2) {
    std::cout << "ERROR: test_project3: there were unexpected edge intersections" << std::endl;
    pass = false;
  }
  else {
    {
      const auto& in = p.edgeIntersections()[0];
      if (in.edge[0] != 1 || in.edge[1] != 1) {
        std::cout << "ERROR: test_project3: edge intersection 0 involves wrong edges "
                  << in.edge[0] << " and " << in.edge[1]
                  << "; expected: 1 and 1" << std::endl;
        pass = false;
      }
      const V expected0{0, 0.894737, 1.12498};
      if (!((in.local[0] - expected0).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: edge intersection 0 at local[0] = "
                  << in.local[0] << "; expected: " << expected0 << std::endl;
        pass = false;
      }
      const V expected1{0, 0.725806, 0.736697};
      if (!((in.local[1] - expected1).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: edge intersection 0 at local[1] = "
                  << in.local[1] << "; expected: " << expected1 << std::endl;
        pass = false;
      }
    }
    {
      const auto& in = p.edgeIntersections()[1];
      if (in.edge[0] != 2 || in.edge[1] != 1) {
        std::cout << "ERROR: test_project3: edge intersection 1 involves wrong edges "
                  << in.edge[0] << " and " << in.edge[1]
                  << "; expected: 2 and 1" << std::endl;
        pass = false;
      }
      const V expected0{0.105263, 0.894737, 1.23024};
      if (!((in.local[0] - expected0).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: edge intersection 1 at local[0] = "
                  << in.local[0] << "; expected: " << expected0 << std::endl;
        pass = false;
      }
      const V expected1{0, 0.874495, 0.761376};
      if (!((in.local[1] - expected1).infinity_norm() < 1e-5)) {
        std::cout << "ERROR: test_project3: edge intersection 1 at local[1] = "
                  << in.local[1] << "; expected: " << expected1 << std::endl;
        pass = false;
      }
    }
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project3.vtk");

  return pass;
}

bool
test_project4()
{
  /* Test with intersting numbers.
   * The result was not checked, but is given to prevent accidential changes... Hopefully it is right.
   */
  using std::get;

  bool pass = true;

  using V = Dune::FieldVector<double, 3>;
  using P = Projection<V>;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{3.3226, 4.41838, 2.5606}, {3.56991, 4.60089, 2.67132},{3.58641, 4.42509, 2.56724}}};
  const Corners c1{{{3.53965, -0.313653, 3.39551},{4.06013, -0.311809, 3.39505},{3.54027, -0.434644, 3.61669}}};

  const auto& corners = std::tie(c0, c1);

  const Normals n0{{{-0.00668383, -0.480515, 0.876961},{-0.00186939, -0.440971, 0.897519},{-0.00727575, -0.477334, 0.878692}}};
  const Normals n1{{{0.00296167, -0.872357, -0.48886},{0.0032529,-0.87389, -0.486112},{0.00375284, -0.877906, -0.478818}}};
  const auto& normals = std::tie(n0, n1);

  P p;
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success = get<0>(p.success());
    if (success.any()) {
      std::cout << "ERROR: test_project4: unexpected Phi(x_i) are in the image triangle" << std::endl;
      pass = false;
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (success.any()) {
      std::cout << "ERROR: test_project4: not all Phi^{-1}(y_i) are outside the preimage triangle" << std::endl;
      pass = false;
    }
  }

  if (p.numberOfEdgeIntersections() != 0) {
    std::cout << "ERROR: test_project4: there were unexpected edge intersections" << std::endl;
    pass = false;
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project4.vtk");

  return pass;
}


bool
test_project5()
{
  using std::get;
  using std::sqrt;

  bool pass = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{1.91, 0.82, 4}, {1.71388, 0.831652, 4.01554}, {1.91, 0.972786, 3.90557}}};
  const Corners c1{{{0, 0, 0}, {0, 0, 3}, {0, 1.5, 3}}};
  const auto& corners = std::tie(c0, c1);

  Normals n0{{{0.00368949, -0.588605, -0.808412}, {-0.18944, -0.577784, -0.793901}, {0.00330402, -0.435137, -0.900358}}};
  for (auto& n : n0)
    n /= n.two_norm();

  Normals n1{{{-0.57735, -0.57735, -0.57735}, {-0.408248, -0.408248, 0.816497}, {-0.707107, 0, 0.707107}}};
  for (auto& n : n1)
    n /= n.two_norm();
  const auto& normals = std::tie(n0, n1);

  P p(0.2,-0.1);
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success0 = get<0>(p.success());
    if (success0.any()) {
      std::cout << "ERROR: test_project5: all Phi(x_i) should be outside of the image triangle" << std::endl;
      pass = false;
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (success.any()) {
      std::cout << "ERROR: test_project5: all Phi^{-1}(y_i) should be outside of the preimage triangle" << std::endl;
      pass = false;
    }
  }

  /* No edge intersections. */
  if (p.numberOfEdgeIntersections() != 0) {
    std::cout << "ERROR: test_project5: there were unexpected edge intersections" << std::endl;
    pass = false;
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project_5.vtk");

  return pass;
}

bool
test_project_overlap()
{
  using std::get;

  bool pass = true;

  typedef Dune::FieldVector<double, 3> V;
  typedef Projection<V> P;
  using Corners = std::array<V, 3>;
  using Normals = Corners;

  const Corners c0{{{0, 0, 1}, {1, 0, 1}, {0, 4, -1}}};
  const Corners c1{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
  const auto& corners = std::tie(c0, c1);

  const Normals n0{{{0, 0, -1}, {0, 0, -1}, {0, 0, -1}}};
  const Normals n1{{{0, 0, 1}, {0, 0, 1}, {0, 0, 1}}};
  const auto& normals = std::tie(n0, n1);

  P p(0.5);
  p.project(corners, normals);

  /* Check image */
  {
    const auto& success = get<0>(p.success());
    if (success != std::bitset<3>(0b011)) {
      std::cout << "ERROR: test_project_overlap: forward projection did not succeed/fail as expected" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, 1}, {0, 4, -1}}};
    const auto& images = get<0>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_overlap: corner " << i
                  << " was projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* Check preimage */
  {
    const auto& success = get<1>(p.success());
    if (!success.all()) {
      std::cout << "ERROR: test_project_overlap: not all Phi^{-1}(y_i) are inside the preimage triangle" << std::endl;
      pass = false;
    }

    const Corners expected{{{0, 0, 1}, {1, 0, 1}, {0, 0.25, 0.5}}};
    const auto& images = get<1>(p.images());
    for (unsigned i = 0; i < images.size(); ++i) {
      if (!((images[i] - expected[i]).infinity_norm() < 1e-8)) {
        std::cout << "ERROR: test_project_overlap: corner " << i
                  << " was inverse-projected on " << images[i]
                  << "; expected: " << expected[i] << std::endl;
        pass = false;
      }
    }
  }

  /* No edge intersections. */
  if (p.numberOfEdgeIntersections() != 0) {
    std::cout << "ERROR: test_project_overlap: there were unexpected edge intersections" << std::endl;
    pass = false;
  }

  if (!pass)
    write(p, corners, normals, "projectiontest_project_overlap.vtk");

  return pass;
}

int main()
{
  bool pass = true;

  pass = pass && test_project_simple();
  pass = pass && test_project_simple2();
  pass = pass && test_project3();
  pass = pass && test_project4();
  pass = pass && test_project5();
  pass = pass && test_project_overlap();

  return pass ? 0 : 1;
}
