#include "Eigen/LU"
#include "Eigen/SVD"

#include "domain.hh"

#include "gb-fit.hh"

SurfaceGeneralizedBezier elevateDegree(const SurfaceGeneralizedBezier &surf) {
  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree() + 1;
  SurfaceGeneralizedBezier result;
  result.initNetwork(n, d);
  size_t l = result.layers();

  // Set the corner control points
  for (size_t i = 0; i < n; ++i)
    result.setControlPoint(i, 0, 0, surf.controlPoint(i, 0, 0));

  // Set the boundary control points
  for (size_t j = 1; j < d; ++j) {
    double eta = (double)j / d;
    for (size_t i = 0; i < n; ++i)
      result.setControlPoint(i, j, 0,
                             surf.controlPoint(i, j - 1, 0) * eta +
                             surf.controlPoint(i, j, 0) * (1 - eta));
  }

  // Set all other control points (except the center)
  // When d is odd, the topmost layer of the original patch
  // is simulated by the control points of the adjacent boundary
  for (size_t j = 1; j <= d / 2; ++j) {
    double eta = (double)j / d;
    for (size_t k = 1; k < l; ++k) {
      double theta = (double)k / d;
      for (size_t i = 0; i < n; ++i) {
        Point3D p1, p2, p3, p4;
        p1 = surf.controlPoint(i, j - 1, k - 1);
        p2 = surf.controlPoint(i, j, k - 1);
        if (d % 2 == 0 || k < d / 2) {
          p3 = surf.controlPoint(i, j - 1, k);
          p4 = surf.controlPoint(i, j, k);
        } else {
          size_t im = (i + n - 1) % n;
          p3 = surf.controlPoint(im, d - k - 1, j - 1);
          if (j == d / 2)
            p4 = surf.centralControlPoint();
          else
            p4 = surf.controlPoint(im, d - k - 1, j);
        }
        result.setControlPoint(i, j, k,
                               p1 * eta * theta + p2 * (1 - eta) * theta +
                               p3 * eta * (1 - theta) + p4 * (1 - eta) * (1 - theta));
      }
    }
  }

  // Set center as the mass center of the innermost control points
  Point3D center(0, 0, 0);
  for (size_t i = 0; i < n; ++i)
    center += result.controlPoint(i, l, l - 1);
  center /= n;
  result.setCentralControlPoint(center);

  result.setupLoop();
  return result;
}

Point2DVector parameterizePoints(const Surface &surf, const PointVector &points) {
  size_t resolution = 15;
  TriMesh mesh = surf.eval(resolution);
  PointVector vertices = mesh.points();
  const Point2DVector &params = surf.domain()->parameters(resolution);
  Point2DVector result; result.reserve(points.size());
  for (const auto &p : points) {
    TriMesh::Triangle tri = mesh.closestTriangle(p);
    const Point3D &a = vertices[tri[0]], &b = vertices[tri[1]], &c = vertices[tri[2]];
    Vector3D n = ((b - a) ^ (c - a)).normalize();
    Point3D q = p + n * ((a - p) * n);
    double x = ((b - q) ^ (c - q)).norm();
    double y = ((a - q) ^ (c - q)).norm();
    double z = ((a - q) ^ (b - q)).norm();
    Point2D uv = (params[tri[0]] * x + params[tri[1]] * y + params[tri[2]] * z) / (x + y + z);
    result.push_back(uv);
  }
  return result;
}

SurfaceGeneralizedBezier fitWithOriginal(const SurfaceGeneralizedBezier &original,
                                         const PointVector &points,
                                         const Point2DVector &params,
                                         double smoothing, size_t fixed_rows) {
  size_t n = original.domain()->vertices().size();
  size_t d = original.degree();
  size_t l = original.layers();
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;                                // # of control points
  size_t bcp = n * (d + 1 - fixed_rows) * fixed_rows; // # of boundary control points
  size_t mcp = cp - bcp;                              // # of movable control points
  size_t m = points.size();                           // # of samples

  // Initialize patch with the fixed boundaries
  SurfaceGeneralizedBezier surf;
  surf.initNetwork(n, d);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < d - 1; ++j)
      for (size_t k = 0; k < 2; ++k)
        surf.setControlPoint(i, j, k, original.controlPoint(i, j, k));
  surf.setupLoop();

  Eigen::MatrixXd A(m + mcp, mcp);
  Eigen::MatrixXd b(m + mcp, 3);
  A.setZero();
  b.setZero();

  // Fill the matrices
  for (size_t j = 0; j < m; ++j) {
    double blend_sum = 0.0;
    for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
      if (col >= d - row) {
        if (++side >= n) {
          side = 0;
          ++row;
        }
        col = row;
      }
      double blend = surf.weight(side, col, row, params[j]);
      if (col < l)
        blend += surf.weight((side + n - 1) % n, d - row, col, params[j]);
      if (col > d - l)
        blend += surf.weight((side + 1) % n, row, d - col, params[j]);
      if (i > bcp)
        A(j, i - bcp) = blend;
      else {
        Point3D p = surf.controlPoint(side, col, row);
        b(j, 0) -= p[0] * blend; b(j, 1) -= p[1] * blend; b(j, 2) -= p[2] * blend;
      }
      blend_sum += blend;
    }
    A(j, 0) = 1.0 - blend_sum;
    b(j, 0) += points[j][0]; b(j, 1) += points[j][1]; b(j, 2) += points[j][2];
  }

  // Smoothing terms
  auto findControlPoint = [n,d,cp](size_t i, size_t j, size_t k) -> size_t {
    for (size_t c = 1, side = 0, col = 0, row = 0; c < cp; ++c, ++col) {
      if (col >= d - row) {
        if (++side >= n) {
          side = 0;
          ++row;
        }
        col = row;
      }
      size_t side_m = (side + n - 1) % n, side_p = (side + 1) % n;
      if ((i == side && j == col && k == row) ||
          (i == side_m && j == d - row && k == col) ||
          (i == side_p && j == row && k == d - col))
        return c;
    }
    return 0;
  };
  auto setWeight = [&A,&b,bcp,surf,findControlPoint]
    (size_t j, size_t side, size_t col, size_t row, double w) {
    size_t i = findControlPoint(side, col, row);
    if (i > bcp)
      A(j, i - bcp) = w;
    else if (i == 0)
      A(j, 0) = w;
    else {
      Point3D p = surf.controlPoint(side, col, row);
      b(j, 0) -= p[0] * w; b(j, 1) -= p[1] * w; b(j, 2) -= p[2] * w;
    }
  };
  for (size_t j = 0; j < n; ++j)
    setWeight(m, j, l, l - 1, -smoothing / n);
  A(m, 0) = smoothing;
  for (size_t i = bcp + 1, side = 0, col = fixed_rows, row = fixed_rows; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    size_t j = m + i - bcp;
    setWeight(j, side, col - 1, row, -smoothing / 4.0);
    setWeight(j, side, col + 1, row, -smoothing / 4.0);
    setWeight(j, side, col, row - 1, -smoothing / 4.0);
    setWeight(j, side, col, row + 1, -smoothing / 4.0);
    A(j, i - bcp) = smoothing;
  }

  // LSQ Fit
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  if (!(A*x).isApprox(b)) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    x = svd.solve(b);
  }

  // Fill control points
  surf.setCentralControlPoint(Point3D(x(0, 0), x(0, 1), x(0, 2)));
  for (size_t i = bcp + 1, side = 0, col = fixed_rows, row = fixed_rows; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    surf.setControlPoint(side, col, row, Point3D(x(i - bcp, 0), x(i - bcp, 1), x(i - bcp, 2)));
  }

  return surf;
}
