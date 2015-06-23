#include "domain.hh"

namespace Transfinite
{

  Domain::
  Domain(Surface *surface) : surface(surface)
  {
  }

  Domain::
  ~Domain()
  {
  }

  // Computes everything from vertices
  // except for parameters, which are computed only when needed
  void Domain::
  invalidate()
  {
    n = vertices.size();
    parameters.clear();
    du.resize(n); dv.resize(n);
    for(size_t i = 0; i < n; ++i) {
      du[i] = vertices[i] - vertices[prev(i)];
      dv[i] = Vector2D(du[i][1], -du[i][0]);
      if((center() - vertices[i]) * dv[i] < 0)
        dv[i] = -dv[i];
    }
    local_vertices.resize(n);
    for(size_t i = 0; i < n; ++i) {
      local_vertices[i].resize(n);
      for(size_t j = 0; j < n; ++j)
        local_vertices[i][j] = toLocal(i, vertices[j]);
    }
  }

  Point2DVector const &Domain::
  verticesGlobal() const
  {
    return vertices;
  }

  Point2DVector const &Domain::
  verticesLocal(size_t i) const
  {
    return local_vertices[i];
  }

  Point2D Domain::
  toLocal(size_t i, Point2D const &p) const
  {
    Vector2D const v = p - vertices[prev(i)];
    if(fabs(dv[i][1]) < epsilon) {
      double const a = (v[1] - v[0] * dv[i][1] / dv[i][0]) /
        (du[i][1] - du[i][0] * dv[i][1] / dv[i][0]);
      return Point2D(v[0] - a * du[i][0] - dv[i][0], a);
    }
    double const a = (v[0] - v[1] * dv[i][0] / dv[i][1]) /
      (du[i][0] - du[i][1] * dv[i][0] / dv[i][1]);
    return Point2D(a, v[1] - a * du[i][1] - dv[i][1]);
  }

  Point2D Domain::
  toGlobal(size_t i, Point2D const &p) const
  {
    return vertices[prev(i)] + du[i] * p[0] + dv[i] * p[1];
  }

  Point2DVector Domain::
  globalParameters(size_t resolution) const
  {
    if(parameters.empty()) {
      parameters.reserve(1 + n * resolution * (resolution + 1) / 2);
      Point2D const c = center();
      parameters.push_back(c);
      for(size_t j = 1; j <= resolution; ++j) {
        double const u = (double)j / resolution;
        for(size_t k = 0; k < n; ++k)
          for(size_t i = 0; i < j; ++i) {
            Point2D const p(Vector2D(c) * (1.0 - u) + (vertices[k] - c) * u);
            parameters.push_back(p);
          }
      }
    }
    return parameters;
  }

  // Builds the topology and also resizes the pointvector of the mesh
  Mesh Domain::
  meshTopology(size_t resolution) const
  {
    Mesh mesh;
    mesh.Points().resize(1 + n * resolution * (resolution + 1) / 2);

    size_t inner_start = 0, outer_vert = 1;
    for(size_t layer = 1; layer <= resolution; ++layer){
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for(size_t side = 0; side < n; ++side) {
        size_t vert = 0;
        while(true) {
          size_t const next_vert = (side == n - 1 && vert == layer - 1
                                    ? outer_start
                                    : outer_vert + 1);
          mesh.AddTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if(++vert == layer)
            break;
          if(layer > 2) {
            size_t const inner_next = (side == n - 1 && vert == layer - 1
                                       ? inner_start
                                       : inner_vert + 1);
            mesh.AddTriangle(inner_vert, next_vert, inner_next);
            inner_vert = inner_next;
          }
        }
      }
      inner_start = outer_start;
    }

    return mesh;
  }

} // Transfinite
