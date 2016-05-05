#include <array>
#include <iostream>

#include "transfinite-kapi.hh"

// Note:
// The API uses double** for point lists,
// but it is much easier to use 2D arrays (double[][]).
// The two template functions below handle the conversion.

template<typename T, size_t M, size_t N>
T **array2pointer(T a[M][N]) {
  T **result;
  result = new T*[M];
  for (size_t i = 0; i < M; ++i) {
    result[i] = new T[N];
    for (size_t j = 0; j < N; ++j)
      result[i][j] = a[i][j];
  }
  return result;
}

template<typename T>
void delete2Dpointer(T **p, size_t m) {
  for (size_t i = 0; i < m; ++i)
    delete[] p[i];
  delete[] p;
}

void print2DCurve(BSplineCurve *c) {
  int degree;
  degree_B_SPLINE_CURVE(c, &degree);
  std::cout << "Degree: " << degree << std::endl;

  double *knots;
  int knots_length;
  knots_B_SPLINE_CURVE(c, &knots, &knots_length);
  std::cout << "Knots:";
  for (int i = 0; i < knots_length; ++i)
    std::cout << ' ' << knots[i];
  std::cout << std::endl;
  delete[] knots;

  int coeffs_length;
  double **coeffs;
  coeffs_B_SPLINE_CURVE(c, &coeffs, &coeffs_length);
  std::cout << "Control points:" << std::endl;
  for (int i = 0; i < coeffs_length; ++i) {
    std::cout << i << ':';
    for (int j = 0; j < 2; ++j)
      std::cout << ' ' << coeffs[i][j];
    std::cout << std::endl;
  }
  delete[] coeffs;
}

int main(int argc, char **argv) {
  double knots1[] = {0,0,0,0,0.267,0.694,1,1,1,1};
  double coeffs1[][3] = {{53.918,-100,      8.475},
                         {56.646, -92.446,  9.017},
                         {59.028, -71.846,  7.977},
                         {44.362,   7.722,  8.555},
                         {12.857,  11.697,  8.492},
                         { 0.019,  11.841,  8.485}};
  double knots2[] = {0,0,0,0,0.538,1,1,1,1};
  double coeffs2[][3] = {{ 0.019,  11.841,  8.485},
                         { 0.846,  -8.867, -2.257},
                         { 0.145, -52.128,-28.903},
                         { 0.602, -89.791,-24.169},
                         { 0.414,-107.179,-23.856}};
  double knots3[] = {0,0,0,0,0.594,1,1,1,1};
  double coeffs3[][3] = {{ 0.414,-107.179,-23.856},
                         {12.999,-107.324,-24.444},
                         {36.495,-101.455,  0.997},
                         {47.675,-100.196,  6.224},
                         {53.918,-100,      8.475}};
  BSplineCurve *curves[3];
  double **coeffs1_p = array2pointer<double, 6, 3>(coeffs1);
  create_B_SPLINE_CURVE(3, knots1, 10, coeffs1_p, 6, &curves[0]);
  delete2Dpointer<double>(coeffs1_p, 6);
  double **coeffs2_p = array2pointer<double, 5, 3>(coeffs2);
  create_B_SPLINE_CURVE(3, knots2, 9, coeffs2_p, 5, &curves[1]);
  delete2Dpointer<double>(coeffs2_p, 5);
  double **coeffs3_p = array2pointer<double, 5, 3>(coeffs3);
  create_B_SPLINE_CURVE(3, knots3, 9, coeffs3_p, 5, &curves[2]);
  delete2Dpointer<double>(coeffs3_p, 5);

  char type[] = "SB", mode[] = "trim";
  BSplineSurface **surface;
  int length;
  transfinite_surface(curves, 3, type, mode, 0.01, &surface, &length);

  delete_B_SPLINE_CURVE(curves[2]);
  delete_B_SPLINE_CURVE(curves[1]);
  delete_B_SPLINE_CURVE(curves[0]);

  int deg_u, deg_v;
  degrees_B_SPLINE_SURFACE(surface[0], &deg_u, &deg_v);
  std::cout << "Degrees: " << deg_u << ", " << deg_v << std::endl;

  double *knots_u, *knots_v;
  int knots_u_length, knots_v_length;
  knots_B_SPLINE_SURFACE(surface[0], &knots_u, &knots_u_length, &knots_v, &knots_v_length);
  std::cout << "U Knots:";
  for (int i = 0; i < knots_u_length; ++i)
    std::cout << ' ' << knots_u[i];
  std::cout << "\nV Knots:";
  for (int i = 0; i < knots_v_length; ++i)
    std::cout << ' ' << knots_v[i];
  std::cout << std::endl;
  delete[] knots_u;
  delete[] knots_v;

  int n_coeffs_u, n_coeffs_v, coeffs_length;
  double **coeffs;
  n_coeffs_B_SPLINE_SURFACE(surface[0], &n_coeffs_u, &n_coeffs_v);
  coeffs_B_SPLINE_SURFACE(surface[0], &coeffs, &coeffs_length);
  std::cout << "Control net:" << std::endl;
  int index = 0;
  for (int i = 0; i < n_coeffs_u; ++i) {
    for (int j = 0; j < n_coeffs_v; ++j) {
      std::cout << i << ',' << j << ':';
      for (int k = 0; k < 3; ++k)
        std::cout << ' ' << coeffs[index][k];
      std::cout << std::endl;
      ++index;
    }
  }
  delete[] coeffs;

  BSplineCurve **trims;
  int trims_length;
  trim_curves_B_SPLINE_SURFACE(surface[0], &trims, &trims_length);
  for (int i = 0; i < trims_length; ++i) {
    std::cout << "Trim curve " << i + 1 << ':' << std::endl;
    print2DCurve(trims[i]);
    delete_B_SPLINE_CURVE(trims[i]);
  }
  delete[] trims;

  delete_B_SPLINE_SURFACE(surface[0]);
  delete[] surface;
  return 0;
}
