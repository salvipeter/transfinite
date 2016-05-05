#pragma once

struct BSplineCurve;
struct BSplineSurface;

extern "C" {

  // Main function

  long int transfinite_surface(BSplineCurve **curves, int curves_length,
                               char *surface_type,
                               char *split_or_trim,
                               double fit_tolerance,
                               BSplineSurface ***surfaces, int *surfaces_length);

  // BSplineCurve members

  long int create_B_SPLINE_CURVE(int degree,
                                 double *knots, int knots_length,
                                 double **coeffs, int coeffs_length,
                                 BSplineCurve **result);

  long int degree_B_SPLINE_CURVE(BSplineCurve *curve,
                                 int *degree);

  long int knots_B_SPLINE_CURVE(BSplineCurve *curve,
                                double **knots, int *knots_length);

  long int coeffs_B_SPLINE_CURVE(BSplineCurve *curve,
                                 double ***coeffs, int *coeffs_length);

  long int delete_B_SPLINE_CURVE(BSplineCurve *curve);

  // BSplineSurface members

  long int degrees_B_SPLINE_SURFACE(BSplineSurface *surface,
                                    int *deg_u, int *deg_v);

  long int knots_B_SPLINE_SURFACE(BSplineSurface *surface,
                                  double **knots_u, int *knots_u_length,
                                  double **knots_v, int *knots_v_length);

  long int n_coeffs_B_SPLINE_SURFACE(BSplineSurface *surface,
                                     int *n_coeffs_u, int *n_coeffs_v);

  long int coeffs_B_SPLINE_SURFACE(BSplineSurface *surface,
                                   double ***coeffs, int *coeffs_length);

  long int trim_curves_B_SPLINE_SURFACE(BSplineSurface *surface,
                                        BSplineCurve ***curves, int *curves_length);

  long int delete_B_SPLINE_SURFACE(BSplineSurface *surface);

}
