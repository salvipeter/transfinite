#include "nelder-mead.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace NelderMead {

  void addmul(Point &x, const Point &d, double alpha) {
    for (size_t i = 0; i < x.size(); ++i)
      x[i] += d[i] * alpha;
  }

  void average(Point &x, const Point &d) {
    for (size_t i = 0; i < x.size(); ++i)
      x[i] = (x[i] + d[i]) / 2.0;
  }

  bool optimize(const Function &f, Point &x, size_t max_iteration, double tolerance,
                double step_length) {
    // Create starting simplex
    size_t n = x.size();
    std::vector<Point> S(n + 1, x);
    for (size_t i = 1; i <= n; ++i)
      S[i][i-1] += step_length;

    double alpha = 1.0, beta = 2.0, gamma = 0.5; // Standard values

    double Delta = std::numeric_limits<double>::infinity();
    std::vector<double> y;
    std::transform(S.begin(), S.end(), std::back_inserter(y), f);

    for (size_t iter = 0; iter < max_iteration && Delta >= tolerance; ++iter) {
      // Find the extreme values
      size_t low = 0, high = 0, second = 0;
      Point centroid(n, 0);     // of the face opposite the high point
      for (size_t i = 0; i <= n; ++i) {
        addmul(centroid, S[i], 1.0 / n);
        if (y[i] < y[low])
          low = i;
        else if (y[i] > y[high]) {
          second = high;
          high = i;
        } else if (y[i] > y[second])
          second = i;
      }
      addmul(centroid, S[high], -1.0 / n);

      // Compute the reflected point
      Point xr = centroid;
      addmul(xr, centroid, alpha);
      addmul(xr, S[high], -alpha);
      double yr = f(xr);

      if (yr < y[low]) {
        // Try to expand further
        Point xe = centroid;
        addmul(xe, xr, beta);
        addmul(xe, centroid, -beta);
        double ye = f(xe);
        // Save the better of the two
        if (ye < yr) {
          S[high] = xe;
          y[high] = ye;
        } else {
          S[high] = xr;
          y[high] = yr;
        }
      } else if (yr > y[second]) {
        // Save the reflected value if it is between the highest and second highest
        if (yr <= y[high]) {
          S[high] = xr;
          y[high] = yr;
        }
        // Try to contract just the highest point
        Point xc = centroid;
        addmul(xc, S[high], gamma);
        addmul(xc, centroid, -gamma);
        double yc = f(xc);
        if (yc > y[high]) {
          // When failed, shrink the whole simplex
          for (size_t i = 0; i <= n; ++i) {
            if (i == low)
              continue;
            average(S[i], S[low]);
            y[i] = f(S[i]);
          }
        } else {
          S[high] = xc;
          y[high] = yc;
        }
      } else {
        // When the reflected value is better than the second highest
        S[high] = xr;
        y[high] = yr;
      }

      // Compute the standard deviation of y
      double mean = std::accumulate(y.begin(), y.end(), 0.0) / (n + 1.0);
      Delta = 0.0;
      for (double yi : y)
        Delta += std::pow(yi - mean, 2);
      Delta /= (n + 1);
    }

    // Find the minimum
    size_t min = 0;
    for (size_t i = 1; i <= n; ++i)
      if (y[i] < y[min])
        min = i;
    x = S[min];

    return Delta < tolerance;
  }

} // namespace
