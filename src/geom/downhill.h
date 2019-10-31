#ifndef DOWNHILL_H
#define DOWNHILL_H

#include <cstddef>
#include <functional>
#include <vector>

class DownhillSimplex
{
public:
  DownhillSimplex() : step(1.0), tolerance(0.0) { }
  void setStepLength(double len) { step = len; }
  void setTolerance(double tol) { tolerance = tol; }
  bool minimize(std::vector<double> &data, size_t iterations,
                std::function<double(const std::vector<double> &)> func);

private:
  static double amotry(std::vector<std::vector<double> > &p, std::vector<double> &y,
		       std::vector<double> &psum,
                       std::function<double(const std::vector<double> &)> funk,
                       size_t ihi, double fac);
  static bool amoeba(std::vector<std::vector<double> > &p, std::vector<double> &y, double tol,
		     size_t nmax, std::function<double(const std::vector<double> &)> funk);

  // Parameters
  static double const TINY;
  double step;
  double tolerance;
};

#endif // DOWNHILL_H
