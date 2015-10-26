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
  void minimize(std::vector<double> &data, size_t iterations,
                std::function<double(const std::vector<double> &)> func);

private:
  static double amotry(std::vector<std::vector<double> > &p, std::vector<double> &y,
		       std::vector<double> &psum,
                       std::function<double(const std::vector<double> &)> funk,
                       int ihi, double fac);
  static void amoeba(std::vector<std::vector<double> > &p, std::vector<double> &y, double tol,
		     unsigned int nmax,
                     std::function<double(const std::vector<double> &)> funk);

  // Parameters
  static double const TINY;
  double step;
  double tolerance;
  size_t iterations;
};

#endif // DOWNHILL_H
