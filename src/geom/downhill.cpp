// based on the algorithm in "Numerical Recipes in C"

#include <cmath>

#include "downhill.h"

double const DownhillSimplex::TINY = 1e-8;

double DownhillSimplex::
amotry(std::vector<std::vector<double> > &p,
       std::vector<double> &y,
       std::vector<double> &psum,
       std::function<double(const std::vector<double> &)> funk,
       size_t ihi, double fac)
{
  size_t const ndim = y.size() - 1;
  double ytry;
  std::vector<double> ptry;

  ptry.reserve(ndim);
  double const fac1 = (1.0 - fac) / (double)ndim;
  double const fac2 = fac1 - fac;
  for(size_t j = 0; j < ndim; ++j)
    ptry.push_back(psum[j] * fac1 - p[ihi][j] * fac2);
  ytry = funk(ptry);

  if(ytry < y[ihi]) {
    y[ihi] = ytry;
    for(size_t j = 0; j < ndim; ++j) {
      psum[j] += ptry[j] - p[ihi][j];
      p[ihi][j] = ptry[j];
    }
  }
    
  return ytry;
}

bool DownhillSimplex::
amoeba(std::vector<std::vector<double> > &p,
       std::vector<double> &y,
       double tol, size_t nmax,
       std::function<double(const std::vector<double> &)> funk)
{
  size_t const mpts = y.size();
  size_t const ndim = mpts - 1;
  size_t ihi, ilo, inhi, nfunk = 0;
  double rtol, ysave, ytry;
  std::vector<double> psum;

  psum.reserve(ndim);
  for(size_t j = 0; j < ndim; ++j) {
    double sum = 0.0;
    for(size_t i = 0; i < mpts; ++i)
      sum += p[i][j];
    psum.push_back(sum);
  }

  for(;;) {
    ilo = 0;
    ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
    for(size_t i = 0; i < mpts; ++i) {
      if(y[i] <= y[ilo])
	ilo = i;
      if(y[i] > y[ihi]) {
	inhi = ihi;
	ihi = i;
      } else if(y[i] > y[inhi] && i != ihi)
	inhi = i;
    }

    rtol = 2.0 * std::abs(y[ihi] - y[ilo]) / 
      (std::abs(y[ihi]) + std::abs(y[ilo]) + TINY);
    if(rtol < tol || nfunk >= nmax) {
      std::swap(y[0], y[ilo]);
      for(size_t i = 0; i < ndim; ++i)
	std::swap(p[0][i], p[ilo][i]);
      return rtol < tol;
    }
    nfunk += 2;

    ytry = amotry(p, y, psum, funk, ihi, -1.0);
    if(ytry <= y[ilo])
      ytry = amotry(p, y, psum, funk, ihi, 2.0);
    else if(ytry >= y[inhi]) {
      ysave = y[ihi];
      ytry = amotry(p, y, psum, funk, ihi, 0.5);
      if(ytry >= ysave) {
	for(size_t i = 0; i < mpts; ++i) {
	  if(i != ilo) {
	    for(size_t j = 0; j < ndim; ++j)
	      p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
	    y[i] = funk(psum);
	  }
	}
	nfunk += ndim;
	for(size_t j = 0; j < ndim; ++j) {
	  double sum = 0.0;
	  for(size_t i = 0; i < mpts; ++i)
	    sum += p[i][j];
	  psum[j] = sum;
	}
      }
    } else
      --nfunk;
  }
}

bool DownhillSimplex::
minimize(std::vector<double> &data, size_t iterations,
         std::function<double(const std::vector<double> &)> func)
{
  size_t n = data.size();

  std::vector<std::vector<double> > p;
  p.resize(n + 1);
  for(size_t i = 0; i <= n; ++i) {
    p[i].reserve(n);
    for(size_t j = 0; j < n; ++j)
      p[i].push_back(data[j]);
    if(i != 0)
      p[i][i-1] += step;
  }

  std::vector<double> y;
  y.reserve(n + 1);
  for(size_t i = 0; i <= n; ++i)
    y.push_back(func(p[i]));

  bool success = amoeba(p, y, tolerance, iterations, func);

  for(size_t i = 0; i < n; ++i)
    data[i] = p[0][i];

  return success;
}
