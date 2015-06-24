#include <algorithm>

#include "geometry.hh"

BSCurve::BSCurve() {}

BSCurve::BSCurve(size_t degree, DoubleVector knots, PointVector cpts)
  : p_(degree), n_(cpts.size() - 1), knots_(knots), cp_(cpts)
{
}

size_t
BSCurve::findSpan(double u) const
{
  if(u == knots_[n_+1])
    return n_;
  return (std::upper_bound(knots_.begin() + p_ + 1, knots_.end(), u) - knots_.begin()) - 1;
}

void
BSCurve::basisFunctions(size_t i, double u, DoubleVector &coeff) const
{
  coeff.clear(); coeff.reserve(p_ + 1);
  coeff.push_back(1.0);
  DoubleVector left(p_ + 1), right(p_ + 1);
  for(size_t j = 1; j <= p_; ++j) {
    left[j]  = u - knots_[i+1-j];
    right[j] = knots_[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
}

void
BSCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2, std::vector<PointVector> &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for(size_t i = 0; i <= r; ++i)
    dcp[0].push_back(cp_[r1+i]);
  for(size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p_ - k + 1;
    for(size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp / (knots_[r1+i+p_+1] - knots_[r1+i+k]));
  }
}

void
BSCurve::basisFunctionsAll(size_t i, double u, std::vector<DoubleVector> &coeff) const
{
  coeff.clear(); coeff.resize(p_ + 1);
  coeff[0].push_back(1.0);
  DoubleVector left(p_ + 1), right(p_ + 1);
  for(size_t j = 1; j <= p_; ++j) {
    coeff[j].reserve(j + 1);
    left[j]  = u - knots_[i+1-j];
    right[j] = knots_[i+j] - u;
    double saved = 0.0;
    for(size_t r = 0; r < j; ++r) {
      double tmp = coeff[j-1][r] / (right[r+1] + left[j-r]);
      coeff[j].push_back(saved + tmp * right[r+1]);
      saved = tmp * left[j-r];
    }
    coeff[j].push_back(saved);
  }
}

Point3D
BSCurve::eval(double u) const {
  double const span = findSpan(u);
  DoubleVector coeff; basisFunctions(span, u, coeff);
  Point3D point(0.0, 0.0, 0.0);
  for(size_t i = 0; i <= p_; ++i)
    point += cp_[span - p_ + i] * coeff[i];
  return point;
}

Point3D
BSCurve::eval(double u, size_t nr_der, VectorVector &der) const {
  size_t const du = std::min(nr_der, p_);
  der.clear();
  size_t span = findSpan(u);
  std::vector<DoubleVector> coeff; basisFunctionsAll(span, u, coeff);
  std::vector<PointVector> dcp; derivativeControlPoints(du, span - p_, span, dcp);
  for(size_t k = 0; k <= du; ++k) {
    der.push_back(Vector3D(0.0, 0.0, 0.0));
    for(size_t j = 0; j <= p_ - k; ++j)
      der[k] += dcp[k][j] * coeff[p_-k][j];
  }
  for(size_t k = p_ + 1; k <= nr_der; ++k)
    der.push_back(Vector3D(0.0, 0.0, 0.0));
  return der[0];
}

void
BSCurve::reverse() {
  size_t k = knots_.size();
  DoubleVector new_knots;
  new_knots.reserve(k);

  double curr = knots_.front();
  for (size_t i = 1, j = k - 1; i < k; ++i, --j) {
    new_knots.push_back(curr);
    curr += knots_[j] - knots_[j-1];
  }
  new_knots.push_back(curr);

  knots_ = new_knots;
  std::reverse(cp_.begin(), cp_.end());
}

void
BSCurve::normalize() {
  size_t k = knots_.size();
  double low = knots_.front(), high = knots_.back(), len = high - low;
  for (size_t i = 0; i < k; ++i) {
    knots_[i] = (knots_[i] - low) / len;
  }
}
