#pragma once

#include <list>

#include "geometry.hh"

#include "BSCF_FN_Curvature.hh"
#include "BSCF_FN_LsqDistance.hh"
#include "BSC_Fitter.hh"

#include "BSSF_FN_Curvature.hh"
#include "BSSF_FN_LsqDistance.hh"
#include "BSSF_FN_Oscillation.hh"
#include "BSSF_KI_LargestSummedDev.hh"
#include "BSS_Fitter.hh"

class CurveFitter::CurveFitterImpl {
public:
  void setTolerance(double tol) { tolerance_ = tol; }
  void setDegree(size_t deg) { f_.Degree() = deg; }
  void setNrControlPoints(size_t nr_cpts) { f_.NrControlPoints() = nr_cpts; }
  void setKnotVector(const DoubleVector &knots) {
    f_.theKnotVector().assign(knots.begin(), knots.end());
  }
  void addControlPoint(size_t i, const Point3D &point) {
    cpts_.push_back(ControlPointConstraint(i, point));
  }
  void addParamPoint(double param, const Point3D &point) {
    ppts_.push_back(ParameterPoint(param, point));
  }
  void fit() {
    if (!ppts_.empty()) {
      BSC_Fitter<Point<3, double>>::ParamPointSetReference point_set = f_.ParamPointSet();
      point_set.push_back(BSC_Fitter<Point<3, double>>::ParamPointGroupType(tolerance_));
      BSC_Fitter<Point<3, double>>::ParamPointGroupReference pg(point_set.front());
      for (const auto &pp : ppts_) {
        Point<1, double> param(pp.u);
        Point<3, double> point(pp.p[0], pp.p[1], pp.p[2]);
        pg.push_back(BSC_Fitter<Point<3, double>>::ParamPointType(point, param));
      }
    }
    if (!cpts_.empty()) {
      size_t nr_ctrl = f_.NrControlPoints();
      BSC_Fitter<Point<3, double>>::ConstraintsType cst(nr_ctrl);
      for (const auto &cp : cpts_)
        cst[cp.i] = Point<3, double>(cp.p[0], cp.p[1], cp.p[2]);
      std::swap(f_.CtrlConstraints(), cst);
    }
    f_.Enlargement() = 0.0;
    f_.AddFunctional(new BSCF_FN_LsqDistance<Point<3, double>>());
    f_.AddFunctional(new BSCF_FN_Curvature<Point<3, double>>());
    f_.OptimizeParameters() = false;
    f_.Fit();
  }
  BSCurve curve() const {
    size_t deg = f_.Curve().Degree();
    const BSplineCurve<Point<3, double>>::KnotVectorType &knots1 = f_.Curve().theKnotVector();
    const BSplineCurve<Point<3, double>>::ControlVectorType &cpts1 = f_.Curve().theControlVector();
    DoubleVector knots(knots1.begin(), knots1.end());
    PointVector cpts;
    size_t n = cpts1.size();
    cpts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      const Point<3, double> &p = cpts1[i];
      cpts.push_back(Point3D(p[0], p[1], p[2]));
    }
    return BSCurve(deg, knots, cpts);
  }
private:
  struct ParameterPoint {
    ParameterPoint(double u, const Point3D &p) : u(u), p(p) {}
    double u;
    Point3D p;
  };
  struct ControlPointConstraint {
    ControlPointConstraint(size_t i, const Point3D &p) : i(i), p(p) {}
    size_t i;
    Point3D p;
  };
  BSC_Fitter<Point<3, double>> f_;
  double tolerance_;
  std::list<ParameterPoint> ppts_;
  std::list<ControlPointConstraint> cpts_;
};

class SurfaceFitter::SurfaceFitterImpl {
public:
  SurfaceFitterImpl() : curvature_weight_(0.0), oscillation_weight_(0.0) {}
  void setTolerance(double tol) { tolerance_ = tol; }
  void setDegreeU(size_t deg_u) { f_.DegreeU() = deg_u; }
  void setDegreeV(size_t deg_v) { f_.DegreeV() = deg_v; }
  void setNrControlPointsU(size_t nr_cpts_u) { f_.NrControlPointsU() = nr_cpts_u; }
  void setNrControlPointsV(size_t nr_cpts_v) { f_.NrControlPointsV() = nr_cpts_v; }
  void setKnotVectorU(const DoubleVector &knots_u) {
    f_.KnotVectorU().assign(knots_u.begin(), knots_u.end());
  }
  void setKnotVectorV(const DoubleVector &knots_v) {
    f_.KnotVectorV().assign(knots_v.begin(), knots_v.end());
  }
  void setMaxNrControlPointsU(size_t max_cpts_u) { f_.MaximalNrControlPointsU() = max_cpts_u; }
  void setMaxNrControlPointsV(size_t max_cpts_v) { f_.MaximalNrControlPointsV() = max_cpts_v; }
  void setCurvatureWeight(double weight) { curvature_weight_ = weight; }
  void setOscillationWeight(double weight) { oscillation_weight_ = weight; }
  void addControlPoint(size_t i, size_t j, const Point3D &point) {
    cpts_.push_back(ControlPointConstraint(i, j, point));
  }
  void addParamPoint(const Point2D &param, const Point3D &point) {
    ppts_.push_back(ParameterPoint(param, point));
  }
  void fit() {
    finalizeSetup();
    f_.Enlargement() = 0.0;
    f_.OptimizeParameters() = false;
    f_.Fit();
  }
  void fitWithCarrierSurface() {
    finalizeSetup();
    f_.LocalOutlierPercentages() = true;
    f_.OptimizeParameters() = true;
    f_.SetKnotInserter(new BSSF_KI_LargestSummedDev<Point<3, double>>());
    f_.CarrierFit();
    f_.Fit();
  }
  BSSurface surface() const {
    BSSurface result;
    result.deg_u_ = f_.Surface().DegreeU();
    result.deg_v_ = f_.Surface().DegreeV();
    const BSplineSurface<Point<3, double>>::KnotVectorType &knots_u = f_.Surface().KnotVectorU();
    result.knots_u_.assign(knots_u.begin(), knots_u.end());
    const BSplineSurface<Point<3, double>>::KnotVectorType &knots_v = f_.Surface().KnotVectorV();
    result.knots_v_.assign(knots_v.begin(), knots_v.end());
    const BSplineSurface<Point<3, double>>::ControlNetType &cnet = f_.Surface().theControlNet();
    size_t n = cnet.sizes()[0], m = cnet.sizes()[1];
    result.cnet_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      result.cnet_[i].reserve(m);
      for (size_t j = 0; j < m; ++j) {
        const Point<3, double> &p = cnet[i][j];
        result.cnet_[i].push_back(Point3D(p[0], p[1], p[2]));
      }
    }
    return result;
  }
protected:
  void finalizeSetup() {
    if (!ppts_.empty()) {
      BSS_Fitter<Point<3, double>>::ParamPointSetReference point_set = f_.ParamPointSet();
      point_set.push_back(BSS_Fitter<Point<3, double>>::ParamPointGroupType(tolerance_));
      BSS_Fitter<Point<3, double>>::ParamPointGroupReference pg(point_set.front());
      for (const auto &pp : ppts_) {
        Point<2, double> param(pp.uv[0], pp.uv[1]);
        Point<3, double> point(pp.p[0], pp.p[1], pp.p[2]);
        pg.push_back(BSS_Fitter<Point<3, double>>::ParamPointType(point, param));
      }
    }
    if (!cpts_.empty()) {
      size_t nr_ctrl[2];
      nr_ctrl[0] = f_.NrControlPointsU(); nr_ctrl[1] = f_.NrControlPointsV();
      BSS_Fitter<Point<3, double>>::ConstraintsType cst(nr_ctrl);
      for (const auto &cp : cpts_)
        cst[cp.i][cp.j] = Point<3, double>(cp.p[0], cp.p[1], cp.p[2]);
      std::swap(f_.CtrlConstraints(), cst);
    }
    f_.AddFunctional(new BSSF_FN_LsqDistance<Point<3, double>>());
    if (curvature_weight_ > 0.0) {
      using CurvatureFn = BSSF_FN_Curvature<Point<3, double>>;
      f_.AddFunctional(new CurvatureFn(CurvatureFn::UV, 1.0, curvature_weight_));
    }
    if (oscillation_weight_ > 0.0)
      f_.AddFunctional(new BSSF_FN_Oscillation<Point<3, double>>(1.0, oscillation_weight_));
  }
private:
  struct ParameterPoint {
    ParameterPoint(const Point2D &uv, const Point3D &p) : uv(uv), p(p) {}
    Point2D uv;
    Point3D p;
  };
  struct ControlPointConstraint {
    ControlPointConstraint(size_t i, size_t j, const Point3D &p) : i(i), j(j), p(p) {}
    size_t i, j;
    Point3D p;
  };
  BSS_Fitter<Point<3, double>> f_;
  double tolerance_, curvature_weight_, oscillation_weight_;
  std::list<ParameterPoint> ppts_;
  std::list<ControlPointConstraint> cpts_;
};
