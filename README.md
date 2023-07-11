# Transfinite library
A library for evaluating transfinite surfaces.

Basic geometry is provided by the [libgeom](http://github.com/salvipeter/libgeom) library (included as a submodule).

Currently the following transfinite surfaces are supported:
- [Side-based](http://doi.org/10.1016/j.cagd.2014.06.006) (SB) patch (a variation of Kato's patch)
- [Corner-based](http://doi.org/10.1016/j.cagd.2014.06.006) (CB) patch (a variation of the Charrot-Gregory patch)
- [Generalized Coons](http://doi.org/10.1016/j.cagd.2014.06.006) (GC) patch
- [Composite Ribbon](http://doi.org/10.1016/j.cagd.2014.06.006) (CR) patch
- [Midpoint](http://arxiv.org/abs/2002.11212) (MP) patch 
- [Midpoint Coons](http://arxiv.org/abs/2002.11212) (MC) patch
- [C0 Coons](http://arxiv.org/abs/2002.11347) patch
- Elastic patch (see M. Sabin: Further transfinite surface developments, 1998)

Some control point based surfaces are also contained in the library:
- [Generalized Bézier](http://doi.org/10.1111/cgf.12833) (GB) patch
- Corner-based variation of the GB patch
- [S-patch](http://doi.org/10.1145/77055.77059)

There are also a couple of discrete surfaces:
- [Harmonic patch](http://doi.org/10.1080/10586458.1993.10504266)
- [Biharmonic patch](http://doi.org/10.1016/j.cag.2016.07.001)

... and also some experimental representations.

Usage example is provided in `src/test/test.cc`, but it follows the pattern

```c++
    CurveVector curves; // std::vector of (shared pointers to) curve objects
    // fill curves ...
    SurfaceMidpointCoons surf;
    surf.setCurves(curves);
    surf.setupLoop();
    surf.update();
    auto p = surf.eval({0.1, 0.3}); // evaluates the surface at a given (2D) domain point
    auto mesh = surf.eval(30);
    // mesh is a triangle mesh with 30 triangles on each sides
 ```

## Complete example

Here is a toy program that defines its own curve type (representing line segments),
and then creates a five-sided patch:

```c++
#include <transfinite/surface-composite-ribbon.hh>

using namespace Transfinite;

class Line : public Curve {
public:
    Line(const Point3D &a, const Point3D &b) : a(a), b(b) { }
    ~Line() { }
    Point3D eval(double u) const { return a * (1 - u) + b * u; }
    Point3D eval(double u, size_t nr_der, VectorVector &der) const {
        // der[0] is the 0th derivative, i.e., the curve point
        der.resize(nr_der + 1);
        der[0] = eval(u);
        if (nr_der > 0)
            der[1] = b - a;
        for (size_t i = 2; i <= nr_der; ++i)
            der[i] = Vector3D(0, 0, 0);
        return der[0]; // returns the curve point
    }
    void reverse() { std::swap(a, b); }
    double arcLength(double from, double to) const {
        return (b - a).norm() * (to - from);
    }
private:
    Point3D a, b;
};

int main(int argc, char **argv) {
    CurveVector curves;

    // Create the boundary curves
    curves.push_back(std::make_shared<Line>(Point3D(0,0,0), Point3D(1,0,0)));
    curves.push_back(std::make_shared<Line>(Point3D(1,0,0), Point3D(2,1,1)));
    curves.push_back(std::make_shared<Line>(Point3D(2,1,1), Point3D(1,2,3)));
    curves.push_back(std::make_shared<Line>(Point3D(1,2,3), Point3D(0,1,0)));
    curves.push_back(std::make_shared<Line>(Point3D(0,1,0), Point3D(0,0,0)));

    SurfaceCompositeRibbon surf;
    surf.setCurves(curves);
    surf.setupLoop();
    surf.update();
    surf.eval(50).writeOBJ("test.obj");
}

```

Note that B-spline curves are implemented out-of-the-box by the `BSplineCurve` class,
which is a wrapper for the `Geometry::BSCurve` class.


