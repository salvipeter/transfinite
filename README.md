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
    CurveVector curves; // std::vector of (shared pointers to) B-spline curves
    // fill curves ...
    SurfaceMidpointCoons surf;
    surf.setCurves(curves);
    surf.setupLoop();
    surf.update();
    auto p = surf.eval({0.1, 0.3}); // evaluates the surface at a given (2D) domain point
    auto mesh = surf.eval(30);
    // mesh is a triangle mesh with 30 triangles on each sides
 ```
