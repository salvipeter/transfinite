#!/usr/bin/env python3

import sys
import pytransurf

def readLOP(filename):
    f = open(filename)
    n = int(f.readline())
    curves = []
    for i in range(n):
        degree = int(f.readline())
        xs = f.readline().split()
        k = int(xs[0])
        if len(xs) != k + 1:
            f.close()
            raise RuntimeError('invalid file')
        knots = [float(x) for x in xs[1:]]
        xs = f.readline().split()
        k = int(xs[0])
        if len(xs) != 3 * k + 1:
            f.close()
            raise RuntimeError('invalid file')
        points = []
        for j in range(k):
            p = [float(x) for x in xs[3*j+1:3*j+4]]
            points.append(p)
        curves.append((degree, knots, points))
    f.close()
    return curves

def writeSurface(file, surface):
    (deg_u, deg_v, knots_u, knots_v, points) = surface
    file.write(str(deg_u) + '\n')                    # degree
    file.write(str(len(knots_u)) + ' 0 ')            # length, cyclic
    file.write(' '.join([str(k) for k in knots_u]))  # knots
    file.write('\n')
    file.write(str(deg_v) + '\n')                    # degree
    file.write(str(len(knots_v)) + ' 0 ')            # length, cyclic
    file.write(' '.join([str(k) for k in knots_v]))  # knots
    file.write('\n')
    n_u = len(knots_u) - deg_u - 1
    n_v = len(knots_v) - deg_v - 1
    file.write(str(n_u) + ' ' + str(n_v) + '\n')     # n m
    for p in points:
        file.write(' '.join([str(x) for x in p]))    # points
        file.write('\n')

def writeTrimmedSurface(file, surface):
    # TODO: trimming
    writeSurface(file, surface[0])

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 5:
        print('Usage: ' + sys.argv[0] + ' filename [surface-type] [fit-type] [tolerance]')
        return
    filename = sys.argv[1]
    if filename.endswith('.lop'):
        filename = filename[:-4]
    curves = readLOP(filename + '.lop')
    surface_type = 'SB'
    if len(sys.argv) > 2:
        surface_type = sys.argv[2]
    fit_type = 'trim'
    if len(sys.argv) > 3:
        fit_type = sys.argv[3]
    tolerance = 0.01
    if len(sys.argv) > 4:
        tolerance = double(sys.argv[4])
    result = pytransurf.transfiniteSurface(curves, surface_type, fit_type, tolerance)
    f = open(filename + '-' + surface_type + ' ' + fit_type + '.bss', 'w')
    if fit_type == 'trim':
        f.write('1\n')
        writeTrimmedSurface(f, result)
    else:
        f.write(str(len(result)) + '\n')
        for s in result:
            writeSurface(f, s)
    f.close()

if __name__ == '__main__':
    main()
