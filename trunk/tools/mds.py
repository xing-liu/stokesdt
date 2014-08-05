#!/usr/bin/env python
import sys, re, os, random
from math import *

# Run as ./mds N PHI

nargin = len(sys.argv) - 1
if (nargin < 1):
    raise Exception('Need to specify N')

N = int(sys.argv[1])
if (N <= 0):
    raise Exception('N must be larger than 0')
PHI = 0.1
SEED = 1234

# Bond lengths
if (nargin > 1):
    PHI = float(sys.argv[2])
    if (PHI <= 0.0):
        raise Exception('PHI must be larger than 0.0')

# Geometry
random.seed(SEED);

L = pow(pi*N/3.0/PHI, 1.0/3.0);
pos = [[random.uniform(0.0, L) for _ in xrange(3)] for _ in xrange(N)]
print 'Generated random monodisperse suspension, N = %d, PHI = %.6f, L = %.6f' %(N, PHI, L)

# Write xyz
xyz_file = './N%d_Phi%.6f.xyz' % (N, PHI)
fh = open(xyz_file, 'w')
fh.write('%d\n' % (N))
fh.write('0 0.0 %.16g %.16g %.16g\n' % (L, L, L))
for par in pos:
    fh.write('N %.16g %.16g %.16g\n' %(par[0], par[1], par[2]))
fh.close()
print 'Writed to ' + xyz_file