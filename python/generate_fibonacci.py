#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
generate_fibonacci.py
Fibonacci grid of directions on a sphere, for the HRTF solvers.

numpy only. The Fibonacci sphere covers the sphere almost uniformly; for
N >= 1000 the quadrature error of the energy mean (diffuse field) is
negligible compared with a Lebedev grid of the same size. Uniform weights
w_i = 1/N.

Usage:
  python generate_fibonacci.py --n 2702 --dist 3.0
  python generate_fibonacci.py --n 2702 --dist 3.0 --out-dir ../run/

Output (in --out-dir):
  pext_fibonacci_NNNN.txt    N on line 1, then N lines x y z (external points file)
  weights_fibonacci_NNNN.txt uniform weights (1/N), for hrtf2hrir.py

Convention: x front, z up, counter-clockwise.
"""

import argparse
import numpy as np
from pathlib import Path


def fibonacci_sphere(n, r=1.0):
    """
    Génère n points quasi-uniformément répartis sur une sphère de rayon r.
    Algorithme : sphère dorée (golden ratio spiral).

    Retourne (pts, weights) :
      pts     : (n, 3) coordonnées cartésiennes
      weights : (n,)   poids uniformes 1/n, normalisés sum=1
    """
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    i = np.arange(n)

    # Angle polaire (de 0 à π) : répartition uniforme en cos(θ)
    theta = np.arccos(1.0 - 2.0 * (i + 0.5) / n)
    # Angle azimutal : spirale dorée
    phi = 2.0 * np.pi * i / golden

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    pts     = np.column_stack([x, y, z])
    weights = np.full(n, 1.0 / n)
    return pts, weights


def main():
    parser = argparse.ArgumentParser(
        description='Fibonacci grid of directions on a sphere, for the HRTF solvers.')
    parser.add_argument('--n',       type=int,   default=2702,
                        help='number of points (default: 2702)')
    parser.add_argument('--dist',    type=float, default=3.0,
                        help='distance of the points, metres (default: 3.0)')
    parser.add_argument('--out-dir', default='.',
                        help='output directory (default: current directory)')
    args = parser.parse_args()

    pts, weights = fibonacci_sphere(args.n, args.dist)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pext_file   = out_dir / f"pext_fibonacci_{args.n:04d}.txt"
    weight_file = out_dir / f"weights_fibonacci_{args.n:04d}.txt"

    # Format Pext Fortran : première ligne = N, puis N lignes x y z
    with open(pext_file, 'w') as f:
        f.write(f"{args.n}\n")
        np.savetxt(f, pts, fmt='%18.10f')

    np.savetxt(weight_file, weights, fmt='%18.10e')

    print("Fibonacci grid written:")
    print(f"  points   : {args.n}")
    print(f"  distance : {args.dist} m")
    print(f"  sum(w)   : {weights.sum():.10f}")
    print(f"  Pext     : {pext_file}")
    print(f"  weights  : {weight_file}")
    print()
    print("Use with hrtf2hrir.py:")
    print(f"  python hrtf2hrir.py config_L.txt \\")
    print(f"    --weights {weight_file}")


if __name__ == '__main__':
    main()
