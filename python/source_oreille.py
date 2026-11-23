#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
source_oreille.py
Source position at the ear for the HRTF solvers, found automatically on the mesh
(gmsh 2.2 ASCII or binary, or STL binary/ASCII; flat triangles). numpy only.

Assumption: the mesh is ALIGNED on the Oy axis, which goes through both ear
canal entrances (right ear at y < 0, left ear at y > 0); python/diag_msh.py
checks it.

  Ray:  source triangle of main_HRTF_Ray = first triangle crossed by Oy from the
        head centre towards the ear (bottom of the blocked canal). Number from 1,
        in the order of the triangles of the file: value 7 of the case file.
  Diff: monopole of main_HRTF_Diff, between the tragus and the concha: the point
        IN THE AIR closest to the canal entrance E whose exact distance to the
        mesh is >= --dmin + --margin mm (3D grid with step --step mm, depth
        0..--depth mm along the axis, +-rlat mm sideways). Final check: exact
        distance >= dmin and winding number ~ 0 (point in the air).
        Coordinates x y z (m): value 7 of the Diff case file.
        For a SERIES of meshes of the same head, compute the Diff source once on
        the reference mesh and keep it; --check x y z verifies that it stays in
        the air and at least dmin from each mesh.

Usage:
  python source_oreille.py head.msh --ear R
  python source_oreille.py head.msh --ear L --dmin 1 --margin 1
  python source_oreille.py head.msh --check -0.004 -0.072 0.0      (given Diff point)
(--oreille D/G, --verifier, --marge, --pas, --prof: former names, still accepted)
"""

import argparse
import numpy as np
from maillage_msh import (lire_msh2, lire_maillage, intersections_droite, distance_point_maillage,
                          est_exterieur, _distance_point_triangles)


def source_ray(P, T, oreille):
    """Premier triangle coupe par Oy du centre vers l'oreille : (indice 0.., point E)."""
    sens = -1.0 if oreille == 'D' else 1.0
    interieur, w = est_exterieur(P, T, np.zeros(3))
    if interieur:
        raise SystemExit("The origin is not INSIDE the head: the mesh is not aligned on Oy "
                         "(axis through both ear canals).")
    idx, t = intersections_droite(P, T, np.zeros(3), np.array([0.0, 1.0, 0.0]))
    cote = sens * t > 0
    if not cote.any():
        raise SystemExit("The Oy axis does not cross the mesh on this side.")
    idx, t = idx[cote], t[cote]
    ordre = np.argsort(np.abs(t))
    return idx[ordre], t[ordre]


def source_diff(P, T, oreille, E, dmin, marge, pas, rlat, prof_max):
    """Monopole : point de l'AIR le plus proche de l'entree du conduit E dont la
    distance exacte au maillage est >= dmin + marge. Grille 3D (pas `pas`) :
    profondeur s le long de l'axe (0..prof_max), decalage lateral x, z (+-rlat).
    Les plans sont parcourus du conduit vers l'exterieur ; on s'arrete des que s
    depasse la meilleure distance a E trouvee (|X - E| >= s)."""
    sens = -1.0 if oreille == 'D' else 1.0
    R = prof_max + rlat + 0.01
    dv = np.linalg.norm(P - E, axis=1)
    sel = np.nonzero((dv[T] < R).any(axis=1))[0]
    A, B, C = P[T[sel, 0]], P[T[sel, 1]], P[T[sel, 2]]
    N = np.cross(B - A, C - A)
    N /= np.linalg.norm(N, axis=1)[:, None]
    # convention d'orientation : signe de (X - Q).n pour un point loin dans l'air
    Xloin = E + np.array([0.0, sens * 0.5, 0.0])
    _, kloin = distance_point_maillage(P, T, Xloin)
    nloin = np.cross(P[T[kloin, 1]] - P[T[kloin, 0]], P[T[kloin, 2]] - P[T[kloin, 0]])
    signe_air = np.sign(np.dot(Xloin - P[T[kloin]].mean(axis=0), nloin))

    g = np.arange(-rlat, rlat + 1e-12, pas)
    gx, gz = np.meshgrid(g, g, indexing='ij')
    gx, gz = gx.ravel(), gz.ravel()
    meilleur, dbest, profils = None, np.inf, []
    for s in np.arange(0.0, prof_max + 1e-12, pas):
        if s >= dbest:
            break
        y = E[1] + sens * s
        X = np.c_[E[0] + gx, np.full_like(gx, y), E[2] + gz]
        # on ne calcule que les points plus proches de E que le meilleur actuel
        X = X[np.linalg.norm(X - E, axis=1) < dbest]
        if len(X) == 0:
            continue
        dpts = np.full(len(X), np.inf)
        kpts = np.zeros(len(X), dtype=np.int64)
        for j0 in range(0, len(X), 256):
            Xs = X[j0:j0 + 256]
            nX = len(Xs)
            d = _distance_point_triangles(np.repeat(Xs, len(A), axis=0), np.tile(A, (nX, 1)),
                                          np.tile(B, (nX, 1)), np.tile(C, (nX, 1))).reshape(nX, -1)
            k = np.argmin(d, axis=1)
            dpts[j0:j0 + nX] = d[np.arange(nX), k]
            kpts[j0:j0 + nX] = k
        Q = (A[kpts] + B[kpts] + C[kpts]) / 3.0
        air = np.sign(np.einsum('ij,ij->i', X - Q, N[kpts])) == signe_air
        ok = air & (dpts >= dmin + marge)
        profils.append((s, np.max(np.where(air, dpts, -np.inf)), None))
        if ok.any():
            dE = np.linalg.norm(X[ok] - E, axis=1)
            j = np.argmin(dE)
            if dE[j] < dbest:
                dbest, meilleur = dE[j], X[ok][j]
    return meilleur, profils


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('maillage')
    ap.add_argument('--ear', '--oreille', dest='oreille', choices=['R', 'L', 'D', 'G'], default='R',
                    help='R: right ear (y < 0), L: left ear (y > 0) (D/G accepted)')
    ap.add_argument('--dmin', type=float, default=1.0, help='minimum distance to the mesh, mm (default 1; '
                    'beyond ~3.5 mm a notch c/(4d) enters the band)')
    ap.add_argument('--margin', '--marge', dest='marge', type=float, default=1.0,
                    help='margin above dmin to choose the plane, mm')
    ap.add_argument('--step', '--pas', dest='pas', type=float, default=0.5, help='search step, mm')
    ap.add_argument('--rlat', type=float, default=10.0, help='half-width of the lateral search, mm')
    ap.add_argument('--depth', '--prof', dest='prof', type=float, default=15.0,
                    help='maximum depth from the canal entrance, mm')
    ap.add_argument('--check', '--verifier', dest='verifier', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                    help="check a GIVEN Diff source point (m) without moving it: distance to the mesh "
                         "and point in the air")
    a = ap.parse_args()
    a.oreille = {'R': 'D', 'L': 'G'}.get(a.oreille, a.oreille)      # D/G internally
    cote = 'right' if a.oreille == 'D' else 'left'

    P, T = lire_maillage(a.maillage)
    print(f"{a.maillage}: {len(T)} triangles, {len(P)} nodes")

    if a.verifier is not None:
        X = np.array(a.verifier)
        d, _ = distance_point_maillage(P, T, X)
        dehors, w = est_exterieur(P, T, X)
        ok = d >= a.dmin * 1e-3 and dehors
        print(f"point {X[0]:.5f} {X[1]:.5f} {X[2]:.5f} (m): distance to the mesh {d*1e3:.2f} mm; "
              f"in the air: {dehors} (winding number {w:.3f}) -> {'OK' if ok else 'FAILED'}")
        raise SystemExit(0 if ok else 1)

    idx, t = source_ray(P, T, a.oreille)
    E = np.array([0.0, t[0], 0.0])
    print(f"\n--- Ray ({cote} ear) ---")
    print(f"source triangle (from 1): {idx[0] + 1}")
    print(f"  canal entrance on Oy: y = {t[0]*1e3:.2f} mm")
    if len(idx) > 1:
        print("  other crossings of Oy on this side (y mm, triangle):",
              ', '.join(f"{tt*1e3:.2f} ({k+1})" for k, tt in zip(idx[1:], t[1:])))

    X, profils = source_diff(P, T, a.oreille, E, a.dmin * 1e-3, a.marge * 1e-3,
                             a.pas * 1e-3, a.rlat * 1e-3, a.prof * 1e-3)
    print(f"\n--- Diff ({cote} ear) ---")
    if X is None:
        print(f"NO point at >= {a.dmin + a.marge} mm from the mesh up to {a.prof} mm from the canal.")
        print("profile (depth mm, largest distance mm):",
              ', '.join(f"{s*1e3:.1f}:{d*1e3:.2f}" for s, d, _ in profils))
        raise SystemExit(1)
    d, _ = distance_point_maillage(P, T, X)
    dehors, w = est_exterieur(P, T, X)
    print(f"monopole x y z (m): {X[0]:.5f} {X[1]:.5f} {X[2]:.5f}")
    print(f"  {abs(X[1]-E[1])*1e3:.1f} mm from the canal along the axis, offset "
          f"x {X[0]*1e3:+.1f} mm, z {X[2]*1e3:+.1f} mm")
    print(f"  exact distance to the mesh: {d*1e3:.2f} mm; in the air: {dehors} (winding number {w:.3f})")
    if d < a.dmin * 1e-3 or not dehors:
        raise SystemExit("Final check FAILED.")


if __name__ == '__main__':
    main()
