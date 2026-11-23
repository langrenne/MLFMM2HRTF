#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""Checks a mesh BEFORE an HRTF computation (the solvers do not detect these errors).

  python diag_msh.py head.msh          gmsh 2.2 ASCII or binary, gmsh 4 ASCII
  python diag_msh.py head.stl          STL, binary or ASCII
  python diag_msh.py tube.msh --tube   + checks of the tube solver (lines, axis, groups)

Reported:
  - element types: the HRTF solvers need flat triangles (T3);
  - closed surface: every edge shared by exactly two triangles;
  - extent: coordinates in metres (an STL in millimetres is converted by the
    solvers, a .msh is not);
  - edge lengths: the octree is limited to 10 levels, avoid edges below ~0.5 mm;
  - triangle quality: minimum angle (very flat triangles degrade the
    Burton-Miller / CHI equations);
  - orientation: signed volume (divergence theorem), V > 0 means normals
    pointing out of the head, into the air, as required;
  - ear axis: the Oy axis must go through both ear canal entrances (right
    ear at y < 0), with the origin inside the head (python/source_oreille.py
    then finds the sources).
Exit code 0 if no problem was found, 1 otherwise.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

TYP = {1: ("line 2n", 2), 2: ("T3", 3), 3: ("Q4", 4),
       8: ("line 3n", 3), 9: ("T6", 6), 15: ("point", 1)}
SURF = (2, 3, 9)        # surface elements
LIGN = (1, 8)           # line elements (tube solver: edges of the tangential Laplacian)


def lire_stl(nom):
    """STL (binary or ASCII) through maillage_msh.lire_stl: merged vertices,
    triangles in the order of the facets, millimetres converted like the solvers."""
    from maillage_msh import lire_stl as lecture
    X, T = lecture(nom)
    P = {i + 1: X[i] for i in range(len(X))}
    elem = [(0, 2, [int(a) + 1, int(b) + 1, int(c) + 1]) for a, b, c in T]
    return 'STL', P, elem, {}


def lire(nom):
    """Returns (format, nodes {id: xyz}, elements [(physical, type, nodes)], physical names)."""
    if nom.lower().endswith('.stl'):
        return lire_stl(nom)
    from maillage_msh import est_msh_binaire, lire_msh2_binaire_brut
    if est_msh_binaire(nom):
        ids, X, el = lire_msh2_binaire_brut(nom)
        P = {int(i): X[k] for k, i in enumerate(ids)}
        elem = [(t[0] if len(t) else 0, typ, v) for typ, t, v in el]
        return 'gmsh 2.2 binary', P, elem, {}
    L = open(nom, errors='replace').read().split('\n')
    i = 0
    while i < len(L) and L[i].strip() != '$MeshFormat':
        i += 1
    if i >= len(L):
        sys.exit(f"  {nom}: no $MeshFormat, not a gmsh file")
    version = float(L[i + 1].split()[0])
    P, elem, phys_nom = {}, [], {}

    def trouve(bal):
        j = 0
        while j < len(L) and L[j].strip() != bal:
            j += 1
        return j + 1 if j < len(L) else None

    j = trouve('$PhysicalNames')
    if j:
        for k in range(int(L[j])):
            v = L[j + 1 + k].split(None, 2)
            phys_nom[int(v[1])] = v[2].strip('"') if len(v) > 2 else ''

    if version < 3.0:
        j = trouve('$Nodes')
        for k in range(int(L[j])):
            v = L[j + 1 + k].split()
            P[int(v[0])] = np.array([float(x) for x in v[1:4]])
        j = trouve('$Elements')
        for k in range(int(L[j])):
            v = [int(float(x)) for x in L[j + 1 + k].split()]
            typ, nt = v[1], v[2]
            elem.append((v[3] if nt >= 1 else 0, typ, v[3 + nt:]))
    else:
        ent = {}                                  # entity -> physical tag (else -entity tag)
        j = trouve('$Entities')
        if j:
            npo, nc, ns, nv = map(int, L[j].split())
            k = j + 1 + npo
            for _ in range(nc):
                v = L[k].split(); nph = int(float(v[7]))
                ent[(1, int(v[0]))] = int(float(v[8])) if nph > 0 else -int(v[0])
                k += 1
            for _ in range(ns):
                v = L[k].split(); nph = int(float(v[7]))
                ent[(2, int(v[0]))] = int(float(v[8])) if nph > 0 else -int(v[0])
                k += 1
        j = trouve('$Nodes')
        nb = int(L[j].split()[0]); k = j + 1
        for _ in range(nb):
            d, t, par, n = map(int, L[k].split()); k += 1
            tg = [int(L[k + m]) for m in range(n)]; k += n
            for m in range(n):
                P[tg[m]] = np.array([float(x) for x in L[k + m].split()])
            k += n
        j = trouve('$Elements')
        nb = int(L[j].split()[0]); k = j + 1
        for _ in range(nb):
            d, t, typ, n = map(int, L[k].split()); k += 1
            ph = ent.get((d, t), 0)
            for m in range(n):
                v = [int(x) for x in L[k + m].split()]
                elem.append((ph, typ, v[1:]))
            k += n
    return f'gmsh {version:g} ASCII', P, elem, phys_nom


def ligne(cle, texte):
    print(f"  {cle:<12}: {texte}")


def main(nom, tube=False):
    fmt, P, elem, phys_nom = lire(nom)
    surf = [e for e in elem if e[1] in SURF]
    lign = [e for e in elem if e[1] in LIGN]
    if not surf:
        sys.exit("  no surface elements")
    problemes = []

    ids = np.array(sorted(P))
    pos = {n: i for i, n in enumerate(ids)}
    X = np.array([P[n] for n in ids])
    T = np.array([[pos[v[0]], pos[v[1]], pos[v[2]]] for _, t, v in surf if t in (2, 9)], dtype=np.int64)

    typs = sorted({e[1] for e in surf})
    print()
    ligne('file', f"{nom} ({fmt})")
    ligne('nodes', f"{len(P)}")
    ligne('surface', f"{len(surf)} elements, types {', '.join(TYP[t][0] for t in typs)}")
    if typs != [2]:
        problemes.append("element types")
        print("               PROBLEM: the HRTF solvers need flat triangles (T3) only")

    # --- closed surface: each edge shared by exactly two triangles ---
    if len(T):
        ar = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
        _, inv, cnt = np.unique(ar, axis=0, return_inverse=True, return_counts=True)
        libres, multiples = int((cnt == 1).sum()), int((cnt > 2).sum())
        if libres == 0 and multiples == 0:
            ligne('closed', "yes (every edge shared by two triangles)")
        else:
            problemes.append("open surface")
            k = np.nonzero(cnt[inv.ravel()] != 2)[0][0] % len(T)
            ligne('closed', f"NO: {libres} free edges, {multiples} edges shared by more than two "
                            f"triangles (e.g. near {np.round(X[T[k]].mean(0), 4)} m)")
            print("               PROBLEM: close the surface (holes, duplicated or T-junction nodes)")

    # --- extent and units ---
    ext = X.max(0) - X.min(0)
    ligne('extent', f"x {ext[0]:.4f}  y {ext[1]:.4f}  z {ext[2]:.4f} m")
    if ext.max() > 2.0:
        problemes.append("units")
        print("               PROBLEM: more than 2 m: coordinates probably in millimetres. The solvers")
        print("               convert an STL automatically, but a .msh must be in metres")

    # --- edge lengths ---
    if len(T):
        h = np.concatenate([np.linalg.norm(X[T[:, a]] - X[T[:, (a + 1) % 3]], axis=1) for a in range(3)])
        ligne('edges', f"{h.min()*1e3:.2f} to {h.max()*1e3:.2f} mm (median {np.median(h)*1e3:.2f} mm)")
        if h.min() < 0.5e-3:
            print("               WARNING: edges below 0.5 mm: the octree is limited to 10 levels")

    # --- triangle quality: minimum angle ---
    if len(T):
        A = [X[T[:, i]] for i in range(3)]
        ang = []
        for i in range(3):
            u = A[(i + 1) % 3] - A[i]; w = A[(i + 2) % 3] - A[i]
            c = (u * w).sum(1) / np.linalg.norm(u, axis=1) / np.linalg.norm(w, axis=1)
            ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
        amin = np.min(ang, axis=0); i0 = int(np.argmin(amin))
        ligne('min angle', f"{amin.min():.1f} deg ({(amin < 20).sum()} triangles below 20 deg, "
                           f"{(amin < 30).sum()} below 30 deg), worst near {np.round(X[T[i0]].mean(0), 4)} m")
        if amin.min() < 10:
            print("               WARNING: very flat triangles degrade the Burton-Miller / CHI equations;")
            print("               remesh this region if possible")

    # --- orientation: signed volume (divergence theorem) ---
    G = X[T].mean(axis=(0, 1))
    V = np.einsum('ij,ij->i', X[T[:, 0]] - G, np.cross(X[T[:, 1]] - G, X[T[:, 2]] - G)).sum() / 6.0
    if V > 0:
        ligne('normals', f"out of the enclosed volume, into the air: OK (volume {V*1e3:.3f} litres)")
    else:
        problemes.append("normals")
        ligne('normals', f"INTO the enclosed volume (volume {V*1e3:.3f} litres)")
        print("               PROBLEM: the HRTF solvers need normals pointing into the air:")
        print("               python/msh_flip.py in.msh out.msh (or .stl)")

    # --- ear axis Oy ---
    if not tube and len(T) and typs == [2]:
        from maillage_msh import intersections_droite, est_exterieur
        dehors, _ = est_exterieur(X, T, np.zeros(3))
        _, t = intersections_droite(X, T, np.zeros(3), np.array([0.0, 1.0, 0.0]))
        yd, yg = t[t < 0], t[t > 0]
        if dehors or len(yd) == 0 or len(yg) == 0:
            problemes.append("ear axis")
            ligne('ear axis', "PROBLEM: the origin must be inside the head and the Oy axis must cross "
                              "the mesh on both sides")
        else:
            ligne('ear axis', f"origin inside; Oy enters the mesh at y = {yd.max()*1e3:.1f} mm (right) "
                              f"and {yg.min()*1e3:+.1f} mm (left)")
            print("               the canal entrances must be there: check, then python/source_oreille.py")

    if tube:
        diag_tube(fmt, P, surf, lign, phys_nom, X)

    print()
    if problemes:
        print(f"  summary     : PROBLEMS: {', '.join(problemes)}")
        return 1
    print("  summary     : no problem found")
    return 0


def diag_tube(fmt, P, surf, lign, phys_nom, X):
    """Checks of the tube solver (main_beltrami_tube): line elements, axis, physical groups."""
    print("\n  --- tube solver ---")
    if lign:
        tl = sorted({e[1] for e in lign})
        ligne('lines', f"{len(lign)} elements, types {[TYP[t][0] for t in tl]}: edge term of the "
                       "tangential Laplacian possible")
    else:
        ligne('lines', "NONE: mesh%NB = 0, the edge term will do nothing")
    if fmt == 'STL':
        return
    ext = X.max(0) - X.min(0)
    ax = int(np.argmax(ext)); lab = "xyz"[ax]
    ligne('main axis', lab)
    if ax != 0:
        print("               WARNING: the tube solver selects the piston by sqrt(x^2+y^2+z^2) <= rayon_P,")
        print(f"               which assumes the axis on x and the piston at x = 0; here the axis is {lab}")
    if all(ph <= 0 for ph, t, v in surf):
        print("  no physical group: lecture_gmsh falls back to geometric ENTITIES (tags below)")
    print(f"\n  {'group':>5} {'faces':>7} {'type':>5}  {'extent on axis':>18}  {'radius':>16}  shape")
    print("  " + "-" * 72)
    grp = {}
    for ph, t, v in surf:
        grp.setdefault(ph, []).append((t, v))
    for ph in sorted(grp, key=lambda x: abs(x)):
        d = grp[ph]
        idx = sorted({n for t, v in d for n in v[:3]})
        C = np.array([P[n] for n in idx])
        a0, a1 = C[:, ax].min(), C[:, ax].max()
        rad = np.delete(C, ax, axis=1)
        r = np.hypot(rad[:, 0], rad[:, 1])
        forme = f"DISC at {lab}={a0:+.4f}" if (a1 - a0) < 1e-9 else "extended surface"
        nom_ph = phys_nom.get(ph, '')
        print(f"  {abs(ph):5d} {len(d):7d} {TYP[d[0][0]][0]:>5}  "
              f"[{a0:+.4f},{a1:+.4f}]  [{r.min():.4f},{r.max():.4f}]  {forme}"
              + (f"  ({nom_ph})" if nom_ph else ""))
    print("\n  The tube solver expects: physical 1 = wall, physical 3 = piston (disc at x = 0).")
    print("  Normals: closed tube (surface enclosing the air) V > 0; open radiating tube")
    print("  (surface enclosing the wall) V < 0. In both cases the normal leaves the fluid.")


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if a != '--tube']
    if len(args) != 1:
        sys.exit(__doc__)
    sys.exit(main(args[0], tube='--tube' in sys.argv))
