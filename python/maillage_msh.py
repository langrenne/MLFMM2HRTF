#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
maillage_msh.py
Outils de geometrie pour les maillages de tete HRTF (gmsh 2.2 ASCII, triangles T3).

Aucune dependance externe (numpy uniquement).

  lire_msh2(fichier)                 -> P (NP,3), T (NE,3) indices 0..NP-1
       Les triangles sont numerotes DANS L'ORDRE DU FICHIER (elements de type 2
       seulement), comme les lit le solveur : l'element k (Fortran, a partir
       de 1) est T[k-1].
  intersections_droite(P, T, O, d)   -> indices des triangles coupes par la
       droite O + t d, et parametres t des points d'intersection
  distance_point_maillage(P, T, X)   -> distance du point X a la surface (exacte,
       point-triangle) et indice du triangle le plus proche
  est_exterieur(P, T, X)             -> True si X est hors du volume ferme par le
       maillage (nombre d'enroulement ~ 0)
"""

import numpy as np


NNOD_GMSH = {1: 2, 2: 3, 3: 4, 4: 4, 8: 3, 9: 6, 15: 1}


def est_msh_binaire(fichier):
    """Vrai si l'en-tete est celui d'un gmsh 2.2 BINAIRE ('2.2 1 8')."""
    with open(fichier, 'rb') as f:
        l1 = f.readline(); l2 = f.readline()
    v = l2.split()
    return l1.strip() == b'$MeshFormat' and len(v) >= 2 and v[1] == b'1'


def lire_msh2_binaire_brut(fichier):
    """gmsh 2.2 binaire (little-endian). Renvoie ids (n), X (n x 3), et la liste des
    elements (type, tags, noeuds) dans l'ordre du fichier."""
    d = open(fichier, 'rb').read()
    pos = 0

    def ligne():
        nonlocal pos
        j = d.index(b'\n', pos); t = d[pos:j]; pos = j + 1
        return t.strip()
    ligne(); ligne()
    if int(np.frombuffer(d[pos:pos + 4], '<i4')[0]) != 1:
        raise SystemExit(f'{fichier}: big-endian binary msh not supported')
    pos += 4; ligne(); ligne()
    while ligne() != b'$Nodes':
        pass
    nn = int(ligne())
    rec = np.frombuffer(d, dtype=np.dtype([('id', '<i4'), ('x', '<f8', 3)]), count=nn, offset=pos)
    pos += 28 * nn; ligne(); ligne(); ligne()
    ne = int(ligne())
    elems, k = [], 0
    while k < ne:
        typ, nb, ntag = np.frombuffer(d[pos:pos + 12], '<i4'); pos += 12
        m = 1 + ntag + NNOD_GMSH[int(typ)]
        blk = np.frombuffer(d, '<i4', count=nb * m, offset=pos).reshape(nb, m); pos += 4 * nb * m
        for r in blk:
            elems.append((int(typ), r[1:1 + ntag].tolist(), r[1 + ntag:].tolist()))
        k += nb
    return rec['id'].astype(np.int64), rec['x'].astype(float), elems


def lire_msh2(fichier):
    """Lit un maillage gmsh 2.2 (ASCII ou binaire) ; ne garde que les triangles (type 2)."""
    if est_msh_binaire(fichier):
        ids, X, elems = lire_msh2_binaire_brut(fichier)
        pos = np.full(ids.max() + 1, -1, dtype=np.int64)
        pos[ids] = np.arange(len(ids))
        T = pos[np.array([e[2] for e in elems if e[0] == 2], dtype=np.int64)]
        return X, T
    with open(fichier) as f:
        lignes = f.read().split('\n')
    i = lignes.index('$Nodes')
    n = int(lignes[i + 1])
    tab = np.array([l.split() for l in lignes[i + 2:i + 2 + n]], dtype=float)
    ids = tab[:, 0].astype(np.int64)
    P = tab[:, 1:4]
    pos = np.full(ids.max() + 1, -1, dtype=np.int64)
    pos[ids] = np.arange(n)
    i = lignes.index('$Elements')
    m = int(lignes[i + 1])
    T = []
    for l in lignes[i + 2:i + 2 + m]:
        v = l.split()
        if v[1] == '2':
            nt = int(v[2])
            T.append([int(x) for x in v[3 + nt:6 + nt]])
    T = pos[np.array(T, dtype=np.int64)]
    return P, T


def lire_stl(fichier):
    """Lit un STL (binaire ou ASCII, detecte par la taille du fichier) ; triangles dans
    l'ordre des facettes (numero Ray = rang de la facette, a partir de 1). Les sommets
    repetes sont fusionnes (coordonnees identiques). Comme le solveur (lecture_stl) :
    si la plus grande etendue depasse 2, fichier suppose en mm et converti en m."""
    import os
    taille = os.path.getsize(fichier)
    with open(fichier, 'rb') as f:
        tete = f.read(84)
    n = int(np.frombuffer(tete[80:84], dtype='<u4')[0]) if len(tete) == 84 else -1
    if taille == 84 + 50 * n:
        rec = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
        t = np.fromfile(fichier, dtype=rec, offset=84, count=n)
        V = t['v'].reshape(-1, 3).astype(float)
    else:
        V = np.array([l.split()[1:4] for l in open(fichier, errors='ignore')
                      if l.lstrip().startswith('vertex')], dtype=float)
    P, inv = np.unique(V, axis=0, return_inverse=True)
    T = inv.reshape(-1, 3).astype(np.int64)
    if (P.max(axis=0) - P.min(axis=0)).max() > 2.0:
        P = P * 1e-3
    return P, T


def lire_maillage(fichier):
    """gmsh 2.2 ASCII (.msh) ou STL (.stl) selon l'extension."""
    if fichier.lower().endswith('.stl'):
        return lire_stl(fichier)
    return lire_msh2(fichier)


def intersections_droite(P, T, O, d, eps=1e-12):
    """Moller-Trumbore vectorise : triangles coupes par la droite O + t d (t reel)."""
    A, B, C = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    e1, e2 = B - A, C - A
    h = np.cross(d, e2)
    a = np.einsum('ij,ij->i', e1, h)
    ok = np.abs(a) > eps
    f = np.zeros_like(a)
    f[ok] = 1.0 / a[ok]
    s = O - A
    u = f * np.einsum('ij,ij->i', s, h)
    q = np.cross(s, e1)
    v = f * (q @ d)
    t = f * np.einsum('ij,ij->i', e2, q)
    coupe = ok & (u >= 0) & (v >= 0) & (u + v <= 1)
    idx = np.nonzero(coupe)[0]
    return idx, t[idx]


def _distance_point_triangles(X, A, B, C):
    """Distance exacte de X aux triangles (A, B, C) (Ericson, Real-Time Collision
    Detection, 5.1.5), vectorisee."""
    ab, ac, ax = B - A, C - A, X - A
    d1 = np.einsum('ij,ij->i', ab, ax)
    d2 = np.einsum('ij,ij->i', ac, ax)
    bx = X - B
    d3 = np.einsum('ij,ij->i', ab, bx)
    d4 = np.einsum('ij,ij->i', ac, bx)
    cx = X - C
    d5 = np.einsum('ij,ij->i', ab, cx)
    d6 = np.einsum('ij,ij->i', ac, cx)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    Q = np.empty_like(A)
    # interieur (par defaut)
    den = va + vb + vc
    den[den == 0] = 1e-300
    v = vb / den
    w = vc / den
    Q[:] = A + ab * v[:, None] + ac * w[:, None]
    # aretes et sommets (ordre inverse de priorite)
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    wbc = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1e-300, (d4 - d3) + (d5 - d6))
    Q[m] = B[m] + (C[m] - B[m]) * wbc[m, None]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    wac = d2 / np.where(d2 - d6 == 0, 1e-300, d2 - d6)
    Q[m] = A[m] + ac[m] * wac[m, None]
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    vab = d1 / np.where(d1 - d3 == 0, 1e-300, d1 - d3)
    Q[m] = A[m] + ab[m] * vab[m, None]
    m = (d6 >= 0) & (d5 <= d6)
    Q[m] = C[m]
    m = (d3 >= 0) & (d4 <= d3)
    Q[m] = B[m]
    m = (d1 <= 0) & (d2 <= 0)
    Q[m] = A[m]
    return np.linalg.norm(X - Q, axis=1)


def distance_point_maillage(P, T, X, rayon=0.03):
    """Distance exacte du point X a la surface. Ne teste que les triangles dont un
    sommet est a moins de `rayon` (m) de X (elargi si aucun)."""
    X = np.asarray(X, dtype=float)
    dv = np.linalg.norm(P - X, axis=1)
    r = rayon
    while True:
        sel = np.nonzero((dv[T] < r).any(axis=1))[0]
        if len(sel) > 0:
            break
        r *= 2
    A, B, C = P[T[sel, 0]], P[T[sel, 1]], P[T[sel, 2]]
    dist = _distance_point_triangles(np.broadcast_to(X, A.shape), A, B, C)
    k = np.argmin(dist)
    # un triangle hors de la boule pourrait etre plus proche que r : verification
    if dist[k] > r:
        return distance_point_maillage(P, T, X, rayon=2 * dist[k])
    return dist[k], sel[k]


def est_exterieur(P, T, X):
    """Nombre d'enroulement (somme des angles solides signes / 4 pi) : ~0 dehors,
    ~ +-1 dedans (signe selon l'orientation des normales)."""
    a = P[T[:, 0]] - X
    b = P[T[:, 1]] - X
    c = P[T[:, 2]] - X
    la, lb, lc = (np.linalg.norm(v, axis=1) for v in (a, b, c))
    num = np.einsum('ij,ij->i', a, np.cross(b, c))
    den = la * lb * lc + np.einsum('ij,ij->i', a, b) * lc + np.einsum('ij,ij->i', a, c) * lb \
        + np.einsum('ij,ij->i', b, c) * la
    w = np.sum(2.0 * np.arctan2(num, den)) / (4.0 * np.pi)
    return abs(w) < 0.5, w
