#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""Reverses the orientation (normals) of all surface triangles of a mesh.

  python msh_flip.py in.msh out.msh      gmsh 2.2 ASCII or binary, gmsh 4 ASCII
  python msh_flip.py in.stl out.stl      STL, binary or ASCII

The HRTF solvers need normals pointing into the air (out of the head);
python/diag_msh.py tells which way they point. The ORDER of the triangles is
kept (the Ray source triangle number stays valid) and everything else in the
file is copied unchanged, except that a binary gmsh 2.2 file is written as
gmsh 2.2 ASCII.

T3 (gmsh type 2):  A B C -> A C B
T6 (gmsh type 9):  A B C mAB mBC mCA -> A C B mCA mBC mAB
STL facet:         vertices 2 and 3 swapped, facet normal reversed
"""
import os
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def retourne(typ, noeuds):
    """Reversed node list of a surface element (None if not a T3 / T6)."""
    if typ == 2:
        a, b, c = noeuds
        return [a, c, b]
    if typ == 9:
        a, b, c, ab, bc, ca = noeuds
        return [a, c, b, ca, bc, ab]
    return None


def flip_msh_ascii(L):
    """In place on the lines of an ASCII gmsh file (2.2 or 4). Returns (n T3, n T6)."""
    i = 0
    while i < len(L) and L[i].strip() != '$MeshFormat':
        i += 1
    version = float(L[i + 1].split()[0])
    while i < len(L) and L[i].strip() != '$Elements':
        i += 1
    if i >= len(L):
        sys.exit('no $Elements section')
    n = {2: 0, 9: 0}
    if version < 3.0:                     # id type ntags tags... nodes
        nb = int(L[i + 1]); i += 2
        for k in range(nb):
            v = L[i + k].split()
            typ, nt = int(v[1]), int(v[2])
            r = retourne(typ, v[3 + nt:])
            if r is not None:
                L[i + k] = ' '.join(v[:3 + nt] + r)
                n[typ] += 1
    else:                                 # blocks: dim tag type n, then id nodes...
        nb = int(L[i + 1].split()[0]); i += 2
        for _ in range(nb):
            typ = int(L[i].split()[2]); ne = int(L[i].split()[3]); i += 1
            for k in range(ne):
                v = L[i + k].split()
                r = retourne(typ, v[1:])
                if r is not None:
                    L[i + k] = ' '.join([v[0]] + r)
                    n[typ] += 1
            i += ne
    return n[2], n[9]


def flip_msh_binaire(fin, fout):
    """gmsh 2.2 binary -> gmsh 2.2 ASCII, same nodes, elements and tags, surfaces reversed."""
    from maillage_msh import lire_msh2_binaire_brut
    ids, X, el = lire_msh2_binaire_brut(fin)
    n = {2: 0, 9: 0}
    with open(fout, 'w') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n')
        f.write(f'{len(ids)}\n')
        for k, i in enumerate(ids):
            f.write(f'{int(i)} {X[k, 0]:.17g} {X[k, 1]:.17g} {X[k, 2]:.17g}\n')
        f.write('$EndNodes\n$Elements\n')
        f.write(f'{len(el)}\n')
        for k, (typ, tags, v) in enumerate(el, start=1):
            r = retourne(typ, list(v))
            if r is not None:
                v = r
                n[typ] += 1
            f.write(' '.join(str(int(x)) for x in [k, typ, len(tags), *tags, *v]) + '\n')
        f.write('$EndElements\n')
    return n[2], n[9]


def flip_stl(fin, fout):
    """STL binary or ASCII: vertices 2 and 3 swapped, facet normal reversed. Returns n."""
    data = open(fin, 'rb').read()
    if len(data) >= 84:
        nt = struct.unpack('<I', data[80:84])[0]
        if len(data) == 84 + 50 * nt:                 # binary
            out = bytearray(data[:84])
            for k in range(nt):
                f = data[84 + 50 * k: 84 + 50 * (k + 1)]
                nx, ny, nz = struct.unpack('<3f', f[0:12])
                out += struct.pack('<3f', -nx, -ny, -nz) + f[12:24] + f[36:48] + f[24:36] + f[48:50]
            open(fout, 'wb').write(out)
            return nt
    L = data.decode(errors='replace').split('\n')     # ASCII
    nt, i = 0, 0
    while i < len(L):
        s = L[i].split()
        if len(s) >= 5 and s[0] == 'facet' and s[1] == 'normal':
            ind = L[i][:len(L[i]) - len(L[i].lstrip())]
            L[i] = f"{ind}facet normal {-float(s[2]):.9g} {-float(s[3]):.9g} {-float(s[4]):.9g}"
            v = [j for j in range(i + 1, min(i + 6, len(L))) if L[j].split()[:1] == ['vertex']]
            if len(v) != 3:
                sys.exit(f'ASCII STL: facet without 3 vertices near line {i + 1}')
            L[v[1]], L[v[2]] = L[v[2]], L[v[1]]
            nt += 1
            i = v[2]
        i += 1
    open(fout, 'w').write('\n'.join(L))
    return nt


def main(fin, fout):
    if fin.lower().endswith('.stl'):
        if not fout.lower().endswith('.stl'):
            sys.exit('an STL input gives an STL output')
        print(f'  {flip_stl(fin, fout)} facets reversed -> {fout}')
        return
    from maillage_msh import est_msh_binaire
    if est_msh_binaire(fin):
        n3, n6 = flip_msh_binaire(fin, fout)
        print(f'  {n3} T3 and {n6} T6 reversed -> {fout} (gmsh 2.2 ASCII)')
        return
    L = open(fin).read().split('\n')
    n3, n6 = flip_msh_ascii(L)
    open(fout, 'w').write('\n'.join(L))
    print(f'  {n3} T3 and {n6} T6 reversed -> {fout}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
